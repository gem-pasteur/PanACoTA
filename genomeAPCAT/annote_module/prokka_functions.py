#!/usr/bin/env python3
# coding: utf-8

"""
Functions to deal with prokka


@author gem
April 2017
"""

import os
import sys
import shutil
import glob
import logging
import subprocess
import shlex
import multiprocessing
import progressbar
import threading

import genomeAPCAT.utils as utils


def run_prokka_all(genomes, threads, force, prok_folder, quiet=False):
    """
    For each genome in genomes, run prokka to annotate the genome.

    Parameters
    ----------
    genomes : dict
        {genome: [name, gpath_cut_gembase, size, nbcont, l90]}
    threads : int
        max number of threads that can be used
    force : bool
        if False, do not override prokka outdir and result dir if they exist. If\
        True, rerun prokka and override existing results, for all genomes.
    prok_folder : str
        folder where prokka results must be written: for each genome,\
        a directory <genome_name>-prokkaRes will be created in this folder, and all the results\
        of prokka for the genome will be written inside
    quiet : bool
        True if nothing must be written to stderr/stdout, False otherwise

    Returns
    -------
    dict
        {genome: boolean} -> with True if prokka ran well, False otherwise.
    """
    main_logger = logging.getLogger("qc_annote.prokka")
    main_logger.info("Annotating all genomes with prokka")
    nbgen = len(genomes)
    bar = None
    if not quiet:
        # Create progressbar
        widgets = ['Annotation: ', progressbar.Bar(marker='█', left='', right=''),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ') - ', progressbar.Timer(), ' - '
                   ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen,
                                      term_width=100).start()
    if threads <= 3:
        cores_prokka = threads
        pool_size = 1
    else:
        # use multiprocessing
        # if there are more threads than genomes, use as many threads as possible per genome
        if len(genomes) <= threads:
            cores_prokka = int(threads / len(genomes))
        # otherwise, use 2 threads per genome (and nb_genome/2 genomes at the same time)
        else:
            cores_prokka = 2
        pool_size = int(threads / cores_prokka)
    pool = multiprocessing.Pool(pool_size)
    # Create a Queue to put logs from processes, and handle them after from a single thread
    m = multiprocessing.Manager()
    q = m.Queue()
    # arguments : (gpath, cores_prokka, name, force, nbcont, q) for each genome
    arguments = [(genomes[g][1], prok_folder, cores_prokka, genomes[g][0],
                  force, genomes[g][3], q)
                 for g in sorted(genomes)]
    try:
        final = pool.map_async(run_prokka, arguments, chunksize=1)
        pool.close()
        # Listen for logs in processes
        lp = threading.Thread(target=utils.logger_thread, args=(q,))
        lp.start()
        if not quiet:
            while True:
                if final.ready():
                    break
                remaining = final._number_left
                bar.update(nbgen - remaining)
            bar.finish()
        pool.join()
        q.put(None)
        lp.join()
        final = final.get()
    # If an error occurs, terminate pool and exit
    except Exception as excp:  # pragma: no cover
        pool.terminate()
        main_logger.error(excp)
        sys.exit(1)
    final = {genome: res for genome, res in zip(sorted(genomes), final)}
    return final


def run_prokka(arguments):
    """
    Run prokka for the given genome.

    Parameters
    ----------
    arguments : tuple
        (gpath, prok_folder, cores_prokka, name, force, nbcont, q) with:

        * gpath: path and filename of genome to annotate
        * prok_folder: path to folder where all prokka folders for all genomes are saved
        * cores_prokka: how many cores can use prokka
        * name: output name of annotated genome
        * force: True if force run (override existing files), False otherwise
        * nbcont: number of contigs in the input genome, to check prokka results
        * q : queue where logs are put

    Returns
    -------
    boolean
        True if eveything went well (all needed output files present,
        corresponding numbers of proteins, genes etc.). False otherwise.
    """
    gpath, prok_folder, threads, name, force, nbcont, q = arguments
    # Set logger for this process
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('run_prokka')
    logger.log(utils.detail_lvl(), "Start annotating {} {}".format(name, gpath))
    # Define prokka directory and logfile, and check their existence
    prok_dir = os.path.join(prok_folder, os.path.basename(gpath) + "-prokkaRes")
    fnull = open(os.devnull, 'w')
    prok_logfile = os.path.join(prok_folder, os.path.basename(gpath) + "-prokka.log")
    if os.path.isdir(prok_dir) and not force:
        logger.warning(("Prokka results folder already exists.").format(prok_dir))
        ok = check_prokka(prok_dir, prok_logfile, name, gpath, nbcont, logger)
        if ok:
            logger.log(utils.detail_lvl(), "Prokka did not run again, "
                                           "formatting step used already generated results of "
                                           "Prokka in {}. If you want to re-run prokka, first "
                                           "remove this result folder, or use '-F' or '--force' "
                                           "option if you want to rerun prokka for all genomes.")
            logger.log(utils.detail_lvl(), "End annotating {} {}".format(name, gpath))
        else:
            logger.warning("Problems in the files contained in your already existing output dir "
                           "({}). Please check it, or remove it to "
                           "re-annotate.".format(prok_dir))
        return ok
    elif os.path.isdir(prok_dir) and force:
        shutil.rmtree(prok_dir)
        logger.debug("Prokka results folder already exists, but removed because --force option "
                     "used")
    cmd = ("prokka --outdir {} --cpus {} "
           "--prefix {} {}").format(prok_dir, threads, name, gpath)
    prokf = open(prok_logfile, "w")
    ok = None
    try:
        subprocess.call(shlex.split(cmd), stdout=fnull, stderr=prokf)
        ok = check_prokka(prok_dir, prok_logfile, name, gpath, nbcont, logger)
        prokf.close()
        if ok:
            logger.log(utils.detail_lvl(), "End annotating {} {}".format(name, gpath))
        return ok
    except Exception as err:  # pragma: no cover
        logger.error("Error while trying to run prokka: {}".format(err))
        prokf.close()
        if ok:
            logger.log(utils.detail_lvl(), "End annotating {} {}".format(name, gpath))
        return False


def check_prokka(outdir, logf, name, gpath, nbcont, logger):
    """
    Prokka writes everything to stderr, and always returns a non-zero return code. So, we
    check if it ran well by checking the content of output directory.

    Parameters
    ----------
    outdir : str
        output directory, where all files are written by prokka
    logf : str
        prokka logfile, containing stderr of prokka (prokka prints everything to stderr)
    name : str
        genome name in gembase format
    gpath : str
        path to fasta file given as input for prokka
    nbcont : int
        number of contigs in fasta file given to prokka
    logger : logging.Logger
        logger object to get logs

    Returns
    -------
    bool
        True if everything went well, False otherwise
    """
    missing_file = False
    problem = False
    if not os.path.isdir(outdir):
        logger.error(("Prokka could not run properly. "
                      "Look at {} for more information.").format(logf))
        missing_file = True
    else:
        oriname = os.path.basename(gpath)
        tblfile = glob.glob(os.path.join(outdir, "*.tbl"))
        faafile = glob.glob(os.path.join(outdir, "*.faa"))
        ffnfile = glob.glob(os.path.join(outdir, "*.ffn"))
        gfffile = glob.glob(os.path.join(outdir, "*.gff"))
        if len(tblfile) == 0:
            logger.error("{} {}: no .tbl file".format(name, oriname))
            missing_file = True
        elif len(tblfile) > 1:
            logger.error("{} {}: several .tbl files".format(name, oriname))
            missing_file = True
        if len(faafile) == 0:
            logger.error("{} {}: no .faa file".format(name, oriname))
            missing_file = True
        elif len(faafile) > 1:
            logger.error("{} {}: several .faa files".format(name, oriname))
            missing_file = True
        if len(ffnfile) == 0:
            logger.error("{} {}: no .ffn file".format(name, oriname))
            missing_file = True
        elif len(ffnfile) > 1:
            logger.error("{} {}: several .ffn files".format(name, oriname))
            missing_file = True
        if len(gfffile) == 0:
            logger.error("{} {}: no .gff file".format(name, oriname))
            missing_file = True
        elif len(gfffile) > 1:
            logger.error("{} {}: several .gff files".format(name, oriname))
            missing_file = True
        if not missing_file:
            tblfile = tblfile[0]
            faafile = faafile[0]
            ffnfile = ffnfile[0]
            fnbcont, tnb_cds, nb_gene, tnb_crispr = count_tbl(tblfile)
            faaprot = count_headers(faafile)
            ffngene = count_headers(ffnfile)
            if nbcont != fnbcont:
                logger.error(("{} {}: no matching number of contigs; "
                              "nbcontig={}; in tbl ={}").format(name, oriname, nbcont, fnbcont))
                problem = True
            if tnb_cds != faaprot:
                logger.error(("{} {}: no matching number of proteins between tbl and faa; "
                              "faa={}; in tbl ={}").format(name, oriname, faaprot, tnb_cds))
                problem = True
            if nb_gene + tnb_crispr != ffngene and nb_gene != ffngene:
                logger.error(("{} {}: no matching number of genes between tbl and ffn; "
                              "ffn={}; in tbl ={}genes {}CRISPR").format(name, oriname,
                                                                         ffngene, nb_gene,
                                                                         tnb_crispr))
                problem = True
    return not problem and not missing_file


def count_tbl(tblfile):
    """
    Count the different features found in the tbl file:

    - number of contigs
    - number of proteins (CDS)
    - number of genes (locus_tag)
    - number of CRISPR arrays (repeat_region)

    Parameters
    ----------
    tblfile : str
        tbl file generated by prokka

    Returns
    -------
    (nbcont, nb_cds, nb_gene, nb_crispr)
        information on features found in the tbl file.
    """
    nbcont = 0
    nb_cds = 0
    nb_gene = 0
    nb_crispr = 0
    with open(tblfile) as tblf:
        for line in tblf:
            if line.startswith(">"):
                nbcont += 1
            if "CDS" in line:
                nb_cds += 1
            if "locus_tag" in line:
                nb_gene += 1
            if "repeat_region" in line:
                nb_crispr += 1
    return nbcont, nb_cds, nb_gene, nb_crispr


def count_headers(seqfile):
    """
    Count how many sequences there are in the given file

    Parameters
    ----------
    seqfile : str
        file containing a sequence in multi-fasta format

    Returns
    -------
    int
        number of contigs in the given multi-fasta file
    """
    nbseq = 0
    with open(seqfile) as faaf:
        for line in faaf:
            if line.startswith(">"):
                nbseq += 1
    return nbseq
