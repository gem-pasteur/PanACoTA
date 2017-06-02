#!/usr/bin/env python3
# coding: utf-8

"""
Functions to deal with prokka


@author gem
April 2017
"""

import os
import shutil
import glob
import logging
import subprocess
import shlex
import multiprocessing
import sys
import progressbar

logger = logging.getLogger()


def run_prokka_all(genomes, threads, force, prok_folder):
    """
    for each genome in genomes, run prokka to annotate the genome.

    * genomes = {genome: [name, gpath_cut_gembase, size, nbcont, l90]}
    * threads: max number of threads that can be used
    * force: if False, do not override prokka outdir and result dir if they exist. If
    True, rerun prokka and override existing results, for all genomes.
    * prok_dir: folder where prokka results must be written: for each genome,
    a directory <genome_name>-prokkaRes will be created in this folder, and all the results
    of prokka for the genome will be written inside

    Returns:
        final: {genome: boolean} -> with True if prokka ran well, False otherwise.
    """
    logger.info("Annotating all genomes with prokka")
    nbgen = len(genomes)
    # Create progressbar
    widgets = ['Annotation: ', progressbar.Bar(marker='█', left='', right='', fill=' '),
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
            cores_prokka = int(threads/len(genomes))
        # otherwise, use 2 threads per genome (and nb_genome/2 genomes at the same time)
        else:
            cores_prokka = 2
        pool_size = int(threads/cores_prokka)
    # arguments : (gpath, cores_prokka, name, force, nbcont) for each genome
    arguments = [(genomes[g][1], prok_folder, cores_prokka, genomes[g][0],
                  force, genomes[g][3])
                 for g in sorted(genomes)]
    pool = multiprocessing.Pool(pool_size)
    try:
        final = pool.map_async(run_prokka, arguments, chunksize=1)
        pool.close()
        while(True):
            if final.ready():
                break
            remaining = final._number_left
            bar.update(nbgen - remaining)
        bar.finish()
        pool.join()
        final = final.get()
    # If an error occurs, terminate pool and exit
    except Exception as excp:  # pragma: no cover
        pool.terminate()
        logger.error(excp)
        sys.exit(1)
    final = {genome: res for genome, res in zip(sorted(genomes), final)}
    return final


def run_prokka(arguments):
    """
    arguments : (gpath, prok_folder, cores_prokka, name, force, nbcont)

    gpath: path and filename of genome to annotate
    prok_dir: path to folder where prokka results must be written
    prok_folder: path to folder where all prokka folders for all genomes are saved
    cores_prokka: how many cores can use prokka
    name: output name of annotated genome
    force: True if force run (override existing files), False otherwise
    nbcont: number of contigs in the input genome, to check prokka results

    returns:
        boolean. True if eveything went well (all needed output files present,
        corresponding numbers of proteins, genes etc.). False otherwise.
    """
    gpath, prok_folder, threads, name, force, nbcont = arguments
    prok_dir = os.path.join(prok_folder, os.path.basename(gpath) + "-prokkaRes")
    FNULL = open(os.devnull, 'w')
    prok_logfile = os.path.join(prok_folder, os.path.basename(gpath) + "-prokka.log")
    if os.path.isdir(prok_dir) and not force:
        logging.warning(("Prokka results folder already exists. Prokka did not run again, "
                         "formatting step used already generated results of Prokka in "
                         "{}. If you want to re-run prokka, first remove this result folder, or "
                         "use '-F' or '--force' option if you want to rerun prokka for "
                         "all genomes.").format(prok_dir))
        ok = check_prokka(prok_dir, prok_logfile, name, gpath, nbcont)
        return ok
    elif os.path.isdir(prok_dir) and force:
        shutil.rmtree(prok_dir)
    cmd = ("prokka --outdir {} --cpus {} "
           "--prefix {} {}").format(prok_dir, threads, name, gpath)
    # logger.debug(cmd)
    prokf = open(prok_logfile, "w")
    try:
        retcode = subprocess.call(shlex.split(cmd), stdout=FNULL, stderr=prokf)
        ok = check_prokka(prok_dir, prok_logfile, name, gpath, nbcont)
        prokf.close()
        return ok
    except Exception as err:  # pragma: no cover
        logging.error("Error while trying to run prokka: {}".format(err))
        prokf.close()
        return False


def check_prokka(outdir, logf, name, gpath, nbcont):
    """
    Prokka writes everything to stderr, and always returns a non-zero return code. So, we
    check if it ran well by checking the content of output directory.

    * outdir: output directory, where all files are written by prokka
    * logf: prokka logfile, containing stderr of prokka (prokka prints everything to stderr)
    * name: genome name in gembase format
    * gpath: path to fasta file given as input for prokka
    * nbcont: number of contigs in fasta file given to prokka
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
        if len(tblfile) == 0:
            logger.error(("{} {}: no .tbl file").format(name, oriname))
            missing_file = True
        elif len(tblfile) > 1:
            logger.error(("{} {}: several .tbl files").format(name, oriname))
            missing_file = True
        if len(faafile) == 0:
            logger.error(("{} {}: no .faa file").format(name, oriname))
            missing_file = True
        elif len(faafile) > 1:
            logger.error(("{} {}: several .faa files").format(name, oriname))
            missing_file = True
        if len(ffnfile) == 0:
            logger.error(("{} {}: no .ffn file").format(name, oriname))
            missing_file = True
        elif len(ffnfile) > 1:
            logger.error(("{} {}: several .ffn files").format(name, oriname))
            missing_file = True
        if not missing_file:
            tblfile = tblfile[0]
            faafile = faafile[0]
            ffnfile = ffnfile[0]
            fnbcont, tnbCDS, tnbGene, tnbCRISPR = count_tbl(tblfile)
            faaprot = count_headers(faafile)
            ffngene = count_headers(ffnfile)
            if nbcont != fnbcont:
                logger.error(("{} {}: no matching number of contigs; "
                              "nbcontig={}; in tbl ={}").format(name, oriname, nbcont, fnbcont))
                problem = True
            if tnbCDS != faaprot:
                logger.error(("{} {}: no matching number of proteins between tbl and faa; "
                              "faa={}; in tbl ={}").format(name, oriname, faaprot, tnbCDS))
                problem = True
            if tnbGene + tnbCRISPR != ffngene and tnbGene != ffngene:
                logger.error(("{} {}: no matching number of genes between tbl and ffn; "
                        "ffn={}; in tbl ={}genes {}CRISPR").format(name, oriname,
                                                                   ffngene, tnbGene, tnbCRISPR))
                problem = True
    return not problem and not missing_file


def count_tbl(tblfile):
    """
    Count the different features found in the tbl file:
    - number of contigs
    - number of proteins (CDS)
    - number of genes (locus_tag)
    - number of CRISPR arrays (repeat_region)
    """
    nbcont = 0
    nbCDS = 0
    nbGene = 0
    nbCRISPR = 0
    with open(tblfile, "r") as tblf:
        for line in tblf:
            if line.startswith(">"):
                nbcont += 1
            if "CDS" in line:
                nbCDS += 1
            if "locus_tag" in line:
                nbGene += 1
            if "repeat_region" in line:
                nbCRISPR += 1
    return nbcont, nbCDS, nbGene, nbCRISPR


def count_headers(seqfile):
    """
    Count how many sequences there are in the given file
    """
    nbseq = 0
    with open(seqfile, "r") as faaf:
        for line in faaf:
            if line.startswith(">"):
                nbseq += 1
    return nbseq
