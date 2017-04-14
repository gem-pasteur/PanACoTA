#!/usr/bin/env python3
# coding: utf-8

"""
Functions to deal with prokka


@author gem
April 2017
"""

import os
import logging
import subprocess
import shlex
import multiprocessing
import sys
import progressbar

logger = logging.getLogger()

def run_prokka_all(genomes, threads, force):
    """
    for each genome in genomes, run prokka to annotate the genome.

    * genomes = {genome: [name, gpath, size, nbcont, l90]}
    * threads: max number of threads that can be used
    * force: if None, do not override outdir (do not run prokka if outdir already exists). If
    '--force', rerun prokka and override existing outdir, for all genomes.

    Returns:
        final: {genome: boolean} -> with True if prokka ran well, False otherwise.
    """
    logger.info("Annotating all genomes with prokka")
    nbgen = len(genomes)
    # Create progressbar
    widgets = ['Annotation: ', progressbar.Bar(marker='█', left='', right='', fill=' '),
               ' ', progressbar.Percentage(), ' - ', progressbar.Timer()
              ]
    bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=100,
                                  redirect_stderr=True, redirect_stdout=True).start()

    if threads <= 3:
        # arguments : (gpath, cores_prokka, name, force, nbcont) for each genome
        arguments = [(genomes[g][1], threads, genomes[g][0], force, genomes[g][3])
                     for g in sorted(genomes)]
        final = []
        for num, arg in enumerate(arguments):
            res = run_prokka(arg)
            final.append(res)
            bar.update(num + 1)
        bar.finish()
    else:
        # use multiprocessing
        # if there are more threads than genomes, use as many threads as possible per genome
        if len(genomes) <= threads:
            cores_prokka = int(threads/len(genomes))
        # otherwise, use 2 threads per genome (and nb_genome/2 genomes at the same time)
        else:
            cores_prokka = 2
        # arguments : (gpath, cores_prokka, name, force, nbcont) for each genome
        arguments = [(genomes[g][1], cores_prokka, genomes[g][0], force, genomes[g][3])
                     for g in sorted(genomes)]
        cores_pool = int(threads/cores_prokka)
        pool = multiprocessing.Pool(cores_pool)
        try:
            final = pool.map_async(run_prokka, arguments)
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
        except AttributeError as excp:
            pool.terminate()
            logger.error(excp.message)
            sys.exit(1)
        except Exception as excp:
            pool.terminate()
            logger.error(excp.message)
            sys.exit(1)
    final = {genome: res for genome, res in zip(sorted(genomes), final)}
    logger.debug(final)
    return final


def run_prokka(arguments):
    """
    arguments : (gpath, cores_prokka, name, force, nbcont)

    returns:
        boolean. True if eveything went well (all needed output files present,
        corresponding numbers of proteins, genes etc.). False otherwise.
    """
    gpath, threads, name, force, nbcont = arguments
    outdir = gpath + "-prokkaRes"
    FNULL = open(os.devnull, 'w')
    prok_logfile = gpath + "-prokka.log"
    if os.path.isdir(outdir) and not force:
        logging.warning(("Prokka results folder already exists. Prokka did not run again, "
                         "formatting step used already generated results of Prokka in "
                         "{}. If you want to re-run prokka, first remove this result folder, or "
                         "use '-F' or '--force' option if you want to rerun prokka for "
                         "all genomes."))
        ok = check_prokka(outdir, prok_logfile, name, gpath, nbcont)
        return ok
    elif os.path.isdir(outdir) and force == "--force":
        prokf = open(prok_logfile, "w")
        cmd = ("prokka --outdir {} --cpus {} {} "
               "--prefix {} {}").format(outdir, threads, force, name, gpath)
    else:
        prokf = open(prok_logfile, "w")
        cmd = ("prokka --outdir {} --cpus {} "
               "--prefix {} {}").format(outdir, threads, name, gpath)
    # logger.debug(cmd)
    # print(cmd)
    try:
        retcode = subprocess.call(shlex.split(cmd), stdout=FNULL, stderr=prokf)
        ok = check_prokka(outdir, prok_logfile, name, gpath, nbcont)
        prokf.close()
        return ok
    except Exception as err:
        logging.error("Error while trying to run prokka: {}".format(err))
        prokf.close()
        return False
        sys.exit(1)


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
        tblfile = os.path.join(outdir, name + ".tbl")
        faafile = os.path.join(outdir, name + ".faa")
        ffnfile = os.path.join(outdir, name + ".ffn")
        if not os.path.isfile(tblfile):
            logger.error(("{} {}: no .tbl file").format(name, oriname))
            missing_file = True
        if not os.path.isfile(faafile):
            logger.error(("{} {}: no .faa file").format(name, oriname))
            missing_file = True
        if not os.path.isfile(ffnfile):
            logger.error(("{} {}: no .ffn file").format(name, oriname))
            missing_file = True
        if not missing_file:
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
            if tnbGene + tnbCRISPR != ffngene:
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
                nbCRISPR == 1
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