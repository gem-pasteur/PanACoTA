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


logger = logging.getLogger()

def run_prokka_all(genomes, threads, force):
    """
    for each genome in genomes, run prokka to annotate the genome.

    """
    if threads <= 3:
        # mettre threads dans argument[1] pour chacun
        arguments = [(gpath, threads, name, force, nbcont)
                     for genome, (name, gpath, _, nbcont, _) in genomes.items()]
        _ = [run_prokka(arg) for arg in arguments]
    else:
        # use multiprocessing
        cores_prokka = 2
        arguments = [(gpath + "-prokkaRes", cores_prokka, name, gpath)
                     for genome, (name, gpath, _, _, _) in genomes.items()]
        cores_pool = int(threads/cores_prokka)
        pool = multiprocessing.Pool(cores_pool)
        try:
            final = pool.map(run_prokka, arguments)
            pool.close()
            pool.join()
        # If an error occurs, terminate pool and exit
        except AttributeError as excp:
            pool.terminate()
            print(excp.message)
            sys.exit(1)
        except Exception as excp:
            pool.terminate()
            print(excp.message)
            sys.exit(1)


def run_prokka(arguments):
    """

    """
    gpath, threads, name, force, nbcont = arguments
    outdir = gpath + "-prokkaRes"
    FNULL = open(os.devnull, 'w')
    prok_logfile = gpath + "-prokka.log"
    prokf = open(prok_logfile, "w")
    if os.path.isdir(outdir) and not force:
        logging.warning(("Prokka results folder already exists. Prokka did not run again, "
                         "formatting step used already generated results of Prokka in "
                         "{}. If you want to re-run prokka, first remove this result folder, or "
                         "use '-F' or '--force' option if you want to rerun prokka for "
                         "all genomes."))
        check_prokka(outdir, prok_logfile, name, gpath, nbcont)
        return
    elif os.path.isdir(outdir) and force == "--force":
        cmd = ("prokka --outdir {} --cpus {} {} "
               "--prefix {} {}").format(outdir, threads, force, name, gpath)
    else:
        cmd = ("prokka --outdir {} --cpus {} "
               "--prefix {} {}").format(outdir, threads, name, gpath)
    logger.debug(cmd)
    try:
        retcode = subprocess.call(shlex.split(cmd), stdout=FNULL, stderr=subprocess.STDOUT)
    except:
        print("error in prokka...")
