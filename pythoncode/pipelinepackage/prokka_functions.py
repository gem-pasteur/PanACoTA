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

def run_prokka_all(genomes, threads):
    """
    for each genome in genomes, run prokka to annotate the genome.

    """
    if threads <= 3:
        # mettre threads dans argument[1] pour chacun
        arguments = [(gpath + "-prokkaRes", threads, name, gpath)
                     for genome, (name, gpath, _, _, _) in genomes.items()]
        _ = [run_prokka(arg) for arg in arguments]
        # run_prokka((outdir, threads, name, gpath))

        logger.debug("not parallel")
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
    outdir, threads, name, gpath = arguments
    cmd = ("prokka --outdir {} --cpus {} "
           "--prefix {} {}").format(outdir, threads, name, gpath)
    logger.debug(cmd)
    FNULL = open(os.devnull, 'w')
    try:
        retcode = subprocess.call(shlex.split(cmd), stdout=FNULL, stderr=subprocess.STDOUT)
    except:
        print("error in prokka...")
