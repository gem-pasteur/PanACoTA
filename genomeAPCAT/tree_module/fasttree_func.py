#!/usr/bin/env python3
# coding: utf-8

"""
Functions to build a bank of all proteins to include in the pangenome

@author gem
June 2017
"""

import os
import logging

from genomeAPCAT import utils

logger = logging.getLogger("tree.fasttree")


def run_tree(alignfile, boot, treefile, quiet, threads, model, *args):
    """
    Run fastme for the given alignment file and options
    """
    define_nb_threads(threads)
    run_fasttree(alignfile, boot, treefile, model, quiet)


def define_nb_threads(threads):
    """
    With fasttree, number of threads to use must be defined before running the
    script, by changing an environment variable.
    """
    os.environ["OMP_NUM_THREADS"] = str(threads)


def run_fasttree(alignfile, boot, treefile, model, quiet):
    """
    Run fasttree on given alignment
    """
    logger.info("Running FasttreeMP...")
    if not boot:
        bootinfo = "-nosupport"
    else:
        bootinfo = "-boot {}".format(boot)
    logfile = alignfile + ".fasttree.log"
    if not treefile:
        treefile = alignfile + ".fasttree_tree.nwk"
    cmd = "FastTreeMP -nt {} -noml -nocat {} -log {} {}".format(model, bootinfo,
                                                                logfile, alignfile)
    if quiet:
        FNULL = open(os.devnull, 'w')
    else:
        FNULL = None
    stdout = open(treefile, "w")
    error = ("Problem while running Fasttree. See log file ({}) for "
             "more information.").format(logfile)
    logger.details(cmd)
    utils.run_cmd(cmd, error, stdout=stdout, eof=True, logger=logger, stderr=FNULL)
