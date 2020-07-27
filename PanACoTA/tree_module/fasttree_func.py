#!/usr/bin/env python3
# coding: utf-8

"""
Function to infer a phylogenetic tree from a persistent genome, using FastTree

@author gem
June 2017
"""

import os
import logging

from PanACoTA import utils

logger = logging.getLogger("tree.fasttree")

def run_tree(alignfile, boot, outdir, quiet, threads, **kwargs):
    """
    Run FastTree for the given alignment file and options

    Parameters
    ----------
    alignfile: str
        path to file containing all persistent families aligned, and grouped by genome
    boot: int or None
        number of bootstraps to calculate, None if no bootstrap asked
    outdir: str or None
        Path to the tree file that must be created
    quiet: bool
        True if nothing must be printed to stderr/stdout, False otherwise
    threads: int
        Maximum number of threads to use
    model: str
        DNA substitution model chosen by user
    kwargs: Object
        Used to be compatible with the 'run_tree' function of other softs like fastME and
        quicktree, which require more arguments
    """
    model = kwargs["model"]
    define_nb_threads(threads)
    run_fasttree(alignfile, boot, outdir, model, quiet)


def define_nb_threads(threads):
    """
    With FastTree, number of threads to use must be defined before running the
    script, by changing an environment variable.

    Parameters
    ----------
    threads: int
        Maximal number of threads to use
    """
    os.environ["OMP_NUM_THREADS"] = str(threads)


def run_fasttree(alignfile, boot, outdir, model, quiet):
    """
    Run FastTree on given alignment

    Parameters
    ----------
    alignfile: str
        Path to file containing all families aligned, grouped by genome
    boot: int or None
        Number of bootstraps to calculate (None if no bootstrap asked)
    treefile: str or None
        Path to the tree file that must be created
    model: str
        DNA substitution model
    quiet: bool
        True if nothing must be printed to stderr/stdout, False otherwise
    """
    logger.info("Running FasttreeMP...")
    if not boot:
        bootinfo = "-nosupport"
    else:
        bootinfo = "-boot {}".format(boot)
    align_name = os.path.basename(alignfile)
    logfile = os.path.join(outdir, align_name + ".fasttree.log")
    treefile = os.path.join(outdir, align_name + ".fasttree_tree.nwk")
    cmd = "FastTreeMP -nt {} -noml -nocat {} -log {} {}".format(model, bootinfo,
                                                                logfile, alignfile)
    logger.info("Fasttree command: " + cmd)
    if quiet:
        fnull = open(os.devnull, 'w')
    else:
        fnull = None
    stdout = open(treefile, "w")
    error = ("Problem while running Fasttree. See log file ({}) for "
             "more information.").format(logfile)
    logger.details(cmd)
    utils.run_cmd(cmd, error, stdout=stdout, eof=True, logger=logger, stderr=fnull)
