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


def define_nb_threads(threads):
    """
    With fasttree, number of threads to use must be defined before running the
    script, by changing an environment variable.
    """
    os.environ["OMP_NUM_THREADS"] = str(threads)


def run_fasttree(alignfile, boot, treefile):
    """
    Run fasttree on given alignment
    """
    if not boot:
        bootinfo = "-nosupport"
    else:
        bootinfo = "-boot {}".format(boot)
    logfile = alignfile + ".fasttree.log"
    if not treefile:
        treefile = alignfile + ".fasttree_tree.nwk"
    cmd = "FastTreeMP -nt -gtr -noml -nocat -log {} {} {}".format(logfile, alignfile, bootinfo)
    stdout = open(treefile, "w")
    error = ("Problem while running Fasttree. See log file ({}) for "
             "more information.").format(logfile)
    utils.run_cmd(cmd, error, stdout=stdout, eof=True, logger=logger)
