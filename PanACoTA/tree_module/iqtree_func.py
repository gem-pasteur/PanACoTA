#!/usr/bin/env python3
# coding: utf-8

"""
Function to infer a phylogenetic tree from a persistent genome, using IQtree

@author gem
July 2020
"""

import os
import logging

from PanACoTA import utils

logger = logging.getLogger("tree.iqtree")

def run_tree(alignfile, boot, treefile, quiet, threads, **kwargs):
    """
    Run IQtree for the given alignment file and options

    Parameters
    ----------
    alignfile: str
        path to file containing all persistent families aligned, and grouped by genome
    boot: int or None
        number of bootstraps to calculate, None if no bootstrap asked
    treefile: str or None
        Path to the tree file that must be created
    quiet: bool
        True if nothing must be printed to stderr/stdout, False otherwise
    threads: int
        Maximum number of threads to use
    kwargs["model"]: str
        DNA substitution model chosen by user
    kwards["wb"]: bool
    	True if all bootstrap pseudo-trees must be saved into a file, False otherwise
    kwargs["mem"]: str
    	Maximal RAM usage in GB | MB | % - Only for iqtree
    kwargs["s"]: str
    	soft to use (iqtree or iqtree2)
    """
    # Get optional arguments
    model = kwargs["model"]
    write_boot = kwargs["wb"]
    memory = kwargs["mem"]
    soft = kwargs["s"]

    logger.info("Running IQtree...")

    # Init non mandatory arguments
    bootinfo = ""
    wb_info = ""
    mem_info = ""
    threadinfo = ""

    # Get info on all options
    if boot:
        bootinfo = f"--B {boot}"
    if write_boot:
    	wb_info = "--boot-trees"
    if memory:
    	mem_info = f"-mem {memory}"
    # Get threads information
    if threads:
    	if soft == "iqtree":
    		threadinfo = f"-nt {threads}"
    	else:
    		threadinfo = f"-T {threads}"

	# get cmd for seqtype
    if soft == "iqtree":
        seqtype = "-st DNA"
    else:
        seqtype = "--seqtype DNA"

    # get prefix cmd:
    if soft == "iqtree":
        prefix = f"-pre {treefile}"
    else:
        prefix = f"--prefix {treefile}"

    if not treefile:
        treefile = alignfile + ".iqtree_tree"
    logfile = alignfile + ".iqtree.log"

    cmd = (f"{soft} -s {alignfile} {threadinfo} -m {model} {mem_info} {bootinfo} {wb_info}"
    	   f"{seqtype} {prefix} -fast")
    logger.info("IQtree command: " + cmd)
    if quiet:
        fnull = open(os.devnull, 'w')
    else:
        fnull = None
    error = (f"Problem while running Fasttree. See log file ({logfile}) for "
             "more information.")
    logger.details(cmd)
    utils.run_cmd(cmd, error, eof=True, logger=logger, stderr=fnull)
