#!/usr/bin/env python3
# coding: utf-8

# ###############################################################################
# This file is part of PanACOTA.                                                #
#                                                                               #
# Authors: Amandine Perrin                                                      #
# Copyright © 2018-2020 Institut Pasteur (Paris).                               #
# See the COPYRIGHT file for details.                                           #
#                                                                               #
# PanACOTA is a software providing tools for large scale bacterial comparative  #
# genomics. From a set of complete and/or draft genomes, you can:               #
#    -  Do a quality control of your strains, to eliminate poor quality         #
# genomes, which would not give any information for the comparative study       #
#    -  Uniformly annotate all genomes                                          #
#    -  Do a Pan-genome                                                         #
#    -  Do a Core or Persistent genome                                          #
#    -  Align all Core/Persistent families                                      #
#    -  Infer a phylogenetic tree from the Core/Persistent families             #
#                                                                               #
# PanACOTA is free software: you can redistribute it and/or modify it under the #
# terms of the Affero GNU General Public License as published by the Free       #
# Software Foundation, either version 3 of the License, or (at your option)     #
# any later version.                                                            #
#                                                                               #
# PanACOTA is distributed in the hope that it will be useful, but WITHOUT ANY   #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS     #
# FOR A PARTICULAR PURPOSE. See the Affero GNU General Public License           #
# for more details.                                                             #
#                                                                               #
# You should have received a copy of the Affero GNU General Public License      #
# along with PanACOTA (COPYING file).                                           #
# If not, see <https://www.gnu.org/licenses/>.                                  #
# ###############################################################################

"""
Function to infer a phylogenetic tree from a persistent genome, using IQtree

@author gem
July 2020
"""

import os
import logging

from PanACoTA import utils

logger = logging.getLogger("tree.iqtree")

def run_tree(alignfile, boot, outdir, quiet, threads, **kwargs):
    """
    Run IQtree for the given alignment file and options

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
    fast = kwargs["f"]
    if not fast:
        fast = ""
    else:
        fast = "-fast"

    logger.info("Running IQtree...")

    # Init non mandatory arguments
    bootinfo = ""
    wb_info = ""
    mem_info = ""
    threadinfo = ""

    # Get info on all options (syntax changes according to IQtree version 1.x or 2.x)
    if boot:
        if soft=="iqtree":
            bootinfo = f"-bb {boot}"
        else:
            bootinfo = f"-B {boot}"
    if write_boot:
        if soft == "iqtree":
            wb_info = "-wbt"
        else:
    	    wb_info = "--boot-trees"
    if memory:
        if soft=="iqtree":
    	    mem_info = f"-mem {memory}"
        else:
            mem_info = f"--mem {memory}"
    # IQtree is always run quietly, but syntax depends on version:
    if soft=="iqtree":
        qu = "-quiet"
    else:
        qu = "--quiet"
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

    # Define treefile name if not given.
    align_name = os.path.basename(alignfile)
    logfile = os.path.join(outdir, align_name + ".iqtree.log")
    treefile = os.path.join(outdir, align_name + ".iqtree_tree")
    # get prefix cmd:
    if soft == "iqtree":
        prefix = f"-pre {treefile}"
    else:
        prefix = f"--prefix {treefile}"
    cmd = (f"{soft} -s {alignfile} {threadinfo} -m {model} {mem_info} {bootinfo} {wb_info} "
    	   f"{seqtype} {prefix} {qu} {fast}")
    logger.details("IQtree command: " + cmd)
    if quiet:
        fnull = open(os.devnull, 'w')
    else:
        fnull = None
    error = (f"Problem while running IQtree. See log file ({logfile}) for "
             "more information.")
    utils.run_cmd(cmd, error, eof=True, logger=logger, stderr=fnull)
