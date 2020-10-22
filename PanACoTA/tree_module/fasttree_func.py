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
    cmd = f"FastTreeMP -nt {model} -noml -nocat {bootinfo} -log {logfile} {alignfile}"
    logger.details("Fasttree command: " + cmd)
    if quiet:
        fnull = open(os.devnull, 'w')
    else:
        fnull = None
    stdout = open(treefile, "w")
    error = ("Problem while running Fasttree. See log file ({}) for "
             "more information.").format(logfile)
    utils.run_cmd(cmd, error, stdout=stdout, eof=True, logger=logger, stderr=fnull)
