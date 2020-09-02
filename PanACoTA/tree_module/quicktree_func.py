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
Functions to infer a phylogenetic tree with quicktree

@author gem
June 2017
"""

from Bio import AlignIO
import os
import logging

from PanACoTA import utils

logger = logging.getLogger("tree.quicktree")


def run_tree(alignfile, boot, outdir, *args, **kwargs):
    """
    Run quicktree for the given alignment file and options

    Parameters
    ----------
    alignfile: str
        Path to file containing alignments of persistent families grouped by genome
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    outdir: str or None
        Path to the tree file that must be created
    args: tuple
        Used to be compatible with the 'run_tree' function of other softs like fastME and
        fastTree which require more arguments like the DNA substitution model, the number of
        threads to use, etc.
    kwargs: dict
        Used to be compatible with the 'run_tree' function of other softs like fastME and
        fastTree which require more arguments like the DNA substitution model, the number of
        threads to use, etc.
    """
    align_name = os.path.basename(alignfile)
    align_stock = os.path.join(outdir, align_name + ".stockholm")
    convert2stockholm(alignfile, align_stock)
    run_quicktree(align_stock, boot, outdir)


def convert2stockholm(infile, outfile):
    """
    Input alignment is in fasta format. Input of quicktree must be in stockholm format.
    Convert it here.

    Parameters
    ----------
    infile: str
        Path to file containing alignments in fasta
    outfile: str
        Path to file which will contain the alignments converted to Stockholm format
    """
    if os.path.isfile(outfile):
        logger.info("Stockholm alignment file already existing.")
        logger.warning(("The Stockholm alignment file {} already exists. The program "
                        "will use it instead of re-converting {}.").format(outfile, infile))
        return
    logger.info("Converting fasta alignment to stockholm format.")
    with open(infile, 'r') as input_handle, open(outfile, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        AlignIO.write(alignments, output_handle, "stockholm")


def run_quicktree(alignfile, boot, outdir):
    """
    Run quicktree on the given alignment.

    Parameters
    ----------
    alignfile: str
        Path to file containing alignments of persistent families grouped by genome,
        in Stockholm format
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    outdir: str or None
        Path to the tree file that must be created
    """
    logger.info("Running Quicktree...")
    bootinfo = ""

    # Get bootstrap information
    if boot:
        bootinfo = f"-boot {boot}"
    # Get output filename and logfile name
    align_name = os.path.basename(alignfile)
    logfile = os.path.join(outdir, align_name + ".quicktree.log")
    treefile = os.path.join(outdir, align_name + ".quicktree_tree.nwk")
    cmd = f"quicktree -in a -out t {bootinfo} {alignfile}"
    outfile = open(treefile, "w")
    logfilef = open(logfile, "w")
    error = (f"Problem while running quicktree. See log file ({logfile}) for "
             "more information.")
    logger.details(cmd)
    utils.run_cmd(cmd, error, stdout=outfile, eof=True, logger=logger, stderr=logfilef)
