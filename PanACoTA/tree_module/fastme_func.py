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
Functions to infer a phylogenetic tree with fastME

@author gem
June 2017
"""

from Bio import AlignIO
import os
import logging

from PanACoTA import utils

logger = logging.getLogger("tree.fastme")

def run_tree(alignfile, boot, outdir, quiet, threads, **kwargs):
    """
    Run fastme for the given alignment file and options

    Parameters
    ----------
    alignfile: str
        Path to file containing alignments of persistent families grouped by genome
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    outdir: str
        output directory to save all results
    quiet: bool
        True if nothing must be printed to stderr/stdout, False otherwise
    threads: int
        Maximum number of threads to use
    kwargs["model"]: str
        DNA substitution model chosen by user
    kwargs["wb"]: bool
        True if all bootstrap pseudo-trees must be saved into a file, False otherwise
    """
    model = kwargs["model"]
    write_boot = kwargs["wb"]
    align_name = os.path.basename(alignfile)
    align_phylip = os.path.join(outdir, align_name + ".phylip")
    convert2phylip(alignfile, align_phylip)
    run_fastme(align_phylip, boot, write_boot, threads, model, outdir, quiet)


def convert2phylip(infile, outfile):
    """
    Input alignment is in fasta format. Input of fastME must be in Phylip-relaxed format.
    Convert it here.

    Parameters
    ----------
    infile: str
        Path to file in fasta format
    outfile: str
        Path to file to generate, in Phylip-relaxed format
    """
    if os.path.isfile(outfile):
        logger.info("Phylip alignment file already existing.")
        logger.warning(("The Phylip alignment file {} already exists. The program "
                        "will use it instead of re-converting {}.").format(outfile, infile))
        return
    logger.info("Converting fasta alignment to PHYLIP-relaxed format.")
    with open(infile, 'r') as input_handle, open(outfile, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        AlignIO.write(alignments, output_handle, "phylip-relaxed")


def run_fastme(alignfile, boot, write_boot, threads, model, outdir, quiet):
    """
    Run fastME on the given alignment.

    Parameters
    ----------
    alignfile: str
        Path to file containing alignments of persistent families grouped by genome
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    write_boot: bool
        True if all bootstrap pseudo-trees must be saved into a file, False otherwise
    threads: int
        Maximum number of threads to use
    model: str or None
        DNA substitution model chosen by user. None if default one
    outdir: str
        output directory to save all results
    quiet: bool
        True if nothing must be printed to stderr/stdout, False otherwise
    """
    logger.info("Running FastME...")
    bootinfo = ""
    threadinfo = ""
    outboot = ""

    # Get bootstrap information
    if boot:
        bootinfo = "-b {}".format(boot)
    # Get threads information
    if threads:
        threadinfo = "-T {}".format(threads)
    # Get output filename
    align_name = os.path.basename(alignfile)
    logfile = os.path.join(outdir, align_name + ".fastme.log")
    treefile = os.path.join(outdir, align_name + ".fastme_tree.nwk")
    # If bootstrap pseudo-trees must be written, define the filename here
    if write_boot:
        outboot = "-B " + os.path.join(outdir, align_name + ".fastme_bootstraps.nwk")
    # Put default model if not given
    if not model:
        model = "T"
    cmd = (f"fastme -i {alignfile} -d{model} -nB -s {threadinfo} {bootinfo} "
           f"-o {treefile} -I {logfile} {outboot}")
    logger.details(cmd)
    if quiet:
        fnull = open(os.devnull, 'w')
    else:
        fnull = None
    error = ("Problem while running FastME. See log file ({}) for "
             "more information.").format(logfile)
    utils.run_cmd(cmd, error, stdout=fnull, eof=True, logger=logger, stderr=fnull)
