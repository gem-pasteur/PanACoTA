#!/usr/bin/env python3
# coding: utf-8

"""
Functions to infer a phylogenetic tree with quicktree

@author gem
June 2017
"""

from Bio import AlignIO
import os
import logging

from genomeAPCAT import utils

logger = logging.getLogger("tree.quicktree")


def run_tree(alignfile, boot, treefile, *args, **kwargs):
    """
    Run quicktree for the given alignment file and options

    Parameters
    ----------
    alignfile: str
        Path to file containing alignments of persistent families grouped by genome
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    treefile: str or None
        Path to file which will contain the tree inferred
    args: tuple
        Used to be compatible with the 'run_tree' function of other softs like fastME and
        fastTree which require more arguments like the DNA substitution model, the number of
        threads to use, etc.
    kwargs: dict
        Used to be compatible with the 'run_tree' function of other softs like fastME and
        fastTree which require more arguments like the DNA substitution model, the number of
        threads to use, etc.
    """
    align_stock = alignfile + ".stockholm"
    convert2stockholm(alignfile, align_stock)
    run_quicktree(align_stock, boot, treefile)


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


def run_quicktree(alignfile, boot, treefile):
    """
    Run quicktree on the given alignment.

    Parameters
    ----------
    alignfile: str
        Path to file containing alignments of persistent families grouped by genome
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    treefile: str or None
        Path to file which will contain the tree inferred
    """
    logger.info("Running Quicktree...")
    bootinfo = ""

    # Get bootstrap information
    if boot:
        bootinfo = "-boot {}".format(boot)
    # Get output filename
    if not treefile:
        treefile = alignfile + ".quicktree_tree.nwk"
    cmd = "quicktree -in a -out t {boot} {infile}".format(boot=bootinfo, infile=alignfile)
    outfile = open(treefile, "w")
    # Define log filename
    logfile = alignfile + ".quicktree.log"
    logfilef = open(logfile, "w")
    error = ("Problem while running quicktree. See log file ({}) for "
             "more information.").format(logfile)
    logger.details(cmd)
    utils.run_cmd(cmd, error, stdout=outfile, eof=True, logger=logger, stderr=logfilef)
