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
    Run fastme for the given alignment file and options
    """
    align_stock = alignfile + ".stockholm"
    convert2stockholm(alignfile, align_stock)
    run_quicktree(align_stock, boot, treefile)


def convert2stockholm(infile, outfile):
    """
    Input alignment is in fasta format. Input of quicktree must be in stockholm format.
    Convert it here.
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
    """
    logger.info("Running Quicktree...")
    bootinfo = ""

    # Get bootstrap information
    if boot:
        bootinfo = "-boot {}".format(boot)
    # Get output filename
    if not treefile:
        treefile = alignfile + ".quicktree_tree.nwk"
    cmd = ("quicktree -in a -out t {boot} {infile}").format(boot=bootinfo, infile=alignfile)
    outfile = open(treefile, "w")
    # Define log filename
    logfile = alignfile + ".quicktree.log"
    logfilef = open(logfile, "w")
    error = ("Problem while running quicktree. See log file ({}) for "
             "more information.").format(logfile)
    logger.details(cmd)
    utils.run_cmd(cmd, error, stdout=outfile, eof=True, logger=logger, stderr=logfilef)


