#!/usr/bin/env python3
# coding: utf-8

"""
Functions to infer a phylogenetic tree with fastME

@author gem
June 2017
"""

from Bio import AlignIO
import os
import logging

from genomeAPCAT import utils

logger = logging.getLogger("tree.fastme")


def run_tree(alignfile, boot, threads, treefile, quiet, model, write_boot):
    """
    Run fastme for the given alignment file and options
    """
    align_phylip = alignfile + ".phylip"
    convert2phylip(alignfile, align_phylip)
    run_fastme(align_phylip, boot, write_boot, threads, model, treefile, quiet)


def convert2phylip(infile, outfile):
    """
    Input alignment is in fasta format. Input of fastME must be in Phylip-relaxed format.
    Convert it here.
    """
    logger.info("Converting fasta alignment to PHYLIP-relaxed format.")
    with open(infile, 'r') as input_handle, open(outfile, 'w') as output_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        AlignIO.write(alignments, output_handle, "phylip-relaxed")


def run_fastme(alignfile, boot, write_boot, threads, model, treefile, quiet):
    """
    Run fastME on the given alignment.
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
    if not treefile:
        treefile = alignfile + ".fastme_tree.nwk"
    # If bootstrap pseudo-trees must be written, define the filename here
    if write_boot:
        outboot = "-B " + alignfile + ".fastme_bootstraps.nwk"
    # Put default model if not given
    if not model:
        model = "T"
    # Define log filename
    logfile = alignfile + ".fastme.log"
    matrix = alignfile + ".fastme_dist-mat.txt"
    cmd = ("fastme -i {align} -d{model} -nB -s {th} {bs} -o {out} -I {log} "
           "-O {mat} {wb}").format(align=alignfile, bs=bootinfo, th=threadinfo, model=model,
                                   out=treefile, log=logfile, mat=matrix, wb=outboot)
    if quiet:
        FNULL = open(os.devnull, 'w')
    else:
        FNULL = None
    error = ("Problem while running FastME. See log file ({}) for "
             "more information.").format(logfile)
    logger.details(cmd)
    utils.run_cmd(cmd, error, stdout=FNULL, eof=True, logger=logger, stderr=FNULL)


