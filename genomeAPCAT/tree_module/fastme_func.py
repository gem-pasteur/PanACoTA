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


def run_tree(alignfile, boot, treefile, quiet, threads, model, write_boot):
    """
    Run fastme for the given alignment file and options

    Parameters
    ----------
    alignfile: str
        Path to file containing alignments of persistent families grouped by genome
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    treefile: str
        Path to file which will contain the tree inferred
    quiet: bool
        True if nothing must be printed to stderr/stdout, False otherwise
    threads: int
        Maximum number of threads to use
    model: str
        DNA substitution model chosen by user
    write_boot: bool
        True if all bootstrap pseudo-trees must be saved into a file, False otherwise
    """
    align_phylip = alignfile + ".phylip"
    convert2phylip(alignfile, align_phylip)
    run_fastme(align_phylip, boot, write_boot, threads, model, treefile, quiet)


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


def run_fastme(alignfile, boot, write_boot, threads, model, treefile, quiet):
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
    treefile: str or None
        Path to file which will contain the tree inferred
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
           "{wb}").format(align=alignfile, bs=bootinfo, th=threadinfo, model=model,
                          out=treefile, log=logfile, wb=outboot)
    if quiet:
        fnull = open(os.devnull, 'w')
    else:
        fnull = None
    error = ("Problem while running FastME. See log file ({}) for "
             "more information.").format(logfile)
    logger.details(cmd)
    utils.run_cmd(cmd, error, stdout=fnull, eof=True, logger=logger, stderr=fnull)
