#!/usr/bin/env python3
# coding: utf-8

"""
Concatenate all alignment files of all families. Then,
group alignments by genome.

@author: GEM, Institut Pasteur
March 2017
"""

import logging
from genomeAPCAT import utils

logger = logging.getLogger("align.post")


def concat_alignments(fam_nums, prefix, quiet):
    """
    Concatenate all family alignment files to a unique file
    """
    logger.info("Concatenating all alignment files")
    list_files = ["{}-mafft-prt2nuc.{}.aln".format(prefix, num_fam) for num_fam in fam_nums]
    output = "{}-complete.cat.aln".format(prefix)
    if quiet:
        utils.cat(list_files, output)
    else:
        utils.cat(list_files, output, title = "Concatenation")

