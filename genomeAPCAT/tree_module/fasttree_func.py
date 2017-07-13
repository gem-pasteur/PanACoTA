#!/usr/bin/env python3
# coding: utf-8

"""
Functions to build a bank of all proteins to include in the pangenome

@author gem
April 2017
"""


import os

from genomeAPCAT import utils


def define_nb_threads(threads):
    """
    With fasttree, number of threads to use must be defined before running the
    script, by changing an environment variable.
    """
    os.environ["OMP_NUM_THREADS"] = str(threads)
    cmd = "FastTreeMP"
    error = "test nb threads"
    utils.run_cmd(cmd, error)
