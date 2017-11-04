#!/usr/bin/env python3

"""
Functional tests for genomeAPCAT pangenome
"""
import os

import pytest

import genomeAPCAT.subcommands.pangenome as pan
from genomeAPCAT import utils


LOGFILE_BASE = "func_test_pangenome"
TEST_FILES = os.path.join("test", "data", "pangenome", "test_files")
DBPATH = os.path.join(TEST_FILES, "example_db", "Proteins")


def setup_module():
    """
    create logger at start of this test module
    """
    utils.init_logger(LOGFILE_BASE, 0, '', verbose=1)


def teardown_module(module):
    """
    Remove log files at the end of this test module
    """
    os.remove(LOGFILE_BASE + ".log")
    os.remove(LOGFILE_BASE + ".log.details")
    os.remove(LOGFILE_BASE + ".log.err")


def test_main(caplog):
    """
    Test
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "testPAN4"
    min_id = 0.8
    outdir = "test_main_pangenome_dir"
    clust_mode = 1
    spe_dir = None
    threads = 1

    pan.main(lstinfo, name, DBPATH, min_id, outdir, clust_mode, spe_dir, threads)
    prtbank = os.path.join(DBPATH, "testPAN4.All.prt")
    assert os.path.isfile(prtbank)
    os.remove(prtbank)
    assert os.path.isdir(outdir)
    #check outdir content (presence of mmseq db, mmseq clust, tmp folder,
    #families found in pangenome.)
