#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for annotate/annotation_functions.py
"""

import pytest
import os
import logging
import shutil

import test.test_unit.utilities_for_tests as tutil
import PanACoTA.annotate_module.annotation_functions as afunc


# Define variables used by several tests
DBDIR = os.path.join("test", "data", "annotate")
GEN_PATH = os.path.join(DBDIR, "genomes")
TMP_PATH = os.path.join(DBDIR, "tmp_files")
EXP_DIR = os.path.join(DBDIR, 'exp_files')
TEST_DIR = os.path.join(DBDIR, 'test_files')
GENEPATH = os.path.join(DBDIR, "generated_by_unit-tests")
logger = logging.getLogger('test_genome_func')

@pytest.fixture(autouse=True)
def setup_teardown_module():
    """
    Remove log files at the end of this test module

    Before each test:
    - init logger
    - create directory to put generated files

    After:
    - remove all log files
    - remove directory with generated results
    """
    os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    print("teardown")


def test_count_headers():
    """
    Count how many sequences there are in the given multi-fasta file
    """
    seqfile = os.path.join(GEN_PATH, "genome4.fasta")
    nb = afunc.count_headers(seqfile)
    assert nb == 5


def test_count_tbl():
    """
    Count the different features found in the tbl file, and return
    nbcont, nbCDS, nbGene, nbCRISPR
    """
    tblfile = os.path.join(TEST_DIR, "original_name.fna-prokkaRes", "prokka_out_for_test.tbl")
    ncont, ncds, ngene, ncris = afunc.count_tbl(tblfile)
    assert ncont == 7
    assert ncds == 13
    assert ngene == 15
    assert ncris == 2
