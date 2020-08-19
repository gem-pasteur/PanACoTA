#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for general_format_functions.py
"""

import os
import logging
import shutil
from io import StringIO
import pytest

import PanACoTA.annotate_module.format_prokka as prokkafunc
import PanACoTA.utils as utils
import test.test_unit.utilities_for_tests as tutil

ANNOTEDIR = os.path.join("test", "data", "annotate")
EXP_ANNOTE = os.path.join(ANNOTEDIR, "exp_files")
TEST_ANNOTE = os.path.join(ANNOTEDIR, "test_files")
GENEPATH = os.path.join(ANNOTEDIR, "generated_by_unit-tests")

LOGFILE_BASE = os.path.join(GENEPATH, "logfile")
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]

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
    utils.init_logger(LOGFILE_BASE, 0, 'test_fastme', verbose=1)
    print("setup")

    yield
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    # shutil.rmtree(GENEPATH)
    print("teardown")


def test_tbl_to_lst_changed_names(caplog):
    """
    Check that generated lstinfo file is as expected, when the genome name is the same as
    it already was in the genome given to prokka.
    The test tblfile contains the following aspects:
    - gene in D strand (start < end)
    - gene in C strand (start > end)
    - CDS features (some with all info = ECnumber, gene name, product etc. ;
    some with missing info)
    - tRNA type
    - repeat_region type (*2)
    - contigs with more than 2 genes
    - contig with only 2 genes (both 'b' loc)
    - contig with 1 gene ('b' loc)
    - contig without gene (should be skipped)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prokka")
    tblfile = os.path.join(TEST_ANNOTE, "prokka_out_tbl_changed-contnames.tbl")
    lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
    contigs = ["test.0417.00002.0001\t50", "test.0417.00002.0002\t50", "test.0417.00002.0003\t50",
               "test.0417.00002.0004\t50", "test.0417.00002.0005\t50", "test.0417.00002.0006\t50",
               "test.0417.00002.0007\t50"]
    name = "test.0417.00002"
    assert prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, logger, changed_name=True)
    exp_lst = os.path.join(EXP_ANNOTE, "res_tbl2lst.lst")
    assert tutil.compare_order_content(exp_lst, lstfile)


def test_tbl_to_lst_not_changed_names(caplog):
    """
    Check that generated lstinfo file is as expected, when the genome name is the same as
    it already was in the genome given to prokka.
    The test tblfile contains the following aspects:
    - gene in D strand (start < end)
    - gene in C strand (start > end)
    - CDS features (some with all info = ECnumber, gene name, product etc. ;
    some with missing info)
    - tRNA type
    - repeat_region type (*2)
    - contigs with more than 2 genes
    - contig with only 2 genes (both 'b' loc)
    - contig with 1 gene ('b' loc)
    - contig without gene (should be skipped)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prokka")
    tblfile = os.path.join(TEST_ANNOTE, "original_name.fna-prokkaRes", "prokka_out_for_test.tbl")
    lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
    contigs = ["test.0417.00002.0001\t50", "test.0417.00002.0002\t50", "test.0417.00002.0003\t50",
               "test.0417.00002.0004\t50", "test.0417.00002.0005\t50", "test.0417.00002.0006\t50",
               "test.0417.00002.0007\t50"]
    name = "test.0417.00002"
    assert prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, logger, changed_name=False)
    exp_lst = os.path.join(EXP_ANNOTE, "res_tbl2lst.lst")
    assert tutil.compare_order_content(exp_lst, lstfile)

