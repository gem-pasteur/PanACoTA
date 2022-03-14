#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the protein_seq_functions in pangenome module
"""

import os
import shutil
import logging
import pytest

import PanACoTA.pangenome_module.protein_seq_functions as psf
import test.test_unit.utilities_for_tests as tutil


# Define variables and functions used by several tests
PANDIR = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PANDIR, "test_files")
PATH_EXP_FILES = os.path.join(PANDIR, "exp_files")
GENEPATH = os.path.join(PANDIR, "generated_by_unit-tests")


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



# Start tests
def test_build_bank(caplog):
    """
    Build a protein bank from a list of genomes, and create it at the same
    place as the database.
    """
    caplog.set_level(logging.DEBUG)
    lstinfo = os.path.join(PATH_TEST_FILES, "list_to_pan.txt")
    dbpath = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    cur_dbpath = os.path.join(GENEPATH, "Proteins")
    shutil.copytree(dbpath, cur_dbpath)
    name = "EXEM"
    spedir = None
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, cur_dbpath, name, spedir, quiet)
    exp_file = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    exp_out = os.path.join(cur_dbpath, name + ".All.prt")

    # Check prt bank filename
    assert outfile == exp_out
    # Check content of bank created
    assert tutil.compare_order_content(exp_file, exp_out)

    # Check logs
    assert ("Building bank with all proteins to test/data/pangenome/"
            "generated_by_unit-tests/Proteins/EXEM.All.prt") in caplog.text


def test_build_bank_spedir_quiet(caplog):
    """
    Build a protein bank from a list of genomes, and create it in a given output directory.
    """
    caplog.set_level(logging.DEBUG)
    lstinfo = os.path.join(PATH_TEST_FILES, "list_to_pan.txt")
    dbpath = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    name = "EXEM"
    spedir = os.path.join(GENEPATH, "test_build_prt", "toto")
    quiet = False
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    exp_file = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    exp_out = os.path.join(spedir, name + ".All.prt")
    assert outfile == exp_out
    assert tutil.compare_order_content(exp_file, exp_out)
    # Check logs
    assert ("Building bank with all proteins to test/data/pangenome/"
            "generated_by_unit-tests/test_build_prt/toto/EXEM.All.prt") in caplog.text


def test_build_bank_exists(caplog):
    """
    Test that when we want to create a bank but the output file already exists, it prints
    a warning, and closes without redoing the bank
    """
    caplog.set_level(logging.DEBUG)
    lstinfo = os.path.join(PATH_TEST_FILES, "list_to_pan.txt")
    dbpath = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    cur_dbpath = os.path.join(GENEPATH, "Proteins")
    shutil.copytree(dbpath, cur_dbpath)
    name = "EXEM"
    spedir = None
    exp_out = os.path.join(cur_dbpath, name + ".All.prt")
    exp_file = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    shutil.copyfile(exp_file, exp_out)
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, cur_dbpath, name, spedir, quiet)
    assert outfile == exp_out
    assert ("Protein bank test/data/pangenome/generated_by_unit-tests/Proteins/"
            "EXEM.All.prt already exists. It will be used by mmseqs.") in caplog.text
    assert caplog.records[0].levelname == "WARNING"
