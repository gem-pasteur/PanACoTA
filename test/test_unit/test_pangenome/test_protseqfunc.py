#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the protein_seq_functions in pangenome module
"""

import os
import shutil
import PanACoTA.pangenome_module.protein_seq_functions as psf
import test.test_unit.utilities_for_tests as utiltest
import logging


# Define variables and functions used by several tests
PATH_TEST_PAN = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PATH_TEST_PAN, "test_files")
PATH_EXP_FILES = os.path.join(PATH_TEST_PAN, "exp_files")


# Start tests
def test_build_bank():
    """
    Build a protein bank from a list of genomes, and create it at the same
    place as the database.
    """
    lstinfo = os.path.join(PATH_TEST_FILES, "list_to_pan.txt")
    dbpath = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    name = "EXEM"
    spedir = None
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    exp_file = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    exp_out = os.path.join(dbpath, name + ".All.prt")

    # Check prt bank filename
    assert outfile == exp_out
    # Check content of bank created
    assert utiltest.compare_order_content(exp_file, exp_out)

    # Remove tmp files
    os.remove(outfile)


def test_build_bank_spedir():
    """
    Build a protein bank from a list of genomes, and create it in a given output directory.
    """
    lstinfo = os.path.join(PATH_TEST_FILES, "list_to_pan.txt")
    dbpath = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    name = "EXEM"
    spedir = os.path.join("test_build_prt", "toto")
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    exp_file = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    exp_out = os.path.join(spedir, name + ".All.prt")
    assert outfile == exp_out
    assert utiltest.compare_order_content(exp_file, exp_out)
    shutil.rmtree("test_build_prt")


def test_build_bank_exists(caplog):
    """
    Test that when we want to create a bank but the output file already exists, it prints
    a warning, and closes without redoing the bank
    """
    caplog.set_level(logging.DEBUG)
    lstinfo = os.path.join(PATH_TEST_FILES, "list_to_pan.txt")
    dbpath = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    name = "EXEM"
    spedir = None
    exp_out = os.path.join(dbpath, name + ".All.prt")
    exp_file = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    shutil.copyfile(exp_file, exp_out)
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    assert outfile == exp_out
    assert ("Protein bank test/data/pangenome/test_files/example_db/Proteins/EXEM.All.prt "
            "already exists. It will be used by mmseqs.") in caplog.text
    assert caplog.records[0].levelname == "WARNING"
    os.remove(outfile)


def test_build_bank_noquiet(caplog):
    """
    Test that when the bank is created without the quiet option, it also works as expected
    """
    caplog.set_level(logging.DEBUG)
    lstinfo = os.path.join(PATH_TEST_FILES, "list_to_pan.txt")
    dbpath = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    name = "EXEM"
    spedir = None
    quiet = False
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    exp_file = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    exp_out = os.path.join(dbpath, name + ".All.prt")
    assert outfile == exp_out
    assert utiltest.compare_order_content(exp_file, outfile)
    assert caplog.records[0].levelname == "INFO"
    assert ("Building bank with all proteins to "
            "test/data/pangenome/test_files/example_db/Proteins/EXEM.All.prt") in caplog.text
    os.remove(outfile)
