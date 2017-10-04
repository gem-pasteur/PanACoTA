#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the protein_seq_functions in pangenome module
"""

import pytest
import os
import shutil
import genomeAPCAT.pangenome_module.protein_seq_functions as psf
import test.test_unit.utilities_for_tests as util_tests


@pytest.fixture(scope="function")
def path_test_pan():
    return os.path.join("test", "data", "pangenome")


@pytest.fixture(scope="function")
def path_test_files():
    return os.path.join(path_test_pan(), "test_files")


@pytest.fixture(scope="function")
def path_exp_files():
    return os.path.join(path_test_pan(), "exp_files")


def test_build_bank(path_test_files, path_exp_files):
    """
    Build a protein bank from a list of genomes, and create it at the same
    place as the database.
    """
    lstinfo = os.path.join(path_test_files, "list_to_pan.txt")
    dbpath = os.path.join(path_test_files, "example_db", "Proteins")
    name = "EXEM"
    spedir = None
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    exp_file = os.path.join(path_exp_files, "exp_EXEM.All.prt")
    exp_out = os.path.join(dbpath, name + ".All.prt")
    assert outfile == exp_out
    with open(outfile, "r") as of, open(exp_file, "r") as ef:
        for l1, l2 in zip(of, ef):
            assert l1 == l2
    os.remove(outfile)


def test_build_bank_spedir(path_test_files, path_exp_files):
    """
    Build a protein bank from a list of genomes, and create it in a given output directory.
    """
    lstinfo = os.path.join(path_test_files, "list_to_pan.txt")
    dbpath = os.path.join(path_test_files, "example_db", "Proteins")
    name = "EXEM"
    spedir = os.path.join("test_build_prt", "toto")
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    exp_file = os.path.join(path_exp_files, "exp_EXEM.All.prt")
    exp_out = os.path.join(spedir, name + ".All.prt")
    assert outfile == exp_out
    with open(outfile, "r") as of, open(exp_file, "r") as ef:
        for l1, l2 in zip(of, ef):
            assert l1 == l2
    shutil.rmtree(spedir)
    shutil.rmtree("test_build_prt")


def test_build_bank_exists(path_test_files, path_exp_files, caplog):
    """
    Test that when we want to create a bank but the output file already exists, it prints
    a warning, and closes without redoing the bank
    """
    lstinfo = os.path.join(path_test_files, "list_to_pan.txt")
    dbpath = os.path.join(path_test_files, "example_db", "Proteins")
    name = "EXEM"
    spedir = None
    exp_out = os.path.join(dbpath, name + ".All.prt")
    exp_file = os.path.join(path_exp_files, "exp_EXEM.All.prt")
    shutil.copyfile(exp_file, exp_out)
    quiet = True
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    assert outfile == exp_out
    assert ("Protein bank test/data/pangenome/test_files/example_db/Proteins/EXEM.All.prt "
            "already exists. It will be used by mmseqs.") in caplog.text
    assert caplog.records[0].levelname == "WARNING"
    os.remove(outfile)


def test_build_bank_noquiet(path_test_files, path_exp_files, caplog):
    """
    Test that when the bank is created without the quiet option, it also works as expected
    """
    lstinfo = os.path.join(path_test_files, "list_to_pan.txt")
    dbpath = os.path.join(path_test_files, "example_db", "Proteins")
    name = "EXEM"
    spedir = None
    quiet = False
    outfile = psf.build_prt_bank(lstinfo, dbpath, name, spedir, quiet)
    exp_file = os.path.join(path_exp_files, "exp_EXEM.All.prt")
    exp_out = os.path.join(dbpath, name + ".All.prt")
    assert outfile == exp_out
    with open(outfile, "r") as of, open(exp_file, "r") as ef:
        for l1, l2 in zip(of, ef):
            assert l1 == l2
    assert ("Building bank with all proteins to EXEM.All.prt") in caplog.text
    assert caplog.records[0].levelname == "INFO"
    os.remove(outfile)
