#!/usr/bin/env python3

"""
Functional tests for genomeAPCAT corepers subcommand
"""
import glob
import os

import genomeAPCAT.subcommands.corepers as corepers


PERS_PATH = os.path.join("test", "data", "persgenome")
TEST_PATH = os.path.join(PERS_PATH, "test_files")
EXP_PATH = os.path.join(PERS_PATH, "exp_files")
PAN = os.path.join(TEST_PATH, "test_pan-for-corepers.txt")


def test_main_default(caplog):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 1
    multi = False
    mixed = False
    corepers.main(PAN, tol, multi, mixed)
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(PAN + ".bin")
    os.remove(PAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    out_pers = os.path.join(TEST_PATH, "PersGenome_test_pan-for-corepers.txt_1.lst")
    exp_pers = os.path.join(EXP_PATH, "exp_coregenome.txt")
    assert os.path.isfile(out_pers)
    with open(out_pers, "r") as outf, open(exp_pers, "r") as expf:
        for line_out, line_exp in zip(outf, expf):
            assert line_out == line_exp
    os.remove(out_pers)
    # Check presence of log files and remove them
    logfile = os.path.join(TEST_PATH, "genomeAPCAT-corepers.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".details")
    assert os.path.isfile(logfile + ".err")
    remove_all(logfile + "*")
    # Check log messages
    assert "Will generate a CoreGenome." in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating Persistent genome of a dataset containing 4 genomes" in caplog.text
    assert ("The persistent genome contains 2 families having at least "
            "4 genomes in each family.") in caplog.text


def remove_all(pattern):
    """
    Remove all files with the given pattern
    """
    files = glob.glob(pattern)
    for f in files:
        os.remove(f)


