#!/usr/bin/env python3

"""
Functional tests for PanACoTA corepers subcommand
"""
import glob
import os
import subprocess
import pytest
import shutil
import logging

import PanACoTA.subcommands.corepers as corepers
import test.test_unit.utilities_for_tests as tutil


PERS_PATH = os.path.join("test", "data", "persgenome")
TEST_PATH = os.path.join(PERS_PATH, "test_files")
EXP_PATH = os.path.join(PERS_PATH, "exp_files")
GENEPATH = os.path.join(PERS_PATH, "generated_by_func-tests")
OPAN = os.path.join(TEST_PATH, "test_pan-for-corepers.txt")
UPAN = os.path.join(GENEPATH, "pangenome.lst")

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
    if not os.path.isdir(GENEPATH):
        os.mkdir(GENEPATH)
    shutil.copyfile(OPAN, UPAN)
    print("setup")

    yield
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")


def test_main_default(capsys):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 1
    multi = False
    mixed = False
    cmd = "cmd"
    out_pers = os.path.join(GENEPATH, "PersGenome_pangenome.lst_1.lst")
    assert corepers.main(cmd, UPAN, tol, multi, mixed, GENEPATH) == out_pers
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(UPAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    exp_pers = os.path.join(EXP_PATH, "exp_coregenome.txt")
    assert os.path.isfile(out_pers)
    assert tutil.compare_order_content(out_pers, exp_pers)
    # Check presence of log files and remove them
    logfile = os.path.join(GENEPATH, "PanACoTA-corepers.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".err")
    # Check log messages
    out, err = capsys.readouterr()
    assert "Will generate a CoreGenome." in out
    assert "Saving all information to a binary file for later use" in out
    assert "Generating Persistent genome of a dataset containing 4 genomes" in out
    assert ("The core genome contains 2 families, each one having exactly "
            "4 members, from the 4 different genomes.") in out


def test_main_pers(capsys):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = False
    mixed = False
    cmd = "cmd"
    out_pers = os.path.join(GENEPATH, "PersGenome_pangenome.lst_0.99.lst")
    assert corepers.main(cmd, UPAN, tol, multi, mixed, GENEPATH) == out_pers
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(UPAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    exp_pers = os.path.join(EXP_PATH, "exp_coregenome.txt")
    assert os.path.isfile(out_pers)
    assert tutil.compare_order_content(out_pers, exp_pers)
    # Check presence of log files and remove them
    logfile = os.path.join(GENEPATH, "PanACoTA-corepers.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".err")
    # Check log messages
    out, err = capsys.readouterr()
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes in each family.") in out
    assert ("To be considered as persistent, a family "
            "must contain exactly 1 member in at least 99.0% of all genomes. "
            "The other genomes are absent from the family.") in out
    assert "Saving all information to a binary file for later use" in out
    assert "Generating Persistent genome of a dataset containing 4 genomes" in out
    assert ("The persistent genome contains 2 families, each one having exactly "
            "1 member from at least 99.0% of the 4 different genomes "
            "(that is 4 genomes). The other genomes are absent from the family.") in out


def test_main_pers_floor_verbose2(capsys):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = False
    mixed = False
    floor = True
    cmd = "cmd"
    floor = True
    out_pers = os.path.join(GENEPATH, "PersGenome_pangenome.lst_F0.99.lst")
    assert corepers.main(cmd, UPAN, tol, multi, mixed, GENEPATH,
                         floor=floor, verbose=2) == out_pers
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(UPAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    exp_pers = os.path.join(EXP_PATH, "exp_pers-floor-strict.txt")
    assert os.path.isfile(out_pers)
    assert tutil.compare_order_content(out_pers, exp_pers)
    # Check presence of log files and remove them
    logfile = os.path.join(GENEPATH, "PanACoTA-corepers.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".err")
    assert os.path.isfile(logfile + ".details")
    # Check log messages
    out, err = capsys.readouterr()
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes in each family.") in out
    assert ("To be considered as persistent, a family must contain exactly 1 "
            "member in at least 99.0% of all genomes. The other genomes are absent from the "
            "family.") in out
    assert "Saving all information to a binary file for later use" in out
    assert "Generating Persistent genome of a dataset containing 4 genomes" in out
    assert ("The persistent genome contains 5 families, each one having exactly "
            "1 member from at least 99.0% of the 4 different genomes (that is 3 genomes). "
            "The other genomes are absent from the family.") in out


def test_main_pers_floor_mixed_debug(capsys):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = False
    mixed = True
    floor = True
    cmd = "cmd"
    out_pers = os.path.join(GENEPATH, "PersGenome_pangenome.lst_F0.99-mixed.lst")
    assert corepers.main(cmd, UPAN, tol, multi, mixed, GENEPATH,
                         floor=floor, verbose = 15) == out_pers
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(UPAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    exp_pers = os.path.join(EXP_PATH, "exp_pers-floor-mixed.txt")
    assert os.path.isfile(out_pers)
    assert tutil.compare_order_content(out_pers, exp_pers)
    # Check presence of log files and remove them
    logfile = os.path.join(GENEPATH, "PanACoTA-corepers.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".details")
    assert os.path.isfile(logfile + ".debug")
    assert os.path.isfile(logfile + ".err")
    # Check log messages
    out, err = capsys.readouterr()
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes in each family.") in out
    assert ("Mixed families are allowed. To be considered as persistent, "
            "a family must have exactly 1 member in 99.0% of the genomes, but in the "
            "remaining 1.0% genomes, there can be 0, 1 or several members.") in out
    assert "Saving all information to a binary file for later use" in out
    assert "Generating Persistent genome of a dataset containing 4 genomes" in out
    assert ("The persistent genome contains 7 families, each one having exactly "
            "1 member from at least 99.0% of the genomes (3 genomes). In the "
            "remaining 1.0% genomes, there can be 0, 1 or several members.") in out


def test_main_pers_floor_multi(capsys):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = True
    mixed = False
    floor = True
    cmd = "cmd"
    outdir = os.path.join(GENEPATH, "outdir")
    out_pers = os.path.join(outdir, "PersGenome_pangenome.lst_F0.99-multi.lst")
    assert corepers.main(cmd, UPAN, tol, multi, mixed, outdir, floor=floor) == out_pers
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(UPAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    exp_pers = os.path.join(EXP_PATH, "exp_pers-floor-multi.txt")
    assert os.path.isfile(out_pers)
    assert tutil.compare_order_content(out_pers, exp_pers)
    # Check presence of log files and remove them
    logfile = os.path.join(outdir, "PanACoTA-corepers.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".err")
    # Check log messages
    out, err = capsys.readouterr()
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes in each family.") in out
    assert ("Multigenic families are allowed (several members in any genome of a family).") in out
    assert "Saving all information to a binary file for later use" in out
    assert "Generating Persistent genome of a dataset containing 4 genomes" in out
    assert ("The persistent genome contains 8 families with members present in "
            "at least 3 different genomes (99.0% of the total number of genomes).") in out


def test_main_from_parse(capsys):
    """
    Test main when we give the output of the parser
    """
    import argparse
    tol = 1
    multi = False
    mixed = False
    floor = False
    args = argparse.Namespace()
    args.pangenome = UPAN
    args.tol = tol
    args.multi = multi
    args.mixed = mixed
    args.floor = floor
    args.outputdir = GENEPATH
    args.verbose = 0
    args.quiet = False
    args.argv = "PanACoTA corepers test_main_from_parse"

    corepers.main_from_parse(args)

    # Check creation of binary file for pangenome, and remove it
    out_pers = os.path.join(GENEPATH, "PersGenome_pangenome.lst_1.lst")
    assert os.path.isfile(UPAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    exp_pers = os.path.join(EXP_PATH, "exp_coregenome.txt")
    assert os.path.isfile(out_pers)
    assert tutil.compare_order_content(out_pers, exp_pers)
    # Check presence of log files and remove them
    logfile = os.path.join(GENEPATH, "PanACoTA-corepers.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".err")
    # Check log messages
    out, err = capsys.readouterr()
    print(out)
    assert "Will generate a CoreGenome." in out
    assert "Saving all information to a binary file for later use" in out
    assert "Generating Persistent genome of a dataset containing 4 genomes" in out
    assert ("The core genome contains 2 families, each one having exactly 4 "
            "members, from the 4 different genomes.") in out
