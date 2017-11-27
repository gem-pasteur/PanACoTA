#!/usr/bin/env python3

"""
Functional tests for genomeAPCAT corepers subcommand
"""
import glob
import os
import subprocess

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


def test_main_default_given_output(caplog):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 1
    multi = False
    mixed = False
    outfile = "test-persistent-genome.txt"
    corepers.main(PAN, tol, multi, mixed, outputfile=outfile)
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(PAN + ".bin")
    os.remove(PAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    out_pers = os.path.join(TEST_PATH, outfile)
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


def test_main_pers(caplog):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = False
    mixed = False
    corepers.main(PAN, tol, multi, mixed)
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(PAN + ".bin")
    os.remove(PAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    out_pers = os.path.join(TEST_PATH, "PersGenome_test_pan-for-corepers.txt_0.99.lst")
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
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes (meaning at least ceil(0.99*nb_strains) genomes) in "
            "each family.\nTo be considered as persistent, a family must contain exactly 1 "
            "member in at least 99.0% of all genomes. The other genomes are absent from the"
            "family.") in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating Persistent genome of a dataset containing 4 genomes" in caplog.text
    assert ("The persistent genome contains 2 families having at least "
            "4 genomes in each family.") in caplog.text


def test_main_pers_floor(caplog):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = False
    mixed = False
    floor = True
    corepers.main(PAN, tol, multi, mixed, floor=floor)
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(PAN + ".bin")
    os.remove(PAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    out_pers = os.path.join(TEST_PATH, "PersGenome_test_pan-for-corepers.txt_F0.99.lst")
    exp_pers = os.path.join(EXP_PATH, "exp_pers-floor-strict.txt")
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
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes (meaning at least floor(0.99*nb_strains) genomes) in "
            "each family.\nTo be considered as persistent, a family must contain exactly 1 "
            "member in at least 99.0% of all genomes. The other genomes are absent from the"
            "family.") in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating Persistent genome of a dataset containing 4 genomes" in caplog.text
    assert ("The persistent genome contains 5 families having at least "
            "3 genomes in each family.") in caplog.text


def test_main_pers_floor_mixed(caplog):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = False
    mixed = True
    floor = True
    corepers.main(PAN, tol, multi, mixed, floor=floor)
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(PAN + ".bin")
    os.remove(PAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    out_pers = os.path.join(TEST_PATH, "PersGenome_test_pan-for-corepers.txt_F0.99-mixed.lst")
    exp_pers = os.path.join(EXP_PATH, "exp_pers-floor-mixed.txt")
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
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes (meaning at least floor(0.99*nb_strains) genomes) in "
            "each family.\nMixed families are allowed. To be considered as persistent, "
            "a family must have exactly 1 member in 99.0% of the genomes, but in the "
            "remaining 1% genomes, there can be 0, 1 or several members.") in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating Persistent genome of a dataset containing 4 genomes" in caplog.text
    assert ("The persistent genome contains 7 families having at least "
            "3 genomes in each family.") in caplog.text


def test_main_pers_floor_multi(caplog):
    """
    Test that with default parameters, it creates the expected core genome.
    """
    tol = 0.99
    multi = True
    mixed = False
    floor = True
    corepers.main(PAN, tol, multi, mixed, floor=floor)
    # Check creation of binary file for pangenome, and remove it
    assert os.path.isfile(PAN + ".bin")
    os.remove(PAN + ".bin")
    # Check presence of persistent genome, and its content, and remove it
    out_pers = os.path.join(TEST_PATH, "PersGenome_test_pan-for-corepers.txt_F0.99-multi.lst")
    exp_pers = os.path.join(EXP_PATH, "exp_pers-floor-multi.txt")
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
    assert ("Will generate a Persistent genome with member(s) in at least 99.0% "
            "of all genomes (meaning at least floor(0.99*nb_strains) genomes) in "
            "each family.\nMultigenic families are allowed (several members "
            "in any genome of a family).") in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating Persistent genome of a dataset containing 4 genomes" in caplog.text
    assert ("The persistent genome contains 8 families having at least "
            "3 genomes in each family.") in caplog.text


def test_main_from_parse(caplog):
    """
    Test main when we give the output of the parser
    """
    import argparse
    tol = 1
    multi = False
    mixed = False
    floor = False
    args = argparse.Namespace()
    args.pangenome = PAN
    args.tol = tol
    args.multi = multi
    args.mixed = mixed
    args.floor = floor
    args.outfile = None
    args.verbose = 0
    args.quiet = False

    corepers.main_from_parse(args)

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


def test_main_all():
    """
    Test when calling corepers from command line, it runs and gives expected output files
    """
    cmd = "genomeAPCAT corepers -p {}".format(PAN)
    ret = subprocess.call(cmd.split())
    assert ret == 0

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
    # check no error
    with open(logfile + ".err", "r") as errf:
        assert errf.readlines() == []
    remove_all(logfile + "*")


def remove_all(pattern):
    """
    Remove all files with the given pattern
    """
    files = glob.glob(pattern)
    for f in files:
        os.remove(f)
