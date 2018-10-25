#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for fastme_func submodule of tree_module
"""

import os
import logging

import genomeAPCAT.tree_module.fastme_func as fme
from genomeAPCAT import utils
from . import utilities


# Define common variables
ALIGN = os.path.join("test", "data", "align", "exp_files", "exp_pers4genomes.grp.aln")
TREEPATH = os.path.join("test", "data", "tree")
EXPPATH = os.path.join(TREEPATH, "exp_files")
LOGFILE_BASE = "log_test_fastme"


def setup_module():
    """
    create logger at start of this test module
    """
    utils.init_logger(LOGFILE_BASE, 0, '', verbose=1)
    print("Createc logger")


def teardown_module():
    """
    Remove log files at the end of this test module
    """
    os.remove(LOGFILE_BASE + ".log")
    os.remove(LOGFILE_BASE + ".log.details")
    os.remove(LOGFILE_BASE + ".log.err")
    print("Remove log files")


def test_convert_phylip(caplog):
    """
    Test that when giving a valid fasta alignment file, it converts it to Stockholm format,
    as expected.
    """
    caplog.set_level(logging.DEBUG)
    outfile = "test_2phylip"
    fme.convert2phylip(ALIGN, outfile)
    exp_stk = os.path.join(EXPPATH, "exp_align_phylip.ph")
    assert os.path.isfile(outfile)
    same_files(outfile, exp_stk)
    os.remove(outfile)
    assert "Converting fasta alignment to PHYLIP-relaxed format" in caplog.text


def test_convert_exists(caplog):
    """
    Test that when asking to convert a file in phylip format, but output file already exists,
    it does not convert again, and writes warning message saying that current file will be used.
    """
    caplog.set_level(logging.DEBUG)
    exp_stk = os.path.join(EXPPATH, "exp_align_phylip.ph")
    assert fme.convert2phylip(ALIGN, exp_stk) is None
    assert 'Phylip alignment file already existing.' in caplog.text
    assert ("The Phylip alignment file test/data/tree/exp_files/exp_align_phylip.ph "
            "already exists. The program will use it instead of re-converting "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln") in caplog.text


def test_run_fme_default(caplog):
    """
    Test that when running fastme without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    boot = None
    write_boot = False
    threads = 1
    model = None
    treefile = "test_tree-fastme-default"
    quiet = False
    fme.run_fastme(align, boot, write_boot, threads, model, treefile, quiet)

    assert "Running FastME..." in caplog.text
    assert ("fastme -i {0} -dT -nB -s -T 1  -o {1} -I "
            "{0}.fastme.log").format(align, treefile) in caplog.text
    assert utilities.is_tree_lengths(treefile)
    assert not utilities.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)
    os.remove(treefile)
    os.remove(logs)


def test_run_fme_boot_j(caplog):
    """
    Test that when running fastme without bootstrap, and with JC69 model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    boot = 105
    write_boot = False
    threads = 1
    model = "J"
    treefile = "test_tree-fastme-boot-JC69"
    quiet = False
    fme.run_fastme(align, boot, write_boot, threads, model, treefile, quiet)

    assert "Running FastME..." in caplog.text
    assert ("fastme -i {0} -dJ -nB -s -T 1 -b 105 -o {1} -I "
            "{0}.fastme.log").format(align, treefile) in caplog.text
    assert not utilities.is_tree_lengths(treefile)
    assert utilities.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)
    bootfile = align + "_fastme_boot.txt"
    assert os.path.isfile(bootfile)
    os.remove(treefile)
    os.remove(logs)
    os.remove(bootfile)


def test_run_fme_boot_write_f84(caplog):
    """
    Test that when running fastme with bootstrap, and with F84 model, it returns a file
    in the expected format (all branches have lengths, + bootstrap value).
    + write bootstrap trees
    """
    caplog.set_level(logging.DEBUG)
    align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    boot = 105
    write_boot = True
    threads = 1
    model = "4"
    treefile = "test_tree-fastme-boot-F84"
    quiet = False
    fme.run_fastme(align, boot, write_boot, threads, model, treefile, quiet)

    assert "Running FastME..." in caplog.text
    assert ("fastme -i {0} -d4 -nB -s -T 1 -b 105 -o {1} -I "
            "{0}.fastme.log -B {0}.fastme_bootstraps.nwk").format(align, treefile) in caplog.text
    assert not utilities.is_tree_lengths(treefile)
    assert utilities.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)
    bootfile = align + ".fastme_bootstraps.nwk"
    assert os.path.isfile(bootfile)
    os.remove(treefile)
    os.remove(logs)
    os.remove(bootfile)


def test_run_fme_notreename_rysym(caplog):
    """
    Test that when running fastme without bootstrap, and with RY-symetric model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    boot = 105
    write_boot = True
    threads = 1
    model = "Y"
    treename = None
    quiet = True
    fme.run_fastme(align, boot, write_boot, threads, model, treename, quiet)
    treefile = align + ".fastme_tree.nwk"
    assert "Running FastME..." in caplog.text
    assert ("fastme -i {0} -dY -nB -s -T 1 -b 105 -o {1} -I "
            "{0}.fastme.log -B {0}.fastme_bootstraps.nwk").format(align, treefile) in caplog.text
    assert not utilities.is_tree_lengths(treefile)
    assert utilities.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)
    bootfile = align + ".fastme_bootstraps.nwk"
    assert os.path.isfile(bootfile)
    os.remove(treefile)
    os.remove(logs)
    os.remove(bootfile)


def test_run_tree(caplog):
    """
    Test generating tree from fasta alignment with fastme (conversion)
    """
    caplog.set_level(logging.DEBUG)
    boot = 110
    treefile = "test_run_tree-fastme"
    quiet = False
    threads = 1
    model = 'T'
    write_boot = False
    fme.run_tree(ALIGN, boot, treefile, quiet, threads, model, write_boot)
    assert "Converting fasta alignment to PHYLIP-relaxed format" in caplog.text
    assert "Running FastME..." in caplog.text
    assert ("fastme -i {0}.phylip -dT -nB -s -T 1 -b 110 -o {1} -I "
            "{0}.phylip.fastme.log").format(ALIGN, treefile) in caplog.text
    assert not utilities.is_tree_lengths(treefile)
    assert utilities.is_tree_bootstrap(treefile)
    logs = ALIGN + ".phylip.fastme.log"
    assert os.path.isfile(logs)
    phylip = ALIGN + ".phylip"
    assert os.path.isfile(phylip)
    bootfile = ALIGN + ".phylip_fastme_boot.txt"
    assert os.path.isfile(bootfile)
    os.remove(treefile)
    os.remove(logs)
    os.remove(bootfile)

    # Redo with phylip alignments already generated
    boot = None
    treefile = "test_run_tree-fastme"
    quiet = False
    threads = 1
    model = 'T'
    write_boot = False
    fme.run_tree(ALIGN, boot, treefile, quiet, threads, model, write_boot)
    assert "Phylip alignment file already existing." in caplog.text
    assert ("The Phylip alignment file {0}.phylip already exists. The program will use it instead "
            "of re-converting {0}.").format(ALIGN) in caplog.text
    assert "Running FastME..." in caplog.text
    assert ("fastme -i {0}.phylip -dT -nB -s -T 1  -o {1} -I "
            "{0}.phylip.fastme.log").format(ALIGN, treefile) in caplog.text
    assert utilities.is_tree_lengths(treefile)
    assert not utilities.is_tree_bootstrap(treefile)
    logs = ALIGN + ".phylip.fastme.log"
    assert os.path.isfile(logs)
    phylip = ALIGN + ".phylip"
    assert os.path.isfile(phylip)
    os.remove(treefile)
    os.remove(logs)
    os.remove(phylip)


def same_files(file_out, file_exp):
    """
    Check that the 2 files have the same content.

    Parameters
    ----------
    file_out : str
        file generated by the test
    file_exp : str
        file containing what should be generated
    """
    with open(file_out, "r") as fo, open(file_exp, "r") as fe:
        lines_out = fo.readlines()
        lines_exp = fe.readlines()
        assert len(lines_exp) == len(lines_out)
        for linout, linexp in zip(lines_out, lines_exp):
            assert linout == linexp
