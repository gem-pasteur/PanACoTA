#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for fastme_func submodule of tree_module
"""

import os
import logging
import pytest
import shutil

import PanACoTA.tree_module.fastme_func as fme
from PanACoTA import utils
import test.test_unit.utilities_for_tests as tutil
from . import utilities as tree_util


# Define common variables
ALPATH = os.path.join("test", "data", "align")
ALIGNMENT = os.path.join(ALPATH, "exp_files", "exp_pers4genomes.grp.aln")
TREEPATH = os.path.join("test", "data", "tree")
EXPPATH = os.path.join(TREEPATH, "exp_files")
GENEPATH = os.path.join(TREEPATH, "generated_by_unit-tests")
LOGFILE_BASE = "log_test_fastme"
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
    utils.init_logger(LOGFILE_BASE, 0, 'test_fastme', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    shutil.rmtree(GENEPATH)
    print("teardown")


def test_convert_phylip(caplog):
    """
    Test that when giving a valid fasta alignment file, it converts it to Stockholm format,
    as expected.
    """
    caplog.set_level(logging.DEBUG)
    outfile = os.path.join(GENEPATH, "test_2phylip")
    fme.convert2phylip(ALIGNMENT, outfile)
    exp_stk = os.path.join(EXPPATH, "exp_align_phylip.ph")
    assert os.path.isfile(outfile)
    tutil.compare_order_content(outfile, exp_stk)
    assert "Converting fasta alignment to PHYLIP-relaxed format" in caplog.text


def test_convert_exists(caplog):
    """
    Test that when asking to convert a file in phylip format, but output file already exists,
    it does not convert again, and writes warning message saying that current file will be used.
    """
    caplog.set_level(logging.DEBUG)
    exp_stk = os.path.join(EXPPATH, "exp_align_phylip.ph")
    fme.convert2phylip(ALIGNMENT, exp_stk)
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
    source_align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    align = os.path.join(GENEPATH, "align_tree.ph")
    # Copy file to output directory
    shutil.copyfile(source_align, align)
    boot = None
    write_boot = False
    threads = 1
    model = None
    treefile = os.path.join(GENEPATH, "test_tree-fastme-default")
    quiet = False
    fme.run_fastme(align, boot, write_boot, threads, model, treefile, quiet)

    assert "Running FastME..." in caplog.text
    assert ("fastme -i test/data/tree/generated_by_unit-tests/align_tree.ph -dT -nB -s -T 1 "
            " -o test/data/tree/generated_by_unit-tests/test_tree-fastme-default "
            "-I test/data/tree/generated_by_unit-tests/align_tree.ph.fastme.log") in caplog.text
    assert tree_util.is_tree_lengths(treefile)
    assert not tree_util.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)


def test_run_fme_boot_j(caplog):
    """
    Test that when running fastme without bootstrap, and with JC69 model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    source_align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    align = os.path.join(GENEPATH, "align_tree.ph")
    # Copy file to output directory
    shutil.copyfile(source_align, align)
    boot = 105
    write_boot = False
    threads = 1
    model = "J"
    treefile = os.path.join(GENEPATH, "test_tree-fastme-boot-JC69.tree")
    quiet = False
    fme.run_fastme(align, boot, write_boot, threads, model, treefile, quiet)

    assert "Running FastME..." in caplog.text
    assert ("fastme -i test/data/tree/generated_by_unit-tests/align_tree.ph -dJ -nB -s -T 1 "
            "-b 105 -o test/data/tree/generated_by_unit-tests/test_tree-fastme-boot-JC69.tree "
            "-I test/data/tree/generated_by_unit-tests/align_tree.ph.fastme.log") in caplog.text
    assert not tree_util.is_tree_lengths(treefile)
    assert tree_util.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)
    bootfile = align + "_fastme_boot.txt"
    assert os.path.isfile(bootfile)


def test_run_fme_boot_write_f84(caplog):
    """
    Test that when running fastme with bootstrap, and with F84 model, it returns a file
    in the expected format (all branches have lengths, + bootstrap value).
    + write bootstrap trees
    """
    caplog.set_level(logging.DEBUG)
    source_align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    align = os.path.join(GENEPATH, "align_tree.ph")
    # Copy file to output directory
    shutil.copyfile(source_align, align)
    boot = 105
    write_boot = True
    threads = 1
    model = "4"
    treefile = os.path.join(GENEPATH, "test_tree-fastme-boot-F84.tree")
    quiet = False
    fme.run_fastme(align, boot, write_boot, threads, model, treefile, quiet)

    assert "Running FastME..." in caplog.text
    assert ("fastme -i test/data/tree/generated_by_unit-tests/align_tree.ph -d4 -nB -s -T 1 "
            "-b 105 "
            "-o test/data/tree/generated_by_unit-tests/test_tree-fastme-boot-F84.tree "
            "-I test/data/tree/generated_by_unit-tests/align_tree.ph.fastme.log "
            "-B test/data/tree/generated_by_unit-tests/"
            "align_tree.ph.fastme_bootstraps.nwk") in caplog.text
    assert not tree_util.is_tree_lengths(treefile)
    assert tree_util.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)
    bootfile = align + ".fastme_bootstraps.nwk"
    assert os.path.isfile(bootfile)


def test_run_fme_notreename_rysym(caplog):
    """
    Test that when running fastme without bootstrap, and with RY-symetric model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    source_align = os.path.join(EXPPATH, "exp_align_phylip.ph")
    align = os.path.join(GENEPATH, "align_tree.ph")
    # Copy file to output directory
    shutil.copyfile(source_align, align)
    boot = 105
    write_boot = True
    threads = 1
    model = "Y"
    treename = None
    quiet = True
    fme.run_fastme(align, boot, write_boot, threads, model, treename, quiet)
    assert "Running FastME..." in caplog.text
    assert ("fastme -i test/data/tree/generated_by_unit-tests/align_tree.ph -dY -nB -s -T 1 "
            "-b 105 "
            "-o test/data/tree/generated_by_unit-tests/align_tree.ph.fastme_tree.nwk "
            "-I test/data/tree/generated_by_unit-tests/align_tree.ph.fastme.log "
            "-B test/data/tree/generated_by_unit-tests/"
            "align_tree.ph.fastme_bootstraps.nwk") in caplog.text
    treefile = align + ".fastme_tree.nwk"
    assert not tree_util.is_tree_lengths(treefile)
    assert tree_util.is_tree_bootstrap(treefile)
    logs = align + ".fastme.log"
    assert os.path.isfile(logs)
    bootfile = align + ".fastme_bootstraps.nwk"
    assert os.path.isfile(bootfile)
    assert os.path.isfile(treefile)


def test_run_tree(caplog):
    """
    Test generating tree from fasta alignment with fastme (conversion)
    """
    caplog.set_level(logging.DEBUG)
    boot = 110
    align = os.path.join(GENEPATH, "align_tree.aln")
    # Copy file to output directory
    shutil.copyfile(ALIGNMENT, align)
    treefile = os.path.join(GENEPATH, "test_run_tree-fastme.tree")
    quiet = False
    threads = 1
    model = 'T'
    write_boot = False
    fme.run_tree(align, boot, treefile, quiet, threads, model, write_boot)
    assert "Converting fasta alignment to PHYLIP-relaxed format" in caplog.text
    assert "Running FastME..." in caplog.text
    assert ("fastme "
            "-i test/data/tree/generated_by_unit-tests/align_tree.aln.phylip -dT -nB -s -T 1 "
            "-b 110 -o test/data/tree/generated_by_unit-tests/test_run_tree-fastme.tree "
            "-I test/data/tree/generated_by_unit-tests/"
            "align_tree.aln.phylip.fastme.log") in caplog.text
    assert not tree_util.is_tree_lengths(treefile)
    assert tree_util.is_tree_bootstrap(treefile)
    logs = align + ".phylip.fastme.log"
    assert os.path.isfile(logs)
    phylip = align + ".phylip"
    assert os.path.isfile(phylip)
    bootfile = align + ".phylip_fastme_boot.txt"
    assert os.path.isfile(bootfile)
    os.remove(treefile)
    os.remove(logs)
    os.remove(bootfile)

    # Redo with phylip alignments already generated
    boot = None
    treefile = os.path.join(GENEPATH, "test_rerun_tree-fastme.tree")
    quiet = False
    threads = 1
    model = 'T'
    write_boot = False
    fme.run_tree(align, boot, treefile, quiet, threads, model, write_boot)
    assert "Phylip alignment file already existing." in caplog.text
    assert ("The Phylip alignment file "
            "test/data/tree/generated_by_unit-tests/align_tree.aln.phylip "
            "already exists. The program will use it instead of re-converting "
            "test/data/tree/generated_by_unit-tests/align_tree.aln") in caplog.text
    assert "Running FastME..." in caplog.text
    assert ("fastme -i test/data/tree/generated_by_unit-tests/align_tree.aln.phylip "
            "-dT -nB -s -T 1 "
            " -o test/data/tree/generated_by_unit-tests/test_rerun_tree-fastme.tree "
            "-I test/data/tree/generated_by_unit-tests/"
            "align_tree.aln.phylip.fastme.log") in caplog.text
    assert tree_util.is_tree_lengths(treefile)
    assert not tree_util.is_tree_bootstrap(treefile)
    logs = align + ".phylip.fastme.log"
    assert os.path.isfile(logs)
    phylip = align + ".phylip"
    assert os.path.isfile(phylip)
    # os.remove(treefile)
    # os.remove(logs)
    # os.remove(phylip)
