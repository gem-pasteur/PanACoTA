#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the iqtree_func submodule in tree module
"""
import subprocess
import os
import logging
import pytest
import shutil

import PanACoTA.tree_module.iqtree_func as ft

# Define common variables
from PanACoTA import utils
from . import utilities as tutil


ALPATH = os.path.join("test", "data", "align")
ALIGNMENT = os.path.join(ALPATH, "exp_files", "exp_pers4genomes.grp.aln")
TREEPATH = os.path.join("test", "data", "tree")
GENEPATH = os.path.join(TREEPATH, "generated_by_unit-tests")
LOGFILE_BASE = "log_test_iqtree"
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
    utils.init_logger(LOGFILE_BASE, 0, 'test_iqtree', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    shutil.rmtree(GENEPATH)
    print("teardown")


def test_run_iqtree_default(caplog):
    """
    Test that when running iqtree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value, treefile
    created with expected name).
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    threads = 1
    quiet = False
    model = "GTR"
    fast = False
    ft.run_tree(ALIGNMENT, boot, GENEPATH, quiet, threads, model=model, wb="", mem="",
                s="iqtree", f=fast)
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-nt 1 -m GTR    -st DNA "
            "-pre test/data/tree/generated_by_unit-tests/exp_pers4genomes.grp.aln.iqtree_tree "
            "-quiet") in caplog.text
    assert tutil.is_tree_lengths(treefile)
    assert not tutil.is_tree_bootstrap(treefile)


def test_run_iqtree2_default(caplog):
    """
    Test that when running iqtree2 without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value, treefile
    created with expected name).
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    threads = 1
    quiet = False
    model = "GTR"
    fast = False
    ft.run_tree(ALIGNMENT, boot, GENEPATH, quiet, threads, model=model, wb="", mem="",
                s="iqtree2", f=fast)
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree2 -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-T 1 -m GTR    --seqtype DNA "
            "--prefix test/data/tree/generated_by_unit-tests/exp_pers4genomes.grp.aln.iqtree_tree "
            "--quiet") in caplog.text
    assert tutil.is_tree_lengths(treefile)
    assert not tutil.is_tree_bootstrap(treefile)


def test_run_iqtree_boot_quiet_TVM(caplog):
    """
    Test that when running iqtree with bootstrap, and with TVM model, it returns a file
    in the expected format (all branches have lengths and bootstrap value).
    Treefile keeps the name given as input.
    """
    caplog.set_level(logging.DEBUG)
    boot = 1000
    threads = 1
    quiet = True
    model = "TVM"
    fast = False
    ft.run_tree(ALIGNMENT, boot, GENEPATH, quiet, threads, model=model, wb="", mem="",
                f=fast, s="iqtree")
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-nt 1 -m TVM  -bb 1000  -st DNA "
            "-pre test/data/tree/generated_by_unit-tests/"
            "exp_pers4genomes.grp.aln.iqtree_tree") in caplog.text
    assert not tutil.is_tree_lengths(treefile)
    assert tutil.is_tree_bootstrap(treefile)


def test_run_iqtree_boot_write_boot(caplog):
    """
    Test that when running iqtree with bootstrap, writing bootstrap trees,
    and with default model, it returns a file
    in the expected format (all branches have lengths and bootstrap value).
    Treefile keeps the name given as input.
    """
    caplog.set_level(logging.DEBUG)
    boot = 1000
    threads = 1
    quiet = True
    model = "GTR"
    fast = False
    ft.run_tree(ALIGNMENT, boot, GENEPATH, quiet, threads, model=model, wb=True, mem="",
                f=fast, s="iqtree")
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-nt 1 -m GTR  -bb 1000 -wbt -st DNA "
            "-pre test/data/tree/generated_by_unit-tests/"
            "exp_pers4genomes.grp.aln.iqtree_tree") in caplog.text
    assert not tutil.is_tree_lengths(treefile)
    assert tutil.is_tree_bootstrap(treefile)


def test_run_iqtree2_boot_write_boot(caplog):
    """
    Test that when running iqtree2 with bootstrap, writing bootstrap trees,
    and with default model, it returns a file
    in the expected format (all branches have lengths and bootstrap value).
    Treefile is named as expected.
    """
    caplog.set_level(logging.DEBUG)
    boot = 1000
    threads = 1
    quiet = False
    model = "GTR"
    treefile = ""
    fast = False
    ft.run_tree(ALIGNMENT, boot, GENEPATH, quiet, threads, model=model, wb=True, mem="",
                f=fast, s="iqtree2")
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree2 -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-T 1 -m GTR  -B 1000 --boot-trees --seqtype DNA "
            "--prefix test/data/tree/generated_by_unit-tests/exp_pers4genomes.grp.aln.iqtree_tree "
            "--quiet") in caplog.text
    assert not tutil.is_tree_lengths(treefile)
    assert tutil.is_tree_bootstrap(treefile)


def test_run_iqtree_fast_mem(caplog):
    """
    Test that when running iqtree without bootstrap, with F81 model,
    with fast option and giving a memory limit, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    Treefile is named as expected
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    threads = 1
    quiet = False
    model = "F81"
    fast = True
    ft.run_tree(ALIGNMENT, boot, GENEPATH, quiet, threads, model=model, wb="", mem="4GB",
                s="iqtree", f=fast)
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-nt 1 -m F81 -mem 4GB   -st DNA "
            "-pre test/data/tree/generated_by_unit-tests/exp_pers4genomes.grp.aln.iqtree_tree "
            "-quiet -fast") in caplog.text
    assert tutil.is_tree_lengths(treefile)
    assert not tutil.is_tree_bootstrap(treefile)


def test_run_iqtree2_fast_mem_quiet(caplog):
    """
    Test that when running iqtree without bootstrap, with default model,
    with fast option and giving a memory limit, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    Treefile is named as expected
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    threads = 1
    quiet = True
    model = "GTR"
    fast = True
    ft.run_tree(ALIGNMENT, boot, GENEPATH, quiet, threads, model=model, wb="", mem="4GB",
                s="iqtree2", f=fast)
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree2 -s test/data/align/exp_files/exp_pers4genomes.grp.aln -T 1 "
            "-m GTR --mem 4GB   --seqtype DNA "
            "--prefix test/data/tree/generated_by_unit-tests/exp_pers4genomes.grp.aln.iqtree_tree "
            "--quiet -fast") in caplog.text
    assert tutil.is_tree_lengths(treefile)
    assert not tutil.is_tree_bootstrap(treefile)

