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
    Test that when running fasttree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    threads = 1
    quiet = False
    model = "GTR"
    treefile = ""
    fast = False
    # Copy input file to the folder used for test results
    aldir = os.path.join(GENEPATH, "aldir")
    os.makedirs(aldir)
    cur_al = os.path.join(aldir, "alignment.grp.aln")
    shutil.copyfile(ALIGNMENT, cur_al)
    ft.run_tree(cur_al, boot, treefile, quiet, threads, model=model, wb="", mem="", 
                s="iqtree", f=fast)
    treefile = cur_al + ".iqtree_tree.treefile"
    assert os.path.isfile(treefile)
    logs = cur_al + ".iqtree_tree.log"
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree -s test/data/tree/generated_by_unit-tests/aldir/alignment.grp.aln "
            "-nt 1 -m GTR    -st DNA "
            "-pre test/data/tree/generated_by_unit-tests/aldir/alignment.grp.aln.iqtree_tree "
            "-quiet") in caplog.text
    assert tutil.is_tree_lengths(treefile)
    assert not tutil.is_tree_bootstrap(treefile)


def test_run_iqtree2_default(caplog):
    """
    Test that when running fasttree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    threads = 1
    quiet = False
    model = "GTR"
    treefile = ""
    fast = False
    # Copy input file to the folder used for test results
    aldir = os.path.join(GENEPATH, "aldir")
    os.makedirs(aldir)
    cur_al = os.path.join(aldir, "alignment.grp.aln")
    shutil.copyfile(ALIGNMENT, cur_al)
    treefile = os.path.join(aldir, "tree_default.iqtree")
    ft.run_tree(cur_al, boot, treefile, quiet, threads, model=model, wb="", mem="", 
                s="iqtree2", f=fast)
    out_treefile = treefile + ".treefile"
    assert os.path.isfile(out_treefile)
    logs = treefile + ".log"
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree2 -s test/data/tree/generated_by_unit-tests/aldir/alignment.grp.aln "
            "-T 1 -m GTR    --seqtype DNA "
            "--prefix test/data/tree/generated_by_unit-tests/aldir/tree_default.iqtree "
            "-quiet") in caplog.text
    assert tutil.is_tree_lengths(out_treefile)
    assert not tutil.is_tree_bootstrap(out_treefile)


def test_run_iqtree_boot(caplog):
    """
    Test that when running fasttree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    boot = 1000
    threads = 1
    quiet = False
    model = "GTR"
    fast = False
    # Copy input file to the folder used for test results
    aldir = os.path.join(GENEPATH, "aldir")
    treefile = os.path.join(aldir, "tree_boot.iqtree")
    os.makedirs(aldir)
    cur_al = os.path.join(aldir, "alignment.grp.aln")
    shutil.copyfile(ALIGNMENT, cur_al)
    ft.run_tree(cur_al, boot, treefile, quiet, threads, model=model, wb="", mem="", 
                f=fast, s="iqtree")
    out_treefile = treefile + ".treefile"
    assert os.path.isfile(out_treefile)
    logs = treefile + ".log"
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree -s test/data/tree/generated_by_unit-tests/aldir/alignment.grp.aln "
            "-nt 1 -m GTR  -bb 1000  -st DNA "
            "-pre test/data/tree/generated_by_unit-tests/aldir/tree_boot.iqtree "
            "-quiet") in caplog.text

    assert not tutil.is_tree_lengths(out_treefile)
    assert tutil.is_tree_bootstrap(out_treefile)


def test_run_iqtree2_boot(caplog):
    """
    Test that when running fasttree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    boot = 1000
    threads = 1
    quiet = False
    model = "GTR"
    treefile = ""
    fast = False
    # Copy input file to the folder used for test results
    aldir = os.path.join(GENEPATH, "aldir")
    os.makedirs(aldir)
    cur_al = os.path.join(aldir, "alignment.grp.aln")
    shutil.copyfile(ALIGNMENT, cur_al)
    ft.run_tree(cur_al, boot, treefile, quiet, threads, model=model, wb="", mem="", 
                f=fast, s="iqtree2")
    treefile = cur_al + ".iqtree_tree.treefile"
    assert os.path.isfile(treefile)
    logs = cur_al + ".iqtree_tree.log"
    assert os.path.isfile(logs)
    assert "Running IQtree..." in caplog.text
    assert ("iqtree2 -s test/data/tree/generated_by_unit-tests/aldir/alignment.grp.aln "
            "-T 1 -m GTR  -B 1000  --seqtype DNA "
            "--prefix test/data/tree/generated_by_unit-tests/aldir/alignment.grp.aln.iqtree_tree "
            "-quiet") in caplog.text
    assert not tutil.is_tree_lengths(treefile)
    assert tutil.is_tree_bootstrap(treefile)


    # TODO: 
    # memory
    # write_bootstrap
    # quiet
    # fast