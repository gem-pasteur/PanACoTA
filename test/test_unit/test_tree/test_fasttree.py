#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the fasttree_func submodule in tree module
"""
import subprocess
import os
import logging
import pytest
import shutil

import PanACoTA.tree_module.fasttree_func as ft

# Define common variables
from PanACoTA import utils
from . import utilities as tutil


ALPATH = os.path.join("test", "data", "align")
ALIGNMENT = os.path.join(ALPATH, "exp_files", "exp_pers4genomes.grp.aln")
TREEPATH = os.path.join("test", "data", "tree")
GENEPATH = os.path.join(TREEPATH, "generated_by_unit-tests")
LOGFILE_BASE = "log_test_fasttree"
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
    utils.init_logger(LOGFILE_BASE, 0, 'test_fasttree', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    shutil.rmtree(GENEPATH)
    print("teardown")


def test_def_nb_threads():
    """
    Test that when changing the number of threads, it is taken into account as expected by
    fasttree binary.
    """
    thread = 10
    ft.define_nb_threads(thread)
    res = subprocess.Popen("FastTreeMP", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = res.communicate()
    assert "OpenMP (10 threads)" in stdout.decode()
    # Put back 1 for next tests, and check that it is taken into account
    ft.define_nb_threads(1)
    res1 = subprocess.Popen("FastTreeMP", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout1, _ = res1.communicate()
    assert "OpenMP (1 threads)" in stdout1.decode()


def test_run_fasttree_default(caplog):
    """
    Test that when running fasttree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    model = "-gtr"
    quiet = False
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree_tree.nwk")
    ft.run_fasttree(ALIGNMENT, boot, GENEPATH, model, quiet)
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree.log")
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport "
            "-log test/data/tree/generated_by_unit-tests/exp_pers4genomes.grp.aln.fasttree.log "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln")in caplog.text
    assert tutil.is_tree_lengths(treefile)
    assert not tutil.is_tree_bootstrap(treefile)


def test_run_fasttree_boot(caplog):
    """
    Test that when running fasttree with 10 bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, and bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    boot = 10
    model = "-gtr"
    quiet = False
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree_tree.nwk")
    ft.run_fasttree(ALIGNMENT, boot, GENEPATH, model, quiet)
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree.log")
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -boot 10 "
            "-log test/data/tree/generated_by_unit-tests/"
            "exp_pers4genomes.grp.aln.fasttree.log "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln") in caplog.text
    assert not tutil.is_tree_lengths(treefile)
    assert tutil.is_tree_bootstrap(treefile)


def test_run_fasttree_boot_noname_jc(caplog):
    """
    Test that when running fasttree with bootstrap, and with JC model, it returns a file
    in the expected format (all branches have lengths, and bootstrap value).
    """
    caplog.set_level(logging.DEBUG)
    boot = 100
    model = ""
    quiet = False
    ft.run_fasttree(ALIGNMENT, boot, GENEPATH, model, quiet)
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree_tree.nwk")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree.log")
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt  -noml -nocat -boot 100 "
            "-log test/data/tree/generated_by_unit-tests/"
            "exp_pers4genomes.grp.aln.fasttree.log "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln") in caplog.text
    assert not tutil.is_tree_lengths(treefile)
    assert tutil.is_tree_bootstrap(treefile)


def test_run_fasttree_default_quiet(caplog):
    """
    Test that when running fasttree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    Check that no bug when launch in quiet mode
    """
    caplog.set_level(logging.DEBUG)
    boot = None
    model = "-gtr"
    quiet = True
    ft.run_fasttree(ALIGNMENT, boot, GENEPATH, model, quiet)
    treefile = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree_tree.nwk")
    assert os.path.isfile(treefile)
    logs = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.fasttree.log")
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport "
            "-log test/data/tree/generated_by_unit-tests/"
            "exp_pers4genomes.grp.aln.fasttree.log "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln") in caplog.text
    assert tutil.is_tree_lengths(treefile)
    assert not tutil.is_tree_bootstrap(treefile)


def test_run_tree(caplog):
    """
    Test that when running run_tree it defines the expected number of threads, and outputs the
    expected tree.
    """
    # without bootstrap
    caplog.set_level(logging.DEBUG)
    boot = None
    quiet = False
    threads = 15
    model = "-gtr"
    ft.run_tree(ALIGNMENT, boot, ".", quiet, threads, model=model)
    # Check it used 15 threads
    res = subprocess.Popen("FastTreeMP", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = res.communicate()
    assert "OpenMP (15 threads)" in stdout.decode()
    # Check outputs
    out_tree = "exp_pers4genomes.grp.aln.fasttree_tree.nwk"
    assert os.path.isfile(out_tree)
    out_log = "exp_pers4genomes.grp.aln.fasttree.log"
    assert os.path.isfile(out_log)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport "
            "-log ./exp_pers4genomes.grp.aln.fasttree.log "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln") in caplog.text
    assert tutil.is_tree_lengths(out_tree)
    assert not tutil.is_tree_bootstrap(out_tree)
    os.remove(out_tree)
    os.remove(out_log)

    # with bootstrap
    boot = 127
    quiet = True
    ft.run_tree(ALIGNMENT, boot, ".", quiet, 1, model=model)
    res = subprocess.Popen("FastTreeMP", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = res.communicate()
    assert "OpenMP (1 threads)" in stdout.decode()
    out_tree = "exp_pers4genomes.grp.aln.fasttree_tree.nwk"
    assert os.path.isfile(out_tree)
    out_log = "exp_pers4genomes.grp.aln.fasttree.log"
    assert os.path.isfile(out_log)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -boot 127 "
            "-log ./exp_pers4genomes.grp.aln.fasttree.log "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln") in caplog.text
    assert not tutil.is_tree_lengths(out_tree)
    assert tutil.is_tree_bootstrap(out_tree)
    os.remove(out_tree)
    os.remove(out_log)
