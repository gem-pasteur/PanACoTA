#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the fasttree_func submodule in tree module
"""
import subprocess

import os

import genomeAPCAT.tree_module.fasttree_func as ft

# Define common variables
from genomeAPCAT import utils
from . import utilities

ALIGN = os.path.join("test", "data", "align", "exp_files", "exp_pers4genomes.grp.aln")
LOGFILE_BASE = "log_test_fasttree"


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
    boot = None
    treefile = "test_tree-fasttree-default"
    model = "-gtr"
    quiet = False
    ft.run_fasttree(ALIGN, boot, treefile, model, quiet)
    assert os.path.isfile(treefile)
    logs = ALIGN + ".fasttree.log"
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport "
            "-log {0}.fasttree.log {0}").format(ALIGN) in caplog.text
    assert utilities.is_tree_lengths(treefile)
    assert not utilities.is_tree_bootstrap(treefile)
    os.remove(treefile)
    os.remove(logs)


def test_run_fasttree_boot(caplog):
    """
    Test that when running fasttree with 10 bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, and bootstrap value).
    """
    boot = 10
    treefile = "test_tree-fasttree-boot"
    model = "-gtr"
    quiet = False
    ft.run_fasttree(ALIGN, boot, treefile, model, quiet)
    assert os.path.isfile(treefile)
    logs = ALIGN + ".fasttree.log"
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -boot 10 "
            "-log {0}.fasttree.log {0}").format(ALIGN) in caplog.text
    assert not utilities.is_tree_lengths(treefile)
    assert utilities.is_tree_bootstrap(treefile)
    os.remove(treefile)
    os.remove(logs)


def test_run_fasttree_boot_noname_jc(caplog):
    """
    Test that when running fasttree with bootstrap, and with JC model, it returns a file
    in the expected format (all branches have lengths, and bootstrap value).
    """
    boot = 100
    model = ""
    quiet = False
    ft.run_fasttree(ALIGN, boot, None, model, quiet)
    treefile = ALIGN + ".fasttree_tree.nwk"
    assert os.path.isfile(treefile)
    logs = ALIGN + ".fasttree.log"
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt  -noml -nocat -boot 100 "
            "-log {0}.fasttree.log {0}").format(ALIGN) in caplog.text
    assert not utilities.is_tree_lengths(treefile)
    assert utilities.is_tree_bootstrap(treefile)
    os.remove(treefile)
    os.remove(logs)


def test_run_fasttree_default_quiet(caplog):
    """
    Test that when running fasttree without bootstrap, and with default model, it returns a file
    in the expected format (all branches have lengths, no bootstrap value).
    Check that no bug when launch in quiet mode
    """
    boot = None
    treefile = "test_tree-fasttree-default"
    model = "-gtr"
    quiet = True
    ft.run_fasttree(ALIGN, boot, treefile, model, quiet)
    logs = ALIGN + ".fasttree.log"
    assert os.path.isfile(logs)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport "
            "-log {0}.fasttree.log {0}").format(ALIGN) in caplog.text
    assert utilities.is_tree_lengths(treefile)
    assert not utilities.is_tree_bootstrap(treefile)
    os.remove(treefile)
    os.remove(logs)


def test_run_tree(caplog):
    """
    Test that when running run_tree it defines the expected number of threads, and outputs the
    expected tree.
    """
    boot = None
    treefile = None
    quiet = False
    threads = 15
    model = "-gtr"
    ft.run_tree(ALIGN, boot, treefile, quiet, threads, model)
    res = subprocess.Popen("FastTreeMP", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = res.communicate()
    assert "OpenMP (15 threads)" in stdout.decode()
    out_tree = ALIGN + ".fasttree_tree.nwk"
    assert os.path.isfile(out_tree)
    out_log = ALIGN + ".fasttree.log"
    assert os.path.isfile(out_log)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport "
            "-log {0}.fasttree.log {0}").format(ALIGN) in caplog.text
    assert utilities.is_tree_lengths(out_tree)
    assert not utilities.is_tree_bootstrap(out_tree)
    os.remove(out_tree)
    os.remove(out_log)

    boot = 127
    quiet = True
    ft.run_tree(ALIGN, boot, treefile, quiet, 1, model)
    res = subprocess.Popen("FastTreeMP", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = res.communicate()
    assert "OpenMP (1 threads)" in stdout.decode()
    assert os.path.isfile(out_tree)
    assert os.path.isfile(out_log)
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -boot 127 "
            "-log {0}.fasttree.log {0}").format(ALIGN) in caplog.text
    assert not utilities.is_tree_lengths(out_tree)
    assert utilities.is_tree_bootstrap(out_tree)
    os.remove(out_tree)
    os.remove(out_log)
