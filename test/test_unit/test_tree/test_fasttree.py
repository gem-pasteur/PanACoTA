#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the fasttree_func submodule in tree module
"""
import subprocess

import os

from Bio import Phylo

import genomeAPCAT.tree_module.fasttree_func as ft

# Define common variables
from genomeAPCAT import utils

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


# def test_run_fasttree_default(caplog):
#     """
#     Test that when running fasttree without bootstrap, and with default model, it returns a file
#     in the expected format (all branches have lengths, no bootstrap value).
#     """
#     boot = None
#     treefile = "test_tree-fasttree-default"
#     model = "-gtr"
#     quiet = False
#     ft.run_fasttree(ALIGN, boot, treefile, model, quiet)
#     assert os.path.isfile(treefile)
#     logs = ALIGN + ".fasttree.log"
#     assert os.path.isfile(logs)
#     assert "Running FasttreeMP..." in caplog.text
#     assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport "
#             "-log {0}.fasttree.log {0}").format(ALIGN) in caplog.text
#     is_tree_lengths(treefile)


# def is_tree_lengths(treefile):
#     """
#     Check that given tree is in newick format with branch lengths only (no bootstrap)

#     Parameters
#     ----------
#     treefile: str
#         Path to file containing tree to check

#     Returns
#     -------
#     bool
#     """
#     mytree = Phylo.read(treefile, "newick")
#     for elem in mytree.find_elements():
#         if isinstance(elem, Phylo.Newick.Clade):
#             if elem.name == ""
#                 # elem.confidence
#         # elem.name = leaf name if leaf, None otherwise
#         # elem.confidence bootstrap value if not leaf and bootstrap asked, None otherwise
#         # elem.branch_length: should always exist
