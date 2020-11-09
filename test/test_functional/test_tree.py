#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for tree subcommand
"""
import os
import subprocess
import shutil
import pytest

from PanACoTA.subcommands import tree
from test.test_unit.test_tree import utilities as tutils

# Define common variables
ALDIR = os.path.join("test", "data", "align")
TREEDIR = os.path.join("test", "data", "tree")
EXPALDIR = os.path.join(ALDIR, "exp_files")
ALIGNMENT = os.path.join(EXPALDIR, "exp_pers4genomes.grp.aln")
GENEPATH = os.path.join(TREEDIR, "generated_by_func_tests")


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
    print("setup")

    yield
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")


def test_main_default(capsys):
    """
    Test that when giving the alignment file, and all default parameters, it runs iqtree,
    and returns expected files
    """
    outdir = GENEPATH
    soft = "iqtree"
    model = "GTR"
    threads = 1
    verbose = 3

    cmd = "cmd test_main_default"
    tree.main(cmd, ALIGNMENT, outdir, soft, model, threads, verbose=verbose)
    # Check output files
    iq_log_file = os.path.join(outdir, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(iq_log_file)
    with open(iq_log_file, "r") as logf:
        iq_lines = logf.readlines()
    assert ("Command: iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln -nt 1 "
            "-m GTR -st DNA -pre test/data/tree/generated_by_func_tests/"
            "exp_pers4genomes.grp.aln.iqtree_tree -quiet") in " ".join(iq_lines)
    assert ("Alignment has 4 sequences with 6438 columns, "
            "81 distinct patterns") in " ".join(iq_lines)
    assert "Analysis results written to:" in " ".join(iq_lines)
    tree_file = os.path.join(outdir, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(tree_file)
    assert tutils.is_tree_lengths(tree_file)
    assert not tutils.is_tree_bootstrap(tree_file)
    logs_base = os.path.join(outdir, "PanACoTA-tree-iqtree.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".err")
    # Check logs
    out, err = capsys.readouterr()
    print(out)
    assert "Running IQtree..." in out
    assert ("IQtree command: iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-nt 1 -m GTR    -st DNA -pre test/data/tree/generated_by_func_tests/"
            "exp_pers4genomes.grp.aln.iqtree_tree -quiet") in out
    assert "END" in out


def test_main_quicktree(capsys):
    """
    Test that when giving the alignment file, running with quicktree, no bootstraps,
    it creates expected files
    """
    outdir = GENEPATH
    soft = "quicktree"
    model = None
    threads = 1
    cmd = "cmd: test_main_quicktree"
    tree.main(cmd, ALIGNMENT, outdir, soft, model, threads, verbose=2)
    # Check output files
    # stockholm alignments
    stockholm = os.path.join(outdir, "exp_pers4genomes.grp.aln.stockholm")
    assert os.path.isfile(stockholm)
    # quicktree logfile
    log_file = stockholm + ".quicktree.log"
    assert os.path.isfile(log_file)
    with open(log_file, "r") as logf:
        assert logf.readlines() == []
    # tree file
    tree_file = stockholm + ".quicktree_tree.nwk"
    assert os.path.isfile(tree_file)
    assert tutils.is_tree_lengths(tree_file)
    assert not tutils.is_tree_bootstrap(tree_file)
    # log files
    logs_base = os.path.join(outdir, "PanACoTA-tree-quicktree.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".details")
    assert os.path.isfile(logs_base + ".err")
    # Check logs
    out, err = capsys.readouterr()
    assert "Converting fasta alignment to stockholm format." in out
    assert "Running Quicktree..." in out
    assert ("quicktree -in a -out t  test/data/tree/generated_by_func_tests/exp_pers4genomes.grp.aln.stockholm") in out
    assert "END" in out


def test_main_quicktree_notverbose(capsys):
    """
    Test that when giving the alignment file, running with quicktree, no bootstraps,
    it creates expected files
    """
    outdir = GENEPATH
    soft = "quicktree"
    model = None
    threads = 1
    cmd = "cmd: test_main_quicktree"
    tree.main(cmd, ALIGNMENT, outdir, soft, model, threads)
    # Check output files
    # stockholm alignments
    stockholm = os.path.join(outdir, "exp_pers4genomes.grp.aln.stockholm")
    assert os.path.isfile(stockholm)
    # quicktree logfile
    log_file = stockholm + ".quicktree.log"
    assert os.path.isfile(log_file)
    with open(log_file, "r") as logf:
        assert logf.readlines() == []
    # tree file
    tree_file = stockholm + ".quicktree_tree.nwk"
    assert os.path.isfile(tree_file)
    assert tutils.is_tree_lengths(tree_file)
    assert not tutils.is_tree_bootstrap(tree_file)
    # log files
    logs_base = os.path.join(outdir, "PanACoTA-tree-quicktree.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".details")
    assert os.path.isfile(logs_base + ".err")
    # Check logs
    out, err = capsys.readouterr()
    assert "Converting fasta alignment to stockholm format." in out
    assert "Running Quicktree..." in out
    assert "END" in out


def test_main_fastme(capsys):
    """
    Test that when giving the alignment file, running with fastme, no bootstraps,
    it creates expected files
    """
    outdir = GENEPATH
    soft = "fastme"
    model = None
    threads = None
    boot = 100
    cmd = "cmd: test_main_quicktree"
    tree.main(cmd, ALIGNMENT, outdir, soft, model, threads, boot=boot, verbose=16, write_boot=True)
    # Check output files
    # phylip alignments
    phylip = os.path.join(outdir, "exp_pers4genomes.grp.aln.phylip")
    assert os.path.isfile(phylip)
    # fastme logfile
    log_file = phylip + ".fastme.log"
    assert os.path.isfile(log_file)
    with open(log_file, "r") as logf:
        fastme_lines = logf.readlines()
    assert "Input data type  DNA" in " ".join(fastme_lines)
    assert "evolutionary model  TN93" in " ".join(fastme_lines)
    # tree file
    tree_file = phylip + ".fastme_tree.nwk"
    assert os.path.isfile(tree_file)
    assert not tutils.is_tree_lengths(tree_file)
    assert tutils.is_tree_bootstrap(tree_file)
    # bootstrap file
    tree_boot_file = phylip + ".fastme_bootstraps.nwk"
    assert os.path.isfile(tree_boot_file)
    with open(tree_boot_file, "r") as tbf:
        lines = tbf.readlines()
    assert len(lines) == boot
    # log files
    logs_base = os.path.join(outdir, "PanACoTA-tree-fastme.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".details")
    assert os.path.isfile(logs_base + ".debug")
    assert os.path.isfile(logs_base + ".err")
    # Check logs
    out, err = capsys.readouterr()
    assert "Converting fasta alignment to PHYLIP-relaxed format." in out
    assert "Running FastME..." in out
    assert ("fastme -i test/data/tree/generated_by_func_tests/exp_pers4genomes.grp.aln.phylip "
            "-dT -nB -s  -b 100 -o test/data/tree/generated_by_func_tests/"
            "exp_pers4genomes.grp.aln.phylip.fastme_tree.nwk "
            "-I test/data/tree/generated_by_func_tests/exp_pers4genomes.grp.aln.phylip.fastme.log "
            "-B test/data/tree/generated_by_func_tests/"
            "exp_pers4genomes.grp.aln.phylip.fastme_bootstraps.nwk") in out
    assert "END" in out


def test_main_fasttree(capsys):
    """
    Test that when giving the alignment file, running with fasttree, no bootstraps,
    it creates expected files
    """
    outdir = GENEPATH
    soft = "fasttree"
    model = "-gtr"
    threads = 1
    boot = 100
    cmd = "cmd: test_main_fasttree"
    tree.main(cmd, ALIGNMENT, outdir, soft, model, threads, boot=boot, verbose=2)
    # Check output files
    # fastme logfile
    log_file = os.path.join(outdir, "exp_pers4genomes.grp.aln.fasttree.log")
    assert os.path.isfile(log_file)
    with open(log_file, "r") as logf:
        fasttree_lines = logf.readlines()
    assert "Read 4 sequences, 6438 positions" in " ".join(fasttree_lines)
    assert "TreeCompleted" in " ".join(fasttree_lines)
    assert ("Nucleotide distances: Jukes-Cantor Joins: balanced Support: "
            "Local boot 100") in " ".join(fasttree_lines)
    # tree file
    tree_file = os.path.join(outdir, "exp_pers4genomes.grp.aln.fasttree_tree.nwk")
    assert os.path.isfile(tree_file)
    assert not tutils.is_tree_lengths(tree_file)
    assert tutils.is_tree_bootstrap(tree_file)
    # log files
    logs_base = os.path.join(outdir, "PanACoTA-tree-fasttree.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".details")
    assert not os.path.isfile(logs_base + ".debug")
    assert os.path.isfile(logs_base + ".err")
    # Check logs
    out, err = capsys.readouterr()
    print(out)
    assert "Running FasttreeMP..." in out
    assert ("Fasttree command: FastTreeMP -nt -gtr -noml -nocat -boot 100 "
            "-log test/data/tree/generated_by_func_tests/exp_pers4genomes.grp.aln.fasttree.log "
            "test/data/align/exp_files/exp_pers4genomes.grp.aln") in out
    assert "END" in out


def test_main_iqtree2_newdir(capsys):
    """
    Test that when giving the alignment file, and all default parameters, it runs iqtree2,
    and returns expected files
    """
    outdir = os.path.join(GENEPATH, "test_iqtree2")
    soft = "iqtree2"
    model = "GTR"
    threads = 1

    cmd = "cmd test_main_default"
    tree.main(cmd, ALIGNMENT, outdir, soft, model, threads, boot=1000, write_boot=True,
              verbose=2)
    # Check output files
    # Check iqtree logfile
    iq_log_file = os.path.join(outdir, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(iq_log_file)
    with open(iq_log_file, "r") as logf:
        iq_lines = logf.readlines()
    assert ("Command: iqtree2 -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-T 1 -m GTR -B 1000 --boot-trees --seqtype DNA "
            "--prefix test/data/tree/generated_by_func_tests/"
            "test_iqtree2/exp_pers4genomes.grp.aln.iqtree_tree "
            "--quiet") in " ".join(iq_lines)
    assert ("Alignment has 4 sequences with 6438 columns, "
            "81 distinct patterns") in " ".join(iq_lines)
    assert "Analysis results written to:" in " ".join(iq_lines)
    # Check treefile
    tree_file = os.path.join(outdir, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(tree_file)
    assert not tutils.is_tree_lengths(tree_file)
    assert tutils.is_tree_bootstrap(tree_file)
    # Check bootstrap tree file
    boot_file = os.path.join(outdir, "exp_pers4genomes.grp.aln.iqtree_tree.ufboot")
    assert os.path.isfile(tree_file)
    with open(boot_file, "r") as tbf:
        lines = tbf.readlines()
    assert len(lines) == 1000
    # Check panacota logfile
    logs_base = os.path.join(outdir, "PanACoTA-tree-iqtree2.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".err")
    assert os.path.isfile(logs_base + ".details")
    # Check logs
    out, err = capsys.readouterr()
    assert "Running IQtree..." in out
    assert ("IQtree command: iqtree2 -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-T 1 -m GTR  -B 1000 --boot-trees --seqtype DNA --prefix test/data/tree/"
            "generated_by_func_tests/test_iqtree2/exp_pers4genomes.grp.aln.iqtree_tree "
            "--quiet") in out
    assert "END" in out


def test_main_from_parse(capsys):
    """
    Test main when we give the output of the parser
    run fastme with K2P model, no bootstrap
    """
    import argparse
    args = argparse.Namespace()
    args.alignment = ALIGNMENT
    args.boot = None
    args.outdir = GENEPATH
    args.soft = "iqtree"
    args.model = "HKY"
    args.write_boot = False
    args.threads = 1
    args.verbose = 2
    args.quiet = False
    args.memory = False
    args.fast = False
    args.argv = "PanACoTA tree test_main_from_parse"
    tree.main_from_parse(args)
    # Check output files
    # iqtree log file
    iq_log_file = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.log")
    assert os.path.isfile(iq_log_file)
    with open(iq_log_file, "r") as logf:
        iq_lines = logf.readlines()
    assert ("Command: iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln -nt 1 "
            "-m HKY -st DNA -pre test/data/tree/generated_by_func_tests/"
            "exp_pers4genomes.grp.aln.iqtree_tree -quiet") in " ".join(iq_lines)
    assert ("Alignment has 4 sequences with 6438 columns, "
            "81 distinct patterns") in " ".join(iq_lines)
    assert "Analysis results written to:" in " ".join(iq_lines)
    # Tree file
    tree_file = os.path.join(GENEPATH, "exp_pers4genomes.grp.aln.iqtree_tree.treefile")
    assert os.path.isfile(tree_file)
    assert tutils.is_tree_lengths(tree_file)
    assert not tutils.is_tree_bootstrap(tree_file)
    # Panacota log file
    logs_base = os.path.join(GENEPATH, "PanACoTA-tree-iqtree.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".err")
    # Check logs
    out, err = capsys.readouterr()
    assert "Running IQtree..." in out
    assert ("IQtree command: iqtree -s test/data/align/exp_files/exp_pers4genomes.grp.aln "
            "-nt 1 -m HKY    -st DNA -pre test/data/tree/generated_by_func_tests/"
            "exp_pers4genomes.grp.aln.iqtree_tree -quiet") in out
    assert "END" in out
