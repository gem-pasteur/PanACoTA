#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of align subcommand
"""
import os
import subprocess

import shutil

import pytest

from genomeAPCAT.subcommands import tree
from ..test_unit.test_tree import utilities

# Define common variables
ALDIR = os.path.join("test", "data", "align", "exp_files")
ALIGN = os.path.join(ALDIR, "exp_pers4genomes.grp.aln")


def test_main_default(caplog):
    """
    Test that when giving the alignment file, and all default parameters, it runs fasttree,
    and returns expected files
    """
    boot = None
    outfile = None
    soft = "fasttree"
    model = "-gtr"
    write_boot = False
    threads = 1
    verbose = 0
    quiet = False
    tree.main(ALIGN, boot, outfile, soft, model, write_boot, threads, verbose, quiet)
    # Check output files
    log_file = ALIGN + ".fasttree.log"
    assert os.path.isfile(log_file)
    with open(log_file, "r") as logf:
        lines = logf.readlines()
    assert ("Command: FastTreeMP -nt -gtr -noml -nocat -nosupport -log {} "
            "{}\n".format(log_file, ALIGN)) in lines
    assert "Read 4 sequences, 6438 positions\n" in lines
    assert "TreeCompleted\n" in lines
    os.remove(log_file)
    tree_file = ALIGN + ".fasttree_tree.nwk"
    assert os.path.isfile(tree_file)
    assert utilities.is_tree_lengths(tree_file)
    assert not utilities.is_tree_bootstrap(tree_file)
    os.remove(tree_file)
    logs_base = os.path.join(ALDIR, "genomeAPCAT-tree-fasttree.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".details")
    assert os.path.isfile(logs_base + ".err")
    os.remove(logs_base)
    os.remove(logs_base + ".details")
    os.remove(logs_base + ".err")
    # Check logs
    assert "Running FasttreeMP..." in caplog.text
    assert ("FastTreeMP -nt -gtr -noml -nocat -nosupport -log {} "
            "{}".format(log_file, ALIGN)) in caplog.text
    assert "END" in caplog.text


def test_main_nofasttree(capsys):
    """
    Test running with default parameters when fasttree is not installed
    """
    # Simulated uninstallation of FastTree
    orig_ft = subprocess.check_output("which FastTreeMP".split()).decode().strip()
    temp_ft = orig_ft + "-orig"
    shutil.move(orig_ft, temp_ft)
    # Define parameters and run
    boot = None
    outfile = None
    soft = "fasttree"
    model = "-gtr"
    write_boot = False
    threads = 1
    verbose = 0
    quiet = False
    with pytest.raises(SystemExit):
        tree.main(ALIGN, boot, outfile, soft, model, write_boot, threads, verbose, quiet)
    out, _ = capsys.readouterr()
    assert out == "FastTreeMP is not installed. 'genomeAPCAT tree' cannot run.\n"
    # 're-install fasttreeMP
    shutil.move(temp_ft, orig_ft)


def test_main_fastme_logdet(caplog):
    """
    Test that when giving the alignment file, running with fastme, LogDet model, and bootstraps,
    it creates expected files
    """
    boot = 100
    outfile = None
    soft = "fastme"
    model = "L"
    write_boot = False
    threads = 1
    verbose = 0
    quiet = False
    tree.main(ALIGN, boot, outfile, soft, model, write_boot, threads, verbose, quiet)
    # Check output files
    # phylip alignments
    phylip = ALIGN + ".phylip"
    assert os.path.isfile(phylip)
    os.remove(phylip)
    # fastme logfile
    log_file = phylip + ".fastme.log"
    assert os.path.isfile(log_file)
    lines_exp = {"Input data type  DNA": False,
                 "evolutionary model  LogDet": False,
                 "Bootstrap: number of replicates  100": False}
    with open(log_file, "r") as logf:
        lines = logf.readlines()
        for line in lines:
            for exp in lines_exp:
                if exp in line:
                    lines_exp[exp] = True
    assert list(lines_exp.values()) == [True] * 3
    os.remove(log_file)
    # tree file
    tree_file = phylip + ".fastme_tree.nwk"
    assert os.path.isfile(tree_file)
    assert not utilities.is_tree_lengths(tree_file)
    assert utilities.is_tree_bootstrap(tree_file)
    os.remove(tree_file)
    # bootstrap file
    boot_file = phylip + "_fastme_boot.txt"
    assert os.path.isfile(boot_file)
    os.remove(boot_file)
    # log files
    logs_base = os.path.join(ALDIR, "genomeAPCAT-tree-fastme.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".details")
    assert os.path.isfile(logs_base + ".err")
    os.remove(logs_base)
    os.remove(logs_base + ".details")
    os.remove(logs_base + ".err")
    # Check logs
    assert "Converting fasta alignment to PHYLIP-relaxed format." in caplog.text
    assert "Running FastME..." in caplog.text
    assert ("fastme -i {0} -dL -nB -s -T 1 -b 100 -o {0}.fastme_tree.nwk "
            "-I {0}.fastme.log".format(phylip)) in caplog.text
    assert "END" in caplog.text


def test_main_nofastme(capsys):
    """
    Test running with default parameters when fastme is not installed
    """
    # Simulated uninstallation of FastTree
    orig_fme = subprocess.check_output("which fastme".split()).decode().strip()
    temp_fme = orig_fme + "-orig"
    shutil.move(orig_fme, temp_fme)
    # Define parameters and run
    boot = 100
    outfile = None
    soft = "fastme"
    model = "L"
    write_boot = False
    threads = 1
    verbose = 0
    quiet = False
    with pytest.raises(SystemExit):
        tree.main(ALIGN, boot, outfile, soft, model, write_boot, threads, verbose, quiet)
    out, _ = capsys.readouterr()
    assert out == "fastme is not installed. 'genomeAPCAT tree' cannot run.\n"
    # 're-install fasttreeMP
    shutil.move(temp_fme, orig_fme)


def test_main_quicktree(caplog):
    """
    Test that when giving the alignment file, running with quicktree, no bootstraps,
    it creates expected files
    """
    boot = None
    outfile = None
    soft = "quicktree"
    model = None
    write_boot = False
    threads = 1
    verbose = 0
    quiet = False
    tree.main(ALIGN, boot, outfile, soft, model, write_boot, threads, verbose, quiet)
    # Check output files
    # stockholm alignments
    stockholm = ALIGN + ".stockholm"
    assert os.path.isfile(stockholm)
    os.remove(stockholm)
    # fastme logfile
    log_file = stockholm + ".quicktree.log"
    assert os.path.isfile(log_file)
    with open(log_file, "r") as logf:
        assert logf.readlines() == []
    os.remove(log_file)
    # tree file
    tree_file = stockholm + ".quicktree_tree.nwk"
    assert os.path.isfile(tree_file)
    assert utilities.is_tree_lengths(tree_file)
    assert not utilities.is_tree_bootstrap(tree_file)
    os.remove(tree_file)
    # log files
    logs_base = os.path.join(ALDIR, "genomeAPCAT-tree-quicktree.log")
    assert os.path.isfile(logs_base)
    assert os.path.isfile(logs_base + ".details")
    assert os.path.isfile(logs_base + ".err")
    os.remove(logs_base)
    os.remove(logs_base + ".details")
    os.remove(logs_base + ".err")
    # Check logs
    assert "Converting fasta alignment to stockholm format." in caplog.text
    assert "Running Quicktree..." in caplog.text
    assert ("quicktree -in a -out t  {}".format(stockholm)) in caplog.text
    assert "END" in caplog.text


def test_main_noquicktree(capsys):
    """
    Test running with default parameters when quicktree is not installed
    """
    # Simulated uninstallation of quicktree
    orig_qui = subprocess.check_output("which quicktree".split()).decode().strip()
    temps_qui = orig_qui + "-orig"
    shutil.move(orig_qui, temps_qui)
    # Define parameters and run
    boot = None
    outfile = None
    soft = "quicktree"
    model = None
    write_boot = False
    threads = 1
    verbose = 0
    quiet = False
    with pytest.raises(SystemExit):
        tree.main(ALIGN, boot, outfile, soft, model, write_boot, threads, verbose, quiet)
    out, _ = capsys.readouterr()
    assert out == "quicktree is not installed. 'genomeAPCAT tree' cannot run.\n"
    # 're-install fasttreeMP
    shutil.move(temps_qui, orig_qui)
