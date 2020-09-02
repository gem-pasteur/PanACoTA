#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing PanACoTA according to already existing dependencies
Here, from basic ubuntu, so not having bioperl, git, etc.
"""
import os

from . import utilities as utils

def teardown_function(function):
    """
    Uninstall PanACoTA and installed dependencies
    """
    cmd = "python3 make uninstall"
    error = "Error uninstall"
    utils.run_cmd(cmd, error)
    os.remove("install.log")
    print("cleaning repo")


def test_install_panacota_base_ubuntu():
    """
    Test that when installing from a computer containing the basic ubuntu, it installs
    PanACoTA, without any dependence (show warning message)
    """
    cmd = "python3 make"
    error = "Error trying to install PanACoTA from ubuntu"
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("prodigal")
    assert not utils.check_installed("PanACoTA")
    assert not utils.check_installed("mafft")
    assert not utils.check_installed("mmseqs")
    assert not utils.check_installed("quicktree")
    assert not utils.check_installed("fastme")
    assert not utils.check_installed("FastTreeMP")
    assert not utils.check_installed("iqtree")
    assert not utils.check_installed("iqtree2")
    assert not utils.check_installed("mash")
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("prodigal")
    assert not utils.check_installed("mafft")
    assert not utils.check_installed("mmseqs")
    assert not utils.check_installed("quicktree")
    assert not utils.check_installed("fastme")
    assert not utils.check_installed("FastTreeMP")
    assert not utils.check_installed("iqtree")
    assert not utils.check_installed("iqtree2")
    assert not utils.check_installed("mash")
    assert utils.check_installed("PanACoTA")
    cmd = "pip3 show PanACoTA"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    with open(stdout, "r") as stdof:
        lines = stdof.readlines()
        found = False
        for line in lines:
            if "Summary: Large scale comparative genomics tools" in line:
                found = True
                break
        assert found is True
    os.remove(stdout)
    logfile = "install.log"
    content = [":: INFO :: Installing PanACoTA...",
               ":: WARNING :: Some dependencies needed "
               "for some subcommands of PanACoTA are not installed. Here is the list of "
               "missing dependencies, and for what they are used. If you plan "
               "to use the subcommands hereafter, first install required dependencies:",
               "prodigal : for annotate subcommand, you at least need prodigal (for syntaxic "
               "annotation only). If you even need functional annotation, also install prokka",
               "- mash (for prepare subcommand, to filter genomes)",
               "- prodigal : for annotate subcommand, you at least need prodigal (for syntaxic ",
               "- prokka (for annotate subcommand, with syntaxic + functional annotation). "
               "If you only need syntaxic annotation, prodigal is enough.",
               "- barrnap. If you use Prokka for functional annotation, it will not predict RNA.",
               "- mmseqs (for pangenome subcommand)",  "* Quicktree", "* IQtree (or IQtree2)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic "
               "tree after)",
               "- One of the 4 following softwares, used to infer a phylogenetic tree:",
               "* FastTree (see README or documentation for more information on how to "
               "install it)", "* FastME", "See more information on how to "
               "download/install those softwares in README or in documentation."]

    # Check output logfile content. Check that all content is present, in any order.
    print("###### LOGFILE :")
    with open(logfile, "r") as logf:
        logf_content = "".join(logf.readlines())
        for linec in content:
            assert linec in logf_content

    # # # Check that needed packages are installed
    assert utils.is_package_installed("argparse")
    assert utils.is_package_installed("progressbar")
    assert utils.is_package_installed("numpy")
    assert utils.is_package_installed("matplotlib")
    assert utils.is_package_installed("Bio")
    assert not os.path.isdir(os.path.join("dependencies"))
    assert not os.path.isdir(os.path.join("binaries"))
    os.remove(logfile)


def test_develop():
    """
    Test installing PanACoTA in developer mode, when barrnap is already installed
    """
    assert not utils.check_installed("PanACoTA")
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("mash")
    cmd = "python3 make develop"
    error = "Error develop"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("PanACoTA")
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("mash")
    assert not utils.check_installed("prokka")
    cmd = "pip3 show PanACoTA"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    # Check that it was not installed
    with open(stdout, "r") as stdof:
        for line in stdof:
            if line.startswith("Location"):
                loc = line.split()[-1]
                assert os.path.isdir(os.path.join(loc, "PanACoTA.egg-info"))
    os.remove(stdout)
    logfile = "install.log"
    content = ["Installing developer packages needed for PanACoTA",
               "Some dependencies needed for some subcommands of PanACoTA are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "- prokka (for annotate subcommand, with syntaxic + functional annotation). "
               "If you only need syntaxic annotation, prodigal is enough.",
               "- prodigal : for annotate subcommand, you at least need prodigal (for syntaxic ",
               "- One of the 4 following softwares, used to infer a phylogenetic tree:",
               "* FastTree (see README or documentation for more information on how to "
               "install it)", "* FastME", "* Quicktree", "* IQtree (or IQtree2)"]
    # Check output logfile content. Check that all content is present, in any order.
    with open(logfile, "r") as logf:
        logf_content = "".join(logf.readlines())
        for linec in content:
            assert linec in logf_content
    # Check that needed packages are installed
    assert utils.is_package_installed("argparse")
    assert utils.is_package_installed("progressbar")
    assert utils.is_package_installed("numpy")
    assert utils.is_package_installed("matplotlib")
    assert utils.is_package_installed("Bio")
    assert utils.is_package_installed("sphinx")
    assert utils.is_package_installed("coverage")
