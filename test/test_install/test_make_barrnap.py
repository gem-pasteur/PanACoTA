#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing PanACoTA according to already existing dependencies
From a computer with ubuntu and barrnap installed only
"""
import os
import glob
import pytest

from . import utilities as utils

@pytest.fixture
def install_panacota():
    print("INSTALLING PANACOTA")
    cmd = "python3 make"
    error = "Error installing"
    utils.run_cmd(cmd, error)


def teardown_function(function):
    """
    Uninstall PanACoTA and installed dependencies
    """
    print("TEARDOWN\n")
    cmd = "python3 make uninstall"
    error = "Error uninstall"
    utils.run_cmd(cmd, error)
    os.remove("install.log")
    print("cleaning repo")


def test_install():
    """
    Test that when installing from a computer containing only barrnap, it installs
    PanACoTA, and returns the list of missing dependencies
    """
    cmd = "python3 make"
    error = "Error trying to install PanACoTA from base"
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("mash")
    assert not utils.check_installed("PanACoTA")
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("prokka")
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("mash")
    assert utils.check_installed("PanACoTA")
    cmd = "pip3 show PanACoTA"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    with open(stdout, "r") as stdof:
        for line in stdof:
            if line.startswith("Location"):
                loc = line.split()[-1]
                assert glob.glob(os.path.join(loc, r'PanACoTA*dist-info'))
    os.remove(stdout)
    logfile = "install.log"
    content = ["Installing PanACoTA",
               "Some dependencies needed for some subcommands of PanACoTA are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mash (for prepare subcommand, to filter genomes)",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "- prokka (for annotate subcommand, with syntaxic + functional annotation)",
               "- prodigal : for annotate subcommand, you at least need prodigal (for syntaxic "
               "annotation only). If you even need functional annotation, also install prokka",
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


def test_upgrade(install_panacota):
    """
    Test upgrading PanACoTA when dependencies are still installed
    # """
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    cmd = "python3 make upgrade"
    error = "Error upgrade"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 2
        assert "Upgrading PanACoTA" in lines[0]
        assert "DONE" in lines[1]


def test_uninstall(install_panacota):
    """
    Test uninstalling PanACoTA
    """
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    cmd = "python3 make uninstall"
    error = "Error uninstalling"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("PanACoTA")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 2
        assert ":: INFO :: Uninstalling PanACoTA" in lines[0]
        assert "DONE" in lines[1]

def test_develop():
    """
    Test installing PanACoTA in developer mode, when barrnap is already installed
    """
    assert not utils.check_installed("PanACoTA")
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    cmd = "python3 make develop"
    error = "Error develop"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    cmd = "pip3 show PanACoTA"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    with open(stdout, "r") as stdof:
        for line in stdof:
            if line.startswith("Location"):
                loc = line.split()[-1]
                assert glob.glob(os.path.join(loc, r'PanACoTA*egg-info'))
    os.remove(stdout)
    logfile = "install.log"
    content = ["Installing developer packages needed for PanACoTA",
               "Some dependencies needed for some subcommands of PanACoTA are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mash (for prepare subcommand, to filter genomes)",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "- prokka (for annotate subcommand, with syntaxic + functional annotation). "
               "If you only need syntaxic annotation, prodigal is enough.",
               "- prodigal : for annotate subcommand, you at least need prodigal (for syntaxic ",
               "- One of the 4 following softwares, used to infer a phylogenetic tree:",
               "* FastTree (see README or documentation for more information on how to "
               "install it)", "* FastME", "* Quicktree", "* IQtree (or IQtree2)", "DONE"]
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


def test_install_user():
    """
    Test that when installing from a computer in user mode, it really installs
    PanACoTA in user mode
    """
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("PanACoTA")
    cmd = "python3 make --user"
    error = "Error trying to install PanACoTA from base"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    # Check logfile content
    logfile = "install.log"
    content = ["Installing PanACoTA in user mode...", "DONE"]
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
