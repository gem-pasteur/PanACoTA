#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing PanACoTA according to already existing dependencies
From a computer with ubuntu, prokka and quicktree installed
"""
import os
import pytest
import glob

from . import utilities as utils

@pytest.fixture
def install_panacota():
    cmd = "python3 make"
    error = "Error installing"
    utils.run_cmd(cmd, error)


def teardown_function(function):
    """
    Uninstall PanACoTA and installed dependencies at the end of each test
    """
    print("TEARDOWN\n")
    cmd = "python3 make uninstall"
    error = "Error uninstall"
    utils.run_cmd(cmd, error)
    os.remove("install.log")
    print("cleaning repo")


def test_install_prokka_quicktree():
    """
    Test that when installing from a computer containing only prokka and quicktree, it installs
    PanACoTA, and returns the list of missing dependencies
    """
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prodigal")
    assert utils.check_installed("quicktree")
    assert not utils.check_installed("iqtree")
    assert utils.check_installed("prokka")
    assert not utils.check_installed("PanACoTA")
    assert utils.check_installed('quicktree')
    cmd = "python3 make"
    error = "Error trying to install PanACoTA from base"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed('quicktree')
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
    content = ["Installing PanACoTA...", "DONE",
               "Some dependencies needed for some subcommands of PanACoTA are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mmseqs (for pangenome subcommand)",
               "- mash (for prepare subcommand, to filter genomes)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "See more information on how to download/install those softwares in README or in "
               "documentation.",
               "- barrnap. If you use Prokka for functional annotation, it will not predict RNA."
               ]
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
    os.remove(logfile)


def test_upgrade(install_panacota):
    """
    Test upgrading PanACoTA when dependencies are still installed
    """
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed('quicktree')
    cmd = "python3 make upgrade"
    error = "Error upgrade"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed('quicktree')
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 2
        assert "Upgrading PanACoTA" in lines[0]
        assert "DONE" in lines[1]


def test_uninstall_withdep(install_panacota):
    """
    Test uninstalling PanACoTA when dependencies are still installed
    """
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("quicktree")
    assert utils.check_installed("PanACoTA")
    cmd = "python3 make uninstall"
    error = "Error uninstalling"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert not utils.check_installed("PanACoTA")
    assert utils.check_installed('quicktree')


def test_develop():
    """
    Test installing PanACoTA in developer mode, when prokka and barrnap are already installed
    """
    assert not utils.check_installed("PanACoTA")
    assert not utils.check_installed("barrnap")
    assert utils.check_installed('quicktree')
    assert utils.check_installed("prokka")
    cmd = "python3 make develop"
    error = "Error develop"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed('quicktree')
    cmd = "pip3 show PanACoTA"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    # Check installation in develop mode
    with open(stdout, "r") as stdof:
        for line in stdof:
            if line.startswith("Location"):
                loc = line.split()[-1]
                assert glob.glob(os.path.join(loc, r'PanACoTA*egg-info'))
    os.remove(stdout)
    logfile = "install.log"
    content = ["Installing PanACoTA...",
               "Installing developer packages needed for PanACoTA",
               "Some dependencies needed for some subcommands of PanACoTA are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "- mash (for prepare subcommand, to filter genomes)",
               "- barrnap. If you use Prokka for functional annotation, it will not predict RNA.",
               "See more information on how to download/install those softwares in README or in "
               "documentation."]

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
    os.remove(logfile)


def test_install_user():
    """
    Test that when installing from a computer containing only prokka, in user mode, it installs
    PanACoTA in /Users and returns list of dependencies
    """
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("mash")
    assert utils.check_installed("prokka")
    assert utils.check_installed("quicktree")
    assert not utils.check_installed("PanACoTA")
    cmd = "python3 make --user"
    error = "Error trying to install PanACoTA from base"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("mash")
    assert utils.check_installed("prokka")
    assert utils.check_installed("quicktree")
    assert utils.check_installed("PanACoTA")
    # Check logfile content
    logfile = "install.log"
    content = ["Installing PanACoTA in user mode...",
               "Some dependencies needed for some subcommands of PanACoTA are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "- mash (for prepare subcommand, to filter genomes)",
               "See more information on how to download/install those softwares in README or in "
               "documentation."]
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
