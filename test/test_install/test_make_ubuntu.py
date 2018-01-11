#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing genomeAPCAT according to already existing dependencies
Here, from ubuntu, so not having bioperl, git, etc.
"""
import os

from . import utilities as utils

def teardown_module():
    """
    Uninstall genomeAPCAT and installed dependencies
    """
    cmd = "python3 make clean"
    error = "Error clean"
    utils.run_cmd(cmd, error)
    cmd = "python3 make uninstall"
    error = "Error uninstall"
    utils.run_cmd(cmd, error)
    os.remove("install.log")
    print("cleaning repo")


def test_build_prokka_only():
    """
    Test that when installing from a computer containing the basic ubuntu, it fails to install
    prokka, but still installs genomeAPCAT, without any dependence (warning message)
    """
    cmd = "apt-get install -y python3-pip"
    error = "Error apt-get install python3 and pip3"
    utils.run_cmd(cmd, error)
    cmd = "python3 make"
    error = "Error trying to install genomeAPCAT from ubuntu"
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("genomeAPCAT")
    assert not utils.check_installed("quicktree")
    assert not utils.check_installed("fastme")
    assert not utils.check_installed("FastTreeMP")
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("mafft")
    assert utils.check_installed("genomeAPCAT")
    cmd = "pip3 show genomeAPCAT"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    with open(stdout, "r") as stdof:
        lines = stdof.readlines()
        assert "/usr/local/lib" in lines[7]
    os.remove(stdout)
    logfile = "install.log"
    content = ["You need wget to install barrnap, the RNA predictor used by prokka.",
               "Installing prokka...",
               "A problem occurred while initializing prokka db. See log above.",
               "Problems while trying to install prokka (see above). While prokka is not "
               "installed, you will not be able to use the 'annotate' subcommand of genomeAPCAT",
               "Finalizing dependencies installation...",
               "Installing genomeAPCAT...",
               "Some dependencies needed for some subcommands of genomeAPCAT are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- prokka (for annotate subcommand)",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "- One of the 3 following softwares, used to infer a phylogenetic tree:",
               "* FastTree (see README or documentation for more information on how to "
               "install it)", "* FastME", "* Quicktree"]
    with open(logfile, "r") as logf:
        for linef, linee in zip(logf, content):
            assert linee in linef
    # Check that needed packages are installed
    assert utils.is_package_installed("argparse")
    assert utils.is_package_installed("progressbar")
    assert utils.is_package_installed("numpy")
    assert utils.is_package_installed("matplotlib")
    assert utils.is_package_installed("Bio")
    assert not os.path.isdir(os.path.join("dependencies", "prokka"))
    os.remove(logfile)


def test_clean():
    """
    Test cleaning dependencies. Should remove 'dependencies' and 'binaries' folders
    """
    assert os.path.isdir("dependencies")
    assert os.path.isdir("binaries")
    cmd = "python3 make clean"
    error = "not clean"
    utils.run_cmd(cmd, error)
    assert not os.path.isdir("dependencies")
    assert not os.path.isdir("binaries")
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("mafft")
    assert utils.check_installed("genomeAPCAT")
