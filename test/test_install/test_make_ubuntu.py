#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing PanACoTA according to already existing dependencies
Here, from ubuntu, so not having bioperl, git, etc.
"""
import os

from . import utilities as utils

def teardown_module():
    """
    Uninstall PanACoTA and installed dependencies
    """
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
    cmd = "python3 make"
    error = "Error trying to install PanACoTA from ubuntu"
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("prodigal")
    assert not utils.check_installed("genomeAPCAT")
    assert not utils.check_installed("mafft")
    assert not utils.check_installed("mmseqs")
    assert not utils.check_installed("quicktree")
    assert not utils.check_installed("fastme")
    assert not utils.check_installed("FastTreeMP")
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert not utils.check_installed("prokka")
    assert not utils.check_installed("prodigal")
    assert not utils.check_installed("mafft")
    assert not utils.check_installed("mmseqs")
    assert not utils.check_installed("quicktree")
    assert not utils.check_installed("fastme")
    assert not utils.check_installed("FastTreeMP")
    assert utils.check_installed("genomeAPCAT")
    cmd = "pip3 show genomeAPCAT"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    with open(stdout, "r") as stdof:
        lines = stdof.readlines()
        found = False
        for line in lines:
            if "/usr/local/lib" in line:
                found = True
                break
        assert found is True
    os.remove(stdout)
    logfile = "install.log"
    content = [":: INFO :: Installing genomeAPCAT...",
               ":: WARNING :: Some dependencies needed for some subcommands of genomeAPCAT are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- prokka (for annotate subcommand, with syntaxic + functional annotation)",
               "- barrnap. If you use Prokka for functional annotation, it will not predict RNA.",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree after)",
               "- mmseqs (for pangenome subcommand)",
               "- One of the 3 following softwares, used to infer a phylogenetic tree:",
               "* FastTree (see README or documentation for more information on how to "
               "install it)", "* FastME", "* Quicktree"]
    print("###### LOGFILE :")
    with open(logfile, "r") as logf:
        for line in logf:
            print(line)
    #     for linef, linee in zip(logf, content):
    #         assert linee in linef
    # # Check that needed packages are installed
    # assert utils.is_package_installed("argparse")
    # assert utils.is_package_installed("progressbar")
    # assert utils.is_package_installed("numpy")
    # assert utils.is_package_installed("matplotlib")
    # assert utils.is_package_installed("Bio")
    # assert not os.path.isdir(os.path.join("dependencies", "prokka"))
    # os.remove(logfile)
