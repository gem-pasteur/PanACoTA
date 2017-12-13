#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing genomeAPCAT according to already existing dependencies
"""
import os

from . import utilities as utils


def teardown_module():
    """
    Uninstall genomeAPCAT and installed dependencies
    """
    # cmd = "python3 make clean"
    # error = "Error clean"
    # utils.run_cmd(cmd, error)
    # cmd = "python3 make uninstall"
    # error = "Error uninstall"
    # utils.run_cmd(cmd, error)
    # os.remove("install.log")
    print("cleaning repo")


def test_build_prokka_quicktree():
    """
    Test that when installing from a computer containing only prokka, it installs
    genomeAPCAT (no barrnap), and returns the list of missing dependencies
    """
    cmd = "python3 make"
    error = "Error trying to install genomeAPCAT from base"
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
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
    content = ["Installing genomeAPCAT...",
               "Some dependencies needed for some subcommands of genomeAPCAT are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "See more information on how to download/install those softwares in README or in "
               "documentation."]
    with open(logfile, "r") as logf:
        lines_f = logf.readlines()
        assert len(lines_f) == len(content)
        for linef, linee in zip(lines_f, content):
            assert linee in linef
    # Check that needed packages are installed
    assert utils.is_package_installed("argparse")
    assert utils.is_package_installed("progressbar")
    assert utils.is_package_installed("numpy")
    assert utils.is_package_installed("matplotlib")
    assert utils.is_package_installed("Bio")
    os.remove(logfile)


def test_upgrade():
    """
    Test upgrading genomeAPCAT when dependencies are still installed
    """
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    cmd = "python3 make upgrade"
    error = "Error upgrade"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 1
        assert "Upgrading genomeAPCAT" in lines[0]
    os.remove(logfile)


def test_uninstall_withdep():
    """
    Test uninstalling genomeAPCAT when dependencies are still installed
    """
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    cmd = "python3 make uninstall"
    error = "Error uninstalling"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert not utils.check_installed("genomeAPCAT")


def test_develop():
    """
    Test installing genomeAPCAT in developer mode, when prokka and barrnap are already installed
    """
    assert not utils.check_installed("genomeAPCAT")
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    cmd = "python3 make develop"
    error = "Error develop"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    cmd = "pip3 show genomeAPCAT"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    # Check that it was not installed
    with open(stdout, "r") as stdof:
        lines = stdof.readlines()
        for line in lines:
            assert "/usr/local/lib" not in line
        assert "/tmp" or "/Users" in lines[7]
    os.remove(stdout)
    logfile = "install.log"
    content = ["Installing genomeAPCAT...",
               "Installing developer packages needed for genomeAPCAT",
               "Some dependencies needed for some subcommands of genomeAPCAT are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "See more information on how to download/install those softwares in README or in "
               "documentation."]
    with open(logfile, "r") as logf:
        lines_f = logf.readlines()
        assert len(lines_f) == len(content)
        for linef, linee in zip(lines_f, content):
            assert linee in linef
    # Check that needed packages are installed
    assert utils.is_package_installed("argparse")
    assert utils.is_package_installed("progressbar")
    assert utils.is_package_installed("numpy")
    assert utils.is_package_installed("matplotlib")
    assert utils.is_package_installed("Bio")
    os.remove(logfile)


def test_clean():
    """
    Test cleaning dependencies installed by make script
    """
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    cmd = "python3 make clean"
    error = "Error trying to clean dependencies"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    assert not os.path.isdir("binaries")
    assert not os.path.isdir("dependencies")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 1
        assert "Cleaning dependencies..." in lines[0]
    os.remove(logfile)


def test_uninstall_prokkadep():
    """
    Test uninstalling genomeAPCAT when dependencies are already removed
    """
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    cmd = "python3 make uninstall"
    error = "Error uninstalling"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("genomeAPCAT")
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")


def test_upgrade_notinstalled_prokkaponly():
    """
    Test upgrading genomeAPCAT when dependencies are not installed (only barrnap),
    and genomeAPCAT is not installed. It just installs genomeAPCAT, without prokka dep
    """
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert not utils.check_installed("genomeAPCAT")
    cmd = "python3 make upgrade"
    error = "Error upgrade"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 1
        assert "Upgrading genomeAPCAT" in lines[0]
    os.remove(logfile)


def test_install_user():
    """
    Test that when installing from a computer containing only prokka, in user mode, it installs
    genomeAPCAT in /Users and returns list of dependencies
    """
    cmd = "python3 make uninstall"
    error = "Error uninstalling"
    utils.run_cmd(cmd, error)
    cmd = "python3 make --user"
    error = "Error trying to install genomeAPCAT from base"
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert not utils.check_installed("genomeAPCAT")
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    cmd = "pip3 show genomeAPCAT"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    with open(stdout, "r") as stdof:
        lines = stdof.readlines()
        assert "/Users" in lines[7]
        assert "/lib/python" in lines[7]
        assert "/site-packages" in lines[7]
    os.remove(stdout)
    logfile = "install.log"
    content = ["Installing genomeAPCAT in user mode...",
               "Some dependencies needed for some subcommands of genomeAPCAT are "
               "not installed. Here is the list of missing dependencies, and for what they are "
               "used. If you plan to use the subcommands hereafter, first install required "
               "dependencies:",
               "- mmseqs (for pangenome subcommand)",
               "- mafft (to align persistent genomes in order to infer a phylogenetic tree "
               "after)",
               "See more information on how to download/install those softwares in README or in "
               "documentation."]
    with open(logfile, "r") as logf:
        for linef, linee in zip(logf, content):
            assert linee in linef
    # Check that needed packages are installed
    assert utils.is_package_installed("argparse")
    assert utils.is_package_installed("progressbar")
    assert utils.is_package_installed("numpy")
    assert utils.is_package_installed("matplotlib")
    assert utils.is_package_installed("Bio")
    os.remove(logfile)