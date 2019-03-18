#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing genomeAPCAT according to already existing dependencies
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
    Uninstall genomeAPCAT and installed dependencies
    """
    print("TEARDOWN\n")
    cmd = "python3 make uninstall"
    error = "Error uninstall"
    utils.run_cmd(cmd, error)
    os.remove("install.log")
    print("cleaning repo")


def test_install_panacota():
    """
    Test that when installing from a computer containing only all dependencies, it returns a message without any warning: everything is ok
    """
    cmd = "python3 make"
    error = "Error trying to install genomeAPCAT from base"
    assert utils.check_installed("barrnap")
    assert not utils.check_installed("genomeAPCAT")
    assert utils.check_installed("prokka")
    # Install panacota
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("genomeAPCAT")
    # Check that panacota is installed (pip3 module exists)
    assert utils.is_package_installed("argparse")
    # Check that it is installed in "final mode"
    cmd = "pip3 show genomeAPCAT"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    with open(stdout, "r") as stdof:
        for line in stdof:
            if line.startswith("Location"):
                loc = line.split()[-1]
                assert glob.glob(os.path.join(loc, r'genomeAPCAT*dist-info'))
                loc = line.split()[-1]
    os.remove(stdout)
    logfile = "install.log"
    content = ["Installing genomeAPCAT...", "DONE"]
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

def test_test():
    genomeAPCAT -h
    print("toto")

    # assert utils.check_installed("genomeAPCAT")
# def test_upgrade(install_panacota):
#     """
#     Test upgrading genomeAPCAT when dependencies are still installed
#     """
#     install_panacota
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert utils.check_installed("genomeAPCAT")
#     cmd = "python3 make upgrade"
#     error = "Error upgrade"
#     utils.run_cmd(cmd, error)
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert utils.check_installed("genomeAPCAT")
#     logfile = "install.log"
#     with open(logfile, "r") as logf:
#         lines = logf.readlines()
#         assert len(lines) == 2
#         assert "Upgrading genomeAPCAT" in lines[0]
#         assert "DONE" in lines[1]


# def test_uninstall(install_panacota):
#     """
#     Test uninstalling genomeAPCAT when dependencies are still installed
#     """
#     install_panacota
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert utils.check_installed("genomeAPCAT")
#     cmd = "python3 make uninstall"
#     error = "Error uninstalling"
#     utils.run_cmd(cmd, error)
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert not utils.check_installed("genomeAPCAT")
#     logfile = "install.log"
#     with open(logfile, "r") as logf:
#         lines = logf.readlines()
#         assert len(lines) == 2
#         assert "Uninstalling genomeAPCAT" in lines[0]
#         assert "DONE" in lines[1]

# def test_develop():
#     """
#     Test installing genomeAPCAT in developer mode, when prokka and barrnap are already installed
#     """
#     assert not utils.check_installed("genomeAPCAT")
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     cmd = "python3 make develop"
#     error = "Error develop"
#     utils.run_cmd(cmd, error)
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert utils.check_installed("genomeAPCAT")
#     cmd = "pip3 show genomeAPCAT"
#     err = "error pip3"
#     stdout = "stdout_pip3show.out"
#     with open(stdout, "w") as stdof:
#         utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
#     # Check installation in develop mode
#     with open(stdout, "r") as stdof:
#         lines = stdof.readlines()
#         for line in lines:
#             if line.startswith("Location"):
#                 loc = line.split()[-1]
#                 assert glob.glob(os.path.join(loc, r'genomeAPCAT*egg-info'))
#     os.remove(stdout)
#     logfile = "install.log"
#     content = ["Installing genomeAPCAT...",
#                "Installing developer packages needed for genomeAPCAT", "DONE"]
#     with open(logfile, "r") as logf:
#         logf_content = "".join(logf.readlines())
#         for linec in content:
#             assert linec in logf_content
#     # Check that needed packages are installed
#     assert utils.is_package_installed("argparse")
#     assert utils.is_package_installed("progressbar")
#     assert utils.is_package_installed("numpy")
#     assert utils.is_package_installed("matplotlib")
#     assert utils.is_package_installed("Bio")
#     assert utils.is_package_installed("sphinx")
#     assert utils.is_package_installed("numpydoc")
#     assert utils.is_package_installed("pytest")
#     assert utils.is_package_installed("pytest-mpl")
#     assert utils.is_package_installed("coverage")


# def test_upgrade_notinstalled():
#     """
#     Test upgrading genomeAPCAT when dependencies are not installed (only barrnap),
#     and genomeAPCAT is not installed. It just installs genomeAPCAT, without prokka dep
#     """
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert not utils.check_installed("genomeAPCAT")
#     cmd = "python3 make upgrade"
#     error = "Error upgrade"
#     utils.run_cmd(cmd, error)
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert utils.check_installed("genomeAPCAT")
#     logfile = "install.log"
#     with open(logfile, "r") as logf:
#         lines = logf.readlines()
#         assert len(lines) == 2
#         assert "Upgrading genomeAPCAT" in lines[0]
#         assert "DONE" in lines[1]
#     os.remove(logfile)


# def test_install_user():
#     """
#     Test that when installing from a computer in user mode, it installs
#     genomeAPCAT in /Users and returns list of dependencies
#     """
#     cmd = "python3 make --user"
#     error = "Error trying to install genomeAPCAT from base"
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert not utils.check_installed("genomeAPCAT")
#     utils.run_cmd(cmd, error)
#     assert utils.check_installed("barrnap")
#     assert utils.check_installed("prokka")
#     assert utils.check_installed("genomeAPCAT")
#     cmd = "pip3 show genomeAPCAT"
#     err = "error pip3"
#     stdout = "stdout_pip3show.out"
#     with open(stdout, "w") as stdof:
#         utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
#     with open(stdout, "r") as stdof:
#         lines = stdof.readlines()
#         for line in lines:
#             assert "/usr" not in line
#             if line.startswith("Location"):
#                 loc = line.split()[-1]
#                 # assert glob.glob(os.path.join(loc, r'genomeAPCAT*dist-info'))
#                 print(loc)
#                 print(os.listdir(loc))
#     os.remove(stdout)
#     logfile = "install.log"
#     content = ["Installing genomeAPCAT in user mode...", "DONE"]
#     # Check output logfile content. Check that all content is present, in any order.
#     with open(logfile, "r") as logf:
#         logf_content = "".join(logf.readlines())
#         for linec in content:
#             assert linec in logf_content
#     # Check that needed packages are installed
#     assert utils.is_package_installed("argparse")
#     assert utils.is_package_installed("progressbar")
#     assert utils.is_package_installed("numpy")
#     assert utils.is_package_installed("matplotlib")
#     assert utils.is_package_installed("Bio")