#!/usr/bin/env python3
# coding: utf-8

"""
Tests for make script, installing PanACoTA according to already existing dependencies
Here, all dependencies are installed
"""
import os
import glob
import pytest
import subprocess
import shlex

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


def test_install_panacota():
    """
    Test that when installing from a computer containing only all dependencies, it returns a message without any warning: everything is ok
    """
    cmd = "python3 make"
    error = "Error trying to install PanACoTA from base"
    assert not utils.check_installed("PanACoTA")
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("mash")
    assert utils.check_installed("FastTreeMP")
    # Install panacota
    utils.run_cmd(cmd, error)
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("mash")
    assert utils.check_installed("FastTreeMP")
    # Check that panacota is installed (pip3 module exists)
    assert utils.is_package_installed("PanACoTA")
    # Check that it is installed in "final mode"
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
    # Check not installed in user mode
    cmd = "pip list --user"
    list_user_packages = str(subprocess.check_output(shlex.split(cmd)))
    assert "PanACoTA" not in list_user_packages
    # Check output logfile content. Check that all content is present, in any order.
    logfile = "install.log"
    content = ["Installing PanACoTA...", "DONE"]
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
    print("test_installed done")


def test_upgrade(install_panacota):
    """
    Test upgrading PanACoTA when dependencies are still installed
    """
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("FastTreeMP")
    assert utils.check_installed("PanACoTA")
    # Upgrade PanACoTA
    cmd = "python3 make upgrade"
    error = "Error upgrade"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed("FastTreeMP")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 2
        assert "Upgrading PanACoTA" in lines[0]
        assert "DONE" in lines[1]
    print("test_upgrade done")


def test_upgrade_notinstalled():
    """
    Test upgrading PanACoTA when dependencies are not installed (only barrnap),
    and PanACoTA is not installed. It just installs PanACoTA, without prokka dep
    """
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("FastTreeMP")
    assert not utils.check_installed("PanACoTA")
    cmd = "python3 make upgrade"
    error = "Error upgrade"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed("FastTreeMP")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 2
        assert "Upgrading PanACoTA" in lines[0]
        assert "DONE" in lines[1]
    os.remove(logfile)
    print("test_upgrade_notinstalled done")


def test_uninstall(install_panacota):
    """
    Test uninstalling PanACoTA when dependencies are still installed
    """
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed("FastTreeMP")
    cmd = "python3 make uninstall"
    error = "Error uninstalling"
    utils.run_cmd(cmd, error)
    assert not utils.check_installed("PanACoTA")
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("mash")
    assert utils.check_installed("FastTreeMP")
    assert utils.check_installed("quicktree")
    assert utils.check_installed("iqtree") or utils.check_installed("iqtree2")
    logfile = "install.log"
    with open(logfile, "r") as logf:
        lines = logf.readlines()
        assert len(lines) == 2
        assert "Uninstalling PanACoTA" in lines[0]
        assert "DONE" in lines[1]
    print("test_uninstall done")


def test_develop():
    """
    Test installing PanACoTA in developer mode, when prokka and barrnap are already installed
    """
    assert not utils.check_installed("PanACoTA")
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    cmd = "python3 make develop"
    error = "Error develop"
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    # Check installed in developper mode (egg-info file is present)
    cmd = "pip3 show PanACoTA"
    err = "error pip3"
    stdout = "stdout_pip3show.out"
    with open(stdout, "w") as stdof:
        utils.run_cmd(cmd, err, stdout=stdof, stderr=stdof)
    # Check installation in develop mode
    with open(stdout, "r") as stdof:
        lines = stdof.readlines()
        for line in lines:
            if line.startswith("Location"):
                loc = line.split()[-1]
                assert glob.glob(os.path.join(loc, r'PanACoTA*egg-info'))
    os.remove(stdout)
    # Check not installed in user mode
    cmd = "pip list --user"
    list_user_packages = str(subprocess.check_output(shlex.split(cmd)))
    assert "PanACoTA" not in list_user_packages
    # check logfile content
    logfile = "install.log"
    content = ["Installing PanACoTA...",
               "Installing developer packages needed for PanACoTA", "DONE"]
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
    assert utils.is_package_installed("numpydoc")
    assert utils.is_package_installed("pytest")
    assert utils.is_package_installed("coverage")
    print("test_develop done")


def test_install_user():
    """
    Test that when installing from a computer in user mode, it really installs
    PanACoTA in user mode (pip showing PanACoTA when asking for user installed packages)
    """
    cmd = "python3 make --user"
    error = "Error trying to install PanACoTA from base"
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert not utils.check_installed("PanACoTA")
    assert utils.check_installed("FastTreeMP")
    utils.run_cmd(cmd, error)
    assert utils.check_installed("barrnap")
    assert utils.check_installed("prokka")
    assert utils.check_installed("PanACoTA")
    assert utils.check_installed("FastTreeMP")
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
