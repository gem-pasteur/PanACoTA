#!/usr/bin/env python3
# coding: utf-8

"""
Functions used for tests
"""
import importlib.util
import logging
import shlex
import subprocess

import sys

import os


def run_cmd(cmd, error, eof=False, **kwargs):
    """
    Run the given command line. If the return code is not 0, print error message.
    if eof (exit on fail) is True, exit program if error code is not 0.

    Parameters
    ----------
    cmd : str
        command to run
    error : str
        error message to print if error while running command
    eof : bool
        True: exit program if command failed, False: do not exit even if command fails
    kwargs : Object
        Can provide a logger, stdout and/or stderr streams

    Returns
    -------
    subprocess.Popen
        returns object of subprocess call (has attributes returncode, pid, communicate etc.)

    """
    if "logger" not in kwargs:
        logger = logging.getLogger("utils.run_cmd")
    else:
        logger = kwargs["logger"]
    if "stdout" not in kwargs:
        kwargs["stdout"] = None
    if "stderr" not in kwargs:
        kwargs["stderr"] = None
    try:
        call = subprocess.Popen(shlex.split(cmd), stdout=kwargs["stdout"],
                                stderr=kwargs["stderr"])
        call.wait()
        retcode = call.returncode
    except OSError:
        logger.error(error + ": " + "{} does not exist".format(cmd))
        if eof:
            sys.exit(-1)
        else:
            return -1
    if retcode != 0:
        logger.error(error)
        if eof:
            sys.exit(retcode)
    return call


def is_package_installed(pkg_name):
    """
    Check if package is installed or not

    Parameters
    ----------
    pkg_name : str
        name of package to check

    Returns
    -------
    bool
        True if package is installed, False otherwise
    """
    spec = ""
    spec = importlib.util.find_spec(pkg_name)
    return spec is not None


def check_installed(cmd):
    """
    Check if the command 'cmd' is in $PATH and can then be executed

    Parameters
    ----------
    cmd : str
        command to run

    Returns
    -------
    bool
        True if installed, False otherwise
    """
    torun = "which " + cmd
    trying = subprocess.Popen(shlex.split(torun), stdout=subprocess.PIPE)
    out, _ = trying.communicate()
    if trying.returncode == 0:
        if os.path.isfile(out.strip()):
            return True
    return False
