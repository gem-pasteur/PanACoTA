#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pipelinepackage.utils as utils
import pytest
import os
import logging
# import util_tests


def test_check_install():
    """
    Try to run prokka, which is installed, and check that there is no problem
    """
    ret = utils.check_installed(["prokka"])
    assert ret == 1


def test_check_install_error(capsys):
    """
    Try to run a command which does not exist, and check that it closes the program
    with exit code 1
    """
    with pytest.raises(SystemExit):
        utils.check_installed(["plop"])
    out, err = capsys.readouterr()
    assert err.startswith("plop failed:")


def test_class_filter():
    """
    Check that for a class LessThanFilter(warning), info and debug are allowed,
    but warning, error and critical are not.
    """
    a = utils.LessThanFilter(logging.WARNING)
    record = logging.LogRecord("root", logging.DEBUG, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.INFO, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.WARNING, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)
    record = logging.LogRecord("root", logging.ERROR, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)
    record = logging.LogRecord("root", logging.CRITICAL, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)

# def test_plot_dist():
#     """
#     Plot a given distribution, and check that output is as expected
#     """
#     values = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 10]
#     limit = 3
#     res_dir = os.path.join("test", "data", "tests_results")
#     os.makedirs(res_dir, exist_ok=True)
#     outfile = os.path.join(res_dir, "distrib.png")
#     reffile = os.path.join("test", "data", "res_plot_distr.png")
#     title = "Distribution test"
#     text = "Max L90 ="
#     utils.plot_distr(values, limit, outfile, title, text)
#     assert util_tests.compare_files(outfile, reffile)
#     os.remove(outfile)