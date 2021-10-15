#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the Clusterisator submodule in pangenome module
"""

import os
import time
import shutil
import glob
import logging
import pytest

from PanACoTA.pangenome_module.ProteinOrtho import ProteinOrtho
import PanACoTA.utils as utils
import test.test_unit.utilities_for_tests as tutil
import glob

LOGFILE_BASE = "logfile_test.txt"
LEVEL = logging.DEBUG
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]
# Define variables shared by several tests
PANDIR = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PANDIR, "test_files")
PATH_EXP_FILES = os.path.join(PANDIR, "exp_files")
GENEPATH = os.path.join(PANDIR, "generated_by_unit-tests")
PRT_PATH = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")

@pytest.hookimpl
def pytest_configure(config):
    logging_plugin = config.pluginmanager.get_plugin("logging-plugin")

    # Change color on existing log level
    logging_plugin.log_cli_handler.formatter.add_color_level(logging.INFO, "cyan")

    # Add color to a custom log level (a custom log level `SPAM` is already set up)
    logging_plugin.log_cli_handler.formatter.add_color_level(logging.SPAM, "blue")


def proteinortho_runner(threads=2, panfile=None, quiet=False, po_mode="diamond", name="EXEM"):
    """

    Parameters
    ----------
    threads : int
        number of threads
    panfile : None or str
        a file to write pangenome. If none, defult is used.
    quiet : bool
        If quiet, the logging is reduced.

    Returns
    -------
    runner : ProteinOrtho
    """
    return ProteinOrtho(po_mode, name, GENEPATH, PRT_PATH, threads, panfile, quiet)

@pytest.fixture
def setup_teardown_module():
    """
    Remove log files at the end of this test module

    Before each test:
    - init logger
    - create directory to put generated files

    After:
    - remove all log files
    - remove directory with generated results
    """
    utils.init_logger(LOGFILE_BASE, logging.DEBUG, 'test_getseq', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    print("teardown")


def test_run(setup_teardown_module):
    runner = proteinortho_runner(threads=2)
    # pytest.skip(runner.tmp_files_cmds[0][0])
    families, panfile = runner.run()
    print(families)


