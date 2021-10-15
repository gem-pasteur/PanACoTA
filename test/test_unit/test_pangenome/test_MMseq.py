#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the MMseq submodule in pangenome module
"""

import os
import sys
import time
import shutil
import glob
import logging
import pytest
import pytest_mock

from PanACoTA.pangenome_module.MMseq import MMseq
import PanACoTA.utils as utils
import PanACoTA.utils_pangenome as utils_pan
import test.test_unit.utilities_for_tests as tutil

LOGFILE_BASE = "logfile_test.txt"
LEVEL = logging.DEBUG
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]
# Define variables shared by several tests
PANDIR = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PANDIR, "test_files")
PATH_EXP_FILES = os.path.join(PANDIR, "exp_files")
GENEPATH = os.path.join(PANDIR, "generated_by_unit-tests")
PRT_BANK = "exp_EXEM.All.prt"
PRT_PATH = os.path.join(PATH_EXP_FILES, PRT_BANK)


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

# Clusters expected for the given databank (ref: [members]
EXP_CLUSTERS = {"GEN2.1017.00001.i0002_00004": ["GEN2.1017.00001.i0002_00004",
                                                "GEN4.1111.00001.i0001_00002",
                                                "GENO.1017.00001.b0002_00003",
                                                "GENO.1216.00002.i0001_00003"],
                "GEN2.1017.00001.b0003_00010": ["GEN2.1017.00001.b0003_00010"],
                "GEN2.1017.00001.b0004_00013": ["GEN2.1017.00001.b0004_00013"],
                "GEN4.1111.00001.b0001_00001": ["GEN2.1017.00001.i0002_00005",
                                                "GEN4.1111.00001.b0001_00001",
                                                "GENO.1017.00001.b0001_00002",
                                                "GENO.1216.00002.b0001_00001",
                                                "GENO.1216.00002.i0001_00002"],
                "GEN4.1111.00001.i0001_00003": ["GEN4.1111.00001.i0001_00003"],
                "GEN4.1111.00001.b0001_00009": ["GEN2.1017.00001.b0004_00011",
                                                "GEN4.1111.00001.b0001_00009",
                                                "GENO.1017.00001.b0002_00011",
                                                "GENO.1216.00002.b0002_00010"],
                "GENO.1017.00001.b0001_00001": ["GEN2.1017.00001.b0002_00006",
                                                "GENO.1017.00001.b0001_00001"],
                "GENO.1017.00001.i0002_00006": ["GEN2.1017.00001.b0003_00007",
                                                "GEN4.1111.00001.i0001_00006",
                                                "GENO.1017.00001.i0002_00006",
                                                "GENO.1017.00001.i0002_00007",
                                                "GENO.1216.00002.i0001_00007"],
                "GENO.1017.00001.i0002_00008": ["GEN2.1017.00001.i0004_00012",
                                                "GENO.1017.00001.i0002_00008"],
                "GENO.1017.00001.i0002_00009": ["GEN2.1017.00001.i0003_00008",
                                                "GEN4.1111.00001.i0001_00007",
                                                "GENO.1017.00001.i0002_00009",
                                                "GENO.1216.00002.b0001_00008"],
                "GENO.1017.00001.i0002_00010": ["GEN2.1017.00001.i0003_00009",
                                                "GEN4.1111.00001.i0001_00008",
                                                "GENO.1017.00001.i0002_00010",
                                                "GENO.1216.00002.b0002_00009"],
                "GENO.1216.00002.i0001_00004": ["GEN2.1017.00001.b0002_00003",
                                                "GENO.1216.00002.i0001_00004"],
                "GENO.1216.00002.i0001_00005": ["GEN2.1017.00001.b0001_00002",
                                                "GEN4.1111.00001.i0001_00004",
                                                "GENO.1017.00001.i0002_00004",
                                                "GENO.1216.00002.i0001_00005"],
                "GENO.1216.00002.i0001_00006": ["GEN2.1017.00001.b0001_00001",
                                                "GEN4.1111.00001.i0001_00005",
                                                "GENO.1017.00001.i0002_00005",
                                                "GENO.1216.00002.i0001_00006"],
                "GENO.1216.00002.b0003_00011": ["GENO.1216.00002.b0003_00011"],
                "GENO.1216.00002.b0003_00012": ["GENO.1216.00002.b0003_00012"]
                }

# protein families expected
FAMILIES4G = [["GEN2.1017.00001.i0002_00004", "GEN4.1111.00001.i0001_00002",
               "GENO.1017.00001.b0002_00003", "GENO.1216.00002.i0001_00003"],
              ["GEN2.1017.00001.b0003_00010"],
              ["GEN2.1017.00001.b0004_00013"],
              ["GEN2.1017.00001.i0002_00005", "GEN4.1111.00001.b0001_00001",
               "GENO.1017.00001.b0001_00002", "GENO.1216.00002.b0001_00001",
               "GENO.1216.00002.i0001_00002"],
              ["GEN4.1111.00001.i0001_00003"],
              ["GEN2.1017.00001.b0004_00011", "GEN4.1111.00001.b0001_00009",
               "GENO.1017.00001.b0002_00011", "GENO.1216.00002.b0002_00010"],
              ["GEN2.1017.00001.b0002_00006", "GENO.1017.00001.b0001_00001"],
              ["GEN2.1017.00001.b0003_00007", "GEN4.1111.00001.i0001_00006",
               "GENO.1017.00001.i0002_00006", "GENO.1017.00001.i0002_00007",
               "GENO.1216.00002.i0001_00007"],
              ["GEN2.1017.00001.i0004_00012", "GENO.1017.00001.i0002_00008"],
              ["GEN2.1017.00001.i0003_00008", "GEN4.1111.00001.i0001_00007",
               "GENO.1017.00001.i0002_00009", "GENO.1216.00002.b0001_00008"],
              ["GEN2.1017.00001.i0003_00009", "GEN4.1111.00001.i0001_00008",
               "GENO.1017.00001.i0002_00010", "GENO.1216.00002.b0002_00009"],
              ["GEN2.1017.00001.b0002_00003", "GENO.1216.00002.i0001_00004"],
              ["GEN2.1017.00001.b0001_00002", "GEN4.1111.00001.i0001_00004",
               "GENO.1017.00001.i0002_00004", "GENO.1216.00002.i0001_00005"],
              ["GEN2.1017.00001.b0001_00001", "GEN4.1111.00001.i0001_00005",
               "GENO.1017.00001.i0002_00005", "GENO.1216.00002.i0001_00006"],
              ["GENO.1216.00002.b0003_00011"],
              ["GENO.1216.00002.b0003_00012"]
              ]

FAMILIES4G_DICT = dict(zip(range(1, len(FAMILIES4G) + 1), FAMILIES4G))

def mmseq_runner(threads=1, panfile=None, quiet=False, min_id=0.8, mode=1):
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
    runner : MMseq
    """
    return MMseq(min_id, mode, GENEPATH, PRT_PATH, threads, panfile, quiet)

# tests
# Clusterisator properties
def test_thread():
    """
    Check if number of threads stored correctly
    """
    assert mmseq_runner().threads == 1

def test_prt_path():
    """
    check if prt_path stored correctly
    """
    assert mmseq_runner().prt_path == PRT_PATH

def test_quiet():
    """
    check if quiet is stored correctly
    """
    assert mmseq_runner(quiet=True).quiet

def test_logger():
    """
    check if logger is not disabled (we are not supposed to test logger here)
    """
    assert not mmseq_runner().logger.disabled

# Clusterisator properties, but values are mmseq-specific

def test_logpath():
    """
    Check if log path is computed correctly
    """
    runner = mmseq_runner()
    assert runner.log_path == os.path.join(GENEPATH, "mmseqs_" + PRT_BANK + "_" + runner.infoname + ".log")

def test_tmpdir():
    """
    Check if tmp dir is computed correctly
    """
    runner = mmseq_runner()
    assert runner.tmpdir == os.path.join(GENEPATH, "tmp_mmseqs_" + PRT_BANK + "_" + runner.infoname)

def test_infoname():
    """
    Check if infoname varies from different algorithm params (it should)
    """
    assert mmseq_runner().infoname != mmseq_runner(threads=2).infoname

def test_mmseqs_parser():
    """
    Test if mmseqs tsv file is correctly parsed to families
    """
    runner = mmseq_runner()
    runner.mmseqstsv = PATH_TEST_FILES + "/mmseq_clust-out.tsv"
    families = runner.parse_to_pangenome()
    assert set(map(lambda a: " ".join(a), families.values())) == set(map(lambda a: " ".join(a), FAMILIES4G))

def test_db_path():
    """
    Test if path to database is computed correctly
    """
    runner = mmseq_runner(threads=2)
    assert runner.mmseqdb == GENEPATH + f"/tmp_mmseqs_{PRT_BANK}_0.8-mode1-th2/{PRT_BANK}-msDB"

def test_clust_path():
    """
    Test if path to clustering results is computed correctly
    """
    runner = mmseq_runner(threads=2)
    assert runner.mmseqclust == GENEPATH + f"/tmp_mmseqs_{PRT_BANK}_0.8-mode1-th2/{PRT_BANK}-clust-0.8-mode1-th2"

# test run results

@pytest.fixture(params=["existing_panfile", "clust_files", "full_db",
                        "clear_folder", "incomplete_db", "clust_files_only"])
def test_cases(setup_teardown_module, request):
    """
    a fixture that returns multiple test cases, for each of which multiple tests are run
    Parameters
    ----------
    setup_teardown_module
    request

    Returns
    -------
    (case : str, runner : Clusterisator)
        case - string that defines the test case
        runner - runner for this case
    """
    tmp_dir = GENEPATH + f"/tmp_mmseqs_{PRT_BANK}_0.8-mode1-th2"

    db_ext = ["", ".index", ".lookup", ".dbtype", "_h", "_h.index", "_h.dbtype"]
    exp_db_base = PATH_TEST_FILES + "/mmseq_db"
    exp_db_files = list(map(lambda ext: exp_db_base + ext, db_ext))

    new_db_base = tmp_dir + "/" + PRT_BANK + "-msDB"
    new_db_files = list(map(lambda ext: new_db_base + ext, db_ext))

    clust_ext = [".0", ".1", ".dbtype", ".tsv"]

    exp_clust_base = PATH_TEST_FILES + "/mmseq_clust-out"
    exp_clust_files = list(map(lambda ext: exp_clust_base + ext, clust_ext))

    new_clust_base = tmp_dir + "/" + PRT_BANK + "-clust-0.8-mode1-th2"
    new_clust_files = list(map(lambda ext: new_clust_base + ext, clust_ext))
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    panfile_out = os.path.join(GENEPATH, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1-th2.lst")


    if request.param == "existing_panfile":
        shutil.copyfile(exp_pan, panfile_out)
        runner = mmseq_runner(threads=2)

    if request.param == "clust_files":
        # in this case databse and clust files do exist
        os.mkdir(tmp_dir)
        for old, new in zip(exp_db_files + exp_clust_files, new_db_files + new_clust_files):
            shutil.copyfile(old, new)

        runner = mmseq_runner(threads=2)

    if request.param == "full_db":
        os.mkdir(tmp_dir)
        for old, new in zip(exp_db_files, new_db_files):
            shutil.copyfile(old, new)

        runner = mmseq_runner(threads=2)
        print(runner.mmseqdb)

    if request.param == "incomplete_db":
        os.mkdir(tmp_dir)
        for old, new in zip(exp_db_files[:2], new_db_files[:2]):
            shutil.copyfile(old, new)

        runner = mmseq_runner(threads=2)

    if request.param == "clust_files_only":
        os.mkdir(tmp_dir)
        for old, new in zip(exp_clust_files, new_clust_files):
            shutil.copyfile(old, new)

        runner = mmseq_runner(threads=2)

    if request.param == "clear_folder":
        runner = mmseq_runner(threads=2)

    return (request.param, runner)


def read_callback(runner, mocker):
    expected_args = {
        "existing_panfile": (
        os.path.join(GENEPATH, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1-th2.lst"), runner.logger),
        "clust_files" : 0,
        "full_db" : 0,
        "clear_folder" : 0,
        "incomplete_db" : 0,
        "clust_files_only" : 0
    }
    return expected_args, mocker.spy(utils_pan, "read_pan_file")

def write_callback(runner, mocker):
    expected_args = {
        "existing_panfile": 0,
        "clust_files" : 1,
        "full_db" : 1,
        "clear_folder" : 1,
        "incomplete_db" : 1,
        "clust_files_only" : 1
    }

    return expected_args, mocker.spy(runner, "write_panfile")

def create_callback(runner, mocker):
    expected_args = {
        "existing_panfile": 0,
        "clust_files" : 1,
        "full_db" : 1,
        "clear_folder" : 1,
        "incomplete_db" : 1,
        "clust_files_only" : 1
    }

    return expected_args, mocker.spy(runner, "create_tmp_files")

def do_callback(runner, mocker):
    expected_args = {
        "existing_panfile": 0,
        "clust_files" : (False,),
        "full_db" : (False,),
        "clear_folder" : (True,),
        "incomplete_db" : (True,),
        "clust_files_only" : (True,)
    }

    return expected_args, mocker.spy(runner, "do_pangenome")

def parse_callback(runner, mocker):
    expected_args = {
        "existing_panfile": 0,
        "clust_files" : 1,
        "full_db" : 1,
        "clear_folder" : 1,
        "incomplete_db" : 1,
        "clust_files_only" : 1
    }
    return expected_args, mocker.spy(runner, "parse_to_pangenome")

def run_cmds_callback(runner, mocker):
    expected_args = {
        "existing_panfile": 0,
        "clust_files" : 0,
        "full_db" : 1,
        "clear_folder" : 2,
        "incomplete_db" : 2,
        "clust_files_only" : 2
    }
    return expected_args, mocker.spy(runner, "run_cmds")

def remove_callback(runner, mocker):
    expected_args = {
        "existing_panfile": 0,
        "clust_files" : 0,
        "full_db" : 0,
        "clear_folder" : 0,
        "incomplete_db" : 2,
        "clust_files_only" : 4
    }

    return expected_args, mocker.spy(utils, "remove")

@pytest.mark.parametrize("spy_producer", [read_callback, write_callback, create_callback, do_callback, parse_callback,
                                          run_cmds_callback, remove_callback])
def test_calls(test_cases, caplog, mocker, spy_producer):
    """
    Test if different functions are called expected number of times in all test cases, or
    are called with expected arguments. If item in dictonary is numeric, it is expected number of calls, if tuple -
    expected arguments (if they are known for sure).
    Parameters
    ----------
    test_cases
    caplog
    mocker
    spy_producer

    Returns
    -------

    """
    case, runner = test_cases
    expected_args, spy = spy_producer(runner, mocker)

    runner.run()

    if case in expected_args.keys():
        e = expected_args[case]

        if isinstance(e, tuple):
            spy.assert_called_once_with(*e)

        if isinstance(e, int):
            assert spy.call_count == e
    else:
        pytest.skip("expected value for this test case not defined")



@pytest.fixture(params=["already exists. PanACoTA will read it to get families.",
                        "Creating temporary files",
                        "Pangenome has 16 families",
                        "temporary files already exists. The program will use them.",
                        "Some, but not all",
                        "Removing existing temporary",
                        "Removing existing clustering and/or pangenome files",
                        "Clustering proteins...",
                        "already exists. The program will now convert it to a pangenome file.",
                        "Parsing "])
def log_phrases(request):
    return request.param

def test_logs(test_cases, caplog, log_phrases):
    """
    check if the runner prints phrases when it should and does not when it shold not
    """
    phrases_dict = {
        "already exists. PanACoTA will read it to get families." : {
            "existing_panfile" : True,
            "clust_files" : False,
            "full_db" : False,
            "clear_folder" : False,
            "incomplete_db" : False,
            "clust_files_only" : False
        },
        "Creating temporary files" : {
            "existing_panfile" : False,
            "clust_files" : False,
            "full_db" : False,
            "clear_folder" : True,
            "incomplete_db" : True,
            "clust_files_only" : True
        },
        "Pangenome has 16 families" : {
            "existing_panfile" : True,
            "clust_files" : True,
            "full_db" : True,
            "clear_folder" : True,
            "incomplete_db" : True,
            "clust_files_only" : True
        },
        "temporary files already exists. The program will use them." : {
            "existing_panfile" : False,
            "clust_files" : True,
            "full_db" : True,
            "clear_folder" : False,
            "incomplete_db" : False,
            "clust_files_only" : False
        },
        "Some, but not all" : {
            "existing_panfile" : False,
            "clust_files" : False,
            "full_db" : False,
            "clear_folder" : False,
            "incomplete_db" : True,
            "clust_files_only" : False
        },
        "Removing existing temporary" : {
            "existing_panfile" : False,
            "clust_files" : False,
            "full_db" : False,
            "clear_folder" : False,
            "incomplete_db" : True,
            "clust_files_only" : False
        },
        "Removing existing clustering and/or pangenome files" : {
            "existing_panfile" : False,
            "clust_files" : False,
            "full_db" : False,
            "clear_folder" : False,
            "incomplete_db" : False,
            "clust_files_only" : True
        },
        "Clustering proteins..." : {
            "existing_panfile" : False,
            "clust_files" : False,
            "full_db" : True,
            "clear_folder" : True,
            "incomplete_db" : True,
            "clust_files_only" : True
        },
        "already exists. The program will now convert it to a pangenome file." : {
            "existing_panfile" : False,
            "clust_files" : True,
            "full_db" : False,
            "clear_folder" : False,
            "incomplete_db" : False,
            "clust_files_only" : False
        },
        "Parsing " : {
            "existing_panfile" : False,
            "clust_files" : True,
            "full_db" : True,
            "clear_folder" : True,
            "incomplete_db" : True,
            "clust_files_only" : True
        }
    }
    case, runner = test_cases

    caplog.set_level(logging.DEBUG)
    runner.run()

    try :
        assert phrases_dict[log_phrases][case] == \
               (log_phrases in caplog.text.replace("\n", ""))
    except KeyError:
        pytest.skip("for this particular test case and phrase result not defined")

