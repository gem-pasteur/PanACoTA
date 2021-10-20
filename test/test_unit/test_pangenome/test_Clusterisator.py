#!/usr/bin/env python3
# coding: utf-8

"""
Functions which are useful for all Clusterisator subclasses testing
"""

import os
import sys
import time
import shutil
import glob
import logging
import pytest
import pytest_mock

import PanACoTA.utils as utils
import PanACoTA.utils_pangenome as utils_pan
import test.test_unit.utilities_for_tests as tutil
from PanACoTA.pangenome_module.MMseq import MMseq
from PanACoTA.pangenome_module.ProteinOrtho import ProteinOrtho

from abc import ABC, abstractmethod, abstractproperty

PANDIR = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PANDIR, "test_files")
PATH_EXP_FILES = os.path.join(PANDIR, "exp_files")
GENEPATH = os.path.join(PANDIR, "generated_by_unit-tests")
LOGFILE_BASE = "logfile_test.txt"
LEVEL = logging.DEBUG
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]

class TestClusterisator(ABC):
    __test__ = False

    @pytest.fixture
    def setup_teardown_module(self):
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
        # shutil.rmtree(GENEPATH)
        for f in LOGFILES:
           if os.path.exists(f):
                os.remove(f)
        print("teardown")

    @pytest.fixture(params=["existing_panfile", "clust_files", "full_db",
                            "clear_folder", "incomplete_db", "clust_files_only"])
    def test_cases(self, setup_teardown_module, request):
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
        if request.param == "existing_panfile":
            shutil.copyfile(*self.exp_new_panfile)
            runner = self.clust_runner(threads=2)

        if request.param == "clust_files":
            # in this case databse and clust files do exist
            os.mkdir(self.TMP_DIR)
            for old, new in self.exp_new_db_files + self.exp_new_clust_files:
                shutil.copyfile(old, new)

            runner = self.clust_runner(threads=2)

        if request.param == "full_db":
            os.mkdir(self.TMP_DIR)
            for old, new in self.exp_new_db_files:
                shutil.copyfile(old, new)

            runner = self.clust_runner(threads=2)

        if request.param == "incomplete_db":
            os.mkdir(self.TMP_DIR)
            for old, new in self.exp_new_db_files[1:]:
                shutil.copyfile(old, new)

            runner = self.clust_runner(threads=2)

        if request.param == "clust_files_only":
            os.mkdir(self.TMP_DIR)
            for old, new in self.exp_new_clust_files:
                shutil.copyfile(old, new)

            runner = self.clust_runner(threads=2)

        if request.param == "clear_folder":
            runner = self.clust_runner(threads=2)

        return (request.param, runner)

    @property
    @abstractmethod
    def TMP_DIR(self):
        pass

    @property
    @abstractmethod
    def exp_new_db_files(self):
        """
        List of pairs (location of expected file, location in temporary folder)
        """
        pass

    @property
    @abstractmethod
    def exp_new_clust_files(self):
        """
        List of pairs (location of expected file, location in temporary folder)
        """
        pass

    @property
    @abstractmethod
    def exp_new_panfile(self):
        """
        (location of expected file, location in temporary folder)
        """
        pass


    @abstractmethod
    def clust_runner(self, threads=1, **kwargs):
        pass

    def test_thread(self):
        """
        Check if number of threads stored correctly
        """
        assert self.clust_runner().threads == 1

    def test_quiet(self):
        """
        check if quiet is stored correctly
        """
        assert self.clust_runner(quiet=True).quiet

    def test_logger(self):
        """
        check if logger is not disabled (we are not supposed to test logger here)
        """
        assert not self.clust_runner().logger.disabled

    def test_infoname(self):
        """
        Check if infoname varies from different algorithm params (it should)
        """
        assert self.clust_runner().infoname != self.clust_runner(threads=2).infoname

    def read_callback(self, runner, mocker):
        expected_args = {
            "existing_panfile": (self.exp_new_panfile[1], runner.logger),
            "clust_files": 0,
            "full_db": 0,
            "clear_folder": 0,
            "incomplete_db": 0,
            "clust_files_only": 0
        }
        return expected_args, mocker.spy(utils_pan, "read_pan_file")

    def write_callback(self, runner, mocker):
        expected_args = {
            "existing_panfile": 0,
            "clust_files": 1,
            "full_db": 1,
            "clear_folder": 1,
            "incomplete_db": 1,
            "clust_files_only": 1
        }

        return expected_args, mocker.spy(runner, "write_panfile")

    def create_callback(self, runner, mocker):
        expected_args = {
            "existing_panfile": 0,
            "clust_files": 1,
            "full_db": 1,
            "clear_folder": 1,
            "incomplete_db": 1,
            "clust_files_only": 1
        }

        return expected_args, mocker.spy(runner, "create_tmp_files")

    def do_callback(self, runner, mocker):
        expected_args = {
            "existing_panfile": 0,
            "clust_files": (False,),
            "full_db": (False,),
            "clear_folder": (True,),
            "incomplete_db": (True,),
            "clust_files_only": (True,)
        }

        return expected_args, mocker.spy(runner, "do_pangenome")

    def parse_callback(self, runner, mocker):
        expected_args = {
            "existing_panfile": 0,
            "clust_files": 1,
            "full_db": 1,
            "clear_folder": 1,
            "incomplete_db": 1,
            "clust_files_only": 1
        }
        return expected_args, mocker.spy(runner, "parse_to_pangenome")

    def run_cmds_callback(self, runner, mocker):
        expected_args = {
            "existing_panfile": 0,
            "clust_files": 0,
            "full_db": 1,
            "clear_folder": 2,
            "incomplete_db": 2,
            "clust_files_only": 2
        }
        return expected_args, mocker.spy(runner, "run_cmds")

    def remove_callback(self, runner, mocker):
        expected_args = {
            "existing_panfile": 0,
            "clust_files": 0,
            "full_db": 0,
            "clear_folder": 0,
            "incomplete_db": len(self.exp_new_db_files) - 1,
            "clust_files_only": 4
        }

        return expected_args, mocker.spy(utils, "remove")

    @pytest.mark.parametrize("spy_producer",
                             [read_callback, write_callback, create_callback, do_callback, parse_callback,
                              run_cmds_callback, remove_callback])
    def test_calls(self, test_cases, caplog, mocker, spy_producer):
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
        expected_args, spy = spy_producer(self, runner, mocker)

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
                            "Pangenome has",
                            "temporary files already exists. The program will use them.",
                            "Some, but not all",
                            "Removing existing temporary",
                            "Removing existing clustering and/or pangenome files",
                            "Clustering proteins...",
                            "already exists. The program will now convert it to a pangenome file.",
                            "Parsing "])
    def log_phrases(self, request):
        return request.param

    def test_logs(self, test_cases, caplog, log_phrases):
        """
        check if the runner prints phrases when it should and does not when it shold not
        """
        phrases_dict = {
            "already exists. PanACoTA will read it to get families.": {
                "existing_panfile": True,
                "clust_files": False,
                "full_db": False,
                "clear_folder": False,
                "incomplete_db": False,
                "clust_files_only": False
            },
            "Creating temporary files": {
                "existing_panfile": False,
                "clust_files": False,
                "full_db": False,
                "clear_folder": True,
                "incomplete_db": True,
                "clust_files_only": True
            },
            "Pangenome has": {
                "existing_panfile": True,
                "clust_files": True,
                "full_db": True,
                "clear_folder": True,
                "incomplete_db": True,
                "clust_files_only": True
            },
            "temporary files already exists. The program will use them.": {
                "existing_panfile": False,
                "clust_files": True,
                "full_db": True,
                "clear_folder": False,
                "incomplete_db": False,
                "clust_files_only": False
            },
            "Some, but not all": {
                "existing_panfile": False,
                "clust_files": False,
                "full_db": False,
                "clear_folder": False,
                "incomplete_db": True,
                "clust_files_only": False
            },
            "Removing existing temporary": {
                "existing_panfile": False,
                "clust_files": False,
                "full_db": False,
                "clear_folder": False,
                "incomplete_db": True,
                "clust_files_only": False
            },
            "Removing existing clustering and/or pangenome files": {
                "existing_panfile": False,
                "clust_files": False,
                "full_db": False,
                "clear_folder": False,
                "incomplete_db": False,
                "clust_files_only": True
            },
            "Clustering proteins...": {
                "existing_panfile": False,
                "clust_files": False,
                "full_db": True,
                "clear_folder": True,
                "incomplete_db": True,
                "clust_files_only": True
            },
            "already exists. The program will now convert it to a pangenome file.": {
                "existing_panfile": False,
                "clust_files": True,
                "full_db": False,
                "clear_folder": False,
                "incomplete_db": False,
                "clust_files_only": False
            },
            "Parsing ": {
                "existing_panfile": False,
                "clust_files": True,
                "full_db": True,
                "clear_folder": True,
                "incomplete_db": True,
                "clust_files_only": True
            }
        }
        case, runner = test_cases

        caplog.set_level(logging.DEBUG)
        runner.run()

        try:
            assert phrases_dict[log_phrases][case] == \
                   (log_phrases in caplog.text.replace("\n", ""))
        except KeyError:
            pytest.skip("for this particular test case and phrase result not defined")

    @abstractmethod
    def test_parse(self):
        pass


class TestMMseq(TestClusterisator):
    __test__ = True
    PRT_BANK = "exp_EXEM.All.prt"
    PRT_PATH = os.path.join(PATH_EXP_FILES, PRT_BANK)

    @property
    def TMP_DIR(self):
        return GENEPATH + f"/tmp_mmseqs_{self.PRT_BANK}_0.8-mode1-th2"

    @property
    def exp_new_db_files(self):
        db_ext = ["", ".index", ".lookup", ".dbtype", "_h", "_h.index", "_h.dbtype"]
        exp_db_base = PATH_TEST_FILES + "/mmseq_db"
        exp_db_files = list(map(lambda ext: exp_db_base + ext, db_ext))

        new_db_base = self.TMP_DIR + "/" + self.PRT_BANK + "-msDB"
        new_db_files = list(map(lambda ext: new_db_base + ext, db_ext))

        return [tuple(a) for a in zip(exp_db_files, new_db_files)]

    @property
    def exp_new_clust_files(self):
        clust_ext = [".0", ".1", ".dbtype", ".tsv"]

        exp_clust_base = PATH_TEST_FILES + "/mmseq_clust-out"
        exp_clust_files = list(map(lambda ext: exp_clust_base + ext, clust_ext))

        new_clust_base = self.TMP_DIR + "/" + self.PRT_BANK + "-clust-0.8-mode1-th2"
        new_clust_files = list(map(lambda ext: new_clust_base + ext, clust_ext))
        return [tuple(a) for a in zip(exp_clust_files, new_clust_files)]

    @property
    def exp_new_panfile(self):
        exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
        panfile_out = os.path.join(GENEPATH, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1-th2.lst")
        return exp_pan, panfile_out

    def clust_runner(self, threads=1, **kwargs):
        min_id = kwargs.get("min_id", 0.8)
        panfile = kwargs.get("panfile", None)
        quiet = kwargs.get("quiet", False)
        mode = kwargs.get("mode", 1)
        return MMseq(min_id, mode, GENEPATH, self.PRT_PATH, threads, panfile, quiet)

    def test_prt_path(self):
        """
        check if prt_path stored correctly
        """
        assert self.clust_runner().prt_path == self.PRT_PATH

    def test_logpath(self):
        """
        Check if log path is computed correctly
        """
        runner = self.clust_runner()
        assert runner.log_path == os.path.join(GENEPATH, "mmseqs_" + self.PRT_BANK + "_" + runner.infoname + ".log")

    def test_tmpdir(self):
        """
        Check if tmp dir is computed correctly
        """
        runner = self.clust_runner()
        assert runner.tmpdir == os.path.join(GENEPATH, "tmp_mmseqs_" + self.PRT_BANK + "_" + runner.infoname)

    def test_clust_path(self):
        """
        Test if path to clustering results is computed correctly
        """
        runner = self.clust_runner(threads=2)
        assert runner.mmseqclust == f"{self.TMP_DIR}/{self.PRT_BANK}-clust-0.8-mode1-th2"

    def test_parse(self):
        runner = self.clust_runner(threads=2)
        runner.mmseqstsv = os.path.join(PATH_TEST_FILES, "mmseq_clust-out.tsv")

        _, exp_fams, _ = utils_pan.read_pan_file(self.exp_new_panfile[0], runner.logger)

        families = runner.parse_to_pangenome()
        assert set(map(lambda a: " ".join(a), families.values())) == set(map(lambda a: " ".join(a), exp_fams.values()))


class TestProteinOrtho(TestClusterisator):
    __test__ = True
    PRT_PATH = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    NAME = "EXEM"

    @property
    def TMP_DIR(self):
        return GENEPATH + f"/tmp_proteinortho_{self.NAME}-All_diamond-search-th2"

    def clust_runner(self, threads=1, **kwargs):
        panfile = kwargs.get("panfile", None)
        quiet = kwargs.get("quiet", False)
        po_mode = kwargs.get("po_mode", "diamond")
        name = kwargs.get("name", self.NAME)
        return ProteinOrtho(po_mode, name, GENEPATH, os.path.join(self.PRT_PATH, f"{self.NAME}-All"),
                            threads, panfile, quiet)

    @property
    def exp_new_panfile(self):
        exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes-proteinortho.lst")
        panfile_out = os.path.join(GENEPATH, "PanGenome-EXEM-clust-diamond-search-th2.lst")
        return exp_pan, panfile_out

    @property
    def exp_new_clust_files(self):
        clust_ext = ["proteinortho-graph", "proteinortho-graph.summary", "proteinortho.html", "proteinortho.tsv"]

        exp_clust_base = PATH_TEST_FILES + "/proteinortho_clust-out"
        exp_clust_files = list(map(lambda ext: exp_clust_base + "." + ext, clust_ext))

        new_clust_base = self.TMP_DIR + "/" + self.NAME
        new_clust_files = list(map(lambda ext: new_clust_base + "." + ext, clust_ext))
        return [tuple(a) for a in zip(exp_clust_files, new_clust_files)]

    @property
    def exp_new_db_files(self):
        db_ext = ["info", "blast-graph"]

        exp_db_base = PATH_TEST_FILES + "/proteinortho_db"
        exp_db_files = list(map(lambda ext: exp_db_base + "." + ext, db_ext))

        new_db_base = self.TMP_DIR + "/" + self.NAME
        new_db_files = list(map(lambda ext: new_db_base + "." + ext, db_ext))
        return [tuple(a) for a in zip(exp_db_files, new_db_files)]

    def test_parse(self):
        runner = self.clust_runner(threads=2, name="proteinortho_clust-out")
        runner.tmpdir = PATH_TEST_FILES

        _, exp_fams, _ = utils_pan.read_pan_file(self.exp_new_panfile[0], runner.logger)
        families = runner.parse_to_pangenome()
        assert set(map(lambda a: " ".join(a), families.values())) == set(map(lambda a: " ".join(a), exp_fams.values()))
