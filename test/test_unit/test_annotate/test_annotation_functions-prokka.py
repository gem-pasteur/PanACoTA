#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for annotate/annotation_functions.py
"""

import pytest
import os
import logging
import shutil
import multiprocessing

import test.test_unit.utilities_for_tests as tutil
import PanACoTA.utils as utils
from PanACoTA.annotate_module import annotation_functions as afunc


# Define variables used by several tests
DBDIR = os.path.join("test", "data", "annotate")
GEN_PATH = os.path.join(DBDIR, "genomes")
TMP_PATH = os.path.join(DBDIR, "tmp_files")
EXP_DIR = os.path.join(DBDIR, 'exp_files')
TEST_DIR = os.path.join(DBDIR, 'test_files')
GENEPATH = os.path.join(DBDIR, "generated_by_unit-tests")
LOGFILE_BASE = os.path.join(GENEPATH, "logfile.log")
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]


@pytest.fixture(autouse=True)
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
    if os.path.isdir(GENEPATH):
        content = os.listdir(GENEPATH)
        for f in content:
            assert f.startswith(".fuse")
    else:
        os.mkdir(GENEPATH)
    print("setup")

    yield
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")


# Create a logger with sublogger inside
def my_logger(name):
    """
    logger given to function called by a subprocess
    """
    m = multiprocessing.Manager()
    q = m.Queue()
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    return q, logging.getLogger(name)


def test_check_prokka_no_outdir():
    """
    Test that prokka returns the right error message when output directory does not exist
    """
    logger = my_logger("test_check_prokka_no_outdir")
    outdir = os.path.join(GENEPATH, "outdir")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = os.path.join(GENEPATH, "toto.fna")
    nbcont = 7
    assert not afunc.check_prokka(outdir, logf, name, gpath, nbcont, logger[1])
    q = logger[0]
    assert q.qsize() == 1
    msg = "Previous annotation could not run properly. Look at prokka.log for more information."
    assert q.get().message == msg


def test_check_prokka_nofna():
    """
    Check that check_prokka returns false when a tbl file is missing, and an error message
    """
    logger = my_logger("test_check_prokka_nofna")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_nofna")
    name = "prokka_out_for_test-missfna"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missfna original_name-error.fna: no .fna file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_sevfna():
    """
    Check that check_prokka returns false when there is more than 1 tbl file,
    and an error message
    """
    logger = my_logger("test_check_prokka_sevfna")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_sevfna")
    name = "prokka_out_for_test-sevfna"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + "2.fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-sevfna original_name-error.fna: several .fna files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_notbl():
    """
    Check that check_prokka returns false when a tbl file is missing, and an error message
    """
    logger = my_logger("test_check_prokka_notbl")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_notbl")
    name = "prokka_out_for_test-misstbl"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-misstbl original_name-error.fna: no .tbl file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_sevtbl():
    """
    Check that check_prokka returns false when there is more than 1 tbl file,
    and an error message
    """
    logger = my_logger("test_check_prokka_sevtbl")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_notbl")
    name = "prokka_out_for_test-misstbl"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + "2.tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-misstbl original_name-error.fna: several .tbl files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_nofaa():
    """
    Check that check_prokka returns false when a faa file is missing, and an error message
    """
    logger = my_logger("test_check_prokka_nofaa")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_nofaa")
    name = "prokka_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missfaa original_name.fna: no .faa file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_sevfaa():
    """
    Check that check_prokka returns false when there is more than 1 faa file,
    and an error message
    """
    logger = my_logger("test_check_prokka_sevfaa")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_nofaa")
    name = "prokka_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + "2.faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missfaa original_name.fna: several .faa files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_noffn():
    """
    Check that check_prokka returns false when a ffn file is missing, and an error message
    """
    logger = my_logger("test_check_prokka_noffn")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prokka_out_for_test-missffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missffn original_name.fna: no .ffn file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_sevffn():
    """
    Check that check_prokka returns false when there is more than 1 ffn file,
    and an error message
    """
    logger = my_logger("test_check_prokka_sevffn")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prokka_out_for_test-missffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + "2.ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missffn original_name.fna: several .ffn files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_nogff():
    """
    Check that check_prokka returns false when a ffn file is missing, and an error message
    """
    logger = my_logger("test_check_prokka_nogff")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prokka_out_for_test-missgff"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missgff original_name.fna: no .gff file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_sevgff():
    """
    Check that check_prokka returns false when there is more than 1 ffn file,
    and an error message
    """
    logger = my_logger("test_check_prokka_sevgff")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prokka_out_for_test-sevgff"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + "2.gff"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-sevgff original_name.fna: several .gff files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_wrong_cont():
    """
    Check that check_prokka returns an error message when the number of contigs in tbl
    file is not as expected
    """
    logger = my_logger("test_check_prokka_wrong_cont")
    outdir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    name = "prokka_out_for_test"
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 10
    assert not afunc.check_prokka(outdir, logf, name, gpath, nbcont, logger[1])
    msg = ("prokka_out_for_test original_name.fna: no matching number of contigs; "
           "nbcontig=10; in tbl =6")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_wrong_tbl_cds():
    """
    Check that check_prokka returns an error message when the number of CDS in tbl
    file is different from the number of headers in faa file
    """
    logger = my_logger("test_check_prokka_wrong_tbl_cds")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "res_checkProkkaWrongTbl")
    os.makedirs(out_dir)
    name = "prokka_out_for_test-wrongCDS"
    tblfile = os.path.join(TEST_DIR, name + ".tbl")
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(out_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))

    shutil.copyfile(tblfile, os.path.join(out_dir, name + ".tbl"))
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg1 = ("prokka_out_for_test-wrongCDS original_name.fna: "
            "no matching number of proteins between tbl and faa; "
            "faa=14; in tbl =12")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg1


def test_check_prokka_ok():
    """
    Check that everything is ok with prokka results (tbl, faa and ffn files exist,
    and number of CDS, CRISPR and genes correspond between them)
    """
    logger = my_logger("test_check_prokka_ok")
    outdir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    name = "prokka_out_for_test"
    logf = os.path.join(GENEPATH, "prokka.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 6
    assert afunc.check_prokka(outdir, logf, name, gpath, nbcont, logger[1])
    q = logger[0]
    assert q.qsize() == 0


def test_run_prokka_out_exists_ok():
    """
    Test that when the output directory already exists, and files inside are OK,
    run_prokka returns True, with a warning message indicating that prokka did not rerun.
    """
    logger = my_logger("test_run_prokka_out_exists_ok")
    utils.init_logger(LOGFILE_BASE, 0, 'prokka_out_exists_ok')
    gpath = "path/to/nogenome/original_name.fna"
    cores_prokka = 1
    name = "prokka_out_for_test"
    force = False
    nbcont = 6
    trn_file = "nofile.trn"
    arguments = (gpath, TEST_DIR, cores_prokka, name, force, nbcont, trn_file, logger[0])
    assert afunc.run_prokka(arguments)

    q = logger[0]
    assert q.qsize() == 4
    # start annotating :
    assert q.get().message.startswith("Start annotating")
    # warning prokka results folder exists:
    assert q.get().message.startswith("Prokka results folder test/data/annotate/"
                                      "test_files/"
                                      "original_name.fna-prokkaRes already exists.")
    # Results in result folder are ok
    assert q.get().message.startswith("Prokka did not run again, formatting step used already "
                                      "generated results of Prokka in "
                                      "test/data/annotate/test_files/original_name.fna-prokkaRes.")
    # End annotation:
    assert q.get().message.startswith("End annotating")


def test_run_prokka_out_exists_error():
    """
    Test that when the output directory already exists, and 1 file is missing,
    run_prokka returns False, and writes the warning message saying that prokka did not
    rerun, + the warning message for the missing file(s).
    """
    logger = my_logger("test_run_prokka_out_exists_error")
    utils.init_logger(LOGFILE_BASE, 0, 'prokka_out_error')
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    new_prok_dir = os.path.join(GENEPATH, "original_name-error-prokkaRes")
    name = "prokka_out_for_test-wrongCDS"
    os.makedirs(new_prok_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".fna"),
                    os.path.join(new_prok_dir, name + ".fna"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(new_prok_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(new_prok_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(new_prok_dir, name + ".gff"))
    gpath = "path/to/nogenome/original_name-error"
    cores_prokka = 1
    force = False
    nbcont = 6
    trn_file = "nofile.trn"
    arguments = (gpath, GENEPATH, cores_prokka, name, force, nbcont, trn_file, logger[0])
    assert not afunc.run_prokka(arguments)
    q = logger[0]
    assert q.qsize() == 4
    # start annotating :
    assert q.get().message.startswith("Start annotating")
    # warning prokka results folder exists:
    assert q.get().message == ("Prokka results folder test/data/annotate/generated_by_unit-tests/"
                               "original_name-error-prokkaRes already exists.")
    # error, no tbl file
    assert q.get().message == "prokka_out_for_test-wrongCDS original_name-error: no .tbl file"
    # warning, files in outdir are not as expected
    assert q.get().message.startswith("Problems in the files contained in your already existing "
                                      "output dir (test/data/annotate/generated_by_unit-tests/"
                                      "original_name-error-prokkaRes)")


def test_run_prokka_out_exists_force():
    """
    Test that when the output directory already exists with wrong files, but force is on,
    prokka is rerun and outputs the right files
    """
    logger = my_logger("test_run_prokka_out_exists_force")
    utils.init_logger(LOGFILE_BASE, 0, 'force')
    gpath = os.path.join(GEN_PATH, "H299_H561.fasta")
    out_prokdir = os.path.join(GENEPATH, "H299_H561.fasta-prokkaRes")
    name = "test_runprokka_H299"
    # Put empty tbl, faa, ffn files in prokka output dir, to check that they are overridden
    os.makedirs(out_prokdir)
    open(os.path.join(out_prokdir, name + ".tbl"), "w").close()
    open(os.path.join(out_prokdir, name + ".faa"), "w").close()
    open(os.path.join(out_prokdir, name + ".ffn"), "w").close()
    cores_prokka = 2
    force = True
    nbcont = 3
    trn_file = "nofile.trn"
    arguments = (gpath, GENEPATH, cores_prokka, name, force, nbcont, trn_file, logger[0])
    assert afunc.run_prokka(arguments)
    # As we used 'force', tbl, faa and ffn files, which were empty, must have been replaced
    # by the prokka output
    exp_dir = os.path.join(EXP_DIR, "H299_H561.fasta-short-contig.fna-prokkaRes",
                           "test_runprokka_H299")
    out_tbl = os.path.join(out_prokdir, name + ".tbl")
    out_faa = os.path.join(out_prokdir, name + ".faa")
    out_ffn = os.path.join(out_prokdir, name + ".ffn")
    assert os.path.isfile(out_tbl)
    # For tbl file, check that, at least, the 3 contigs were considered,
    # and that the number of CDS is as expected.
    # Before, we checked that the output
    # was exactly as expected. But it changes with the different versions of prokka, so
    # we cannot compare the whole file.
    with open(out_tbl, "r") as outt:
        lines = [line.strip() for line in outt.readlines()]
        # Check that there are 3 contigs
        feature = 0
        for line in lines:
            if 'Feature' in line:
                feature += 1
        assert feature == 3
        # Check that there are 16 CDS
        CDS = 0
        for line in lines:
            if "CDS" in line:
                CDS += 1
        assert CDS == 16
    # Check that faa and ffn files are as expected
    assert os.path.isfile(out_faa)
    assert tutil.compare_order_content(exp_dir + ".faa", out_faa)
    assert os.path.isfile(out_ffn)
    assert tutil.compare_order_content(exp_dir + ".ffn", out_ffn)
    q = logger[0]
    assert q.qsize() == 4
    assert q.get() .message.startswith("Start annotating test_runprokka_H299 from test/data/"
                                       "annotate/genomes/H299_H561.fasta with Prokka")
    assert q.get() .message == ("Prokka results folder already exists, but removed because "
                                "--force option used")
    assert q.get().message == ("Prokka command: prokka "
                               "--outdir test/data/annotate/generated_by_unit-tests/"
                               "H299_H561.fasta-prokkaRes --cpus 2 --prefix test_runprokka_H299 "
                               "--centre prokka test/data/annotate/genomes/H299_H561.fasta")
    assert q.get() .message.startswith("End annotating test_runprokka_H299 "
                                       "from test/data/annotate/genomes/H299_H561.fasta")


def test_run_prokka_out_doesnt_exist_ok():
    """
    Test that when the output directory does not exist, it creates it, and runs prokka
    with all expected outfiles
    """
    logger = my_logger("test_run_prokka_out_doesnt_exist")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_prokka_out_doesnt_exist')
    gpath = os.path.join(GEN_PATH, "H299_H561.fasta")
    out_dir = os.path.join(GENEPATH, "H299_H561.fasta-prokkaRes")
    cores_prokka = 2
    name = "test_runprokka_H299"
    force = False
    nbcont = 3
    trn_file = "nofile.trn"
    arguments = (gpath, GENEPATH, cores_prokka, name, force, nbcont, trn_file, logger[0])
    assert afunc.run_prokka(arguments)
    # Check content of tbl, ffn and faa files
    exp_dir = os.path.join(EXP_DIR, "H299_H561.fasta-short-contig.fna-prokkaRes",
                           "test_runprokka_H299")
    out_tbl = os.path.join(out_dir, name + ".tbl")
    out_faa = os.path.join(out_dir, name + ".faa")
    out_ffn = os.path.join(out_dir, name + ".ffn")
    out_gff = os.path.join(out_dir, name + ".gff")
    assert os.path.isfile(out_tbl)
    # For tbl file, check that, at least, the 3 contigs were considered,
    # and that the number of CDS is as expected.
    # Before, we checked that the output
    # was exactly as expected. But it changes with the different versions of prokka, so
    # we cannot compare the whole file.
    with open(out_tbl, "r") as outt:
        lines = [line.strip() for line in outt.readlines()]
        # Check that there are 3 contigs
        feature = 0
        for line in lines:
            if 'Feature' in line:
                feature += 1
        assert feature == 3
        # Check that there are 16 CDS
        CDS = 0
        for line in lines:
            if "CDS" in line:
                CDS += 1
        assert CDS == 16
    # Check that faa and ffn files are as expected
    assert os.path.isfile(out_faa)
    assert tutil.compare_order_content(exp_dir + ".faa", out_faa)
    assert os.path.isfile(out_ffn)
    assert tutil.compare_order_content(exp_dir + ".ffn", out_ffn)
    q = logger[0]
    assert q.qsize() == 3
    assert q.get().message.startswith("Start annotating")
    assert q.get().message == ("Prokka command: prokka "
                               "--outdir test/data/annotate/generated_by_unit-tests/"
                               "H299_H561.fasta-prokkaRes --cpus 2 --prefix test_runprokka_H299 "
                               "--centre prokka test/data/annotate/genomes/H299_H561.fasta")
    assert q.get().message.startswith("End annotating")


def test_run_prokka_out_problem_running():
    """
    Check that when a problem occurs while trying to run prokka, run_prokka returns False,
    and the error message indicating to read in the log why it couldn't run
    """
    logger = my_logger("test_run_prokka_out_problem_running")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_prokka_out_problem_running')
    gpath = os.path.join(GEN_PATH, "H299_H561bis.fasta")
    cores_prokka = 2
    name = "test_runprokka_H299-error"
    force = False
    nbcont = 3
    logf = os.path.join(GENEPATH, "H299_H561.fasta-prokka.log")
    trn_file = "nofile.trn"
    arguments = (gpath, GENEPATH, cores_prokka, name, force, nbcont, trn_file, logger[0])
    assert not afunc.run_prokka(arguments)
    q = logger[0]
    assert q.qsize() == 3
    assert q.get().message.startswith("Start annotating")
    assert q.get().message == ("Prokka command: prokka "
                               "--outdir test/data/annotate/generated_by_unit-tests/"
                               "H299_H561bis.fasta-prokkaRes --cpus 2 "
                               "--prefix test_runprokka_H299-error "
                               "--centre prokka test/data/annotate/genomes/H299_H561bis.fasta")
    assert q.get().message == ("Error while trying to run prokka on test_runprokka_H299-error "
                               "from test/data/annotate/genomes/H299_H561bis.fasta")
