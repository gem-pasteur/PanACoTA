#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for annotate/annotation_functions.py
"""

import pytest
import os
import logging
import shutil

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
    import multiprocessing
    m = multiprocessing.Manager()
    q = m.Queue()
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    return q, logging.getLogger(name)


def test_prodigal_train(caplog):
    """
    Check prodigal training on a genome
    """
    caplog.set_level(logging.DEBUG)
    train_gpath = os.path.join(GEN_PATH, "A_H738.fasta")
    gtrain = afunc.prodigal_train(train_gpath, GENEPATH)
    assert ("prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
            "-t test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn") in caplog.text
    assert ("End annotating A_H738.fasta (from test/data/annotate/genomes/A_H738.fasta)")
    assert gtrain == ("test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn")


def test_prodigal_train_error(caplog):
    """
    Check prodigal training on a genome too small
    """
    caplog.set_level(logging.DEBUG)
    train_gpath = os.path.join(GEN_PATH, "H299_H561.fasta")
    gtrain = afunc.prodigal_train(train_gpath, GENEPATH)
    assert ("prodigal command: prodigal -i test/data/annotate/genomes/H299_H561.fasta "
            "-t test/data/annotate/generated_by_unit-tests/H299_H561.fasta.trn") in caplog.text
    assert ("Error while trying to train prodigal on H299_H561.fasta") in caplog.text
    assert gtrain == ""


def test_prodigal_train_exists(caplog):
    """
    Check prodigal training when the training file already exists
    """
    caplog.set_level(logging.DEBUG)
    train_gpath = os.path.join(GEN_PATH, "H299_H561.fasta")
    # Create (empty) trn file
    gtrain = os.path.join(GENEPATH, "H299_H561.fasta.trn")
    open(gtrain, "w").close()
    gtrain_found = afunc.prodigal_train(train_gpath, GENEPATH)
    assert gtrain_found == gtrain
    assert ("A training file already exists (test/data/annotate/generated_by_unit-tests/"
            "H299_H561.fasta.trn). It will be used to annotate all genomes.") in caplog.text


def test_check_prodigal_nofaa():
    """
    Check that check_prodigal returns false when a faa file is missing, and an error message
    """
    logger = my_logger("test_check_prodigal_nofaa")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join(GENEPATH, "out_test_nofaa")
    name = "prodigal_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prodigal.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "prodigal_out_for_test-missfaa original_name.fna: no or several .faa file(s)"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prodigal_sevfaa():
    """
    Check that check_prodigal returns false when there is more than 1 faa file,
    and an error message
    """
    logger = my_logger("test_check_prodigal_sevfaa")
    ori_prod_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_nofaa")
    name = "prodigal_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prod_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prod_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prod_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + "2.faa"))
    shutil.copyfile(os.path.join(ori_prod_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prodigal.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "prodigal_out_for_test-missfaa original_name.fna: no or several .faa file(s)"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prodigal_noffn():
    """
    Check that check_prodigal returns false when a ffn file is missing, and an error message
    """
    logger = my_logger("test_check_prodigal_noffn")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test-missffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prodigal.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "prodigal_out_for_test-missffn original_name.fna: no or several .ffn file"
    q = logger[0]
    assert q.qsize() == 1
    assert msg in q.get().message


def test_check_prodigal_sevffn():
    """
    Check that check_prodigal returns false when there is more than 1 ffn file,
    and an error message
    """
    logger = my_logger("test_check_prodigal_sevffn")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test-sevffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + "2.ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prodigal.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "prodigal_out_for_test-sevffn original_name.fna: no or several .ffn file"
    q = logger[0]
    assert q.qsize() == 1
    assert msg in q.get().message


def test_check_prodigal_nogff():
    """
    Check that check_prodigal returns false when a ffn file is missing, and an error message
    """
    logger = my_logger("test_check_prodigal_nogff")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test-missgff"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    logf = os.path.join(GENEPATH, "prodigal.log")
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "prodigal_out_for_test-missgff original_name.fna: no or several .gff file"
    q = logger[0]
    assert q.qsize() == 1
    assert msg in q.get().message


def test_check_prodigal_sevgff():
    """
    Check that check_prodigal returns false when there is more than 1 ffn file,
    and an error message
    """
    logger = my_logger("test_check_prodigal_sevgff")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test-sevgff"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + "2.gff"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = os.path.join(GENEPATH, "prodigal.log")
    gpath = "path/to/nogenome/original_name.fna"
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "prodigal_out_for_test-sevgff original_name.fna: no or several .gff file"
    q = logger[0]
    assert q.qsize() == 1
    assert msg in q.get().message


def test_check_prodigal_emptyfaa():
    """
    Check that check_prodigal returns false when there are all expected files, but faa
    file is empty
    """
    logger = my_logger("test_check_prodigal_ok")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test"
    logf = os.path.join(GENEPATH, "prodigal.log")
    os.makedirs(out_dir)
    open(os.path.join(out_dir, name + ".faa"), "w").close()
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    gpath = "path/to/nogenome/original_name.fna"
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "Genome prodigal_out_for_test (from original_name.fna): At least one of your Prodigal result file is empty."
    q = logger[0]
    assert q.qsize() == 1
    assert msg in q.get().message


def test_check_prodigal_emptyffn():
    """
    Check that check_prodigal returns false when there are all expected files, but ffn
    file is empty
    """
    logger = my_logger("test_check_prodigal_ok")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test"
    logf = os.path.join(GENEPATH, "prodigal.log")
    os.makedirs(out_dir)
    open(os.path.join(out_dir, name + ".ffn"), "w").close()
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    gpath = "path/to/nogenome/original_name.fna"
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "Genome prodigal_out_for_test (from original_name.fna): At least one of your Prodigal result file is empty."
    q = logger[0]
    assert q.qsize() == 1
    assert msg in q.get().message


def test_check_prodigal_emptygff():
    """
    Check that check_prodigal returns false when there are all expected files, but gff
    file is empty
    """
    logger = my_logger("test_check_prodigal_ok")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test"
    logf = os.path.join(GENEPATH, "prodigal.log")
    os.makedirs(out_dir)
    open(os.path.join(out_dir, name + ".gff"), "w").close()
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    gpath = "path/to/nogenome/original_name.fna"
    assert not afunc.check_prodigal(gpath, name, out_dir, logger[1])
    msg = "Genome prodigal_out_for_test (from original_name.fna): At least one of your Prodigal result file is empty."
    q = logger[0]
    assert q.qsize() == 1
    assert msg in q.get().message


def test_check_prodigal_ok():
    """
    Check that everything is ok with prodigal results (tbl, faa and ffn files exist,
    and number of CDS, CRISPR and genes correspond between them)
    """
    logger = my_logger("test_check_prodigal_ok")
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    out_dir = os.path.join(GENEPATH, "out_test_noffn")
    name = "prodigal_out_for_test"
    logf = os.path.join(GENEPATH, "prodigal.log")
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    gpath = "path/to/nogenome/original_name.fna"
    assert afunc.check_prodigal(gpath, name, out_dir, logger[1])


def test_run_prodigal_out_exists_ok():
    """
    Test that when the output directory already exists, and files inside are OK,
    run_prodigal returns True, with a warning message indicating that prodigal did not rerun.
    """
    logger = my_logger("test_run_prodigal_out_exists_ok")
    utils.init_logger(LOGFILE_BASE, 0, 'prodigal_out_exists_ok')
    gpath = "path/to/nogenome/original_name.fna"
    cores_prodigal = 1
    name = "prodigal.outtest.ok"
    force = False
    nbcont = 7
    trn_file = os.path.join(TEST_DIR, "A_H738-and-B2_A3_5.fna.trn")
    arguments = (gpath, TEST_DIR, cores_prodigal, name, force, nbcont, trn_file, logger[0])
    assert afunc.run_prodigal(arguments)

    q = logger[0]
    assert q.qsize() == 4
    # start annotating :
    assert q.get().message.startswith("Start annotating prodigal.outtest.ok (from "
                                      "path/to/nogenome/original_name.fna sequence) with Prodigal")
    # # warning prodigal results folder exists:
    assert q.get().message.startswith("Prodigal results folder test/data/annotate/test_files/"
                                      "original_name.fna-prodigalRes already exists.")
    # Results in result folder are ok
    assert q.get().message.startswith("Prodigal did not run again. Formatting step will use "
                                      "already generated results of Prodigal in "
                                      "test/data/annotate/test_files/"
                                      "original_name.fna-prodigalRes.")
    # End annotation:
    assert q.get().message.startswith("End annotating")


def test_run_prodigal_out_exists_error():
    """
    Test that when the output directory already exists, and 1 file is missing,
    run_prodigal returns False, and writes the warning message saying that prodigal did not
    rerun, + the warning message for the missing file(s).
    """
    logger = my_logger("test_run_prodigal_out_exists_error")
    utils.init_logger(LOGFILE_BASE, 0, 'prodigal_out_error')
    ori_prok_dir = os.path.join(TEST_DIR, "original_name.fna-prodigalRes")
    ori_name = "prodigal.outtest.ok"
    new_prok_dir = os.path.join(GENEPATH, "original_name-error-prodigalRes")
    name = "prodigal_out_for_test-wrongCDS"
    os.makedirs(new_prok_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(new_prok_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(new_prok_dir, name + ".faa"))
    open(os.path.join(new_prok_dir, name + ".gff"), "w").close()
    gpath = "path/to/nogenome/original_name-error"
    cores_prodigal = 1
    force = False
    trn_file = os.path.join(TEST_DIR, "A_H738-and-B2_A3_5.fna.trn")
    nbcont = 7
    arguments = (gpath, GENEPATH, cores_prodigal, name, force, nbcont, trn_file, logger[0])
    assert not afunc.run_prodigal(arguments)
    q = logger[0]
    assert q.qsize() == 4
    # start annotating :
    assert q.get().message.startswith("Start annotating")
    # warning prodigal results folder exists:
    assert q.get().message == ("Prodigal results folder test/data/annotate/"
                               "generated_by_unit-tests/"
                               "original_name-error-prodigalRes already exists.")
    # error, empty gff
    msg = ("Genome prodigal_out_for_test-wrongCDS (from original_name-error): "
           "At least one of your Prodigal result file is empty.")
    assert q.get().message == msg
    # warning, files in outdir are not as expected
    assert q.get().message.startswith("Problems in the files contained in your already existing "
                                      "output dir (test/data/annotate/generated_by_unit-tests/"
                                      "original_name-error-prodigalRes")


def test_run_prodigal_out_exists_force():
    """
    Test that when the output directory already exists with wrong files, but force is on,
    prodigal is rerun and outputs the right files
    """
    logger = my_logger("test_run_prodigal_out_exists_force")
    utils.init_logger(LOGFILE_BASE, 0, 'force')
    gpath = os.path.join(GEN_PATH, "H299_H561.fasta")
    out_prokdir = os.path.join(GENEPATH, "H299_H561.fasta-prodigalRes")
    name = "test_runprodigal_H299"
    # Put empty tbl, faa, ffn files in prodigal output dir, to check that they are overridden
    os.makedirs(out_prokdir)
    open(os.path.join(out_prokdir, name + ".gff"), "w").close()
    open(os.path.join(out_prokdir, name + ".faa"), "w").close()
    open(os.path.join(out_prokdir, name + ".ffn"), "w").close()
    cores_prodigal = 2
    force = True
    nbcont = 3
    trn_file = os.path.join(TEST_DIR, "A_H738-and-B2_A3_5.fna.trn")
    arguments = (gpath, GENEPATH, cores_prodigal, name, force, nbcont, trn_file, logger[0])
    assert afunc.run_prodigal(arguments)
    # As we used 'force', tbl, faa and ffn files, which were empty, must have been replaced
    # by the prodigal output
    exp_dir = os.path.join(EXP_DIR, "H299_H561.fasta-prodigalRes",
                           "ESCO.1015.00001")
    out_gff = os.path.join(out_prokdir, name + ".gff")
    out_faa = os.path.join(out_prokdir, name + ".faa")
    out_ffn = os.path.join(out_prokdir, name + ".ffn")
    # Check that faa and ffn files are as expected
    assert os.path.isfile(out_faa)
    assert tutil.compare_order_content(exp_dir + ".faa", out_faa)
    assert os.path.isfile(out_ffn)
    assert tutil.compare_order_content(exp_dir + ".ffn", out_ffn)
    q = logger[0]
    assert q.qsize() == 4
    assert q.get() .message.startswith("Prodigal results folder already exists, but is "
                                       "removed because --force option was used")
    assert q.get() .message.startswith("Start annotating test_runprodigal_H299 (from test/data/"
                                       "annotate/genomes/H299_H561.fasta sequence) "
                                       "with Prodigal")
    assert q.get().message.startswith("Prodigal command: prodigal -i test/data/annotate/genomes/"
                                      "H299_H561.fasta -d test/data/annotate/"
                                      "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                                      "test_runprodigal_H299.ffn -a test/data/annotate/"
                                      "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                                      "test_runprodigal_H299.faa -f gff -o test/data/annotate/"
                                      "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                                      "test_runprodigal_H299.gff -t "
                                      "test/data/annotate/test_files/A_H738-and-B2_A3_5.fna.trn "
                                      "-q")
    assert q.get() .message.startswith("End annotating test_runprodigal_H299 "
                                       "(from test/data/annotate/genomes/H299_H561.fasta)")


def test_run_prodigal_out_doesnt_exist():
    """
    Test that when the output directory does not exist, it creates it, and runs prodigal
    with all expected outfiles
    """
    logger = my_logger("test_run_prodigal_out_doesnt_exist")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_prodigal_out_doesnt_exist')
    gpath = os.path.join(GEN_PATH, "H299_H561.fasta")
    out_dir = os.path.join(GENEPATH, "H299_H561.fasta-prodigalRes")
    cores_prodigal = 2
    name = "test_runprodigal_H299"
    force = False
    trn_file = os.path.join(TEST_DIR, "A_H738-and-B2_A3_5.fna.trn")
    nbcont = 3
    arguments = (gpath, GENEPATH, cores_prodigal, name, force, nbcont, trn_file, logger[0])
    assert afunc.run_prodigal(arguments)
    # Check content of tbl, ffn and faa files
    exp_dir = os.path.join(EXP_DIR, "H299_H561.fasta-prodigalRes",
                           "ESCO.1015.00001")
    out_faa = os.path.join(out_dir, name + ".faa")
    out_ffn = os.path.join(out_dir, name + ".ffn")
    out_gff = os.path.join(out_dir, name + ".gff")
    # Check that faa and ffn files are as expected
    assert os.path.isfile(out_faa)
    assert tutil.compare_order_content(exp_dir + ".faa", out_faa)
    assert os.path.isfile(out_ffn)
    assert tutil.compare_order_content(exp_dir + ".ffn", out_ffn)
    assert os.path.isfile(out_ffn)
    assert tutil.compare_order_content(exp_dir + ".gff", out_gff)
    q = logger[0]
    assert q.qsize() == 3
    assert q.get().message.startswith("Start annotating")
    assert q.get().message == ("Prodigal command: prodigal -i test/data/annotate/genomes/"
                               "H299_H561.fasta -d test/data/annotate/"
                               "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                               "test_runprodigal_H299.ffn -a test/data/annotate/"
                               "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                               "test_runprodigal_H299.faa -f gff -o test/data/annotate/"
                               "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                               "test_runprodigal_H299.gff -t "
                               "test/data/annotate/test_files/A_H738-and-B2_A3_5.fna.trn "
                               "-q")
    assert q.get().message.startswith("End annotating")


def test_run_prodigal_small():
    """
    Test that when the output directory does not exist, it creates it, and runs prodigal
    with all expected outfiles. Here, we run prodigal with --small option (on a small genome)
    """
    logger = my_logger("test_run_prodigal_small")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_prodigal_small')
    gpath = os.path.join(GEN_PATH, "H299_H561.fasta")
    out_dir = os.path.join(GENEPATH, "H299_H561.fasta-prodigalRes")
    cores_prodigal = 2
    name = "test_runprodigal_small_H299"
    force = False
    trn_file = "small option"
    nbcont = 3
    arguments = (gpath, GENEPATH, cores_prodigal, name, force, nbcont, trn_file, logger[0])
    assert afunc.run_prodigal(arguments)

    # Check content of tbl, ffn and faa files
    exp_dir = os.path.join(EXP_DIR, "H299_H561.fasta_small-prodigalRes",
                           "test_runprodigal_small_H299")
    out_faa = os.path.join(out_dir, name + ".faa")
    out_ffn = os.path.join(out_dir, name + ".ffn")
    out_gff = os.path.join(out_dir, name + ".gff")
    # Check that faa and ffn files are as expected
    assert os.path.isfile(out_faa)
    assert tutil.compare_order_content(exp_dir + ".faa", out_faa)
    assert os.path.isfile(out_ffn)
    assert tutil.compare_order_content(exp_dir + ".ffn", out_ffn)
    assert os.path.isfile(out_ffn)
    assert tutil.compare_order_content(exp_dir + ".gff", out_gff)
    # Check logs
    q = logger[0]
    assert q.qsize() == 3
    assert q.get().message.startswith("Start annotating")
    prodigal_cmd = q.get().message
    assert ("Prodigal command: prodigal -i test/data/annotate/genomes/"
            "H299_H561.fasta -d test/data/annotate/"
            "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
            "test_runprodigal_small_H299.ffn -a test/data/annotate/"
            "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
            "test_runprodigal_small_H299.faa -f gff -o test/data/annotate/"
            "generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
            "test_runprodigal_small_H299.gff -p meta -q")in prodigal_cmd
    assert q.get().message.startswith("End annotating")


def test_run_prodigal_out_problem_running():
    """
    Check that when a problem occurs while trying to run prodigal, run_prodigal returns False,
    and the error message indicating to read in the log why it couldn't run
    """
    logger = my_logger("test_run_prodigal_out_problem_running")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_prodigal_out_problem_running')
    gpath = os.path.join(GEN_PATH, "H299_H561bis.fasta")
    cores_prodigal = 2
    name = "test_runprodigal_H299-error"
    force = False
    nbcont = 3
    trn_file = os.path.join(TEST_DIR, "A_H738-and-B2_A3_5.fna.trn")
    logf = os.path.join(GENEPATH, "H299_H561bis.fasta-prodigal.log")
    arguments = (gpath, GENEPATH, cores_prodigal, name, force, nbcont, trn_file, logger[0])
    assert not afunc.run_prodigal(arguments)
    # Check that output directory is empty
    outdir = os.path.join(GENEPATH, "H299_H561bis.fasta-prodigalRes")
    assert os.listdir(outdir) == []
    # Check logs
    q = logger[0]
    assert q.qsize() == 3
    assert q.get().message.startswith("Start annotating")
    assert q.get().message.startswith("Prodigal command: prodigal -i test/data/annotate/genomes/"
                               "H299_H561bis.fasta -d test/data/annotate/"
                               "generated_by_unit-tests/H299_H561bis.fasta-prodigalRes/"
                               "test_runprodigal_H299-error.ffn -a test/data/annotate/"
                               "generated_by_unit-tests/H299_H561bis.fasta-prodigalRes/"
                               "test_runprodigal_H299-error.faa -f gff -o test/data/annotate/"
                               "generated_by_unit-tests/H299_H561bis.fasta-prodigalRes/"
                               "test_runprodigal_H299-error.gff -t "
                               "test/data/annotate/test_files/A_H738-and-B2_A3_5.fna.trn "
                               "-q")
    assert q.get().message.startswith("Error while trying to run prodigal. See test/data/"
                                      "annotate/generated_by_unit-tests/"
                                      "H299_H561bis.fasta-prodigal.log.err.")


def test_run_prodigal_noout_notrain():
    """
    Prodigal result directory does not exist (not already run)
    training file does not exist (probably, problem while trying to train)
    -> return  False
    """
    logger = my_logger("test_run_prodigal_out_exists_error")
    utils.init_logger(LOGFILE_BASE, 0, 'prodigal_out_error')
    gpath = "path/to/nogenome/original_name-error"
    cores_prodigal = 1
    name = "prodigal_out_for_test-wrongCDS"
    force = False
    nbcont = 7
    trn_file = "ghost_trn_file"
    arguments = (gpath, GENEPATH, cores_prodigal, name, force, nbcont, trn_file, logger[0])
    assert not afunc.run_prodigal(arguments)
    q = logger[0]
    assert q.qsize() == 0
