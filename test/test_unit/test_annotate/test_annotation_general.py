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
import PanACoTA.annotate_module.annotation_functions as afunc


# Define variables used by several tests
DBDIR = os.path.join("test", "data", "annotate")
GEN_PATH = os.path.join(DBDIR, "genomes")
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


def test_count_headers():
    """
    Count how many sequences there are in the given multi-fasta file
    """
    seqfile = os.path.join(GEN_PATH, "genome4.fasta")
    nb = afunc.count_headers(seqfile)
    assert nb == 5


def test_count_tbl():
    """
    Count the different features found in the tbl file, and return
    nbcont, nbCDS, nbGene, nbCRISPR
    """
    tblfile = os.path.join(TEST_DIR, "original_name.fna-prokkaRes", "prokka_out_for_test.tbl")
    ncont, ncds, ngene = afunc.count_tbl(tblfile)
    assert ncont == 6
    assert ncds == 14
    assert ngene == 16


def test_run_all_prodigal():
    """
    Check that there is no problem when running prodigal on all genomes
    Start and end are not necessarily in the same order (ex: start1, start2, end2, end1)
    """
    logger = my_logger("test_run_all_parallel_more_threads")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_parallel_more_threads')
    # genomes = {genome: [name, gpath, annot_path, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join(GEN_PATH, genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join(GEN_PATH, genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    trn_gname = genome2
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_gname,
                                     prodigal_only=True, quiet=True)
    assert final[genome1]
    assert final[genome2]
    q = logger[0]
    assert q.qsize() == 10
    assert q.get().message == "Annotating all genomes with prodigal"
    assert q.get().message == "Prodigal will train using test/data/annotate/genomes/A_H738.fasta"
    assert q.get().message == ("prodigal command: prodigal -i "
                               "test/data/annotate/genomes/A_H738.fasta -t "
                               "test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn")
    assert q.get().message == "End training on test/data/annotate/genomes/A_H738.fasta"
    messages = []
    for i in range(6):
        a = q.get().message
        messages.append(a)
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "(from test/data/annotate/genomes/H299_H561.fasta sequence) "
                            "with Prodigal")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "(from test/data/annotate/genomes/A_H738.fasta sequence) "
                            "with Prodigal")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    assert message_start_annot1 in messages
    assert message_start_annot2 in messages
    # Prodigal cmd
    message_cmd1 = ("Prodigal command: prodigal -i test/data/annotate/genomes/H299_H561.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                    "test_runall_1by1_1.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "H299_H561.fasta-prodigalRes/test_runall_1by1_1.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/"
                    "H299_H561.fasta-prodigalRes/test_runall_1by1_1.gff -t "
                    "test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn -q")
    message_cmd2 = ("Prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "A_H738.fasta-prodigalRes/test_runall_1by1_2.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.gff -t "
                    "test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn -q")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    message_end_annot1 = ("End annotating test_runall_1by1_1 (from test/data/annotate/genomes/"
                          "H299_H561.fasta)")
    message_end_annot2 = ("End annotating test_runall_1by1_2 (from test/data/annotate/genomes/"
                          "A_H738.fasta)")
    assert message_end_annot1 in messages
    assert message_end_annot2 in messages


def test_run_all_prodigal_small():
    """
    Check that there is no problem when running prodigal with --small option on all genomes
    Start and end are not necessarily in the same order (ex: start1, start2, end2, end1)
    """
    logger = my_logger("test_run_all_prodigal_small")
    # genomes = {genome: [name, gpath, annot_path, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join(GEN_PATH, genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join(GEN_PATH, genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    trn_gname = genome2
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_gname,
                                     prodigal_only=True, quiet=True, small=True)
    assert final[genome1]
    assert final[genome2]
    q = logger[0]
    # assert q.qsize() == 10
    assert q.get().message == "Annotating all genomes with prodigal"
    # assert q.get().message == "Prodigal will train using test/data/annotate/genomes/A_H738.fasta"
    # assert q.get().message == ("prodigal command: prodigal -i "
    #                            "test/data/annotate/genomes/A_H738.fasta -t "
    #                            "test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn")
    # assert q.get().message == "End training on test/data/annotate/genomes/A_H738.fasta"
    messages = []
    for i in range(6):
        a = q.get().message
        messages.append(a)
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "(from test/data/annotate/genomes/H299_H561.fasta sequence) "
                            "with Prodigal")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "(from test/data/annotate/genomes/A_H738.fasta sequence) "
                            "with Prodigal")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    assert message_start_annot1 in messages
    assert message_start_annot2 in messages
    # Prodigal cmd
    message_cmd1 = ("Prodigal command: prodigal -i test/data/annotate/genomes/H299_H561.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/H299_H561.fasta-prodigalRes/"
                    "test_runall_1by1_1.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "H299_H561.fasta-prodigalRes/test_runall_1by1_1.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/"
                    "H299_H561.fasta-prodigalRes/test_runall_1by1_1.gff -p meta -q")
    message_cmd2 = ("Prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "A_H738.fasta-prodigalRes/test_runall_1by1_2.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.gff -p meta -q")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    message_end_annot1 = ("End annotating test_runall_1by1_1 (from test/data/annotate/genomes/"
                          "H299_H561.fasta)")
    message_end_annot2 = ("End annotating test_runall_1by1_2 (from test/data/annotate/genomes/"
                          "A_H738.fasta)")
    assert message_end_annot1 in messages
    assert message_end_annot2 in messages


def test_run_all_prodigal_error_train():
    """
    Check that when we want to train on a genome but it fails, it returns False for all genomes
    Here, it fails because genome to train on is too small
    """
    logger = my_logger("test_run_all_parallel_more_threads")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_parallel_more_threads')
    # genomes = {genome: [name, gpath, annot_path, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join(GEN_PATH, genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join(GEN_PATH, genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    trn_gname = genome1
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_gname,
                                     prodigal_only=True, quiet=True)
    assert not final[genome1]
    assert not final[genome2]
    q = logger[0]
    assert q.qsize() == 4
    assert q.get().message == "Annotating all genomes with prodigal"
    assert q.get().message == ("Prodigal will train using "
                               "test/data/annotate/genomes/H299_H561.fasta")
    assert q.get().message == ("prodigal command: prodigal -i "
                               "test/data/annotate/genomes/H299_H561.fasta -t "
                               "test/data/annotate/generated_by_unit-tests/H299_H561.fasta.trn")
    assert q.get().message == ("Error while trying to train prodigal on H299_H561.fasta. See "
                               "test/data/annotate/generated_by_unit-tests/"
                               "H299_H561.fasta.trn-prodigal-train.log.err.")


def test_run_all_prodigal_train_exists_error():
    """
    Check that when we give a wrong training file, it does not train again, but
    fails while annotating
    """
    logger = my_logger("test_run_prodigal_train_exist_error")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_prodigal_train_exist_error')
    # genomes = {genome: [name, gpath, annot_path, size, nbcont, l90]}
    genome1 = "toto.fasta"
    gpath1 = os.path.join(GEN_PATH, genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join(GEN_PATH, genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    trn_gname = genome1
    # Create empty trn file
    trn_file = os.path.join(GENEPATH, "toto.fasta.trn")
    open(trn_file, "w").close()
    # Run annotation all
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_gname,
                                     prodigal_only=True, quiet=False)
    assert not final[genome1]
    assert not final[genome2]
    q = logger[0]
    assert q.qsize() == 9
    assert q.get().message == "Annotating all genomes with prodigal"
    assert q.get().message == ("Prodigal will train using "
                               "test/data/annotate/genomes/toto.fasta")
    assert q.get().message == ("A training file already exists (test/data/annotate/"
                               "generated_by_unit-tests/toto.fasta.trn). It will be used "
                               "to annotate all genomes.")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    messages = []
    for i in range(6):
        a = q.get().message
        messages.append(a)
    # Check start annotation messages
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "(from test/data/annotate/genomes/toto.fasta sequence) "
                            "with Prodigal")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "(from test/data/annotate/genomes/A_H738.fasta sequence) "
                            "with Prodigal")
    assert message_start_annot1 in messages
    assert message_start_annot2 in messages
    # Prodigal cmd
    message_cmd1 = ("Prodigal command: prodigal -i test/data/annotate/genomes/toto.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/toto.fasta-prodigalRes/"
                    "test_runall_1by1_1.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "toto.fasta-prodigalRes/test_runall_1by1_1.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/"
                    "toto.fasta-prodigalRes/test_runall_1by1_1.gff -t "
                    "test/data/annotate/generated_by_unit-tests/toto.fasta.trn -q")
    message_cmd2 = ("Prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "A_H738.fasta-prodigalRes/test_runall_1by1_2.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.gff -t "
                    "test/data/annotate/generated_by_unit-tests/toto.fasta.trn -q")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    message_end_annot1 = ("Error while trying to run prodigal. See "
                          "test/data/annotate/generated_by_unit-tests/"
                          "toto.fasta-prodigal.log.err.")
    message_end_annot2 = ("Error while trying to run prodigal. See "
                          "test/data/annotate/generated_by_unit-tests/"
                          "A_H738.fasta-prodigal.log.err.")
    assert message_end_annot1 in messages
    assert message_end_annot2 in messages


def test_run_all_prodigal_train_exists_ok():
    """
    Check that when we want to train on a genome but it fails, it returns False for all genomes
    Here, it fails because genome to train on is too small
    """
    logger = my_logger("test_run_prodigal_train_exist_error")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_prodigal_train_exist_error')
    # genomes = {genome: [name, gpath, annot_path, size, nbcont, l90]}
    genome1 = "toto.fasta"
    gpath1 = os.path.join(GEN_PATH, genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join(GEN_PATH, genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    trn_gname = genome1
    # Copy trn file to outdir, so that panacota detects that it already exists
    orig_trn_file = os.path.join(TEST_DIR, "A_H738-and-B2_A3_5.fna.trn")
    trn_file = os.path.join(GENEPATH, "toto.fasta.trn")
    shutil.copyfile(orig_trn_file, trn_file)
    # Run annotation all
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_gname,
                                     prodigal_only=True, quiet=False)
    assert not final[genome1]
    assert final[genome2]
    q = logger[0]
    assert q.qsize() == 9
    assert q.get().message == "Annotating all genomes with prodigal"
    assert q.get().message == ("Prodigal will train using "
                               "test/data/annotate/genomes/toto.fasta")
    assert q.get().message == ("A training file already exists (test/data/annotate/"
                               "generated_by_unit-tests/toto.fasta.trn). It will be used "
                               "to annotate all genomes.")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    messages = []
    for i in range(6):
        a = q.get().message
        messages.append(a)
    # Check start annotation messages
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "(from test/data/annotate/genomes/toto.fasta sequence) "
                            "with Prodigal")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "(from test/data/annotate/genomes/A_H738.fasta sequence) "
                            "with Prodigal")
    assert message_start_annot1 in messages
    assert message_start_annot2 in messages
    # Prodigal cmd
    message_cmd1 = ("Prodigal command: prodigal -i test/data/annotate/genomes/toto.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/toto.fasta-prodigalRes/"
                    "test_runall_1by1_1.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "toto.fasta-prodigalRes/test_runall_1by1_1.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/"
                    "toto.fasta-prodigalRes/test_runall_1by1_1.gff -t "
                    "test/data/annotate/generated_by_unit-tests/toto.fasta.trn -q")
    message_cmd2 = ("Prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "A_H738.fasta-prodigalRes/test_runall_1by1_2.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.gff -t "
                    "test/data/annotate/generated_by_unit-tests/toto.fasta.trn -q")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    message_end_annot1 = ("Error while trying to run prodigal. See "
                          "test/data/annotate/generated_by_unit-tests/"
                          "toto.fasta-prodigal.log.err.")
    message_end_annot2 = ("End annotating test_runall_1by1_2 (from test/data/annotate/genomes/"
                          "A_H738.fasta)")
    assert message_end_annot1 in messages
    assert message_end_annot2 in messages


def test_run_all_prodigal_error_annotate():
    """
    Check running prodigal on 2 genomes:
    - prodigal train ok
    - running on genome1 ok
    - running on genome2 error
    """
    logger = my_logger("test_run_all_parallel_more_threads")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_parallel_more_threads')
    # genomes = {genome: [name, gpath, annot_path, size, nbcont, l90]}
    genome1 = "toto.fasta"
    gpath1 = os.path.join(GEN_PATH, genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join(GEN_PATH, genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    trn_gname = genome2
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_gname,
                                     prodigal_only=True, quiet=False)
    assert not final[genome1]
    assert final[genome2]
    q = logger[0]
    assert q.qsize() == 10
    assert q.get().message == "Annotating all genomes with prodigal"
    assert q.get().message == "Prodigal will train using test/data/annotate/genomes/A_H738.fasta"
    assert q.get().message == ("prodigal command: prodigal -i "
                               "test/data/annotate/genomes/A_H738.fasta -t "
                               "test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn")
    assert q.get().message == "End training on test/data/annotate/genomes/A_H738.fasta"
    messages = []
    for i in range(6):
        a = q.get().message
        messages.append(a)
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "(from test/data/annotate/genomes/toto.fasta sequence) "
                            "with Prodigal")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "(from test/data/annotate/genomes/A_H738.fasta sequence) "
                            "with Prodigal")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    assert message_start_annot1 in messages
    assert message_start_annot2 in messages
    # Prodigal cmd
    message_cmd1 = ("Prodigal command: prodigal -i test/data/annotate/genomes/toto.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/toto.fasta-prodigalRes/"
                    "test_runall_1by1_1.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "toto.fasta-prodigalRes/test_runall_1by1_1.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/"
                    "toto.fasta-prodigalRes/test_runall_1by1_1.gff -t "
                    "test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn -q")
    message_cmd2 = ("Prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "A_H738.fasta-prodigalRes/test_runall_1by1_2.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.gff -t "
                    "test/data/annotate/generated_by_unit-tests/A_H738.fasta.trn -q")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    message_end_annot1 = ("Error while trying to run prodigal. See "
                          "test/data/annotate/generated_by_unit-tests/"
                          "toto.fasta-prodigal.log.err.")
    message_end_annot2 = ("End annotating test_runall_1by1_2 (from test/data/annotate/genomes/"
                          "A_H738.fasta)")
    assert message_end_annot1 in messages
    assert message_end_annot2 in messages


def test_run_all_prodigal_outexists_error():
    """
    trn file already exists, and output folder too. No force option. Output folder is empty
    -> error message while checking prodigal
    """
    logger = my_logger("test_run_all_parallel_more_threads")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_parallel_more_threads')
    # genomes = {genome: [name, gpath, annot_path, size, nbcont, l90]}
    genome1 = "toto.fasta"
    genome2 = "A_H738.fasta"
    genomes = {genome1: ["test_runall_1by1_1", genome1, genome1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", genome2, genome2, 456464645, 1, 465]}
    # Create prodigal result directories
    prodigaldir_g1 = os.path.join(GENEPATH, "A_H738.fasta-prodigalRes")
    prodigaldir_g2 = os.path.join(GENEPATH, "toto.fasta-prodigalRes")
    os.makedirs(prodigaldir_g1)
    os.makedirs(prodigaldir_g2)
    # Other parameters
    threads = 1
    force = False
    # Add existing training file
    orig_trn_file = os.path.join(TEST_DIR, "A_H738-and-B2_A3_5.fna.trn")
    trn_file = os.path.join(GENEPATH, "toto.fasta.trn")
    shutil.copyfile(orig_trn_file, trn_file)
    trn_gname = genome1
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_gname,
                                     prodigal_only=True, quiet=False)
    assert not final[genome1]
    assert not final[genome2]
    q = logger[0]
    assert q.qsize() == 15
    assert q.get().message == "Annotating all genomes with prodigal"
    assert q.get().message == "Prodigal will train using toto.fasta"
    assert q.get().message == ("A training file already exists (test/data/annotate/"
                               "generated_by_unit-tests/toto.fasta.trn). It will "
                               "be used to annotate all genomes.")
    messages = []
    for i in range(12):
        a = q.get().message
        messages.append(a)
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "(from toto.fasta sequence) with Prodigal")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    assert message_start_annot1 in messages
    # Prodigal cmd
    message_exists1 = ("Prodigal results folder test/data/annotate/generated_by_unit-tests/"
                    "toto.fasta-prodigalRes already exists.")
    message_errorfaa = ("test_runall_1by1_1 toto.fasta: no or several .faa file(s)")
    message_errorffn = ("test_runall_1by1_1 toto.fasta: no or several .ffn file(s)")
    message_errorgff = ("test_runall_1by1_1 toto.fasta: no or several .gff file(s)")
    message_error1 = ("Problems in the files contained in your already existing output dir "
                    "(test/data/annotate/generated_by_unit-tests/toto.fasta-prodigalRes). "
                    "Please check it, or remove it to re-annotate.")
    assert message_exists1 in messages
    assert message_errorfaa in messages
    assert message_errorffn in messages
    assert message_errorgff in messages
    assert message_error1 in messages
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "(from A_H738.fasta sequence) with Prodigal")
    assert message_start_annot2 in messages
    message_error_annot2 = ("Problems in the files contained in your already existing output dir "
                          "(test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes). "
                          "Please check it, or remove it to re-annotate.")
    assert message_error_annot2 in messages


def test_run_all_1by1_prokka():
    """
    Check that when running with 3 threads (not parallel), prokka runs as expected,
    and returns True for each genome
    -> Runs 1 by 1, with prokka using 3 cpus
    Start and end must be ordered: (start1, end1, start2, end2) or (start2, end2, start1, end1)
    """
    logger = my_logger("test_runall_1by1_1")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_1by1')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join(GEN_PATH, genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join(GEN_PATH, genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, gpath2, 456464645, 1, 465]}
    threads = 3
    force = False
    trn_file = "nofile.trn"
    annot_folder = os.path.join(GENEPATH, "annot-folder")
    os.makedirs(annot_folder)
    final = afunc.run_annotation_all(genomes, threads, force, annot_folder, trn_file)
    assert final[genome1]
    assert final[genome2]
    q = logger[0]
    assert q.qsize() == 7
    assert q.get().message == 'Annotating all genomes with prokka'
    # Messages for start and end annotation of the different genomes
    message_start_annot1 = ("Start annotating test_runall_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_cmd1 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "annot-folder/H299_H561.fasta-prokkaRes --cpus 3")
    message_end_annot1 = ("End annotating test_runall_1by1_1 from test/data/annotate/genomes/"
                            "H299_H561.fasta.")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    message_cmd2 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "annot-folder/A_H738.fasta-prokkaRes --cpus 3")
    message_end_annot2 = ("End annotating test_runall_1by1_2 from test/data/annotate/genomes/"
                            "A_H738.fasta.")
    qget = q.get().message
    # Check logs. Given that it is executed in parallel, we cannot know in which order messages
    # will appear
    assert qget == message_start_annot1 or message_start_annot2
    if qget == message_start_annot1:
        # Ending annotation of first genome (same genome as started because running 1by1)
        assert q.get().message.startswith(message_cmd1)
        assert q.get().message == message_end_annot1
    else:
        assert q.get().message.startswith(message_cmd2)
        assert q.get().message == message_end_annot2
    qget2 = q.get().message
    assert qget2 == message_start_annot1 or message_start_annot2
    if qget2 == message_start_annot2:
        # Ending annotation of first genome (same genome as started because running 1by1)
        assert q.get().message.startswith(message_cmd2)
        assert q.get().message == message_end_annot2
    else:
        assert q.get().message.startswith(message_cmd1)
        assert q.get().message == message_end_annot1


def test_run_all_prokka_parallel_less_threads():
    """
    Check that there is no problem when running with less threads than genomes (each genomes
    uses 2 threads)
    Genomes H299 and A_H738 should run well, but genomes genome* have problems (no CDS found),
    so check_prokka should return false.
    """
    logger = my_logger("test_run_all_parallel_more_threads")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_4threads')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    gnames = ["H299_H561.fasta", "A_H738.fasta", "genome1.fasta", "genome2.fasta", "genome3.fasta"]
    gpaths = [os.path.join(GEN_PATH, name) for name in gnames]
    genomes = {gnames[0]: ["test_runall_1by1_1", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["test_runall_1by1_2", gpaths[1], gpaths[1], 456464645, 1, 1],
               gnames[2]: ["test_runall_1by1_3", gpaths[2], gpaths[2], 456464645, 4, 1],
               gnames[3]: ["test_runall_1by1_4", gpaths[3], gpaths[3], 456464645, 3, 1],
               gnames[4]: ["test_runall_1by1_5", gpaths[4], gpaths[4], 456464645, 1, 1]
               }
    threads = 4
    force = False
    trn_file = "nofile.trn"
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_file)
    assert final[gnames[0]]
    assert final[gnames[1]]
    assert not final[gnames[2]]
    assert not final[gnames[3]]
    assert not final[gnames[4]]
    q = logger[0]
    # Check size of logs
    # -> starting log -> 1 log
    # -> for each genome ok (2 first ones): start annotate, prokka cmd, end annotate -> 6 logs
    # -> for each genome not ok (3 others):
    #           start annotate, prokka cmd, problem, end annotate -> 12 logs
    assert q.qsize() == 19
    assert q.get().message == "Annotating all genomes with prokka"
    # messages start annotation
    messages = []
    for i in range(18):
        a = q.get().message
        messages.append(a)
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "from test/data/annotate/genomes/H299_H561.fasta "
                            "with Prokka")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "from test/data/annotate/genomes/A_H738.fasta "
                            "with Prokka")
    message_start_annot3 = ("Start annotating test_runall_1by1_4 "
                            "from test/data/annotate/genomes/genome2.fasta "
                            "with Prokka")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    assert message_start_annot1 in messages
    assert message_start_annot2 in messages
    assert message_start_annot3 in messages
    # messages Prokka cmd
    message_cmd1 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "H299_H561.fasta-prokkaRes --cpus 2 --prefix test_runall_1by1_1 "
                    "--centre prokka test/data/annotate/genomes/H299_H561.fasta")
    message_cmd2 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "A_H738.fasta-prokkaRes --cpus 2 --prefix test_runall_1by1_2 "
                    "--centre prokka test/data/annotate/genomes/A_H738.fasta")
    message_cmd3 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "genome1.fasta-prokkaRes --cpus 2 --prefix test_runall_1by1_3 "
                    "--centre prokka test/data/annotate/genomes/genome1.fasta")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    assert message_cmd3 in messages
    # Messages end annotation cmd
    message_end1 = ("End annotating test_runall_1by1_1 from "
                    "test/data/annotate/genomes/H299_H561.fasta.")
    message_end2 = ("End annotating test_runall_1by1_3 from "
                    "test/data/annotate/genomes/genome1.fasta.")
    message_end3 = ("End annotating test_runall_1by1_5 from "
                    "test/data/annotate/genomes/genome3.fasta.")
    assert message_end1 in messages
    assert message_end2 in messages
    assert message_end3 in messages
    # Messages error annotation cmd
    message_err1 = "test_runall_1by1_3 genome1.fasta: several .faa files"
    message_err2 = "test_runall_1by1_4 genome2.fasta: several .faa files"
    message_err3 = "test_runall_1by1_5 genome3.fasta: several .faa files"
    assert message_err1 in messages
    assert message_err2 in messages
    assert message_err3 in messages


def test_run_all_parallel_prokka_more_threads():
    """
    Check that there is no problem when running with more threads than genomes
    (6 threads and 2 genome: each genome uses 3 threads)
    Genomes H299 should run well but genome1.fasta should get an error
    """
    logger = my_logger("test_run_all_parallel_more_threads")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_4threads')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    gnames = ["H299_H561.fasta", "genome1.fasta"]
    gpaths = [os.path.join(GEN_PATH, name) for name in gnames]
    genomes = {gnames[0]: ["test_runall_1by1_1", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["test_runall_1by1_2", gpaths[1], gpaths[1], 456464645, 4, 1],
               }
    threads = 6
    force = False
    trn_file = "nofile.trn"
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH, trn_file)
    assert final[gnames[0]]
    assert not final[gnames[1]]
    q = logger[0]
    # Check size of logs
    # -> starting log -> 1 log
    # -> for genome ok : start annotate, prokka cmd, end annotate -> 3 logs
    # -> for genome not ok : start annotate, prokka cmd, problem, end annotate -> 4 logs
    assert q.qsize() == 8
    assert q.get().message == "Annotating all genomes with prokka"
    # messages start annotation
    messages = []
    for i in range(7):
        a = q.get().message
        messages.append(a)
    message_start_annot1 = ("Start annotating test_runall_1by1_1 "
                            "from test/data/annotate/genomes/H299_H561.fasta "
                            "with Prokka")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 "
                            "from test/data/annotate/genomes/genome1.fasta "
                            "with Prokka")
    # Check that all messages exist. We cannot know in which order,
    # as 'genomes' is a dict, hence unordered, and as computation is done in parallel
    assert message_start_annot1 in messages
    assert message_start_annot2 in messages
    # messages Prokka cmd
    message_cmd1 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "H299_H561.fasta-prokkaRes --cpus 3 --prefix test_runall_1by1_1 "
                    "--centre prokka test/data/annotate/genomes/H299_H561.fasta")
    message_cmd2 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "genome1.fasta-prokkaRes --cpus 3 --prefix test_runall_1by1_2 "
                    "--centre prokka test/data/annotate/genomes/genome1.fasta")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    # Messages end annotation cmd
    message_end1 = ("End annotating test_runall_1by1_1 from "
                    "test/data/annotate/genomes/H299_H561.fasta.")
    message_end2 = ("End annotating test_runall_1by1_2 from "
                    "test/data/annotate/genomes/genome1.fasta.")
    assert message_end1 in messages
    assert message_end2 in messages
    # Messages error annotation cmd
    message_err1 = "test_runall_1by1_2 genome1.fasta: several .faa files"
    assert message_err1 in messages

