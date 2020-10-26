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


def test_run_all_1by1():
    """
    Check that when running with 3 threads (not parallel), prokka runs as expected,
    and returns True for each genome
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
    annot_folder = os.path.join(GENEPATH, "annot-folder")
    os.makedirs(annot_folder)
    final = afunc.run_annotation_all(genomes, threads, force, annot_folder)
    assert final[genome1]
    assert final[genome2]
    q = logger[0]
    assert q.qsize() == 7
    assert q.get().message == 'Annotating all genomes with prokka'
    # Messages for start and end annotation of the different genomes
    message_start_annot1 = ("Start annotating test_runall_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_cmd1 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "annot-folder/H299_H561.fasta-prokkaRes")
    message_end_annot1 = ("End annotating test_runall_1by1_1 from test/data/annotate/genomes/"
                            "H299_H561.fasta.")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    message_cmd2 = ("Prokka command: prokka --outdir test/data/annotate/generated_by_unit-tests/"
                    "annot-folder/A_H738.fasta-prokkaRes")
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


def test_run_all_parallel_more_threads():
    """
    Check that there is no problem when running with more threads than genomes (each genome
    uses nb_threads/nb_genome threads)
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
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH,
                                     prodigal_only=True, small=True, quiet=True)
    assert final[genome1]
    assert final[genome2]
    q = logger[0]
    assert q.qsize() == 7
    assert q.get().message == 'Annotating all genomes with prodigal'
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
                    "H299_H561.fasta-prodigalRes/test_runall_1by1_1.gff -q -p meta")
    message_cmd2 = ("Prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
                    "-d test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.ffn -a test/data/annotate/generated_by_unit-tests/"
                    "A_H738.fasta-prodigalRes/test_runall_1by1_2.faa -f gff "
                    "-o test/data/annotate/generated_by_unit-tests/A_H738.fasta-prodigalRes/"
                    "test_runall_1by1_2.gff -q -p meta")
    assert message_cmd1 in messages
    assert message_cmd2 in messages
    message_end_annot1 = ("End annotating test_runall_1by1_1 (from test/data/annotate/genomes/"
                          "H299_H561.fasta)")
    message_end_annot2 = ("End annotating test_runall_1by1_2 (from test/data/annotate/genomes/"
                          "A_H738.fasta)")
    assert message_end_annot1 in messages
    assert message_end_annot2 in messages


def test_run_all_parallel_less_threads():
    """
    Check that there is no problem when running with less threads than genomes (each genomes
    uses 2 threads)
    Genomes H299 and A_H738 should run well, but genomes genome* have problems (no CDS found),
    so check_prokka should return false.
    """
    logger = my_logger("test_run_all_parallel_more_threads")
    utils.init_logger(LOGFILE_BASE, 0, 'test_run_all_parallel_more_threads')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    gnames = ["H299_H561.fasta", "A_H738.fasta", "genome1.fasta", "genome2.fasta", "genome3.fasta"]
    gpaths = [os.path.join(GEN_PATH, name) for name in gnames]
    genomes = {gnames[0]: ["test_runall_1by1_1", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["test_runall_1by1_2", gpaths[1], gpaths[1], 456464645, 1, 1],
               gnames[2]: ["test_runall_1by1_2", gpaths[2], gpaths[2], 456464645, 4, 1],
               gnames[3]: ["test_runall_1by1_2", gpaths[3], gpaths[3], 456464645, 3, 1],
               gnames[4]: ["test_runall_1by1_2", gpaths[4], gpaths[4], 456464645, 1, 1]
               }
    threads = 4
    force = False
    final = afunc.run_annotation_all(genomes, threads, force, GENEPATH)
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
    # Check at least 1st log
    assert q.get().message == "Annotating all genomes with prokka"


