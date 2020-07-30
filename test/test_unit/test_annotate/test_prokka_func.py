#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import os
import shutil
import logging
import PanACoTA.annotate_module.prokka_functions as pfunc
import PanACoTA.utils as utils

logfile_base = "panacota"

# Define methods and variables shared by several tests
def my_logger():
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
    return q, logging.getLogger('process')


def teardown_function(function):
    """
    Remove logfiles when test is done
    """
    if os.path.isfile(logfile_base + ".log"):
        os.remove(logfile_base + ".log")
    if os.path.isfile(logfile_base + ".log.err"):
        os.remove(logfile_base + ".log.err")
    if os.path.isfile(logfile_base + ".log.details"):
        os.remove(logfile_base + ".log.details")











def test_run_all_1by1():
    """
    Check that when running with 3 threads (not parallel), prokka runs as expected,
    and returns True for each genome
    Start and end must be ordered: (start1, end1, start2, end2) or (start2, end2, start1, end1)
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join("test", "data", "annotate", "genomes", genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join("test", "data", "annotate", "genomes", genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, 456464645, 1, 465]}
    threads = 3
    force = False
    prok_folder = os.path.join("test", "data", "annotate")
    res_dir = os.path.join("test", "data", "annotate")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
    #  arguments = [(genomes[g][1], prok_folder, cores_prokka, genomes[g][0],
                  # force, genomes[g][3], q)
    assert final[genome1]
    assert final[genome2]
    shutil.rmtree(os.path.join(prok_folder, genome1 + "-prokkaRes"))
    shutil.rmtree(os.path.join(prok_folder, genome2 + "-prokkaRes"))
    os.remove(os.path.join(res_dir, genome1 + "-prokka.log"))
    os.remove(os.path.join(res_dir, genome2 + "-prokka.log"))
    q = logger[0]
    assert q.qsize() == 5
    assert q.get().message == 'Annotating all genomes with prokka'
    # Messages for start and end annotation of the different genomes
    message_start_annot1 = ("Start annotating test_runall_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    message_end_annot1 = ("End annotating test_runall_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_end_annot2 = ("End annotating test_runall_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    qget = q.get().message
    assert qget == message_start_annot1 or message_start_annot2
    if qget == message_start_annot1:
        # Ending annotation of first genome (same genome as started because running 1by1)
        assert q.get().message == message_end_annot1
    else:
        assert q.get().message == message_end_annot2
    qget2 = q.get().message
    assert qget2 == message_start_annot1 or message_start_annot2
    if qget2 == message_start_annot2:
        # Ending annotation of first genome (same genome as started because running 1by1)
        assert q.get().message == message_end_annot2
    else:
        assert q.get().message == message_end_annot1


def test_run_all_parallel_more_threads():
    """
    Check that there is no problem when running with more threads than genomes (each genome
    uses nb_threads/nb_genome threads)
    Start and end are not necessarily in the same order (ex: start1, start2, end2, end1)
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join("test", "data", "annotate", "genomes", genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join("test", "data", "annotate", "genomes", genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    prok_folder = os.path.join("test", "data", "annotate")
    res_dir = os.path.join("test", "data", "annotate")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
    assert final[genome1]
    assert final[genome2]
    shutil.rmtree(os.path.join(prok_folder, genome1 + "-prokkaRes"))
    shutil.rmtree(os.path.join(prok_folder, genome2 + "-prokkaRes"))
    os.remove(os.path.join(res_dir, genome1 + "-prokka.log"))
    os.remove(os.path.join(res_dir, genome2 + "-prokka.log"))
    q = logger[0]
    assert q.qsize() == 5
    assert q.get().message == 'Annotating all genomes with prokka'
    message_start_annot1 = ("Start annotating test_run_all_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_start_annot2 = ("Start annotating test_run_all_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    # Starting annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_start_annot1 or message_start_annot2
    message_end_annot1 = ("End annotating test_run_all_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_end_annot2 = ("End annotating test_run_all_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    # Ending annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_end_annot1 or message_end_annot2
    # Starting annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_start_annot1 or message_start_annot2
    # Ending annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_end_annot1 or message_end_annot2


def test_run_all_parallel_less_threads():
    """
    Check that there is no problem when running with less threads than genomes (each genomes
    uses 2 threads)
    Genomes H299 and A_H738 should run well, but genomes genome* have problems (no CDS found),
    so check_prokka should return false.
    """
    utils.init_logger(logfile_base, 0, '')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    gnames = ["H299_H561.fasta", "A_H738.fasta", "genome1.fasta", "genome2.fasta", "genome3.fasta"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["test_runall_1by1_1", gpaths[0], 12656, 3, 1],
               gnames[1]: ["test_runall_1by1_2", gpaths[1], 456464645, 1, 1],
               gnames[2]: ["test_runall_1by1_2", gpaths[2], 456464645, 4, 1],
               gnames[3]: ["test_runall_1by1_2", gpaths[3], 456464645, 3, 1],
               gnames[4]: ["test_runall_1by1_2", gpaths[4], 456464645, 1, 1]
               }
    threads = 4
    force = False
    prok_folder = os.path.join("test", "data", "annotate")
    res_dir = os.path.join("test", "data", "annotate")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
    assert final[gnames[0]]
    assert final[gnames[1]]
    assert not final[gnames[2]]
    assert not final[gnames[3]]
    assert not final[gnames[4]]
    for name in gnames:
        shutil.rmtree(os.path.join(prok_folder, name + "-prokkaRes"))
        os.remove(os.path.join(res_dir, name + "-prokka.log"))

