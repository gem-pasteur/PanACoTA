#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for annote_pipeline.py
"""

import pytest
import os
import glob
import subprocess
import shutil
import time
import logging
import argparse

from genomeAPCAT.subcommands import qc_and_annote as annot
import genomeAPCAT.utils as utils


@pytest.fixture(scope="module")
def logger():
    print("--- Init logger ---")
    logfile_base = "test_main_from_parse"
    utils.init_logger(logfile_base, 0, '', verbose=1)
    yield logging.getLogger()
    print("--- Clean logs ---")
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_main_from_parse(logger):
    """
    Test that when a tmp folder is given by user, tmp files are saved in it,
    and prokka files too.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcfromparse")
    tmpdir = os.path.join("test", "data", "tmp_funcGivenTmp")
    name = "ESCO"
    l90 = 1
    date = "0417"
    args = argparse.Namespace()
    args.list_file = list_file
    args.db_path = dbpath
    args.res_path = resdir
    args.name = name
    args.date = date
    args.l90 = l90
    args.nbcont = 999
    args.cutn = 0
    args.threads = 1
    args.force = False
    args.qc_only = False
    args.tmpdir = tmpdir
    args.prokkadir = None
    args.verbose=False
    args.quiet=False
    annot.main_from_parse(args)
    # Check that tmp files exist in the right folder
    assert os.path.isfile(os.path.join(tmpdir, "A_H738.fasta-all.fna-short-contig.fna"))
    # Test that prokka folder is in the right directory
    assert os.path.isdir(os.path.join(tmpdir, "A_H738.fasta-all.fna-short-contig.fna-prokkaRes"))
    shutil.rmtree(tmpdir, ignore_errors=True)
    shutil.rmtree(resdir, ignore_errors=True)
    os.remove(os.path.join(dbpath, "A_H738.fasta-all.fna"))
    os.remove(os.path.join(dbpath, "H299_H561.fasta-all.fna"))


def test_main_noValidGenome(capsys):
    """
    Test that when, in the list file, all genomes are wrong (do not correspond to
    filenames in the given dbpath), it closes the program with an error message.
    """
    list_file = os.path.join("test", "data", "test_files", "list_no_genome.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_NoGenome")
    name = "ESCO"
    date = "0417"
    with pytest.raises(SystemExit):
        annot.main(list_file, dbpath, resdir, name, date)
    _, err = capsys.readouterr()
    assert ("We did not find any genome listed in test/data/test_files/list_no_genome.txt "
            "in the folder test/data/genomes. Please check your list to give valid "
            "genome names.") in err
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_given_tmp(logger):
    """
    Test that when a tmp folder is given by user, tmp files are saved in it,
    and prokka files too.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcGivenTmp")
    tmpdir = os.path.join("test", "data", "tmp_funcGivenTmp")
    name = "ESCO"
    l90 = 1
    date = "0417"
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, cutn=0,
                                         tmp_dir=tmpdir)
    assert skip == []
    assert skipf == []
    expk = {"A_H738.fasta-all.fna":
                ["ESCO.1015.00001",
                 'test/data/tmp_funcGivenTmp/A_H738.fasta-all.fna-short-contig.fna',
                 20031, 5, 1]
           }
    assert kept == expk
    expg = {"A_H738.fasta-all.fna":
                ["ESCO.1015.00001",
                 'test/data/tmp_funcGivenTmp/A_H738.fasta-all.fna-short-contig.fna',
                 20031, 5, 1],
            'H299_H561.fasta-all.fna':
                ['ESCO.1015',
                 'test/data/tmp_funcGivenTmp/H299_H561.fasta-all.fna-short-contig.fna',
                 13259, 7, 3],
            'B2_A3_5.fasta-changeName.fna':
                ['ESCO.1116',
                 'test/data/tmp_funcGivenTmp/B2_A3_5.fasta-changeName.fna-short-contig.fna',
                 120529, 5, 4]
           }
    assert allg == expg
    # Check that tmp files exist in the right folder
    assert os.path.isfile(os.path.join(tmpdir, "A_H738.fasta-all.fna-short-contig.fna"))
    # Test that prokka folder is in the right directory
    assert os.path.isdir(os.path.join(tmpdir, "A_H738.fasta-all.fna-short-contig.fna-prokkaRes"))
    shutil.rmtree(resdir, ignore_errors=True)
    shutil.rmtree(tmpdir, ignore_errors=True)
    os.remove(os.path.join(dbpath, "A_H738.fasta-all.fna"))
    os.remove(os.path.join(dbpath, "H299_H561.fasta-all.fna"))


def test_main_given_prokka(logger):
    """
    Test that when a prokka folder is given by user, tmp files are saved in result/tmp_files,
    and prokka files are saved in the given prokka folder.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcGivenTmp")
    prokdir = os.path.join("test", "data", "prokka_funcGivenTmp")
    name = "ESCO"
    l90 = 1
    date = "0417"
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, cutn=0,
                                         prok_dir=prokdir)
    assert skip == []
    assert skipf == []
    expk = {"A_H738.fasta-all.fna":
                ["ESCO.1015.00001",
                 'test/data/res_test_funcGivenTmp/tmp_files/A_H738.fasta-all.fna-short-contig.fna',
                 20031, 5, 1]
           }
    assert kept == expk
    expg = {"A_H738.fasta-all.fna":
                ["ESCO.1015.00001",
                 'test/data/res_test_funcGivenTmp/tmp_files/A_H738.fasta-all.fna-short-contig.fna',
                 20031, 5, 1],
            'H299_H561.fasta-all.fna':
                ['ESCO.1015',
                 ('test/data/res_test_funcGivenTmp/tmp_files/'
                  'H299_H561.fasta-all.fna-short-contig.fna'),
                 13259, 7, 3],
            'B2_A3_5.fasta-changeName.fna':
                ['ESCO.1116',
                 ('test/data/res_test_funcGivenTmp/tmp_files/'
                  'B2_A3_5.fasta-changeName.fna-short-contig.fna'),
                 120529, 5, 4]
           }
    assert allg == expg
    # Check that tmp files exist in the right folder (result/tmp_files)
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "A_H738.fasta-all.fna-short-contig.fna"))
    # Test that prokka folder is in the right directory (given)
    assert os.path.isdir(os.path.join(prokdir, "A_H738.fasta-all.fna-short-contig.fna-prokkaRes"))
    shutil.rmtree(resdir, ignore_errors=True)
    shutil.rmtree(prokdir, ignore_errors=True)
    os.remove(os.path.join(dbpath, "A_H738.fasta-all.fna"))
    os.remove(os.path.join(dbpath, "H299_H561.fasta-all.fna"))


def test_main_given_tmpAndprokka(logger):
    """
    Test that when a tmp folder and a prokka folder are given by user, tmp files are saved in
    tmp folder, and prokka files in prokka folder.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcGivenTmp")
    prokdir = os.path.join("test", "data", "prokka_funcGivenTmp")
    tmpdir = os.path.join("test", "data", "tmp_funcGivenTmp")
    name = "ESCO"
    l90 = 1
    date = "0417"
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, cutn=0,
                                         prok_dir=prokdir, tmp_dir=tmpdir)
    assert skip == []
    assert skipf == []
    expk = {"A_H738.fasta-all.fna":
                ["ESCO.1015.00001",
                 'test/data/tmp_funcGivenTmp/A_H738.fasta-all.fna-short-contig.fna',
                 20031, 5, 1]
           }
    assert kept == expk
    expg = {"A_H738.fasta-all.fna":
                ["ESCO.1015.00001",
                 'test/data/tmp_funcGivenTmp/A_H738.fasta-all.fna-short-contig.fna',
                 20031, 5, 1],
            'H299_H561.fasta-all.fna':
                ['ESCO.1015',
                 'test/data/tmp_funcGivenTmp/H299_H561.fasta-all.fna-short-contig.fna',
                 13259, 7, 3],
            'B2_A3_5.fasta-changeName.fna':
                ['ESCO.1116',
                 'test/data/tmp_funcGivenTmp/B2_A3_5.fasta-changeName.fna-short-contig.fna',
                 120529, 5, 4]
           }
    assert allg == expg
    # Check that tmp files exist in the right folder (result/tmp_files)
    assert os.path.isfile(os.path.join(tmpdir, "A_H738.fasta-all.fna-short-contig.fna"))
    # Test that prokka folder is in the right directory (given)
    assert os.path.isdir(os.path.join(prokdir, "A_H738.fasta-all.fna-short-contig.fna-prokkaRes"))
    assert os.path.isfile(os.path.join(prokdir, "A_H738.fasta-all.fna-short-contig.fna-prokka.log"))
    shutil.rmtree(resdir, ignore_errors=True)
    shutil.rmtree(prokdir, ignore_errors=True)
    shutil.rmtree(tmpdir, ignore_errors=True)
    os.remove(os.path.join(dbpath, "A_H738.fasta-all.fna"))
    os.remove(os.path.join(dbpath, "H299_H561.fasta-all.fna"))


def test_main_allDiscardNbcont(logger):
    """
    Test that when the genomes given in list file have high nbcontigs compared
    to given threshold, they appear in the list of genomes, but are not kept for analysis.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcAllDiscardNbcont")
    name = "ESCO"
    l90 = 100
    nbcont = 1
    cutn = 0
    threads = 1
    date = "0417"
    force = False
    qc_only = False
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
                                         cutn, threads, force, qc_only)
    assert skip == []
    assert skipf == []
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna':
               ['ESCO.0417', ('test/data/res_test_funcAllDiscardNbcont/tmp_files/'
                              'B2_A3_5.fasta-changeName.fna-short-contig.fna'),
                120529, 5, 4],
            'H299_H561.fasta':
               ['ESCO.0417', ('test/data/res_test_funcAllDiscardNbcont/tmp_files/'
                              'H299_H561.fasta-short-contig.fna'),
                13143, 3, 3],
            'genome1.fasta':
               ['ESCO.0417', ('test/data/res_test_funcAllDiscardNbcont/tmp_files/'
                              'genome1.fasta-short-contig.fna'),
                51, 4, 2]
            }
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_allDiscardL90(logger):
    """
    Test that when the genomes given in list file have high L90 compared
    to given threshold, they appear in the list of genomes, but are not kept for analysis.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcAllDiscardL90")
    name = "ESCO"
    l90 = 1
    nbcont = 999
    cutn = 0
    threads = 1
    date = "0417"
    force = False
    qc_only = False
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
                                         cutn, threads, force, qc_only)
    assert skip == []
    assert skipf == []
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna':
                ['ESCO.0417', ('test/data/res_test_funcAllDiscardL90/tmp_files/'
                               'B2_A3_5.fasta-changeName.fna-short-contig.fna'),
                 120529, 5, 4],
            'H299_H561.fasta':
                ['ESCO.0417', ('test/data/res_test_funcAllDiscardL90/tmp_files/'
                               'H299_H561.fasta-short-contig.fna'),
                 13143, 3, 3],
            'genome1.fasta':
                ['ESCO.0417', ('test/data/res_test_funcAllDiscardL90/tmp_files/'
                               'genome1.fasta-short-contig.fna'),
                 51, 4, 2]
            }
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_QC():
    """
    Test that when only QC is run, it returns only the list of all genomes + list of genomes
    that would be kept for an analysis. It does not return the list of genomes for
    which annotation or format had problems, as this step did not run.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcQC")
    name = "ESCO"
    l90 = 1
    nbcont = 999
    cutn = 0
    threads = 1
    date = "0417"
    force = False
    qc_only = True
    allg, kept = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
                                         cutn, threads, force, qc_only)
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna':
                ['ESCO.0417', ('test/data/res_test_funcQC/tmp_files/'
                               'B2_A3_5.fasta-changeName.fna-short-contig.fna'),
                 120529, 5, 4],
            'H299_H561.fasta':
                ['ESCO.0417', ('test/data/res_test_funcQC/tmp_files/'
                               'H299_H561.fasta-short-contig.fna'),
                 13143, 3, 3],
            'genome1.fasta':
                ['ESCO.0417', ('test/data/res_test_funcQC/tmp_files/'
                               'genome1.fasta-short-contig.fna'),
                 51, 4, 2]
            }
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_wrongSeqNames(logger):
    """
    Test that when some genome names given in the list file do not exist in the
    db path, they are not considered in the list of genomes to annotate/format. The others
    appear in the list of genomes.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-wrongNames.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcAllDiscardL90")
    name = "ESCO"
    l90 = 1
    nbcont = 999
    cutn = 0
    threads = 1
    date = "0417"
    force = False
    qc_only = False
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
                                         cutn, threads, force, qc_only)
    assert skip == []
    assert skipf == []
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna':
                ['ESCO.0417', ('test/data/res_test_funcAllDiscardL90/tmp_files/'
                               'B2_A3_5.fasta-changeName.fna-short-contig.fna'),
                 120529, 5, 4],
            'H299_H561.fasta':
                ['ESCO.0417', ('test/data/res_test_funcAllDiscardL90/tmp_files/'
                               'H299_H561.fasta-short-contig.fna'),
                 13143, 3, 3]
            }
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_onExistingProkkaDir(logger):
    """
    Test that, when the pipeline is run with a given prokka dir, where prokka results already
    exist, and are ok, all runs well
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-exist_dir.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcExistingProkka")
    prokdir = os.path.join("test", "data", "exp_files")
    name = "ESCO"
    date = "0417"
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
                                         prok_dir=prokdir)
    assert skip == []
    assert skipf == []
    expg = {'H299_H561.fasta':
                ['ESCO.1015.00001',
                 ('test/data/res_test_funcExistingProkka/tmp_files/'
                  'H299_H561.fasta-short-contig.fna'),
                 13143, 3, 3],
            'B2_A3_5.fasta-changeName.fna':
                ['ESCO.1116.00002',
                 ('test/data/res_test_funcExistingProkka/tmp_files/'
                  'B2_A3_5.fasta-changeName.fna-short-contig.fna'),
                 120529, 5, 4]
           }
    assert allg == expg
    assert kept == expg
    # Check that tmp files exist in the right folder (result/tmp_files)
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "B2_A3_5.fasta-changeName.fna-short-contig.fna"))
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "H299_H561.fasta-short-contig.fna"))
    # Test that prokka folder was not recreated
    prok_fold = os.path.join(resdir, "tmp_files",
                             "B2_A3_5.fasta-changeName.fna-short-contig.fna-prokkaRes")
    assert not os.path.isdir(prok_fold)
    assert not os.path.isdir(os.path.join(resdir, "tmp_files",
                                          "H299_H561.fasta-short-contig.fna-prokkaRes"))
    # Test that result files are in result dir
    assert os.path.isfile(os.path.join(resdir, "LSTINFO-list_genomes-func-test-exist_dir.lst"))
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_onExistingProkkaDirErrorProkk(logger, capsys):
    """
    Test that, when the pipeline is run with a given prokka dir, where prokka results have
    problems (no tbl file), it returns an error message and the genome with problems
    is in skipped.
    """
    list_file = os.path.join("test", "data", "test_files",
                             "list_genomes-func-test-exist-dir-err.txt")
    dbpath = os.path.join("test", "data", "genomes")
    genome_ori = os.path.join(dbpath, "B2_A3_5.fasta-changeName.fna")
    genome_here = os.path.join(dbpath, "B2_A3_5.fasta-problems.fna")
    shutil.copyfile(genome_ori, genome_here)
    resdir = os.path.join("test", "data", "res_test_ProkErr")
    prokdir = os.path.join("test", "data", "exp_files")
    name = "ESCO"
    date = "0417"
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
                                         prok_dir=prokdir, verbose=True)
    assert skip == ['B2_A3_5.fasta-problems.fna']
    assert skipf == []
    expg = {'H299_H561.fasta':
                ['ESCO.1015.00001',
                 'test/data/res_test_ProkErr/tmp_files/H299_H561.fasta-short-contig.fna',
                 13143, 3, 3],
            'B2_A3_5.fasta-problems.fna':
                ['ESCO.1116.00002',
                 'test/data/res_test_ProkErr/tmp_files/B2_A3_5.fasta-problems.fna-short-contig.fna',
                 120529, 5, 4]
           }
    assert allg == expg
    assert kept == expg
    # Check that tmp files exist in the right folder (result/tmp_files)
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "B2_A3_5.fasta-problems.fna-short-contig.fna"))
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "H299_H561.fasta-short-contig.fna"))
    # Test that prokka folder was not recreated
    assert not os.path.isdir(os.path.join(resdir, "tmp_files",
                                          "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes"))
    assert not os.path.isdir(os.path.join(resdir, "tmp_files",
                                          "H299_H561.fasta-short-contig.fna-prokkaRes"))
    # Test that result files are in result dir
    assert os.path.isfile(os.path.join(resdir, "LSTINFO-list_genomes-func-test-exist-dir-err.lst"))
    # Check error messages
    _, err = capsys.readouterr()
    assert ("Prokka had problems while annotating some genomes, or did not find any gene. "
            "Hence, they are not formatted, and absent from your output database. Please look "
            "at their Prokka logs (<output_directory>/tmp_files/<genome_name>-prokka.log) and "
            "to the current error log (<output_directory>/<input_filename>.log.err) to get "
            "more information, and run again to annotate and format them. Here are the "
            "genomes (problem with prokka or no gene found):") in err
    assert "- B2_A3_5.fasta-problems.fna" in err
    shutil.rmtree(resdir, ignore_errors=True)
    os.remove(genome_here)


def test_main_onExistingProkkaDirErrorForm(logger, capsys):
    """
    Test that, when the pipeline is run with a given prokka dir, where prokka results are ok
    (good files), but have problems inside (wrong header format), it returns an error
    message and the genome with problems is in skipped_format.
    """
    list_file = os.path.join("test", "data", "test_files",
                             "list_genomes-func-test-exist-dir-err.txt")
    dbpath = os.path.join("test", "data", "genomes")
    genome_ori = os.path.join(dbpath, "B2_A3_5.fasta-changeName.fna")
    genome_here = os.path.join(dbpath, "B2_A3_5.fasta-problems.fna")
    shutil.copyfile(genome_ori, genome_here)
    resdir = os.path.join("test", "data", "res_test_ProkErr")
    prokdir = os.path.join("test", "data", "exp_files")
    tblInit = os.path.join(prokdir, "B2_A3_5.fasta-changeName.fna-short-contig.fna-prokkaRes",
                           "test.0417.00002.tbl")
    tblHere = os.path.join(prokdir, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                           "test.0417.00002.tbl")
    shutil.copyfile(tblInit, tblHere)
    name = "ESCO"
    date = "0417"
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
                                         prok_dir=prokdir, verbose=True)
    assert skip == []
    assert skipf == ['B2_A3_5.fasta-problems.fna']
    expg = {'H299_H561.fasta':
                ['ESCO.1015.00001',
                 'test/data/res_test_ProkErr/tmp_files/H299_H561.fasta-short-contig.fna',
                 13143, 3, 3],
            'B2_A3_5.fasta-problems.fna':
                ['ESCO.1116.00002',
                 ('test/data/res_test_ProkErr/tmp_files/'
                  'B2_A3_5.fasta-problems.fna-short-contig.fna'),
                 120529, 5, 4]
           }
    assert allg == expg
    assert kept == expg
    # Check that tmp files exist in the right folder (result/tmp_files)
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "B2_A3_5.fasta-problems.fna-short-contig.fna"))
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "H299_H561.fasta-short-contig.fna"))
    # Test that prokka folder was not recreated
    assert not os.path.isdir(os.path.join(resdir, "tmp_files",
                                          "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes"))
    assert not os.path.isdir(os.path.join(resdir, "tmp_files",
                                          "H299_H561.fasta-short-contig.fna-prokkaRes"))
    # Test that result files are in result dir
    assert os.path.isfile(os.path.join(resdir, "LSTINFO-list_genomes-func-test-exist-dir-err.lst"))
    # Check error messages
    _, err = capsys.readouterr()
    assert ("Some genomes were annotated by prokka, but could not be formatted, "
            "and are hence absent from your output database. Please look at log "
            "files to get more information about why they could not be "
            "formatted.") in err
    assert "- B2_A3_5.fasta-problems.fna" in err
    shutil.rmtree(resdir, ignore_errors=True)
    os.remove(genome_here)
    os.remove(tblHere)


def test_run_exist_resdir(capsys):
    """
    Test that when the pipeline is called, with a given resdir which already contains
    results, the program ends, with an error message.
    """
    resdir = os.path.join("test", "data", "test_func_resdir")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "Proteins"))
    open(os.path.join(resdir, "Proteins", "toto.prt"), "w").close()
    with pytest.raises(SystemExit):
        annot.main("list_file.lst", "path/db", resdir, "toto", "0123")
    out, err = capsys.readouterr()
    print("out", out)
    print("err", err)
    assert ("ERROR: Your output directory already has .prt files in the "
            "Proteins folder. Provide another result directory, or remove the "
            "files in this one.") in err
    shutil.rmtree(resdir)


def test_main_onExistResDirForce(logger):
    """
    Test that, when the pipeline is run on an existing result directory, but force option is on,
    it removes the result folders and runs again.
    """
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
    dbpath = os.path.join("test", "data", "genomes")
    resdir = os.path.join("test", "data", "res_test_funcExistResdirForce")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "Proteins"))
    open(os.path.join(resdir, "Proteins", "toto.prt"), "w").close()
    assert os.path.isfile(os.path.join(resdir, "Proteins", "toto.prt"))
    name = "ESCO"
    date = "0417"
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
                                         force=True, threads=8)
    assert skip == ["genome1.fasta"]
    assert skipf == []
    # 2 solutions as genome2 and H299 have the same L90 and nbcontig: randomly choose the
    # order between the 2.
    expg = {'H299_H561.fasta':
                ['ESCO.0417.00002',
                 ('test/data/res_test_funcExistResdirForce/tmp_files/'
                  'H299_H561.fasta-short-contig.fna'),
                 13143, 3, 3],
            'B2_A3_5.fasta-changeName.fna':
                ['ESCO.0417.00003',
                 ('test/data/res_test_funcExistResdirForce/tmp_files/'
                  'B2_A3_5.fasta-changeName.fna-short-contig.fna'),
                 120529, 5, 4],
            'genome1.fasta':
                ['ESCO.0417.00001',
                 ('test/data/res_test_funcExistResdirForce/tmp_files/'
                  'genome1.fasta-short-contig.fna'),
                 51, 4, 2]
           }

    assert allg == expg
    assert kept == expg
    # Check that tmp files exist in the right folder (result/tmp_files)
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "B2_A3_5.fasta-changeName.fna-short-contig.fna"))
    assert os.path.isfile(os.path.join(resdir, "tmp_files",
                                       "H299_H561.fasta-short-contig.fna"))
    # Test that prokka files exist in the right folder (resdir/tmp_files)
    assert os.path.isdir(os.path.join(resdir, "tmp_files",
                                      "B2_A3_5.fasta-changeName.fna-short-contig.fna-prokkaRes"))
    assert os.path.isdir(os.path.join(resdir, "tmp_files",
                                      "H299_H561.fasta-short-contig.fna-prokkaRes"))
    # Test that result folders were removed before run: toto.prt should not be in Proteins anymore
    assert not os.path.isfile(os.path.join(resdir, "Proteins", "toto.prt"))
    shutil.rmtree(resdir, ignore_errors=True)


def test_annote_all():
    """
    Test that when we call the pipeline with all default parameters, all expected output files
    are created. Check the content of result files (LSTINFO, Genes, Proteins, Replicons),
    lstinfo, discarded and log files.
    Test that prokka files are in the right folder (result/tmp_files), and that tmp files
    are also in results/tmp_files (no tmp or prokka dir given, so default is used)
    """
    date = time.strftime("%m%y")
    fulldate = time.strftime("%Y-%m-%d")
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
    dbpath = os.path.join("test", "data", "genomes")
    respath = os.path.join("test", "data", "res_test_funcDefault")
    name = "GENO"
    # check that log files do not already exist. If they do, remove them
    log_files = [os.path.join(respath, "genomeAPCAT-annotate_list_genomes-func-test-default" + ext)
                 for ext in [".log", ".log.err", ".log.details"]]
    for lf in log_files:
        if os.path.isfile(lf):
            os.remove(lf)
    cmd = "genomeAPCAT annotate {} -d {} -r {} -n {}".format(list_file, dbpath, respath, name)
    ret = subprocess.call(cmd.split())
    assert ret == 0
    # Get output files
    lstfile = [os.path.join(respath, "LSTINFO-list_genomes-func-test-default.lst")]
    discfile = [os.path.join(respath, "discarded-list_genomes-func-test-default.lst")]
    QCfiles = [os.path.join(respath, qc + "list_genomes-func-test-default.png")
               for qc in ["QC_L90-", "QC_nb-contigs-"]]
    genomes = {"B2_A3_5.fasta-changeName.fna": "ESCO.1116.00002",
               "H299_H561.fasta-all.fna": "ESCO.1015.00001",
               "A_H738.fasta-all.fna": "GENO.1015.00001"}
    split5Nfiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna") for g in genomes]
    proklogfiles = [os.path.join(respath, "tmp_files",
                                 g + "-split5N.fna-prokka.log")
                    for g in genomes]
    prokka_files = [os.path.join(respath, "tmp_files",
                                 g + "-split5N.fna-prokkaRes",
                                 genomes[g])
                    for g in genomes]
    lstinffiles = [os.path.join(respath, "LSTINFO", name + ".lst")
                   for name in list(genomes.values())]
    prtfiles = [os.path.join(respath, "Proteins", name + ".prt")
                for name in list(genomes.values())]
    genfiles = [os.path.join(respath, "Genes", name + ".gen")
                for name in list(genomes.values())]
    repfiles = [os.path.join(respath, "Replicons", name + ".fna")
                for name in list(genomes.values())]

    # Check all output files exist
    for f in (log_files + lstfile + discfile + QCfiles + split5Nfiles +
              proklogfiles + lstinffiles + prtfiles + genfiles + repfiles):
        assert os.path.isfile(f)
    for f in prokka_files:
        assert os.path.isfile(f + ".tbl")
        assert os.path.isfile(f + ".faa")
        assert os.path.isfile(f + ".ffn")

    # Check content of result database files
    exp_dir = os.path.join("test", "data", "exp_files", "results_test_func-default")
    for f in (lstinffiles + prtfiles + genfiles + repfiles):
        exp_file = os.sep.join(f.split(os.sep)[-2:])
        exp_file_path = os.path.join(exp_dir, exp_file)
        with open(exp_file_path, "r") as expf, open(f, "r") as outf:
            for line_exp, line_out in zip(expf, outf):
                assert line_exp == line_out

    # Check that log files (error, details and log) contain expected information
    with open(log_files[0], "r") as logf:
        infos = []
        for line in logf:
            assert line.startswith("[" + fulldate)
            assert line.count("::") == 2
            infos.append(line.split("::")[-1].strip())
        assert "Reading genomes" in infos
        assert ("Cutting genomes at each stretch of at least 5 'N', and then, calculating "
                "genome size, number of contigs and L90.") in infos
        assert ("Annotating all genomes with prokka") in infos
        assert ("Formatting all genomes") in infos
    with open(log_files[1], "r") as logf:
        infos = []
        for line in logf:
            assert line.startswith("[" + fulldate)
            assert line.count("::") == 2
            infos.append(line.split("::")[-1].strip())
        assert ("toto.fst genome file does not exist. Its file will be ignored when "
                "concatenating ['A_H738.fasta', 'toto.fst', 'genome1.fasta', "
                "'genome.fst']") in infos
        assert ("genome.fst genome file does not exist. Its file will be ignored when "
                "concatenating ['A_H738.fasta', 'toto.fst', 'genome1.fasta', "
                "'genome.fst']") in infos
    with open(log_files[2], "r") as logd:
        infos = []
        for line in logd:
            assert line.startswith("[" + fulldate)
            assert line.count("::") == 2
            infos.append(line.split("::")[-1].strip())
        print(infos)
        assert ("Start annotating GENO.1015.00001 {}/A_H738.fasta-all.fna"
                "-split5N.fna").format(os.path.join(respath, "tmp_files")) in infos
        assert ("Start annotating ESCO.1015.00001 {}/H299_H561.fasta-all.fna"
                "-split5N.fna").format(os.path.join(respath, "tmp_files")) in infos
        assert ("Start annotating ESCO.1116.00002 {}/B2_A3_5.fasta-changeName.fna"
                "-split5N.fna").format(os.path.join(respath, "tmp_files")) in infos
        assert ("End annotating GENO.1015.00001 {}/A_H738.fasta-all.fna"
                "-split5N.fna").format(os.path.join(respath, "tmp_files")) in infos
        assert ("End annotating ESCO.1015.00001 {}/H299_H561.fasta-all.fna"
                "-split5N.fna").format(os.path.join(respath, "tmp_files")) in infos
        assert ("End annotating ESCO.1116.00002 {}/B2_A3_5.fasta-changeName.fna"
                "-split5N.fna").format(os.path.join(respath, "tmp_files")) in infos
    # Check file content of genomes "-all" files
    # H299_H561.fasta and genome6.fasta
    exp_Hand6 = os.path.join(dbpath, "H299_H561-and-genome6.fna")
    Handg6 = os.path.join(dbpath, "H299_H561.fasta-all.fna")
    with open(exp_Hand6, "r") as expf, open(Handg6, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    os.remove(Handg6)
    # A_H738.fasta and genome1.fasta
    exp_Aand1 = os.path.join(dbpath, "A_H738-and-genome1.fna")
    Aand1 = os.path.join(dbpath, "A_H738.fasta-all.fna")
    with open(exp_Aand1, "r") as expf, open(Aand1, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    os.remove(Aand1)
    # Check lstinfo file content
    explst = os.path.join(exp_dir, "LSTINFO-list_genomes-func-test-default.lst")
    with open(lstfile[0], "r") as expf, open(explst, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    shutil.rmtree(respath)

