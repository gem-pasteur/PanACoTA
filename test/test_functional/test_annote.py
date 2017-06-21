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


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert ("-d DB_PATH -r RES_PATH [-n NAME] [-Q] [--l90 L90]") in err
    assert "[--nbcont NBCONT] [--cutN CUTN] [--date DATE] [--tmp TMPDIR]" in err
    assert "[--prok PROKKADIR] [-F] [--threads THREADS] [-v] [-q] [-h]" in err
    assert "list_file" in err
    assert "the following arguments are required: list_file, -d, -r" in err


def test_parser_noname(capsys):
    """
    Test that when the script is called without any name for the genomes not -Q option,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath".split())
    _, err = capsys.readouterr()
    assert ("You must specify your genomes dataset name in 4 characters with "
            "'-n name' option (type -h for more information). Or, if you do not want "
            "to annotate and format your genomes but just to run quality control, use "
            "option '-Q") in err


def test_parser_wrongname(capsys):
    """
    Test that when the script is called with a genome name with more than 4 characters,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n genome".split())
    _, err = capsys.readouterr()
    assert ("The genome name must contain 4 characters. For example, this name can "
            "correspond to the 2 first letters of genus, and 2 first letters of "
            "species, e.g. ESCO for Escherichia Coli.") in err


def test_parser_negativeCont(capsys):
    """
    Test that when the script is called with a limit of contig number higher than 9999,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --nbcont -5".split())
    _, err = capsys.readouterr()
    assert ("The maximum number of contigs allowed must be a positive number.") in err


def test_parser_highCont(capsys):
    """
    Test that when the script is called with a negative limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --nbcont 10005".split())
    _, err = capsys.readouterr()
    assert ("We do not support genomes with more than 9999 contigs.") in err


def test_parser_wrongCont(capsys):
    """
    Test that when the script is called with a non integer limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --nbcont 10.5".split())
    _, err = capsys.readouterr()
    assert ("argument --nbcont: invalid int value: 10.5") in err


def test_parser_wrongl90(capsys):
    """
    Test that when the user does not give an int for the l90 limit, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --l90 l90".split())
    _, err = capsys.readouterr()
    assert ("argument --l90: invalid int value: 'l90'") in err


def test_parser_wrongCut(capsys):
    """
    Test that when the user does not give an int for the cutN value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --cutN 10.5".split())
    _, err = capsys.readouterr()
    assert ("argument --cutN: invalid int value: '10.5'") in err


def test_parser_wrongThread(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --threads 10.5".split())
    _, err = capsys.readouterr()
    assert ("argument --threads: invalid int value: '10.5'") in err


def test_parser_wrongDate(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --date 417".split())
    _, err = capsys.readouterr()
    assert (("The date must contain 4 characters. Usually, it contains 4 digits, "
             "corresponding to the month (2 digits) and year (2 digits).")) in err


def test_parser_qAndv(capsys):
    """
    Test that when the user wants both quiet and verbose option, it gives an error message
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 -q -v".split())
    _, err = capsys.readouterr()
    assert (("Choose between a verbose output (-v) or quiet output (-q)."
             " You cannot have both...")) in err


def test_parser_default():
    """
    Test that when run with the minimum required arguments, all default values are
    as expected.
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    OPTIONS = annot.parse(parser, "list_file -d dbpath -r respath -n g123".split())
    assert OPTIONS.list_file == "list_file"
    assert OPTIONS.db_path == "dbpath"
    assert OPTIONS.res_path == "respath"
    assert OPTIONS.name == "g123"
    assert OPTIONS.l90 == 100
    assert OPTIONS.nbcont == 999
    assert OPTIONS.cutn == 5
    assert OPTIONS.threads == 1
    assert OPTIONS.date == time.strftime("%m%y")
    assert OPTIONS.force == False
    assert OPTIONS.qc_only == False


def test_parser_values():
    """
    Test that values for L90, nbcontig, cutn, threads, date are taken into account
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    OPTIONS = annot.parse(parser, ("list_file -d dbpath -r respath -n g123 --l90 2 "
                                   "--nbcont 10 --cutN 0 --threads 8 --date toto").split())
    assert OPTIONS.list_file == "list_file"
    assert OPTIONS.db_path == "dbpath"
    assert OPTIONS.res_path == "respath"
    assert OPTIONS.name == "g123"
    assert OPTIONS.l90 == 2
    assert OPTIONS.nbcont == 10
    assert OPTIONS.cutn == 0
    assert OPTIONS.threads == 8
    assert OPTIONS.date == "toto"
    assert OPTIONS.force == False
    assert OPTIONS.qc_only == False


def test_parser_force():
    """
    Test that when run with '-F' option, force is initialized to "--force".
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    OPTIONS = annot.parse(parser, "list_file -d dbpath -r respath -n g123 -F".split())
    assert OPTIONS.list_file == "list_file"
    assert OPTIONS.db_path == "dbpath"
    assert OPTIONS.res_path == "respath"
    assert OPTIONS.name == "g123"
    assert OPTIONS.l90 == 100
    assert OPTIONS.nbcont == 999
    assert OPTIONS.cutn == 5
    assert OPTIONS.threads == 1
    assert OPTIONS.date == time.strftime("%m%y")
    assert OPTIONS.force == True
    assert OPTIONS.qc_only == False


def test_parser_wrongforce(capsys):
    """
    Test that when run with '-F' option + a value, it returns an error message.
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 -F 10".split())
    _, err = capsys.readouterr()
    assert "unrecognized arguments: 10" in err


def test_parser_qc():
    """
    Test that when run with '-Q' option (for QC only) and no name given for the genome, it
    is set to "NONE"
    """
    parser = argparse.ArgumentParser(description=("Annotate all genomes"), add_help=False)
    annot.build_parser(parser)
    OPTIONS = annot.parse(parser, "list_file -d dbpath -r respath -Q".split())
    assert OPTIONS.list_file == "list_file"
    assert OPTIONS.db_path == "dbpath"
    assert OPTIONS.res_path == "respath"
    assert OPTIONS.name == "NONE"
    assert OPTIONS.l90 == 100
    assert OPTIONS.nbcont == 999
    assert OPTIONS.cutn == 5
    assert OPTIONS.threads == 1
    assert OPTIONS.date == time.strftime("%m%y")
    assert OPTIONS.force == False
    assert OPTIONS.qc_only == True


def test_main_from_parse():
    """
    Test that when a tmp folder is given by user, tmp files are saved in it,
    and prokka files too.
    """
    logfile_base = "test_main_from_parse"
    utils.init_logger(logfile_base, 0, verbose=1)
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
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_main_given_tmp():
    """
    Test that when a tmp folder is given by user, tmp files are saved in it,
    and prokka files too.
    """
    logfile_base = "test_main_from_parse"
    utils.init_logger(logfile_base, 0, verbose=1)
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
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_main_given_prokka():
    """
    Test that when a prokka folder is given by user, tmp files are saved in result/tmp_files,
    and prokka files are saved in the given prokka folder.
    """
    logfile_base = "test_main_from_parse"
    utils.init_logger(logfile_base, 0, verbose=1)
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
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


# def test_main_given_tmpAndprokka():
#     """
#     Test that when a tmp folder and a prokka folder are given by user, tmp files are saved in
#     tmp folder, and prokka files in prokka folder.
#     """
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     resdir = os.path.join("test", "data", "res_test_funcGivenTmp")
#     prokdir = os.path.join("test", "data", "prokka_funcGivenTmp")
#     tmpdir = os.path.join("test", "data", "tmp_funcGivenTmp")
#     name = "ESCO"
#     l90 = 1
#     date = "0417"
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, cutn=0,
#                                          prok_dir=prokdir, tmp_dir=tmpdir)
#     assert skip == []
#     assert skipf == []
#     expk = {"A_H738.fasta-all.fna":
#                 ["ESCO.1015.00001",
#                  'test/data/tmp_funcGivenTmp/A_H738.fasta-all.fna-gembase.fna',
#                  20031, 5, 1]
#            }
#     assert kept == expk
#     expg = {"A_H738.fasta-all.fna":
#                 ["ESCO.1015.00001",
#                  'test/data/tmp_funcGivenTmp/A_H738.fasta-all.fna-gembase.fna',
#                  20031, 5, 1],
#             'H299_H561.fasta-all.fna':
#                 ['ESCO.1015',
#                  'test/data/genomes/H299_H561.fasta-all.fna',
#                  13259, 7, 3],
#             'B2_A3_5.fasta-changeName.fna':
#                 ['ESCO.1116',
#                  'test/data/genomes/B2_A3_5.fasta-changeName.fna',
#                  120529, 5, 4]
#            }
#     assert allg == expg
#     # Check that tmp files exist in the right folder (result/tmp_files)
#     assert os.path.isfile(os.path.join(tmpdir, "A_H738.fasta-all.fna-gembase.fna"))
#     # Test that prokka folder is in the right directory (given)
#     assert os.path.isdir(os.path.join(prokdir, "A_H738.fasta-all.fna-gembase.fna-prokkaRes"))
#     assert os.path.isfile(os.path.join(prokdir, "A_H738.fasta-all.fna-gembase.fna-prokka.log"))
#     shutil.rmtree(resdir, ignore_errors=True)
#     shutil.rmtree(prokdir, ignore_errors=True)
#     shutil.rmtree(tmpdir, ignore_errors=True)
#     os.remove(os.path.join(dbpath, "A_H738.fasta-all.fna"))
#     os.remove(os.path.join(dbpath, "H299_H561.fasta-all.fna"))


# def test_main_allDiscardNbcont():
#     """
#     Test that when the genomes given in list file have high nbcontigs compared
#     to given threshold, they appear in the list of genomes, but are not kept for analysis.
#     """
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     resdir = os.path.join("test", "data", "res_test_funcAllDiscardNbcont")
#     name = "ESCO"
#     l90 = 100
#     nbcont = 1
#     cutn = 0
#     threads = 1
#     date = "0417"
#     force = False
#     qc_only = False
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
#                                          cutn, threads, force, qc_only)
#     assert skip == []
#     assert skipf == []
#     assert kept == {}
#     expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
#                                              'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
#             'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3],
#             'genome1.fasta': ['ESCO.0417', 'test/data/genomes/genome1.fasta', 51, 4, 2]}
#     assert allg == expg
#     shutil.rmtree(resdir, ignore_errors=True)


# def test_main_allDiscardL90():
#     """
#     Test that when the genomes given in list file have high L90 compared
#     to given threshold, they appear in the list of genomes, but are not kept for analysis.
#     """
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     resdir = os.path.join("test", "data", "res_test_funcAllDiscardL90")
#     name = "ESCO"
#     l90 = 1
#     nbcont = 999
#     cutn = 0
#     threads = 1
#     date = "0417"
#     force = False
#     qc_only = False
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
#                                          cutn, threads, force, qc_only)
#     assert skip == []
#     assert skipf == []
#     assert kept == {}
#     expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
#                                              'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
#             'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3],
#             'genome1.fasta': ['ESCO.0417', 'test/data/genomes/genome1.fasta', 51, 4, 2]}
#     assert allg == expg
#     shutil.rmtree(resdir, ignore_errors=True)


# def test_main_QC():
#     """
#     Test that when only QC is run, it returns only the list of all genomes + list of genomes
#     that would be kept for an analysis. It does not return the list of genomes for
#     which annotation or format had problems, as this step did not run.
#     """
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     resdir = os.path.join("test", "data", "res_test_funcQC")
#     name = "ESCO"
#     l90 = 1
#     nbcont = 999
#     cutn = 0
#     threads = 1
#     date = "0417"
#     force = False
#     qc_only = True
#     allg, kept = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
#                                          cutn, threads, force, qc_only)
#     assert kept == {}
#     expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
#                                              'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
#             'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3],
#             'genome1.fasta': ['ESCO.0417', 'test/data/genomes/genome1.fasta', 51, 4, 2]}
#     assert allg == expg
#     shutil.rmtree(resdir, ignore_errors=True)


# def test_main_wrongSeqNames():
#     """
#     Test that when some genome names given in the list file do not exist in the
#     db path, they are not considered in the list of genomes to annotate/format. The others
#     appear in the list of genomes.
#     """
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-wrongNames.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     resdir = os.path.join("test", "data", "res_test_funcAllDiscardL90")
#     name = "ESCO"
#     l90 = 1
#     nbcont = 999
#     cutn = 0
#     threads = 1
#     date = "0417"
#     force = False
#     qc_only = False
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, l90, nbcont,
#                                          cutn, threads, force, qc_only)
#     assert skip == []
#     assert skipf == []
#     assert kept == {}
#     expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
#                                              'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
#             'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3]}
#     assert allg == expg
#     shutil.rmtree(resdir, ignore_errors=True)


# def test_main_onExistingProkkaDir():
#     """
#     Test that, when the pipeline is run with a given prokka dir, where prokka results already
#     exist, and are ok, all runs well
#     """
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-exist_dir.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     resdir = os.path.join("test", "data", "res_test_funcExistingProkka")
#     prokdir = os.path.join("test", "data", "exp_files")
#     name = "ESCO"
#     date = "0417"
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
#                                          prok_dir=prokdir)
#     assert skip == []
#     assert skipf == []
#     expg = {'H299_H561.fasta':
#                 ['ESCO.1015.00001',
#                  'test/data/res_test_funcExistingProkka/tmp_files/H299_H561.fasta-gembase.fna',
#                  13143, 3, 3],
#             'B2_A3_5.fasta-changeName.fna':
#                 ['ESCO.1116.00002',
#                  'test/data/res_test_funcExistingProkka/tmp_files/B2_A3_5.fasta-changeName.fna-gembase.fna',
#                  120529, 5, 4]
#            }
#     assert allg == expg
#     assert kept == expg
#     # Check that tmp files exist in the right folder (result/tmp_files)
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "B2_A3_5.fasta-changeName.fna-gembase.fna"))
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "H299_H561.fasta-gembase.fna"))
#     # Test that prokka folder was not recreated
#     assert not os.path.isdir(os.path.join(resdir, "tmp_files",
#                                           "B2_A3_5.fasta-changeName.fna-gembase.fna-prokkaRes"))
#     assert not os.path.isdir(os.path.join(resdir, "tmp_files",
#                                           "H299_H561.fasta-gembase.fna-prokkaRes"))
#     # Test that result files are in result dir
#     assert os.path.isfile(os.path.join(resdir, "LSTINFO-list_genomes-func-test-exist_dir.lst"))
#     shutil.rmtree(resdir, ignore_errors=True)


# def test_main_onExistingProkkaDirErrorProkk(capsys):
#     """
#     Test that, when the pipeline is run with a given prokka dir, where prokka results have
#     problems (no tbl file), it returns an error message and the genome with problems
#     is in skipped.
#     """
#     list_file = os.path.join("test", "data", "test_files",
#                              "list_genomes-func-test-exist-dir-err.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     genome_ori = os.path.join(dbpath, "B2_A3_5.fasta-changeName.fna")
#     genome_here = os.path.join(dbpath, "B2_A3_5.fasta-problems.fna")
#     shutil.copyfile(genome_ori, genome_here)
#     resdir = os.path.join("test", "data", "res_test_ProkErr")
#     prokdir = os.path.join("test", "data", "exp_files")
#     name = "ESCO"
#     date = "0417"
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
#                                          prok_dir=prokdir, verbose=True)
#     assert skip == ['B2_A3_5.fasta-problems.fna']
#     assert skipf == []
#     expg = {'H299_H561.fasta':
#                 ['ESCO.1015.00001',
#                  'test/data/res_test_ProkErr/tmp_files/H299_H561.fasta-gembase.fna',
#                  13143, 3, 3],
#             'B2_A3_5.fasta-problems.fna':
#                 ['ESCO.1116.00002',
#                  'test/data/res_test_ProkErr/tmp_files/B2_A3_5.fasta-problems.fna-gembase.fna',
#                  120529, 5, 4]
#            }
#     assert allg == expg
#     assert kept == expg
#     # Check that tmp files exist in the right folder (result/tmp_files)
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "B2_A3_5.fasta-problems.fna-gembase.fna"))
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "H299_H561.fasta-gembase.fna"))
#     # Test that prokka folder was not recreated
#     assert not os.path.isdir(os.path.join(resdir, "tmp_files",
#                                           "B2_A3_5.fasta-problems.fna-gembase.fna-prokkaRes"))
#     assert not os.path.isdir(os.path.join(resdir, "tmp_files",
#                                           "H299_H561.fasta-gembase.fna-prokkaRes"))
#     # Test that result files are in result dir
#     assert os.path.isfile(os.path.join(resdir, "LSTINFO-list_genomes-func-test-exist-dir-err.lst"))
#     # Check error messages
#     _, err = capsys.readouterr()
#     assert ("Prokka had problems while annotating some genomes, or did not find any gene. "
#             "Hence, they are not formatted, and absent from your output database. Please look "
#             "at their Prokka logs (<output_directory>/tmp_files/<genome_name>-prokka.log) and "
#             "to the current error log (<output_directory>/<input_filename>.log.err) to get "
#             "more information, and run again to annotate and format them. Here are the "
#             "genomes (problem with prokka or no gene found):") in err
#     assert "- B2_A3_5.fasta-problems.fna" in err
#     shutil.rmtree(resdir, ignore_errors=True)
#     os.remove(genome_here)


# def test_main_onExistingProkkaDirErrorForm(capsys):
#     """
#     Test that, when the pipeline is run with a given prokka dir, where prokka results are ok
#     (good files), but have problems inside (wrong header format), it returns an error
#     message and the genome with problems is in skipped_format.
#     """
#     list_file = os.path.join("test", "data", "test_files",
#                              "list_genomes-func-test-exist-dir-err.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     genome_ori = os.path.join(dbpath, "B2_A3_5.fasta-changeName.fna")
#     genome_here = os.path.join(dbpath, "B2_A3_5.fasta-problems.fna")
#     shutil.copyfile(genome_ori, genome_here)
#     resdir = os.path.join("test", "data", "res_test_ProkErr")
#     prokdir = os.path.join("test", "data", "exp_files")
#     tblInit = os.path.join(prokdir, "B2_A3_5.fasta-changeName.fna-gembase.fna-prokkaRes",
#                            "test.0417.00002.tbl")
#     tblHere = os.path.join(prokdir, "B2_A3_5.fasta-problems.fna-gembase.fna-prokkaRes",
#                            "test.0417.00002.tbl")
#     shutil.copyfile(tblInit, tblHere)
#     name = "ESCO"
#     date = "0417"
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
#                                          prok_dir=prokdir, verbose=True)
#     assert skip == []
#     assert skipf == ['B2_A3_5.fasta-problems.fna']
#     expg = {'H299_H561.fasta':
#                 ['ESCO.1015.00001',
#                  'test/data/res_test_ProkErr/tmp_files/H299_H561.fasta-gembase.fna',
#                  13143, 3, 3],
#             'B2_A3_5.fasta-problems.fna':
#                 ['ESCO.1116.00002',
#                  'test/data/res_test_ProkErr/tmp_files/B2_A3_5.fasta-problems.fna-gembase.fna',
#                  120529, 5, 4]
#            }
#     assert allg == expg
#     assert kept == expg
#     # Check that tmp files exist in the right folder (result/tmp_files)
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "B2_A3_5.fasta-problems.fna-gembase.fna"))
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "H299_H561.fasta-gembase.fna"))
#     # Test that prokka folder was not recreated
#     assert not os.path.isdir(os.path.join(resdir, "tmp_files",
#                                           "B2_A3_5.fasta-problems.fna-gembase.fna-prokkaRes"))
#     assert not os.path.isdir(os.path.join(resdir, "tmp_files",
#                                           "H299_H561.fasta-gembase.fna-prokkaRes"))
#     # Test that result files are in result dir
#     assert os.path.isfile(os.path.join(resdir, "LSTINFO-list_genomes-func-test-exist-dir-err.lst"))
#     # Check error messages
#     _, err = capsys.readouterr()
#     assert ("Some genomes were annotated by prokka, but could not be formatted, "
#             "and are hence absent from your output database. Please look at log "
#             "files to get more information about why they could not be "
#             "formatted.") in err
#     assert "- B2_A3_5.fasta-problems.fna" in err
#     shutil.rmtree(resdir, ignore_errors=True)
#     os.remove(genome_here)
#     os.remove(tblHere)


# def test_run_exist_resdir(capsys):
#     """
#     Test that when the pipeline is called, with a given resdir which already contains
#     results, the program ends, with an error message.
#     """
#     resdir = os.path.join("test", "data", "test_func_resdir")
#     # Create output directory with a lst file in LSTINFO
#     os.makedirs(os.path.join(resdir, "Proteins"))
#     open(os.path.join(resdir, "Proteins", "toto.prt"), "w").close()
#     with pytest.raises(SystemExit):
#         annot.main("list_file.lst", "path/db", resdir, "toto", "0123")
#     out, err = capsys.readouterr()
#     assert ("ERROR: Your output directory already has .prt files in the "
#             "Proteins folder. Provide another result directory, or remove the "
#             "files in this one.") in err
#     shutil.rmtree(resdir)


# def test_main_onExistResDirForce():
#     """
#     Test that, when the pipeline is run on an existing result directory, but force option is on,
#     it removes the result folders and runs again.
#     """
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-allDiscard.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     resdir = os.path.join("test", "data", "res_test_funcExistResdirForce")
#     # Create output directory with a lst file in LSTINFO
#     os.makedirs(os.path.join(resdir, "Proteins"))
#     open(os.path.join(resdir, "Proteins", "toto.prt"), "w").close()
#     assert os.path.isfile(os.path.join(resdir, "Proteins", "toto.prt"))
#     name = "ESCO"
#     date = "0417"
#     allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, date, cutn=0,
#                                          force=True, threads=8)
#     assert skip == ["genome1.fasta"]
#     assert skipf == []
#     # 2 solutions as genome2 and H299 have the same L90 and nbcontig: randomly choose the
#     # order between the 2.
#     expg = {'H299_H561.fasta':
#                 ['ESCO.0417.00002',
#                  'test/data/res_test_funcExistResdirForce/tmp_files/H299_H561.fasta-gembase.fna',
#                  13143, 3, 3],
#             'B2_A3_5.fasta-changeName.fna':
#                 ['ESCO.0417.00003',
#                  'test/data/res_test_funcExistResdirForce/tmp_files/B2_A3_5.fasta-changeName.fna-gembase.fna',
#                  120529, 5, 4],
#             'genome1.fasta':
#                 ['ESCO.0417.00001',
#                  'test/data/res_test_funcExistResdirForce/tmp_files/genome1.fasta-gembase.fna',
#                  51, 4, 2]
#            }

#     assert allg == expg
#     assert kept == expg
#     # Check that tmp files exist in the right folder (result/tmp_files)
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "B2_A3_5.fasta-changeName.fna-gembase.fna"))
#     assert os.path.isfile(os.path.join(resdir, "tmp_files",
#                                        "H299_H561.fasta-gembase.fna"))
#     # Test that prokka files exist in the right folder (resdir/tmp_files)
#     assert os.path.isdir(os.path.join(resdir, "tmp_files",
#                                       "B2_A3_5.fasta-changeName.fna-gembase.fna-prokkaRes"))
#     assert os.path.isdir(os.path.join(resdir, "tmp_files",
#                                       "H299_H561.fasta-gembase.fna-prokkaRes"))
#     # Test that result folders were removed before run: toto.prt should not be in Proteins anymore
#     assert not os.path.isfile(os.path.join(resdir, "Proteins", "toto.prt"))
#     shutil.rmtree(resdir, ignore_errors=True)


# def test_annote_all():
#     """
#     Test that when we call the pipeline with all default parameters, all expected output files
#     are created. Check the content of result files (LSTINFO, Genes, Proteins, Replicons),
#     lstinfo, discarded and log files.
#     Test that prokka files are in the right folder (result/tmp_files), and that tmp files
#     are also in results/tmp_files (no tmp or prokka dir given, so default is used)
#     """
#     date = time.strftime("%m%y")
#     fulldate = time.strftime("%Y-%m-%d")
#     list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
#     dbpath = os.path.join("test", "data", "genomes")
#     respath = os.path.join("test", "data", "res_test_funcDefault")
#     name = "GENO"
#     cmd = "genomeAPCAT annotate {} -d {} -r {} -n {}".format(list_file, dbpath, respath, name)
#     ret = subprocess.call(cmd.split())
#     assert ret == 0
#     # Get output files
#     log_files = [os.path.join(respath, "genomeAPCAT-annotate_list_genomes-func-test-default" + ext)
#                  for ext in [".log", ".log.err"]]
#     lstfile = [os.path.join(respath, "LSTINFO-list_genomes-func-test-default.lst")]
#     discfile = [os.path.join(respath, "discarded-list_genomes-func-test-default.lst")]
#     QCfiles = [os.path.join(respath, qc + "list_genomes-func-test-default.png")
#                for qc in ["QC_L90-", "QC_nb-contigs-"]]
#     genomes = {"B2_A3_5.fasta-changeName.fna": "ESCO.1116.00002",
#                "H299_H561.fasta-all.fna": "ESCO.1015.00001",
#                "A_H738.fasta-all.fna": "GENO.1015.00001"}
#     split5Nfiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna") for g in genomes]
#     gembasefiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna")
#                     for g in genomes]
#     proklogfiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna-prokka.log")
#                     for g in genomes]
#     prokka_files = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna-prokkaRes",
#                                  genomes[g])
#                     for g in genomes]
#     lstinffiles = [os.path.join(respath, "LSTINFO", name + ".lst")
#                    for name in list(genomes.values())]
#     prtfiles = [os.path.join(respath, "Proteins", name + ".prt")
#                 for name in list(genomes.values())]
#     genfiles = [os.path.join(respath, "Genes", name + ".gen")
#                 for name in list(genomes.values())]
#     repfiles = [os.path.join(respath, "Replicons", name + ".fna")
#                 for name in list(genomes.values())]

#     # Check all output files exist
#     for f in (log_files + lstfile + discfile + QCfiles + split5Nfiles + gembasefiles +
#               proklogfiles + lstinffiles + prtfiles + genfiles + repfiles):
#         assert os.path.isfile(f)
#     for f in prokka_files:
#         assert os.path.isfile(f + ".tbl")
#         assert os.path.isfile(f + ".faa")
#         assert os.path.isfile(f + ".ffn")

#     # TODO: check next steps, update with changed lstfile.
#     # Check content of result database files
#     exp_dir = os.path.join("test", "data", "exp_files", "results_test_func-default")
#     for f in (lstinffiles + prtfiles + genfiles + repfiles):
#         exp_file = os.sep.join(f.split(os.sep)[-2:])
#         exp_file_path = os.path.join(exp_dir, exp_file)
#         with open(exp_file_path, "r") as expf, open(f, "r") as outf:
#             for line_exp, line_out in zip(expf, outf):
#                 assert line_exp == line_out

#     # Check that log files (error and log) contain expected information
#     with open(log_files[0], "r") as logf:
#         infos = []
#         for line in logf:
#             assert line.startswith("[" + fulldate)
#             assert line.count("::") == 2
#             infos.append(line.strip().split("::")[-1].strip())
#         assert "Reading genomes" in infos
#         assert ("Cutting genomes at each stretch of at least 5 'N', and then, calculating "
#                 "genome size, number of contigs and L90.") in infos
#         assert ("Annotating all genomes with prokka") in infos
#         assert ("Formatting all genomes") in infos
#     with open(log_files[1], "r") as logf:
#         infos = []
#         for line in logf:
#             assert line.startswith("[" + fulldate)
#             assert line.count("::") == 2
#             infos.append(line.strip().split("::")[-1].strip())
#         assert ("toto.fst genome file does not exist. Its file will be ignored when "
#                 "concatenating ['A_H738.fasta', 'toto.fst', 'genome1.fasta', "
#                 "'genome.fst']") in infos
#         assert ("genome.fst genome file does not exist. Its file will be ignored when "
#                 "concatenating ['A_H738.fasta', 'toto.fst', 'genome1.fasta', "
#                 "'genome.fst']") in infos
#     # Check file content of genomes "-all" files
#     # H299_H561.fasta and genome6.fasta
#     exp_Hand6 = os.path.join(dbpath, "H299_H561-and-genome6.fna")
#     Handg6 = os.path.join(dbpath, "H299_H561.fasta-all.fna")
#     with open(exp_Hand6, "r") as expf, open(Handg6, "r") as outf:
#         for line_exp, line_out in zip(expf, outf):
#             assert line_exp == line_out
#     os.remove(Handg6)
#     # A_H738.fasta and genome1.fasta
#     exp_Aand1 = os.path.join(dbpath, "A_H738-and-genome1.fna")
#     Aand1 = os.path.join(dbpath, "A_H738.fasta-all.fna")
#     with open(exp_Aand1, "r") as expf, open(Aand1, "r") as outf:
#         for line_exp, line_out in zip(expf, outf):
#             assert line_exp == line_out
#     os.remove(Aand1)
#     # Check lstinfo file content
#     explst = os.path.join(exp_dir, "LSTINFO-list_genomes-func-test-default.lst")
#     with open(lstfile[0], "r") as expf, open(explst, "r") as outf:
#         for line_exp, line_out in zip(expf, outf):
#             assert line_exp == line_out
#     shutil.rmtree(respath)

