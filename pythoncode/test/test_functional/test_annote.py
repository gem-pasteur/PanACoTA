#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for annote_pipeline.py
"""

import pytest
import os
import subprocess
import shutil
import time
import logging

from pipelinepackage import annote_pipeline as annot


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    with pytest.raises(SystemExit):
        annot.parse("".split())
    _, err = capsys.readouterr()
    assert ("[-h] -d DB_PATH -r RES_PATH [-n NAME] [--l90 L90]") in err
    assert "[--nbcont NBCONT] [--cutN CUTN] [--threads THREADS]" in err
    assert "[--date DATE] [-F] [-Q]" in err
    assert "list_file" in err
    assert "the following arguments are required: list_file, -d, -r" in err


def test_parser_noname(capsys):
    """
    Test that when the script is called without any name for the genomes not -Q option,
    it returns an error message
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath".split())
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
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n genome".split())
    _, err = capsys.readouterr()
    assert ("The genome name must contain 4 characters. For example, this name can "
            " correspond to the 2 first letters of genus, and 2 first letters of "
            "species, e.g. ESCO for Escherichia Coli.") in err

def test_parser_negativeCont(capsys):
    """
    Test that when the script is called with a limit of contig number higher than 9999,
    it returns an error message
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 --nbcont -5".split())
    _, err = capsys.readouterr()
    assert ("The maximum number of contigs allowed must be a positive number.") in err


def test_parser_highCont(capsys):
    """
    Test that when the script is called with a negative limit of contig number,
    it returns an error message
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 --nbcont 10005".split())
    _, err = capsys.readouterr()
    assert ("We do not support genomes with more than 9999 contigs.") in err


def test_parser_wrongCont(capsys):
    """
    Test that when the script is called with a non integer limit of contig number,
    it returns an error message
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 --nbcont 10.5".split())
    _, err = capsys.readouterr()
    assert ("argument --nbcont: invalid int value: 10.5") in err


def test_parser_wrongl90(capsys):
    """
    Test that when the user does not give an int for the l90 limit, it returns an
    error message.
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 --l90 l90".split())
    _, err = capsys.readouterr()
    assert ("argument --l90: invalid int value: 'l90'") in err


def test_parser_wrongCut(capsys):
    """
    Test that when the user does not give an int for the cutN value, it returns an
    error message.
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 --cutN 10.5".split())
    _, err = capsys.readouterr()
    assert ("argument --cutN: invalid int value: '10.5'") in err


def test_parser_wrongThread(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 --threads 10.5".split())
    _, err = capsys.readouterr()
    assert ("argument --threads: invalid int value: '10.5'") in err


def test_parser_wrongDate(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 --date 417".split())
    _, err = capsys.readouterr()
    assert (("The date must contain 4 characters. Usually, it contains 4 digits, "
             "corresponding to the month (2 digits) and year (2 digits).")) in err


def test_parser_default():
    """
    Test that when run with the minimum required arguments, all default values are
    as expected.
    """
    OPTIONS = annot.parse("list_file -d dbpath -r respath -n g123".split())
    assert OPTIONS.list_file == "list_file"
    assert OPTIONS.db_path == "dbpath"
    assert OPTIONS.res_path == "respath"
    assert OPTIONS.name == "g123"
    assert OPTIONS.l90 == 100
    assert OPTIONS.nbcont == 999
    assert OPTIONS.cutn == 5
    assert OPTIONS.threads == 1
    assert OPTIONS.date == time.strftime("%m%y")
    assert OPTIONS.force == None
    assert OPTIONS.qc_only == False


def test_parser_values():
    """
    Test that values for L90, nbcontig, cutn, threads, date are taken into account
    """
    OPTIONS = annot.parse(("list_file -d dbpath -r respath -n g123 --l90 2 "
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
    assert OPTIONS.force == None
    assert OPTIONS.qc_only == False


def test_parser_force():
    """
    Test that when run with '-F' option, force is initialized to "--force".
    """
    OPTIONS = annot.parse("list_file -d dbpath -r respath -n g123 -F".split())
    assert OPTIONS.list_file == "list_file"
    assert OPTIONS.db_path == "dbpath"
    assert OPTIONS.res_path == "respath"
    assert OPTIONS.name == "g123"
    assert OPTIONS.l90 == 100
    assert OPTIONS.nbcont == 999
    assert OPTIONS.cutn == 5
    assert OPTIONS.threads == 1
    assert OPTIONS.date == time.strftime("%m%y")
    assert OPTIONS.force == "--force"
    assert OPTIONS.qc_only == False


def test_parser_wrongforce(capsys):
    """
    Test that when run with '-F' option + a value, it returns an error message.
    """
    with pytest.raises(SystemExit):
        annot.parse("list_file -d dbpath -r respath -n g123 -F 10".split())
    _, err = capsys.readouterr()
    assert "unrecognized arguments: 10" in err


def test_parser_qc():
    """
    Test that when run with '-Q' option (for QC only) and no name given for the genome, it
    is set to "NONE"
    """
    OPTIONS = annot.parse("list_file -d dbpath -r respath -Q".split())
    assert OPTIONS.list_file == "list_file"
    assert OPTIONS.db_path == "dbpath"
    assert OPTIONS.res_path == "respath"
    assert OPTIONS.name == "NONE"
    assert OPTIONS.l90 == 100
    assert OPTIONS.nbcont == 999
    assert OPTIONS.cutn == 5
    assert OPTIONS.threads == 1
    assert OPTIONS.date == time.strftime("%m%y")
    assert OPTIONS.force == None
    assert OPTIONS.qc_only == True


def test_logger_default(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = "logfile_test.txt"
    level = logging.DEBUG
    annot.init_logger(logfile, level, "default")
    logger = logging.getLogger("default")
    logger.debug("info debug")
    logger.info("info info")
    logger.warning("info warning")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert out == "  * info debug\n  * info info\n"
    assert err == "  * info warning\n  * info critical\n"
    with open(logfile, "r") as logf:
        assert logf.readline().endswith("] :: DEBUG :: info debug\n")
        assert logf.readline().endswith("] :: INFO :: info info\n")
        assert logf.readline().endswith("] :: WARNING :: info warning\n")
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    with open(logfile + ".err", "r") as logf:
        assert logf.readline().endswith("] :: WARNING :: info warning\n")
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    os.remove(logfile)
    os.remove(logfile + ".err")


def test_logger_info(capsys):
    """
    Test that when logger is initialized with "INFO" level, it does not return DEBUG info.
    """
    logfile = "logfile_test.txt"
    level = logging.INFO
    annot.init_logger(logfile, level, "info")
    logger = logging.getLogger("info")
    logger.debug("info debug")
    logger.info("info info")
    logger.warning("info warning")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert out == "  * info info\n"
    assert err == "  * info warning\n  * info critical\n"
    with open(logfile, "r") as logf:
        assert logf.readline().endswith("] :: INFO :: info info\n")
        assert logf.readline().endswith("] :: WARNING :: info warning\n")
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    with open(logfile + ".err", "r") as logf:
        assert logf.readline().endswith("] :: WARNING :: info warning\n")
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    os.remove(logfile)
    os.remove(logfile + ".err")


def test_logger_warning(capsys):
    """
    Test that when logger is initialized with "WARNING" level, it does not return
    anything in stdout, as DEBUG and INFO are not returned.
    """
    logfile = "logfile_test.txt"
    level = logging.WARNING
    annot.init_logger(logfile, level, "warn")
    logger = logging.getLogger("warn")
    logger.debug("info debug")
    logger.info("info info")
    logger.warning("info warning")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert out == ""
    assert err == "  * info warning\n  * info critical\n"
    with open(logfile, "r") as logf:
        assert logf.readline().endswith("] :: WARNING :: info warning\n")
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    with open(logfile + ".err", "r") as logf:
        assert logf.readline().endswith("] :: WARNING :: info warning\n")
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    os.remove(logfile)
    os.remove(logfile + ".err")


def test_logger_critical(capsys):
    """
    Test that when logger is initialized with "CRITICAL" level, it only returns
    CRITICAL information.

    """
    logfile = "logfile_test.txt"
    level = logging.CRITICAL
    annot.init_logger(logfile, level, "crit")
    logger = logging.getLogger("crit")
    logger.debug("info debug")
    logger.info("info info")
    logger.warning("info warning")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert out == ""
    assert err == "  * info critical\n"
    with open(logfile, "r") as logf:
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    with open(logfile + ".err", "r") as logf:
        assert logf.readline().endswith("] :: CRITICAL :: info critical\n")
    os.remove(logfile)
    os.remove(logfile + ".err")


def test_main_allDiscardNbcont():
    """

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
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, l90, nbcont,
                                         cutn, threads, date, force, qc_only)
    assert skip == []
    assert skipf == []
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
                                             'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
            'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3],
            'genome1.fasta': ['ESCO.0417', 'test/data/genomes/genome1.fasta', 51, 4, 2],
            'genome2.fasta': ['ESCO.0417', 'test/data/genomes/genome2.fasta', 67, 3, 3]}
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_allDiscardL90():
    """

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
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, l90, nbcont,
                                         cutn, threads, date, force, qc_only)
    assert skip == []
    assert skipf == []
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
                                             'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
            'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3],
            'genome1.fasta': ['ESCO.0417', 'test/data/genomes/genome1.fasta', 51, 4, 2],
            'genome2.fasta': ['ESCO.0417', 'test/data/genomes/genome2.fasta', 67, 3, 3]}
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_QC():
    """

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
    allg, kept = annot.main(list_file, dbpath, resdir, name, l90, nbcont,
                                         cutn, threads, date, force, qc_only)
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
                                             'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
            'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3],
            'genome1.fasta': ['ESCO.0417', 'test/data/genomes/genome1.fasta', 51, 4, 2],
            'genome2.fasta': ['ESCO.0417', 'test/data/genomes/genome2.fasta', 67, 3, 3]}
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_main_wrongSeqNames():
    """

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
    allg, kept, skip, skipf = annot.main(list_file, dbpath, resdir, name, l90, nbcont,
                                         cutn, threads, date, force, qc_only)
    assert skip == []
    assert skipf == []
    assert kept == {}
    expg = {'B2_A3_5.fasta-changeName.fna': ['ESCO.0417', ('test/data/genomes/'
                                             'B2_A3_5.fasta-changeName.fna'), 120529, 5, 4],
            'H299_H561.fasta': ['ESCO.0417', 'test/data/genomes/H299_H561.fasta', 13143, 3, 3]}
    assert allg == expg
    shutil.rmtree(resdir, ignore_errors=True)


def test_annote_allDefault():
    """
    Test that when we call the pipeline with all default parameters, all expected output files
    are created. Check the content of result files (LSTINFO, Genes, Proteins, Replicons),
    lstinfo, discarded and log files.
    """
    date = time.strftime("%m%y")
    fulldate = time.strftime("%Y-%m-%d")
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
    dbpath = os.path.join("test", "data", "genomes")
    respath = os.path.join("test", "data", "res_test_funcDefault")
    name = "GENO"
    cmd = "annote_pipeline.py {} -d {} -r {} -n {}".format(list_file, dbpath, respath, name)
    ret = subprocess.call(cmd.split())
    assert ret == 0
    # Get output files
    log_files = [os.path.join(respath, "annote-genomes-list_genomes-func-test-default" + ext)
                 for ext in [".log", ".log.err"]]
    lstfile = [os.path.join(respath, "LSTINFO-list_genomes-func-test-default.lst")]
    discfile = [os.path.join(respath, "discarded-list_genomes-func-test-default.lst")]
    QCfiles = [os.path.join(respath, qc + "list_genomes-func-test-default.png")
               for qc in ["QC_L90-", "QC_nb-contigs-"]]
    genomes = {"B2_A3_5.fasta-changeName.fna": "GENO.{}.00003".format(date),
               "H299_H561.fasta": "GENO.{}.00002".format(date),
               "A_H738.fasta": "GENO.{}.00001".format(date)}
    split5Nfiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna") for g in genomes]
    gembasefiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna")
                    for g in genomes]
    proklogfiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna-prokka.log")
                    for g in genomes]
    prokka_files = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna-prokkaRes",
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
    for f in (log_files + lstfile + discfile + QCfiles + split5Nfiles + gembasefiles +
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
        exp_path = os.path.join(exp_dir, exp_file)
        with open(exp_path, "r") as expf, open(f, "r") as outf:
            for line_exp, line_out in zip(expf, outf):
                assert line_exp == line_out
    # Check that err file is empty
    with open(log_files[1], "r") as errf:
        lines = errf.readlines()
        assert lines == []
    # Check that log file contains expected information (info level)
    with open(log_files[0], "r") as logf:
        infos = []
        for line in logf:
            assert line.startswith("[" + fulldate)
            assert line.count("::") == 2
            infos.append(line.strip().split("::")[-1].strip())
        assert "Reading genomes" in infos
        assert ("Cutting genomes at each stretch of at least 5 'N', and then, calculating "
                "genome size, number of contigs and L90.") in infos
        assert ("Annotating all genomes with prokka") in infos
        assert ("Formatting all genomes") in infos
    # Check lstinfo file content
    explst = os.path.join(exp_dir, "LSTINFO-list_genomes-func-test-default.lst")
    with open(lstfile[0], "r") as expf, open(explst, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    shutil.rmtree(respath)