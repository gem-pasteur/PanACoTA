#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for 'annotate' subcommand
"""

from PanACoTA.subcommands import annotate as annot
import test.test_unit.utilities_for_tests as tutil

import pytest
import os
import subprocess
import shutil
import time
import argparse
import matplotlib
import logging
import glob
matplotlib.use('AGG')


# LOGFILE_BASE = "test_main_from_parse"
# Define variables used by several tests
DBDIR = os.path.join("test", "data", "annotate")
GEN_PATH = os.path.join(DBDIR, "genomes")
EXP_DIR = os.path.join(DBDIR, 'exp_files')
TEST_DIR = os.path.join(DBDIR, 'test_files')
GENEPATH = os.path.join(DBDIR, "generated_by_func-tests")


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
    if not os.path.isdir(GENEPATH):
        print("setup")
        os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")


def test_main_from_parse():
    """
    Test that when a tmp folder is given by user, tmp files are saved in it,
    and prokka files too.
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-default.txt")
    tmpdir = os.path.join(GENEPATH, "tmp_funcGivenTmp")
    name = "ESCO"
    l90 = 100
    date = "0417"
    args = argparse.Namespace()
    args.list_file = list_file
    args.db_path = GEN_PATH
    args.res_path = GENEPATH
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
    args.verbose = False
    args.quiet = False
    args.from_info = False
    args.prodigal_only = False
    args.small = False
    args.annotdir = False
    args.argv = ["annotate", "test_annote.py", "test_main_from_parse"]
    args.prodigal_only = False
    annot.main_from_parse(args)
    # Check that only 2 fna files (over 3) were created in tmp files. The 3rd
    # one is annotated from the original seq, nothing to merge
    assert len(glob.glob(os.path.join(tmpdir, '*.fna'))) == 2
    # Check that tmp files exist in the right folder
    assert os.path.isfile(os.path.join(tmpdir, "A_H738.fasta-all.fna"))
    assert os.path.isfile(os.path.join(tmpdir, "H299_H561.fasta-all.fna"))
    # Test that prokka folder is in the right directory
    assert os.path.isdir(os.path.join(tmpdir, "A_H738.fasta-all.fna-prokkaRes"))
    assert os.path.isdir(os.path.join(tmpdir, "B2_A3_5.fasta-changeName.fna-prokkaRes"))
    assert os.path.isdir(os.path.join(tmpdir, "H299_H561.fasta-all.fna-prokkaRes"))


def test_main_novalid_genome(capsys):
    """
    Test that when, in the list file, all genomes are wrong (do not correspond to
    filenames in the given dbpath), it closes the program with an error message.
    """
    list_file = os.path.join(TEST_DIR, "list_no_genome.txt")
    name = "ESCO"
    date = "0417"
    with pytest.raises(SystemExit):
        annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date)
    _, err = capsys.readouterr()
    assert ("We did not find any genome listed in test/data/annotate/test_files/"
            "list_no_genome.txt "
            "in the folder test/data/annotate/genomes. Please check your list to give valid "
            "genome names.") in err


def test_main_given_tmp_verbose3(capsys):
    """
    Test that when a tmp folder is given by user, tmp files are saved in it,
    and prokka files too.
    + check that, with verbose=3, warning and details are written to stdout

    Giving 4 genomes in list_files
    - for 1 genome, toto.fst does not exist, and will not be in the concatenated file
    - 2 concatenated files
    - 4 files to annotate
    - 4 prokkaRes
    - 1 genome with problems: no CDS found
    - 3 genomes in result dirs
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-default.txt")
    tmpdir = os.path.join(GENEPATH, "tmp_funcGivenTmp")
    name = "ESCO"
    l90 = 10
    date = "0417"
    verbose = 3
    info_file = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-default.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, l90,
                      cutn=3, tmp_dir=tmpdir, verbose=verbose) == (info_file, 3)
    out, err = capsys.readouterr()
    # Check that warnings are written to stderr
    assert "WARNING" in err
    assert ("toto.fst genome file does not exist. Its file will be ignored when "
            "concatenating ['A_H738.fasta', 'genome1.fasta', 'toto.fst']") in err
    # Check that tmp files exist in the right folder
    # -> 2 fna files created (concatenations)
    # -> + 3 files created (split 5N)
    assert os.path.isfile(os.path.join(tmpdir, "A_H738.fasta-all.fna"))
    assert os.path.isfile(os.path.join(tmpdir, "H299_H561.fasta-all.fna"))
    assert len(glob.glob(os.path.join(tmpdir, '*.fna'))) == 6
    assert len(glob.glob(os.path.join(tmpdir, '*split3N.fna'))) == 4
    # Check that split contigs were renamed with unique ID at the begining of the header
    res_file = os.path.join(tmpdir, "A_H738.fasta-all.fna_prokka-split3N.fna")
    exp_file = os.path.join(EXP_DIR, "exp_A_H738.fasta-all.fna_prokka-split3N.fna")
    assert tutil.compare_order_content(exp_file, res_file)
    # Check that even for complete genome, contig was renamed with ID
    res_file = os.path.join(tmpdir, "complete_genome.fna_prokka-split3N.fna")
    exp_file = os.path.join(EXP_DIR, "exp_complete_genome.fna_prokka-split3N.fna")
    assert tutil.compare_order_content(exp_file, res_file)
    # Test that prokka folder is in the right directory
    # Only 1 genome annotated by  prokka (the 2 others do not have appropriate L90/nbcont)
    assert os.path.isdir(os.path.join(tmpdir, "A_H738.fasta-all.fna_prokka-split3N.fna-prokkaRes"))
    assert not os.path.isdir(os.path.join(tmpdir, "H299_H561.fasta-all.fna-prokkaRes"))


def test_main_all_discard_nbcont(capsys):
    """
    Test that when the genomes given in list file have high nbcontigs compared
    to given threshold, error message as there are no genome to annotate
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-default.txt")
    name = "ESCO"
    nbcont = 0
    cutn = 0
    date = "0417"
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, nbcont=nbcont,
                      cutn=cutn) == ("", 0)
    # check that there are the 2 concatenated genomes in tmppath.
    # The third genome is listfile is composed of only 1 file, so no need to concatenate, nor
    # to change the file as we do not use cutn
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files", "H299_H561.fasta-all.fna"))
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files", "A_H738.fasta-all.fna"))
    out, err = capsys.readouterr()
    assert 'No genome kept for annotation' in out


def test_main_qc():
    """
    Test that when only QC is run, it writes:
    - the list of all genomes with their characteristics
    - the list of genomes that would be discarded for annotation
    - the 2 png files
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-default.txt")
    name = "ESCO"
    cutn = 0
    threads = 1
    l90 = 1
    date = "0417"
    force = False
    qc_only = True
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, l90=l90,
                      cutn=cutn, qc_only=qc_only) == ("", 0)
    # Check files are here
    lstfile = os.path.join(GENEPATH, "ALL-GENOMES-info-list_genomes-func-test-default.lst")
    exp_lstfile = os.path.join(EXP_DIR, "exp_ALL-GENOMES-QC.lst")
    discardedfile = os.path.join(GENEPATH, "discarded-list_genomes-func-test-default.lst")
    exp_discarded = os.path.join(EXP_DIR, "exp_discarded_QC.lst")
    assert os.path.isfile(lstfile)
    assert os.path.isfile(discardedfile)
    assert os.path.isfile(os.path.join(GENEPATH,
                          "QC_L90-list_genomes-func-test-default.png"))
    assert os.path.isfile(os.path.join(GENEPATH,
                          "QC_nb-contigs-list_genomes-func-test-default.png"))
    # Check content of discarded genomes
    assert tutil.compare_file_content(lstfile, exp_lstfile)
    assert tutil.compare_file_content(discardedfile, exp_discarded)


def test_main_existresdirforce(capsys):
    """
    Test that, when the pipeline is run on an existing result directory, but force option is on,
    it removes the result folders and runs again.
    Result folders contain expected files, the ones put before are removed
    Giving 4 genomes
    - 4 are OK
    - trained on complete_genome_big.fna
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-force.txt")
    # Create output directory with a prt file in Proteins folder
    protdir = os.path.join(GENEPATH, "Proteins")
    os.makedirs(protdir)
    open(os.path.join(protdir, "toto.prt"), "w").close()
    assert os.path.isfile(os.path.join(protdir, "toto.prt"))
    name = "ESCO"
    date = "0417"
    l90 = 5
    cutn = 3
    info_file = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-force.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, force=True, l90=l90,
                      prodigal_only=True, cutn = cutn) == (info_file, 4)
    out, err = capsys.readouterr()

    # Check that tmp files exist in the right folder
    # -> 2 fna files created (concatenations)
    # -> + 4 files created (split 3N)
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files", "A_H738.fasta-all.fna"))
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files", "H299_H561.fasta-all.fna"))
    assert len(glob.glob(os.path.join(GENEPATH, "tmp_files", '*.fna'))) == 6
    assert len(glob.glob(os.path.join(GENEPATH, "tmp_files", '*split3N.fna'))) == 4
    # Check that split contigs were renamed with unique ID at the begining of the header
    res_file = os.path.join(GENEPATH, "tmp_files", "A_H738.fasta-all.fna_prodigal-split3N.fna")
    exp_file = os.path.join(EXP_DIR, "exp_A_H738.fasta-all.fna_prokka-split3N.fna")
    assert tutil.compare_order_content(exp_file, res_file)
    # Check that even for complete genome, contig was renamed with ID
    res_file = os.path.join(GENEPATH, "tmp_files", "complete_genome_big.fna_prodigal-split3N.fna")
    exp_file = os.path.join(EXP_DIR, "exp_complete_genome_big.fna-split3N.fna")
    assert tutil.compare_order_content(exp_file, res_file)
    # Check that it trained on expected genome, and training file is ok
    trn_file = os.path.join(GENEPATH, "tmp_files", "complete_genome_big.fna"
                                                   "_prodigal-split3N.fna.trn")
    exp_trn_file = os.path.join(EXP_DIR, "exp_complete_genome_big.fna.trn")
    assert tutil.compare_files_bin(trn_file, exp_trn_file)
    # Check that tmp files exist in the right folder (result/tmp_files)
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files",
                                       "B2_A3_5.fasta-changeName.fna_prodigal-split3N.fna"))
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files",
                                       "H299_H561.fasta-all.fna_prodigal-split3N.fna"))
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files",
                                       "H299_H561.fasta-all.fna"))
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files",
                                       "A_H738.fasta-all.fna_prodigal-split3N.fna"))
    assert os.path.isfile(os.path.join(GENEPATH, "tmp_files",
                                       "A_H738.fasta-all.fna"))
    # Test all result folders are empty (in particular Proteins) as no genome is annotated
    assert os.path.isdir(protdir)
    assert len(os.listdir(protdir)) == 4
    assert not os.path.isfile(os.path.join(protdir, "toto.prt"))
    assert os.path.isfile(os.path.join(protdir, "ESCO.0417.00001.prt"))
    assert os.path.isfile(os.path.join(protdir, "ESCO.1015.00002.prt"))
    assert os.path.isfile(os.path.join(protdir, "ESCO.1015.00003.prt"))
    assert os.path.isfile(os.path.join(protdir, "ESCO.1116.00004.prt"))


def test_run_exist_resdir(caplog):
    """
    Test that when the pipeline is called, with a given resdir which already contains
    results, the program ends, with an error message.
    """
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(GENEPATH, "Proteins"))
    open(os.path.join(GENEPATH, "Proteins", "toto.prt"), "w").close()
    with pytest.raises(SystemExit):
        annot.main("cmd", "list_file.lst", "path/db", GENEPATH, "toto", "0123")
    assert ("ERROR: Your output directory already has .prt files in the "
            "Proteins folder. Provide another result directory, or remove the "
            "files in this one.") in caplog.text
    # File was not removed
    assert os.path.isfile(os.path.join(GENEPATH, "Proteins", "toto.prt"))


def test_main_onexistingprokkadir(capsys):
    """
    Test that, when the pipeline is run with a given prokka dir, where prokka results already
    exist, and are ok, all runs well, no re-annotation, just format


    main function arguments:
    cmd, list_file, db_path, res_dir, name, date, l90=100, nbcont=999, cutn=5,
    threads=1, force=False, qc_only=False, from_info=None, tmp_dir=None, res_annot_dir=None,
    verbose=0, quiet=False, prodigal_only=False):

    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")
    name = "ESCO"
    date = "0417"
    lstout = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-exist_dir.lst")
    lstexp = os.path.join(EXP_DIR, "exp_LSTINFO-func-annot_exists-prokkadir.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, cutn=0,
                      res_annot_dir=EXP_DIR, verbose=3) == (lstout, 2)
    out, err = capsys.readouterr()
    # Check that tmp files folder is empty (prokka res are somewhere else)
    assert len(os.listdir(os.path.join(GENEPATH, "tmp_files"))) == 0
    # Test that result files are in result dir
    assert os.path.isfile(lstout)
    assert tutil.compare_order_content(lstout, lstexp)
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Prokka results folder "
            "test/data/annotate/exp_files/B2_A3_5.fasta-changeName.fna-prokkaRes "
            "already exists") in " ".join(log_content)
    assert ("Prokka did not run again, formatting step used already generated results of Prokka "
            "in test/data/annotate/exp_files/H299_H561.fasta-prokkaRes. "
            "If you want to re-run prokka, first remove this result folder, or use '-F' or "
            "'--force' option if you want to rerun prokka "
            "for all genomes.") in ' '.join(log_content)


def test_main_onexistingprodigaldir(capsys):
    """
    Test that, when the pipeline is run with a given prodigal dir, where prodigal results already
    exist, and are ok, all runs well, no re-annotation, just format

    - trains
    - no re-annotation
    - format

    main function arguments:
    cmd, list_file, db_path, res_dir, name, date, l90=100, nbcont=999, cutn=5,
         threads=1, force=False, qc_only=False, from_info=None, tmp_dir=None, res_annot_dir=None,
         verbose=0, quiet=False, prodigal_only=False

    """
    # FOLDER with all results
    # Create result folder, with existing prodigal folders (which are OK)
    res_folder = os.path.join(GENEPATH, "results-prodigal")
    os.makedirs(res_folder)
    # copy prodigalRes folders
    B2_A3_5_folder = os.path.join(EXP_DIR, "B2_A3_5.fasta-changeName.fna-prodigalRes")
    H299_folder = os.path.join(EXP_DIR, "H299_H561.fasta-prodigalRes")
    res_B2_A3_5_folder = os.path.join(res_folder, "B2_A3_5.fasta-changeName.fna-prodigalRes")
    res_H299_folder = os.path.join(res_folder, "H299_H561.fasta-prodigalRes")
    shutil.copytree(B2_A3_5_folder, res_B2_A3_5_folder)
    shutil.copytree(H299_folder, res_H299_folder)

    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")
    name = "ESCO"
    date = "0417"
    lstout = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-exist_dir.lst")
    lstexp = os.path.join(EXP_DIR, "exp_LSTINFO-func-annot_exists-prokkadir.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, cutn=0,
                      res_annot_dir=res_folder, verbose=3, prodigal_only=True) == (lstout, 2)
    out, err = capsys.readouterr()
    # Check that tmp files folder is empty (prokka res are somewhere else)
    assert len(os.listdir(os.path.join(GENEPATH, "tmp_files"))) == 0
    # Test that result files are in result dir
    assert os.path.isfile(lstout)
    assert tutil.compare_order_content(lstout, lstexp)
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Prodigal will train using "
            "test/data/annotate/genomes/H299_H561.fasta") in " ".join(log_content)
    assert ("prodigal command: prodigal -i test/data/annotate/genomes/H299_H561.fasta "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/H299_H561.fasta.trn")
    assert("Error while trying to train prodigal on H299_H561.fasta. "
           "See test/data/annotate/generated_by_func-tests/results-prodigal/"
           "H299_H561.fasta.trn-prodigal-train.log.err") in " ".join(log_content)
    assert ("Prodigal results folder test/data/annotate/generated_by_func-tests/"
            "results-prodigal/H299_H561.fasta-prodigalRes "
            "already exists") in " ".join(log_content)
    assert ("Prodigal did not run again. Formatting step will use already generated results of "
            "Prodigal in test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes. "
            "If you want to re-run Prodigal, first remove this result folder, or use '-F' or "
            "'--force' option.") in ' '.join(log_content)
    assert "Formatting all genomes" in " ".join(log_content)
    assert "Annotation step done" in " ".join(log_content)


def test_main_onexistingprodigaldir_train_exists(capsys):
    """
    Test that, when the pipeline is run with a given prodigal dir, where prodigal results already
    exist, and are ok, all runs well, no re-annotation, just format

    - no train
    - no reannote
    - format

    2 genomes in list file: B2_A3_5.fasta-changeName.fna and H299_H561.fasta
    """
    # FOLDER with all results
    # Create result folder, with existing prodigal folders (which are OK)
    res_folder = os.path.join(GENEPATH, "results-prodigal")
    os.makedirs(res_folder)
    # copy prodigalRes folders
    B2_A3_5_folder = os.path.join(EXP_DIR, "B2_A3_5.fasta-changeName.fna-prodigalRes")
    H299_folder = os.path.join(EXP_DIR, "H299_H561.fasta-prodigalRes")
    res_B2_A3_5_folder = os.path.join(res_folder, "B2_A3_5.fasta-changeName.fna-prodigalRes")
    res_H299_folder = os.path.join(res_folder, "H299_H561.fasta-prodigalRes")
    shutil.copytree(B2_A3_5_folder, res_B2_A3_5_folder)
    shutil.copytree(H299_folder, res_H299_folder)
    # Add a training file in result folder
    trn_file = os.path.join(res_folder, "H299_H561.fasta.trn")
    open(trn_file, "w").close()

    # Function arguments
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")
    name = "ESCO"
    date = "0417"
    lstout = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-exist_dir.lst")
    lstexp = os.path.join(EXP_DIR, "exp_LSTINFO-func-annot_exists-prokkadir.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, cutn=0,
                      res_annot_dir=res_folder, verbose=3, prodigal_only=True) == (lstout, 2)
    out, err = capsys.readouterr()
    # Check that tmp files folder is empty (prokka res are somewhere else)
    assert len(os.listdir(os.path.join(GENEPATH, "tmp_files"))) == 0
    # Test that result files are in result dir
    assert os.path.isfile(lstout)
    assert tutil.compare_order_content(lstout, lstexp)
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("A training file already exists (test/data/annotate/generated_by_func-tests/"
            "results-prodigal/H299_H561.fasta.trn). It will be used to annotate "
            "all genomes.") in " ".join(log_content)
    assert ("Prodigal results folder test/data/annotate/generated_by_func-tests/"
            "results-prodigal/H299_H561.fasta-prodigalRes "
            "already exists") in " ".join(log_content)
    assert ("Prodigal results folder test/data/annotate/generated_by_func-tests/"
            "results-prodigal/B2_A3_5.fasta-changeName.fna-prodigalRes "
            "already exists") in " ".join(log_content)
    assert ("Prodigal did not run again. Formatting step will use already generated results of "
            "Prodigal in test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes. "
            "If you want to re-run Prodigal, first remove this result folder, or use '-F' or "
            "'--force' option.") in ' '.join(log_content)
    assert "Formatting all genomes" in " ".join(log_content)
    assert "Annotation step done" in " ".join(log_content)


def test_main_prodigal_train_empty(capsys):
    """
    Test that, when the pipeline is run with a given prodigal dir, where prodigal results do
    not exist, and the given trn file is empty
    -> error, with prodigal command

    - no train
    - try reannote but fails -> exits

    2 genomes in list file: B2_A3_5.fasta-changeName.fna and H299_H561.fasta
    """
    # FOLDER with all results
    # Create result folder, with empty trn file
    res_folder = os.path.join(GENEPATH, "results-prodigal")
    os.makedirs(res_folder)
    trn_file = os.path.join(res_folder, "H299_H561.fasta.trn")
    open(trn_file, "w").close()

    # Function arguments
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")
    name = "ESCO"
    date = "0417"
    lstout = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-exist_dir.lst")
    lstexp = os.path.join(EXP_DIR, "exp_LSTINFO-func-annot_exists-prokkadir.lst")
    with pytest.raises(SystemExit):
        annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, cutn=0,
                   res_annot_dir=res_folder, verbose=3, prodigal_only=True)
    out, err = capsys.readouterr()
    # Test that result files are in result dir
    assert os.path.isfile(lstout)
    assert tutil.compare_order_content(lstout, lstexp)
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    # Check logs
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("A training file already exists (test/data/annotate/generated_by_func-tests/"
            "results-prodigal/H299_H561.fasta.trn). It will be used to annotate "
            "all genomes.") in " ".join(log_content)
    assert ("Start annotating ESCO.1116.00002 (from test/data/annotate/genomes/"
            "B2_A3_5.fasta-changeName.fna sequence) with Prodigal") in " ".join(log_content)
    assert ("Start annotating ESCO.1015.00001 (from test/data/annotate/genomes/"
            "H299_H561.fasta sequence) with Prodigal") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/B2_A3_5.fasta-changeName.fna "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "B2_A3_5.fasta-changeName.fna-prodigalRes/ESCO.1116.00002.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "B2_A3_5.fasta-changeName.fna-prodigalRes/ESCO.1116.00002.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "B2_A3_5.fasta-changeName.fna-prodigalRes/ESCO.1116.00002.gff "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/H299_H561.fasta.trn "
            "-q") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/H299_H561.fasta "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/ESCO.1015.00001.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/ESCO.1015.00001.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/ESCO.1015.00001.gff "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/H299_H561.fasta.trn "
            "-q") in " ".join(log_content)
    assert ("Error while trying to run prodigal. See test/data/annotate/generated_by_func-tests/"
            "results-prodigal/H299_H561.fasta-prodigal.log.err.") in " ".join(log_content)
    assert ("Error while trying to run prodigal. See test/data/annotate/generated_by_func-tests/"
            "results-prodigal/"
            "B2_A3_5.fasta-changeName.fna-prodigal.log.err") in " ".join(log_content)
    assert ("No genome was correctly annotated, no need to format them.") in " ".join(log_content)


def test_main_prodigal_train_ok(capsys):
    """
    Test that, when the pipeline is run with a given prodigal dir, where prodigal train exists and is ok:

    - no train
    - reannotate
    - format

    2 genomes in list file: B2_A3_5.fasta-changeName.fna and H299_H561.fasta
    """
    # FOLDER with all results
    # Create result folder, with existing prodigal folders (which are OK)
    res_folder = os.path.join(GENEPATH, "results-prodigal")
    os.makedirs(res_folder)
    # Add a valid training file in result folder
    orig_trn_file = os.path.join(EXP_DIR, "exp_complete_genome_big.fna.trn")
    trn_file = os.path.join(res_folder, "H299_H561.fasta.trn")
    shutil.copyfile(orig_trn_file, trn_file)

    # Function arguments
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")
    name = "ESCO"
    date = "0417"
    lstout = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-exist_dir.lst")
    lstexp = os.path.join(EXP_DIR, "exp_LSTINFO-func-annot_exists-prokkadir.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, cutn=0,
                      res_annot_dir=res_folder, verbose=3, prodigal_only=True) == (lstout, 2)
    out, err = capsys.readouterr()
    # Check that tmp files folder is empty (prokka res are somewhere else)
    assert len(os.listdir(os.path.join(GENEPATH, "tmp_files"))) == 0
    # Test that result files are in result dir
    assert os.path.isfile(lstout)
    assert tutil.compare_order_content(lstout, lstexp)
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("A training file already exists (test/data/annotate/generated_by_func-tests/"
            "results-prodigal/H299_H561.fasta.trn). It will be used to annotate "
            "all genomes.") in " ".join(log_content)
    assert ("Start annotating ESCO.1116.00002 (from test/data/annotate/genomes/"
            "B2_A3_5.fasta-changeName.fna sequence) with Prodigal") in " ".join(log_content)
    assert ("Start annotating ESCO.1015.00001 (from test/data/annotate/genomes/"
            "H299_H561.fasta sequence) with Prodigal") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/B2_A3_5.fasta-changeName.fna "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "B2_A3_5.fasta-changeName.fna-prodigalRes/ESCO.1116.00002.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "B2_A3_5.fasta-changeName.fna-prodigalRes/ESCO.1116.00002.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "B2_A3_5.fasta-changeName.fna-prodigalRes/ESCO.1116.00002.gff "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/H299_H561.fasta.trn "
            "-q") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/H299_H561.fasta "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/ESCO.1015.00001.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/ESCO.1015.00001.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/ESCO.1015.00001.gff "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/H299_H561.fasta.trn "
            "-q") in " ".join(log_content)
    assert ("End annotating ESCO.1015.00001 (from test/data/annotate/genomes/"
            "H299_H561.fasta)") in " ".join(log_content)
    assert "Formatting all genomes" in " ".join(log_content)
    assert "Annotation step done" in " ".join(log_content)


def test_main_prodigal_ok(capsys):
    """
    Test that, when the pipeline is run with a given empty prodigal dir, it does everything

    - train
    - reannotate
    - format

    2 genomes in list file: B2_A3_5.fasta-changeName.fna and H299_H561.fasta
    """
    # FOLDER with all results
    # Create result folder, with existing prodigal folders (which are OK)
    res_folder = os.path.join(GENEPATH, "results-prodigal")
    os.makedirs(res_folder)

    # Function arguments
    list_file = os.path.join(GENEPATH, "list_genomes.txt")
    with open(list_file, "w") as lf:
        lf.write("A_H738.fasta \n")
        lf.write("H299_H561.fasta::TOTO")
    name = "ESCO"
    date = "0417"
    lstout = os.path.join(GENEPATH, "LSTINFO-list_genomes.lst")
    lstexp = os.path.join(EXP_DIR, "exp_LSTINFO-func-annot_exists-prokkadir.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, cutn=0,
                      res_annot_dir=res_folder, verbose=3, prodigal_only=True) == (lstout, 2)
    out, err = capsys.readouterr()
    # Test that result files are in result dir
    assert os.path.isfile(lstout)
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Prodigal will train using "
            "test/data/annotate/genomes/A_H738.fasta") in " ".join(log_content)
    assert ("prodigal command: prodigal -i test/data/annotate/genomes/A_H738.fasta "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/"
            "A_H738.fasta.trn") in " ".join(log_content)
    assert ("Start annotating ESCO.0417.00001 (from test/data/annotate/genomes/"
            "A_H738.fasta sequence) with Prodigal") in " ".join(log_content)
    assert ("Start annotating TOTO.0417.00001 (from test/data/annotate/genomes/"
            "H299_H561.fasta sequence) with Prodigal") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/A_H738.fasta "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "A_H738.fasta-prodigalRes/ESCO.0417.00001.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "A_H738.fasta-prodigalRes/ESCO.0417.00001.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "A_H738.fasta-prodigalRes/ESCO.0417.00001.gff "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/A_H738.fasta.trn "
            "-q") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/H299_H561.fasta "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/TOTO.0417.00001.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/TOTO.0417.00001.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/TOTO.0417.00001.gff "
            "-t test/data/annotate/generated_by_func-tests/results-prodigal/A_H738.fasta.trn "
            "-q") in " ".join(log_content)
    assert ("End annotating ESCO.0417.00001 (from test/data/annotate/genomes/"
            "A_H738.fasta)") in " ".join(log_content)
    assert ("End annotating TOTO.0417.00001 (from test/data/annotate/genomes/"
            "H299_H561.fasta)") in " ".join(log_content)
    assert "Formatting all genomes" in " ".join(log_content)
    assert "Annotation step done" in " ".join(log_content)


def test_main_prodigal_small_ok(capsys):
    """
    Test that, when the pipeline is run with a given prodigal dir, and --small option, it does:

    - no train
    - reannotate
    - format

    2 genomes in list file: B2_A3_5.fasta-changeName.fna and H299_H561.fasta
    """
    # FOLDER with all results
    # Create result folder, with existing prodigal folders (which are OK)
    res_folder = os.path.join(GENEPATH, "results-prodigal")
    os.makedirs(res_folder)

    # Function arguments
    list_file = os.path.join(GENEPATH, "list_genomes_prodigal_small.txt")
    with open(list_file, "w") as lf:
        lf.write("A_H738.fasta \n")
        lf.write("H299_H561.fasta::TOTO")
    name = "ESCO"
    date = "0417"
    lstout = os.path.join(GENEPATH, "LSTINFO-list_genomes_prodigal_small.lst")
    # lstexp = os.path.join(EXP_DIR, "exp_LSTINFO-func-annot_exists-prokkadir.lst")
    assert annot.main("cmd", list_file, GEN_PATH, GENEPATH, name, date, cutn=0,
                      res_annot_dir=res_folder, verbose=3, prodigal_only=True, small=True) == (lstout, 2)
    out, err = capsys.readouterr()
    # Test that result files are in result dir
    assert os.path.isfile(lstout)
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes_prodigal_small.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Start annotating ESCO.0417.00001 (from test/data/annotate/genomes/"
            "A_H738.fasta sequence) with Prodigal") in " ".join(log_content)
    assert ("Start annotating TOTO.0417.00001 (from test/data/annotate/genomes/"
            "H299_H561.fasta sequence) with Prodigal") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/A_H738.fasta "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "A_H738.fasta-prodigalRes/ESCO.0417.00001.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "A_H738.fasta-prodigalRes/ESCO.0417.00001.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "A_H738.fasta-prodigalRes/ESCO.0417.00001.gff "
            "-p meta -q") in " ".join(log_content)
    assert ("Prodigal command: "
            "prodigal -i test/data/annotate/genomes/H299_H561.fasta "
            "-d test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/TOTO.0417.00001.ffn "
            "-a test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/TOTO.0417.00001.faa "
            "-f gff -o test/data/annotate/generated_by_func-tests/results-prodigal/"
            "H299_H561.fasta-prodigalRes/TOTO.0417.00001.gff "
            "-p meta -q") in " ".join(log_content)
    assert ("End annotating ESCO.0417.00001 (from test/data/annotate/genomes/"
            "A_H738.fasta)") in " ".join(log_content)
    assert ("End annotating TOTO.0417.00001 (from test/data/annotate/genomes/"
            "H299_H561.fasta)") in " ".join(log_content)
    assert "Formatting all genomes" in " ".join(log_content)
    assert "Annotation step done" in " ".join(log_content)


def test_main_existing_prokkadir_errorannot():
    """
    Test that, when the pipeline is run with a given prokka dir, where prokka results have
    problems (no tbl file and no gff file), it returns an error message and the genome with
    problems is in skipped.
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")

    # Create directory with all files needed for the test
    genome_path_used = os.path.join(GENEPATH, "genomes")
    os.makedirs(genome_path_used)
    ori_genome1 = os.path.join(GEN_PATH, "B2_A3_5.fasta-changeName.fna")
    ori_prok_g1 = os.path.join(EXP_DIR, "B2_A3_5.fasta-changeName.fna-prokkaRes")
    used_genome1 = os.path.join(genome_path_used, "B2_A3_5.fasta-changeName.fna")
    used_prok_g1 = used_genome1 + "-prokkaRes"
    # Copy original fasta file to genepath/genomes
    shutil.copyfile(ori_genome1, used_genome1)
    # Copy prokka results to genepath/genomes/gname-prokkaRes
    shutil.copytree(ori_prok_g1, used_prok_g1)
    # Same think for 2nd genome
    ori_genome2 = os.path.join(GEN_PATH, "H299_H561.fasta")
    ori_prok_g2 = os.path.join(EXP_DIR, "H299_H561.fasta-prokkaRes")
    used_genome2 = os.path.join(genome_path_used, "H299_H561.fasta")
    used_prok_g2 = used_genome2 + "-prokkaRes"
    shutil.copyfile(ori_genome2, used_genome2)
    shutil.copytree(ori_prok_g2, used_prok_g2)

    # Remove tbl file for genome1, so that check_prokka returns an error
    os.remove(os.path.join(used_prok_g1, "ESCO.1116.00002.tbl"))

    # # Run annotation
    name = "ESCO"
    date = "0417"
    info_file = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-exist_dir.lst")
    assert annot.main("cmd", list_file, genome_path_used, GENEPATH, name, date, cutn=0,
                      res_annot_dir=genome_path_used) == (info_file, 1)

    # Check that only 1 genome was formated (the other one had problems with prokka)
    prot_dir = os.path.join(GENEPATH, "Proteins")
    assert len(os.listdir(prot_dir)) == 1
    rep_dir = os.path.join(GENEPATH, "Replicons")
    assert len(os.listdir(rep_dir)) == 1

    # Check that the genome formated was not re-annotated
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    with open(logfile, "r") as lf:
        log_content = lf.readlines()
    assert ("Prokka results folder "
            "test/data/annotate/generated_by_func-tests/genomes/H299_H561.fasta-prokkaRes "
            "already exists") in " ".join(log_content)
    # Check that genome not formated because error in prokka res
    assert ("ESCO.1116.00002 B2_A3_5.fasta-changeName.fna: no .tbl file") in " ".join(log_content)
    assert ("Problems in the files contained in your already existing output dir "
            "(test/data/annotate/generated_by_func-tests/genomes/"
            "B2_A3_5.fasta-changeName.fna-prokkaRes). Please check it, "
            "or remove it to re-annotate.") in ' '.join(log_content)


def test_main_existing_prodigaldir_errorannot(capsys):
    """
    Test that, when the pipeline is run with a given prodigal dir, where prodigal results have
    problems for both genomes, it returns an error message and the genomes with
    problems are in skipped. Error message with no genome to format.
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")

    # Create directory with all files needed for the test
    genome_path_used = os.path.join(GENEPATH, "genomes")
    os.makedirs(genome_path_used)
    ori_genome1 = os.path.join(GEN_PATH, "B2_A3_5.fasta-changeName.fna")
    ori_prok_g1 = os.path.join(EXP_DIR, "B2_A3_5.fasta-changeName.fna-prodigalRes")
    used_genome1 = os.path.join(genome_path_used, "B2_A3_5.fasta-changeName.fna")
    used_prok_g1 = used_genome1 + "-prodigalRes"
    # Copy original fasta file to genepath/genomes
    shutil.copyfile(ori_genome1, used_genome1)
    # Copy prokka results to genepath/genomes/gname-prokkaRes
    shutil.copytree(ori_prok_g1, used_prok_g1)
    # Same thing for 2nd genome
    ori_genome2 = os.path.join(GEN_PATH, "H299_H561.fasta")
    ori_prok_g2 = os.path.join(EXP_DIR, "H299_H561.fasta-prodigalRes")
    used_genome2 = os.path.join(genome_path_used, "H299_H561.fasta")
    used_prok_g2 = used_genome2 + "-prodigalRes"
    shutil.copyfile(ori_genome2, used_genome2)
    shutil.copytree(ori_prok_g2, used_prok_g2)

    # Remove faa file for genome1, so that check_prodigal returns an error
    os.remove(os.path.join(used_prok_g1, "ESCO.1116.00002.faa"))
    # Remove gff file for genome1, so that check_prodigal returns an error
    os.remove(os.path.join(used_prok_g2, "ESCO.1015.00001.gff"))

    # Run annotation
    name = "ESCO"
    date = "0417"
    with pytest.raises(SystemExit):
        annot.main("cmd", list_file, genome_path_used, GENEPATH, name, date, cutn=0,
                   res_annot_dir=genome_path_used, prodigal_only=True, verbose=15)

    # Check that Replicons & co folders are not created
    prot_dir = os.path.join(GENEPATH, "Proteins")
    assert not os.path.isdir(prot_dir)
    rep_dir = os.path.join(GENEPATH, "Replicons")
    assert not os.path.isdir(rep_dir)

    # Check that not formatted because exists + error
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Error: No genome was correctly annotated, "
            "no need to format them") in ' '.join(log_content)
    assert ("Prodigal results folder test/data/annotate/generated_by_func-tests/genomes/"
            "H299_H561.fasta-prodigalRes already exists.") in ' '.join(log_content)
    assert ("ESCO.1116.00002 B2_A3_5.fasta-changeName.fna: "
            "no or several .faa file(s)") in ' '.join(log_content)
    assert ("ESCO.1015.00001 H299_H561.fasta: "
            "no or several .gff file(s)") in ' '.join(log_content)


def test_main_existing_prokkadir_errorformat():
    """
    Test that, when the pipeline is run with a given prokka dir, where prokka results are ok
    (all expected files), but have problems inside (wrong header format), it returns an error
    message and the genome with problems is in skipped_format.
    """
    list_file = os.path.join(TEST_DIR, "list_genomes-func-test-exist_dir.txt")

    # Create directory with all files needed for the test
    genome_path_used = os.path.join(GENEPATH, "genomes")
    os.makedirs(genome_path_used)
    ori_genome1 = os.path.join(GEN_PATH, "B2_A3_5.fasta-changeName.fna")
    # orig prokka dir has a tbl file with wrong format
    ori_prok_g1 = os.path.join(EXP_DIR, "B2_A3_5.fasta-changeName.fna-short-contig.fna-prokkaRes")
    used_genome1 = os.path.join(genome_path_used, "B2_A3_5.fasta-changeName.fna")
    used_prok_g1 = used_genome1 + "-prokkaRes"
    # Copy original fasta file to genepath/genomes
    shutil.copyfile(ori_genome1, used_genome1)
    # Copy prokka results to genepath/genomes/gname-prokkaRes
    shutil.copytree(ori_prok_g1, used_prok_g1)
    # and add .fna file to prokka-dir
    used_fna1 = os.path.join(used_prok_g1, "B2_A3_5.fasta-changeName.fna")
    shutil.copyfile(ori_genome1, used_fna1)

    # Same thing for 2nd genome
    ori_genome2 = os.path.join(GEN_PATH, "H299_H561.fasta")
    ori_prok_g2 = os.path.join(EXP_DIR, "H299_H561.fasta-prokkaRes")
    used_genome2 = os.path.join(genome_path_used, "H299_H561.fasta")
    used_prok_g2 = used_genome2 + "-prokkaRes"
    # Copy original fasta file to tmp resdir
    shutil.copyfile(ori_genome2, used_genome2)
    # Copy folder with prokka result files to genepath result path
    shutil.copytree(ori_prok_g2, used_prok_g2)
    # and add .fna file to prokka-dir
    used_fna = os.path.join(used_prok_g2, "H299_H561.fasta")
    shutil.copyfile(ori_genome2, used_fna)

    # Run annotation
    name = "ESCO"
    date = "0417"
    info_file = os.path.join(GENEPATH, "LSTINFO-list_genomes-func-test-exist_dir.lst")
    assert annot.main("cmd", list_file, genome_path_used, GENEPATH, name, date, cutn=0,
               res_annot_dir=genome_path_used) == (info_file, 1)

    # Check that genome 1 is not formatted, while no error with prokka
    logfile = os.path.join(GENEPATH,
                           "PanACoTA-annotate_list_genomes-func-test-exist_dir.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Prokka results folder test/data/annotate/generated_by_func-tests/genomes/"
            "B2_A3_5.fasta-changeName.fna-prokkaRes already exists.") in ' '.join(log_content)
    assert ("Prokka did not run again, formatting step used already generated "
            "results of Prokka in test/data/annotate/generated_by_func-tests/genomes/"
            "B2_A3_5.fasta-changeName.fna-prokkaRes.") in ' '.join(log_content)
    # Error while trying to format:
    assert ("'changesHead.0417.00010.0005' found in "
            "test/data/annotate/generated_by_func-tests/genomes/"
            "B2_A3_5.fasta-changeName.fna-prokkaRes/test.0417.00002.tbl "
            "does not exist in test/data/annotate/generated_by_func-tests/genomes/"
            "B2_A3_5.fasta-changeName.fna") in ' '.join(log_content)
    assert ("Problems while generating LSTINFO file for ESCO.1116.00002") in ' '.join(log_content)


def test_main_frominfo(capsys):
    """
    test that it runs well when giving an info file instead of list file + db etc.
    It does not re-calculate L90 and nbcont
    """
    listfile = None
    dbpath = None
    name = "TOTO"
    date = "1205"
    infofile = os.path.join(TEST_DIR, "lstinfo.lst")
    out_infofile = os.path.join(GENEPATH, "LSTINFO-lstinfo.lst")
    assert annot.main("cmd", listfile, dbpath, GENEPATH, name, date, from_info=infofile,
                      prodigal_only=True) == (out_infofile, 3)
    out, err = capsys.readouterr()
    # Check logs
    assert ("Generating distribution of L90 and #contigs graphs.") in out

    # Check output files present
    protdir = os.path.join(GENEPATH, "Proteins")
    assert len(os.listdir(protdir)) == 3
    gffdir = os.path.join(GENEPATH, "gff3")
    assert len(os.listdir(gffdir)) == 3
    lstdir = os.path.join(GENEPATH, "LSTINFO")
    assert len(os.listdir(lstdir)) == 3

    # Check genomes are renamed as expected, and with expected L90/nbcont values
    exp_lstinfo = os.path.join(EXP_DIR, "exp_LSTINFO-test-main-frominfo.lst")
    res_lstinfo = os.path.join(GENEPATH, "LSTINFO-lstinfo.lst")
    assert tutil.compare_order_content(exp_lstinfo, res_lstinfo)


def test_main_novalid_genome_frominfo(capsys):
    """
    Test that when, in the list file, all genomes are wrong (do not correspond to
    filenames in the given dbpath), it closes the program with an error message.
    """
    listfile = None
    dbpath = None
    name = "TOTO"
    date = "1205"
    infofile = os.path.join(TEST_DIR, "lstinfo-no-genome.lst")
    with pytest.raises(SystemExit):
        annot.main("cmd", listfile, dbpath, GENEPATH, name, date, from_info=infofile,
                   prodigal_only=True)
    out, err = capsys.readouterr()
    # Check logs
    assert ("No genome listed in test/data/annotate/test_files/lstinfo-no-genome.lst "
            "was found.") in err

    # Check output folders not created
    protdir = os.path.join(GENEPATH, "Proteins")
    assert not os.path.isdir(protdir)
    gffdir = os.path.join(GENEPATH, "gff3")
    assert not os.path.isdir(gffdir)
    lstdir = os.path.join(GENEPATH, "LSTINFO")
    assert not os.path.isdir(lstdir)

    # Check tmp_files is empty
    tmp = os.path.join(GENEPATH, "tmp_files")
    assert len(os.listdir(tmp)) == 0
