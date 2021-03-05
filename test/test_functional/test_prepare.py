#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for PanACoTA 'prepare' subcommand
"""

from PanACoTA.subcommands import prepare
import test.test_unit.utilities_for_tests as tutil

import pytest
import os
import subprocess
import shutil
import time
import argparse
import logging
import glob


# LOGFILE_BASE = "test_main_from_parse"
# Define variables used by several tests
DBDIR = os.path.join("test", "data", "prepare")
GEN_PATH = os.path.join(DBDIR, "genomes")
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
    Run
    """
    args = argparse.Namespace()
    args.argv = ["prepare", "test_func_prepare"]
    args.ncbi_species_name = "Acetobacter orleanensis"
    args.ncbi_species_taxid = "104099"
    args.ncbi_taxid = ""
    args.ncbi_section = "refseq"
    args.outdir = GENEPATH
    args.tmp_dir = ""
    args.parallel = 1
    args.norefseq = False
    args.db_dir = ""
    args.only_mash = False
    args.info_file = ""
    args.l90 = 100
    args.nbcont = 999
    args.cutn = 0
    args.min_dist = 1e-4
    args.max_dist = 0.06
    args.verbose = 0
    args.quiet = False
    args.levels = ""

    prepare.main_from_parse(args)

    # Check output files
    summary =  os.path.join(GENEPATH, "assembly_summary-Acetobacter_orleanensis.txt")
    assert os.path.isfile(summary)
    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(GENEPATH, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) >= 4
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty as we do not split
    tmp_folder = os.listdir(os.path.join(GENEPATH, "tmp_files"))
    assert len(tmp_folder) == 0
    # Check Database_init folder created, with at list 4 ".fna" genomes
    fna_files = glob.glob(os.path.join(GENEPATH, "Database_init", "*.fna"))
    assert len(fna_files) >= 4


def test_main_from_parse_longspeciesname():
    """
    Run
    """
    args = argparse.Namespace()
    args.argv = ["prepare", "test_func_prepare"]
    args.ncbi_species_name = "Salmonella enterica subsp. enterica serovar Paratyphi C"
    args.ncbi_species_taxid = ""
    args.ncbi_taxid = ""
    args.ncbi_section = "refseq"
    args.outdir = GENEPATH
    args.tmp_dir = ""
    args.parallel = 1
    args.norefseq = False
    args.db_dir = ""
    args.only_mash = False
    args.info_file = ""
    args.l90 = 100
    args.nbcont = 999
    args.cutn = 0
    args.min_dist = 1e-4
    args.max_dist = 0.06
    args.verbose = 0
    args.quiet = False
    args.levels = ""

    prepare.main_from_parse(args)

    # Check output files
    summary =  os.path.join(GENEPATH, "assembly_summary-Salmonella_enterica_subsp._enterica_serovar_Paratyphi_C.txt")
    assert os.path.isfile(summary)
    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(GENEPATH, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) >= 1
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty as we do not split
    tmp_folder = os.listdir(os.path.join(GENEPATH, "tmp_files"))
    assert len(tmp_folder) == 0
    # Check Database_init folder created, with at list 4 ".fna" genomes
    fna_files = glob.glob(os.path.join(GENEPATH, "Database_init", "*.fna"))
    assert len(fna_files) >= 1


def test_main_not_only_mash_infoexists():
    """
    We run without option only_mash, but still provide a lstinfo file
    -> will change its name to .back to save it when the new file will be created
    """
    NCBI_species_name = ""
    NCBI_species_taxid = "104099"
    NCBI_taxid = ""
    NCBI_section = "refseq"
    levels = ""
    outdir = GENEPATH
    tmp_dir = os.path.join(outdir, "temporary_directory")
    threads = 1
    norefseq = False
    db_dir = ""
    only_mash = False
    info_file = os.path.join(outdir, "LSTINFO-existing.lst")
    open(info_file, "w").close()  #create empty info file, to check it is renamed
    l90 = 100
    nbcont = 999
    cutn = 5
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 2
    quiet = False
    out_info_file = os.path.join(outdir, "LSTINFO-104099-filtered-0.0001_0.06.txt")
    assert prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, NCBI_section, outdir, tmp_dir,
                        threads, norefseq, db_dir, only_mash, info_file, l90, nbcont,
                        cutn, min_dist, max_dist, verbose, quiet) == out_info_file

    # Check output files
    summary =  os.path.join(GENEPATH, "assembly_summary-104099.txt")
    assert os.path.isfile(summary)
    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(GENEPATH, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) >= 4
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty as we do not split
    tmp_files = glob.glob(os.path.join(tmp_dir, "*.fna_prepare-split5N.fna"))
    assert len(tmp_files) >= 4
    # Check Database_init folder created, with at list 4 ".fna" genomes
    fna_files = glob.glob(os.path.join(GENEPATH, "Database_init", "*.fna"))
    assert len(fna_files) >= 4
    # Check that LSTINFO file existing was renamed and still empty
    # And new LSTINFO file created
    assert os.path.isfile(info_file + ".back")
    assert os.stat(info_file + ".back").st_size == 0


def test_main_wrong_taxid(capsys):
    """
    We run without option only_mash, but still provide a lstinfo file
    -> will change its name to .back to save it when the new file will be created
    """
    NCBI_species_name = ""
    NCBI_taxid = "123"
    NCBI_species_taxid = ""
    NCBI_section = "genbank"
    levels = ""
    outdir = ""
    tmp_dir = os.path.join("123", "temporary_directory")
    threads = 1
    norefseq = False
    info_file = ""
    db_dir = ""
    only_mash = False
    l90 = 100
    nbcont = 999
    cutn = 5
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 2
    quiet = False
    res_outdir = "123"
    with pytest.raises(SystemExit):
        prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, NCBI_section, 
                     outdir, tmp_dir, threads, norefseq,
                     db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist, max_dist,
                     verbose, quiet)
    _, err = capsys.readouterr()
    assert ("Could not download genomes. Check that you gave valid NCBI taxid and/or "
            "NCBI species name. If you gave both, check that given taxID and name really "
            "correspond to the same species.") in err
    # Check output files
    summary =  os.path.join(res_outdir, "assembly_summary-123.txt")
    assert not os.path.isfile(summary)
    ngd_outdir = os.path.join(res_outdir, "genbank", "bacteria")
    assert not os.path.isdir(ngd_outdir)
    # # Check logfiles are here
    log_files = glob.glob(os.path.join(res_outdir, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty asnothing is downloaded
    assert len(os.listdir(tmp_dir)) == 0
    # Check Database_init folder created, with at list 4 ".fna" genomes
    assert not os.path.isdir(os.path.join(res_outdir, "Database_init"))

    # Remove output directory
    shutil.rmtree(res_outdir, ignore_errors=True)


def test_main_norefseq_wrongdbpath(capsys):
    """
    We run with option norefseq, but given db_dir does not exist.
    -> error message
    """
    NCBI_species_name = ""
    NCBI_species_taxid = ""
    NCBI_taxid = ""
    NCBI_section = "refseq"
    levels = ""
    outdir = GENEPATH
    tmp_dir = os.path.join(outdir, "temporary_directory")
    threads = 1
    norefseq = True
    db_dir = "dbdir"
    only_mash = False
    l90 = 100
    nbcont = 999
    cutn = 5
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 15
    quiet = False
    info_file = ""
    with pytest.raises(SystemExit):
        prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, NCBI_section,
                     outdir, tmp_dir, threads, norefseq,
                     db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist, max_dist,
                     verbose, quiet)
    _, err = capsys.readouterr()
    assert ("You asked to skip refseq downloads") in err
    assert ("Database folder dbdir supposed to contain fasta sequences does not exist. Please "
            "give a valid folder, or leave the default directory (no '-d' option)") in err
    # Check output files
    summary =  os.path.join(GENEPATH, "assembly_summary-123.txt")
    assert not os.path.isfile(summary)
    ngd_outdir = os.path.join(GENEPATH, "refseq", "bacteria")
    assert not os.path.isdir(ngd_outdir)
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 4  #.log.debug as we put verbose = 15
    # Check tmp files folder created, but empty asnothing is downloaded
    assert len(os.listdir(tmp_dir)) == 0
    # Check Database_init folder created, with at list 4 ".fna" genomes
    assert not os.path.isdir(os.path.join(GENEPATH, "Database_init"))


def test_main_norefseq_nodefault_dbdir_nor_refseq(capsys):
    """
    We run with option norefseq, but given db_dir does not exist.
    -> error message
    """
    NCBI_species_name = ""
    NCBI_species_taxid = ""
    NCBI_taxid = ""
    NCBI_section = "genbank"
    levels = ""
    outdir = GENEPATH
    tmp_dir = ""
    threads = 1
    norefseq = True
    db_dir = ""
    only_mash = False
    l90 = 100
    nbcont = 999
    cutn = 5
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 2
    quiet = False
    info_file = ""
    with pytest.raises(SystemExit):
        prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, 
                     NCBI_section, outdir, tmp_dir, threads, norefseq,
                     db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist, max_dist,
                     verbose, quiet)
    _, err = capsys.readouterr()
    assert ("You asked to skip genbank downloads") in err
    assert ("Database folder test/data/prepare/generated_by_func-tests/Database_init supposed "
            "to contain fasta sequences does not exist. We will check if the download folder "
            "(with compressed sequences) exists.") in err
    assert ("Folder test/data/prepare/generated_by_func-tests/genbank/bacteria "
            "does not exist. You do not have any genome to analyse. Possible reasons:\n") in err
    assert ("- if you want to rerun analysis in the same folder as "
            "sequences were downloaded (my_outdir/Database_init or "
            "my_outdir/genbank), make sure you have '-o my_outdir' option\n") in err
    assert ("- if you want to rerun analysis and save them in a new "
            "output folder called 'new_outdir', make sure you have '-o new_outdir' option, "
            "and you specified where the uncompressed sequences to use are "
            "('-d sequence_database_path'") in err
    # # Check output files
    summary =  os.path.join(GENEPATH, "assembly_summary-123.txt")
    assert not os.path.isfile(summary)
    ngd_outdir = os.path.join(GENEPATH, "genbank", "bacteria")
    assert not os.path.isdir(ngd_outdir)
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty asnothing is downloaded
    assert len(os.listdir(os.path.join(GENEPATH, "tmp_files"))) == 0
    # Check Database_init folder created, with at list 4 ".fna" genomes
    assert not os.path.isdir(os.path.join(GENEPATH, "Database_init"))


def test_main_norefseq_nodefault_dbdir_but_refseq(capsys):
    """
    We run with option norefseq, but given db_dir does not exist.
    -> error message
    """
    NCBI_species_name = ""
    NCBI_species_taxid = "123"
    NCBI_taxid = ""
    NCBI_section = "genbank"
    levels = ""
    # Copy refseq/bacteria and content into outdirectory
    outdir = GENEPATH
    tmp_dir = ""
    threads = 1
    norefseq = True
    orig_dbdir = os.path.join(GEN_PATH, "refseq")
    refseq_db_dir = os.path.join(GENEPATH, "genbank")
    shutil.copytree(orig_dbdir, refseq_db_dir)
    db_dir = ""
    only_mash = False
    l90 = 100
    nbcont = 999
    cutn = 0
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 2
    quiet = False
    info_file = ""
    out_info_file = os.path.join(outdir, f"LSTINFO-123-filtered-0.0001_0.06.txt")
    assert prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, 
                        NCBI_section, outdir, tmp_dir, threads,
                        norefseq, db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist,
                        max_dist, verbose, quiet) == out_info_file
    out, err = capsys.readouterr()
    assert ("You asked to skip genbank downloads") in err
    assert ("Database folder test/data/prepare/generated_by_func-tests/"
            "Database_init supposed "
            "to contain fasta sequences does not exist. We will check if the download folder "
            "(with compressed sequences) exists.") in err
    assert ("Uncompressing genome files") in out
    assert ("Total number of genomes for 123: 3") in out
    assert ("Computing pairwise distances between all genomes") in out
    assert ("Final number of genomes in dataset: 1") in out
    # Check output files
    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(GENEPATH, "genbank", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) == 3
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty as we do not split
    tmp_folder = os.listdir(os.path.join(GENEPATH, "tmp_files"))
    assert len(tmp_folder) == 0
    # Check Database_init folder created, with the 3 ".fna" genomes
    fna_files = glob.glob(os.path.join(GENEPATH, "Database_init", "*.fna"))
    assert len(fna_files) == 3


def test_main_norefseq_defaultdbdir(capsys):
    """
    We run with option norefseq, but given db_dir does not exist.
    -> error message
    """
    NCBI_species_name = ""
    NCBI_species_taxid = ""
    NCBI_taxid = ""
    NCBI_section = "refseq"
    levels = ""
    # Copy refseq/bacteria and content into outdirectory
    outdir = GENEPATH
    tmp_dir = ""
    threads = 1
    norefseq = True
    orig_dbdir = os.path.join(GEN_PATH, "genomes_comparison")
    refseq_db_dir = os.path.join(GENEPATH, "Database_init")
    shutil.copytree(orig_dbdir, refseq_db_dir)
    db_dir = ""
    only_mash = False
    l90 = 100
    nbcont = 999
    cutn = 0
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 2
    quiet = False
    info_file = ""
    out_info_file = os.path.join(outdir, "LSTINFO-NA-filtered-0.0001_0.06.txt")
    assert prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, 
                        NCBI_section, outdir, tmp_dir, threads,
                        norefseq, db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist,
                        max_dist, verbose, quiet) == out_info_file
    out, err = capsys.readouterr()
    assert ("You asked to skip refseq downloads") in err
    assert ("Total number of genomes for NA: 5") in out
    assert ("Computing pairwise distances between all genomes") in out
    assert ("Final number of genomes in dataset: 1") in out
    # Check output files
    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(GENEPATH, "refseq", "bacteria")
    assert not os.path.isdir(ngd_outdir)
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty as we do not split
    tmp_folder = os.listdir(os.path.join(GENEPATH, "tmp_files"))
    assert len(tmp_folder) == 0
    # Check Database_init folder created, with the 3 ".fna" genomes
    fna_files = glob.glob(os.path.join(GENEPATH, "Database_init", "*.fna"))
    assert len(fna_files) == 5


def test_main_norefseq_givendbdir(capsys):
    """
    We run with option norefseq, but given db_dir does not exist.
    -> error message
    """
    NCBI_species_name = ""
    NCBI_species_taxid = ""
    NCBI_taxid = ""
    NCBI_section = "refseq"
    levels = ""
    # Copy refseq/bacteria and content into outdirectory
    outdir = GENEPATH
    tmp_dir = ""
    threads = 1
    norefseq = True
    orig_dbdir = os.path.join(GEN_PATH, "genomes_comparison")
    refseq_db_dir = os.path.join(GENEPATH, "genomes_comparison")
    shutil.copytree(orig_dbdir, refseq_db_dir)
    db_dir = refseq_db_dir
    only_mash = False
    l90 = 100
    nbcont = 999
    cutn = 2
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 2
    quiet = False
    info_file = ""
    out_info_file = os.path.join(outdir, "LSTINFO-NA-filtered-0.0001_0.06.txt")
    assert prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, 
                        NCBI_section, outdir, tmp_dir, threads,
                        norefseq, db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist,
                        max_dist, verbose, quiet) == out_info_file
    out, err = capsys.readouterr()
    assert ("You asked to skip refseq downloads") in err
    assert ("Total number of genomes for NA: 5") in out
    assert ("Computing pairwise distances between all genomes") in out
    assert ("Final number of genomes in dataset: 1") in out
    # Check output files
    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(GENEPATH, "refseq", "bacteria")
    assert not os.path.isdir(ngd_outdir)
    # Check logfiles are here
    log_files = glob.glob(os.path.join(GENEPATH, "*log*"))
    assert len(log_files) == 3
    # Check tmp files folder created, but empty as we do not split
    tmp_files = glob.glob(os.path.join(GENEPATH, "tmp_files", "*.fna_prepare-split2N.fna"))
    assert len(tmp_files) == 5


def test_only_mash(capsys):
    """
    Running only mash step (giving genomes and corresponding LSTINFO file)
    """
    NCBI_species_name = ""
    NCBI_species_taxid = ""
    NCBI_taxid = ""
    NCBI_section = "refseq"
    levels = ""
    outdir = GENEPATH
    tmp_dir = ""
    threads = 1
    norefseq = False
    db_dir = ""
    only_mash = True
    info_file = os.path.join(TEST_DIR, "test_lstinfo_onlymash.lst")
    l90 = 100
    nbcont = 999
    cutn = 5
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 1
    quiet = False
    out_info_file = os.path.join(outdir, "LSTINFO-NA-filtered-0.0001_0.06.txt")
    assert prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels, 
                        NCBI_section, outdir, tmp_dir, threads,
                        norefseq, db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist,
                        max_dist, verbose, quiet) == out_info_file
    out, err = capsys.readouterr()
    assert ("You asked to run only mash steps") in err
    assert ("You want to run only mash steps. Getting information from "
            "test/data/prepare/test_files/test_lstinfo_onlymash.lst") in out
    assert ("Found 5 genomes in total") in out
    assert ("Computing pairwise distances between all genomes") in out
    assert ("Sorting all 5 genomes by quality") in out
    assert ("Final number of genomes in dataset: 1") in out

    # Check output files
    assert len(os.listdir(os.path.join(outdir, "tmp_files"))) == 0
    # Check logfiles are here
    log_files = glob.glob(os.path.join(outdir, "*log*"))
    assert len(log_files) == 3
    # Check content of output lstinfo file
    out_lst = os.path.join(outdir, "LSTINFO-NA-filtered-0.0001_0.06.txt")
    exp_lst = os.path.join(DBDIR, "exp_files", "exp_lstinfo_run_only-mash.lst")
    assert tutil.compare_order_content(out_lst, exp_lst)


def test_only_mash_empty_lstinfo(capsys):
    """
    Running only mash step giving an empty lstinfo file -> error, no genome found
    """
    NCBI_species_name = ""
    NCBI_species_taxid = ""
    NCBI_taxid = ""
    NCBI_section = "refseq"
    levels = ""
    outdir = GENEPATH
    tmp_dir = ""
    threads = 1
    norefseq = False
    db_dir = ""
    only_mash = True
    # Create empty lstinfo file
    info_file = os.path.join(GENEPATH, "LSTINFO-empty.lst")
    open(info_file, "w").close()
    l90 = 100
    nbcont = 999
    cutn = 5
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 1
    quiet = False
    with pytest.raises(SystemExit):
        prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels,
                     NCBI_section, outdir, tmp_dir, threads, norefseq,
                     db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist, max_dist,
                     verbose, quiet)
    out, err = capsys.readouterr()
    assert ("You asked to run only mash steps") in err
    assert ("You want to run only mash steps. Getting information from "
            "test/data/prepare/generated_by_func-tests/LSTINFO-empty.lst") in out
    assert ("No genome listed in test/data/prepare/generated_by_func-tests/LSTINFO-empty.lst "
            "was found.") in err

    # Check output files
    assert len(os.listdir(os.path.join(outdir, "tmp_files"))) == 0
    # Check logfiles are here
    log_files = glob.glob(os.path.join(outdir, "*log*"))
    assert len(log_files) == 3
    # Check lstinfo file is still here and still empty
    assert os.stat(info_file).st_size == 0


def test_only_mash_no_lstinfo(capsys):
    """
    Running only mash step giving an info file which does not exist -> error missing infofile
    """
    NCBI_species_name = ""
    NCBI_species_taxid = ""
    NCBI_taxid = ""
    NCBI_section = "refseq"
    levels = ""
    outdir = GENEPATH
    tmp_dir = ""
    threads = 1
    norefseq = False
    db_dir = ""
    only_mash = True
    # Create empty lstinfo file
    info_file = "info_file.lst"
    l90 = 100
    nbcont = 999
    cutn = 5
    min_dist = 1e-4
    max_dist = 0.06
    verbose = 1
    quiet = False
    with pytest.raises(SystemExit):
        prepare.main("cmd", NCBI_species_name, NCBI_species_taxid, NCBI_taxid, levels,
                     NCBI_section, outdir, tmp_dir, threads, norefseq,
                     db_dir, only_mash, info_file, l90, nbcont, cutn, min_dist, max_dist,
                     verbose, quiet)
    out, err = capsys.readouterr()
    assert ("You asked to run only mash steps") in err
    assert ("Your info file info_file.lst does not exist. Please provide the  "
            "right name/path, or remove the '--mash-only option to rerun "
            "quality control.") in err

    # Check output files
    assert len(os.listdir(os.path.join(outdir, "tmp_files"))) == 0
    # Check logfiles are here
    log_files = glob.glob(os.path.join(outdir, "*log*"))
    assert len(log_files) == 3
    # Check that outdir contains only 4 elements: 3 logs + tmp_files repo
    files = os.listdir(outdir)
    files = [f for f in files if "fuse" not in f]
    assert len(files) == 4
