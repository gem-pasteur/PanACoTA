#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the download_genomes_func submodule in prepare module
"""
import os
import logging
import glob
import shutil
import pytest

import PanACoTA.prepare_module.download_genomes_func as downg


DATA_TEST_DIR = os.path.join("test", "data", "prepare")

def test_to_database():
    """
    Test that all fna.gz files are uncompressed and moved to a created Database_init folder
    """
    out_dir = os.path.join(DATA_TEST_DIR, "genomes")
    nb_gen, db_init_dir = downg.to_database(out_dir)
    db_dir = os.path.join(DATA_TEST_DIR, "genomes", "Database_init")
    assert os.path.isdir(db_dir)
    files_all = glob.glob(os.path.join(db_dir, "*"))
    files_fna = glob.glob(os.path.join(db_dir, "*.fna"))
    # Check that there are only 3 files in result database
    assert len(files_all) == len(files_fna)
    # And that those files are .fna files
    assert len(files_fna) == 3
    # Check that we have as many genomes as expected, and that the output database has the
    # expected name
    assert nb_gen == 3
    assert db_init_dir == db_dir
    assert os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))

    # Remove database created
    shutil.rmtree(db_dir)


def test_to_database_nofolder_refseq(caplog):
    """
    Test behavior when the folder that should contain refseq downloaded genomes does not exist
    -> should exit with error message
    """
    caplog.set_level(logging.DEBUG)
    with pytest.raises(SystemExit):
        downg.to_database(DATA_TEST_DIR)

    assert "ERROR" in caplog.text
    assert ("The folder containing genomes downloaded from NCBI refseq "
            "(test/data/prepare/refseq/bacteria) does not exist.") in caplog.text
    assert ("Check that you really downloaded sequences (fna.gz) and that they are "
            "in this folder") in caplog.text


def test_to_database_nofolder_per_genome(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, but there are no folders inside
    -> should exit with error message
    """
    empty_dir = os.path.join(DATA_TEST_DIR, "refseq", "bacteria")
    os.makedirs(empty_dir)
    caplog.set_level(logging.DEBUG)
    with pytest.raises(SystemExit):
        downg.to_database(DATA_TEST_DIR)
    # Check error message is as expected
    assert "ERROR" in caplog.text
    assert ("The folder supposed to contain genomes downloaded from NCBI refseq "
            "(test/data/prepare/refseq/bacteria) exists but is empty") in caplog.text
    assert ("Check that you really downloaded sequences (fna.gz)") in caplog.text

    # Remove files/folders specific to test
    shutil.rmtree(os.path.join(DATA_TEST_DIR, "refseq"))


def test_to_database_1empty_genome_folder(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, there are subfolders inside,
    but 1 of them is empty: warning message informing that this genome will be ignored
    """
    out_dir = os.path.join(DATA_TEST_DIR, "genomes")
    gz_genomes_folder = os.path.join(out_dir, "refseq", "bacteria")

    # Empty 1 directory: move its file to 'out_dir'
    to_move_filename = "ACOR002.0519.fna.gz" # File that must be moved
    to_empty_dir = "ACOR002" # Directory containing file to move
    to_move_file = os.path.join(gz_genomes_folder, to_empty_dir, to_move_filename)
    shutil.move(to_move_file, os.path.join(out_dir, to_move_filename))

    # Run to_database
    nb_gen, db_dir = downg.to_database(out_dir)
    assert nb_gen == 2
    assert db_dir == os.path.join(out_dir, "Database_init")

    # Check that a warning message was raised, indicating that genome is ignored
    caplog.set_level(logging.DEBUG)
    assert "WARNING" in caplog.text
    assert ("Problem with genome in ACOR002: no compressed fasta file downloaded. "
            "This genome will be ignored.") in caplog.text
    assert not os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR001.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))

    # Remove files/folders specific to test
    shutil.move(os.path.join(out_dir, to_move_filename), to_move_file)
    shutil.rmtree(db_dir)


def test_to_database_several_genomes(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, there are subfolders inside,
    but 1 of them contains more than 1 genome: warning message informing that this
    genome will be ignored
    """
    out_dir = os.path.join(DATA_TEST_DIR, "genomes")
    gz_genomes_folder = os.path.join(out_dir, "refseq", "bacteria")

    # Create a new gz file in one of the genome directories
    to_create_filename = "ACOR002.0519.bis.fna.gz" # File that must be moved
    to_fill_dir = "ACOR002" # Directory containing file to move
    to_create_path = os.path.join(gz_genomes_folder, to_fill_dir, to_create_filename)
    # Create empty gz file
    open(to_create_path, "w").close()

    # Run to_database, and check that only 2 genomes were considered
    nb_gen, db_dir = downg.to_database(out_dir)
    assert nb_gen == 2
    assert db_dir == os.path.join(out_dir, "Database_init")

    # Check that a warning message was raised, indicating that genome is ignored
    caplog.set_level(logging.DEBUG)
    assert "WARNING" in caplog.text
    assert ("Problem with genome in ACOR002: several compressed fasta files found. "
            "This genome will be ignored.") in caplog.text
    assert not os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR001.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))

    # Remove test files/folders
    os.remove(to_create_path)
    shutil.rmtree(db_dir)


def test_to_database_1genome_wrong_format(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, there is 1 genome per subfolder,
    but 1 genome cannot be unzipped
    """
    out_dir = os.path.join(DATA_TEST_DIR, "genomes")
    gz_genomes_folder = os.path.join(out_dir, "refseq", "bacteria")

    # Name of directory directly containing the original gz file
    to_corrupt_dir = "ACOR001"
    # Name of original gz file that must be moved to be saved
    to_empty_filename = "ACOR001.0519.fna.gz"
    # Complete path to this original gz file
    to_empty_path = os.path.join(gz_genomes_folder, to_corrupt_dir, to_empty_filename) #
    # copy real gz genome file to outdir to save it, and create a fake one in place of it
    shutil.copy(to_empty_path, os.path.join(out_dir, to_empty_filename))
    # Create fake gz file (txt file)
    false_gz = open(to_empty_path, "w")
    false_gz.write("This is not a gz file")
    false_gz.close()

    # Run to_database
    nb_gen, db_dir = downg.to_database(out_dir)
    assert nb_gen == 2
    assert db_dir == os.path.join(out_dir, "Database_init")

    # Check that a error message was raised, indicating that genome is ignored
    caplog.set_level(logging.DEBUG)
    assert "ERROR" in caplog.text
    assert ("Error while trying to uncompress "
            "test/data/prepare/genomes/Database_init/ACOR001.0519.fna.gz. "
            "This genome will be ignored") in caplog.text
    # Check that there are only 2 files in the database, and that they correspond
    # to uncompressed gz files
    list_db = os.listdir(db_dir)
    assert len(list_db) == 2
    assert not os.path.isfile(os.path.join(db_dir, to_empty_filename))
    assert os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))

    # Remove test files/Folders
    shutil.move(os.path.join(out_dir, to_empty_filename), to_empty_path)
    shutil.rmtree(db_dir)
