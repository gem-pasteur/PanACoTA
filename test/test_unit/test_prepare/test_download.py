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
    # print("CAPLOG = ")
    # print(caplog.text)
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
    shutil.rmtree(os.path.join(DATA_TEST_DIR, "refseq"))


test folders for genomes, but 1 empty
test folders for genomes, but 1 with several genome
# # def test_to_database_nofasta():
#     """
#     Test behavior when there are no fasta files in the given database
#     """
#     wrong_path = os.path.join(DATA_TEST_DIR)
#     os.makedirs(wrong_path)
#     nbgen, db_init_path = downg.to_database(wrong_path)
