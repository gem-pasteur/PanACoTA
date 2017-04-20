#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pytest
import os
import pipelinepackage.prokka_functions as pfunc

def test_count_tbl():
    """
    Count the different features found in the tbl file, and return
    nbcont, nbCDS, nbGene, nbCRISPR
    """
    tblfile = os.path.join("test", "data", "test_files", "prokka_out_for_test.tbl")
    ncont, ncds, ngene, ncris = pfunc.count_tbl(tblfile)
    assert ncont == 7
    assert ncds == 14
    assert ngene == 16
    assert ncris == 2


def test_count_headers():
    """
    Count how many sequences there are in the given multi-fasta file
    """
    seqfile = os.path.join("test", "data", "genomes", "genome4.fasta")
    nb = pfunc.count_headers(seqfile)
    assert nb == 5

def test_check_prokka_no_outdir():
    """

    """
    check log

# def test_check_prokka_notbl():
#     """
#     Check that check_prokka returns false when a tbl file is missing
#     """
#     outdir = os.path.join("test", "data", "test_files")
#     name = "prokka_out_for_test-misstbl"
#     logf = "prokka.log"
#     gpath = "path/to/nogenome/original_name.fna"
#     nbcont = 7
#     assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)


# def test_check_prokka_nofaa():
#     """
#     Check that check_prokka returns false when a tbl file is missing
#     """
#     outdir = os.path.join("test", "data", "test_files")
#     name = "prokka_out_for_test-missfaa"
#     logf = "prokka.log"
#     gpath = "path/to/nogenome/original_name.fna"
#     nbcont = 7
#     assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)


def test_check_prokka_noffn():
    """
    Check that check_prokka returns false when a tbl file is missing
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-missffn"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)

def test_check_prokka_wrong_cont():
    """

    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 10
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    + check log

check wrong cds/prot
check wrong gene


# def test_check_prokka_ok():
#     """
#     Check that everything is ok with prokka results (tbl, faa and ffn files exist,
#     and number of CDS, CRISPR and genes correspond between them)
#     """
#     outdir = os.path.join("test", "data", "test_files")
#     name = "prokka_out_for_test"
#     logf = "prokka.log"
#     gpath = "path/to/nogenome/original_name.fna"
#     nbcont = 7
#     assert pfunc.check_prokka(outdir, logf, name, gpath, nbcont)


