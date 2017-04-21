#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pytest
import os
import shutil
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


def test_check_prokka_no_outdir(capsys):
    """
    Test that prokka returns the right error message when output directory does not exist
    """
    outdir = "toto"
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "Prokka could not run properly. Look at prokka.log for more information.\n"


def test_check_prokka_notbl(capsys):
    """
    Check that check_prokka returns false when a tbl file is missing, and an error message
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-misstbl"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-misstbl original_name.fna: no .tbl file\n"
    os.remove(os.path.join(ori_dir, name + ".faa"))
    os.remove(os.path.join(ori_dir, name + ".ffn"))


def test_check_prokka_nofaa(capsys):
    """
    Check that check_prokka returns false when a faa file is missing, and an error message
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-missfaa"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".tbl"), os.path.join(ori_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missfaa original_name.fna: no .faa file\n"
    os.remove(os.path.join(ori_dir, name + ".tbl"))
    os.remove(os.path.join(ori_dir, name + ".ffn"))


def test_check_prokka_noffn(capsys):
    """
    Check that check_prokka returns false when a ffn file is missing, and an error message
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-missffn"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".tbl"), os.path.join(ori_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missffn original_name.fna: no .ffn file\n"
    os.remove(os.path.join(ori_dir, name + ".tbl"))
    os.remove(os.path.join(ori_dir, name + ".faa"))


def test_check_prokka_wrong_cont(capsys):
    """
    Check that check_prokka returns an error message when the number of contigs in tbl
    file is not as expected
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 10
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == ("prokka_out_for_test original_name.fna: no matching number of contigs; "
                   "nbcontig=10; in tbl =7\n")


def test_check_prokka_wrong_tblCDS(capsys):
    """
    Check that check_prokka returns an error message when the number of CDS in tbl
    file is different from the number of headers in faa file
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-wrongCDS"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err.split("\n")[0] == ("prokka_out_for_test-wrongCDS original_name.fna: "
                                  "no matching number of proteins between tbl and faa; "
                                  "faa=14; in tbl =13")
    assert err.split("\n")[1] == ("prokka_out_for_test-wrongCDS original_name.fna: "
                                  "no matching number of genes between tbl and ffn; "
                                  "ffn=18; in tbl =15genes 2CRISPR")
    os.remove(os.path.join(ori_dir, name + ".ffn"))
    os.remove(os.path.join(ori_dir, name + ".faa"))


def test_check_prokka_wrong_tblCRISPR(capsys):
    """
    Check that check_prokka returns an error message when the number of headers in ffn
    file is different from the number of CDS + CRISPR in tbl file (1CRISPR in tbl, 2 in ffn)
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-wrongtblCRISP"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    print(err)
    assert err == ("prokka_out_for_test-wrongtblCRISP original_name.fna: "
                   "no matching number of genes between tbl and ffn; "
                   "ffn=18; in tbl =16genes 1CRISPR\n")
    os.remove(os.path.join(ori_dir, name + ".ffn"))
    os.remove(os.path.join(ori_dir, name + ".faa"))


def test_check_prokka_ok():
    """
    Check that everything is ok with prokka results (tbl, faa and ffn files exist,
    and number of CDS, CRISPR and genes correspond between them)
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert pfunc.check_prokka(outdir, logf, name, gpath, nbcont)


