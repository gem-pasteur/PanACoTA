#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the download_genomes_func submodule in prepare module
"""
import os
import logging
import shutil

import PanACoTA.prepare_module.filter_genomes as filterg

DATA_TEST_DIR = os.path.join("test", "data", "prepare")

def test_write_output(caplog):
    """
    Check that the files with kept genomes and discarded genomes are created
    """
    corresp_genomes = {"ACOR001": "ACOR001.0519.fna.gz", "ACOR002": "ACOR002.0519.fna.gz",
                       "ACOR003": "ACOR003.0519.fna.gz"}
    sorted_genomes = [os.path.join(DATA_TEST_DIR, "genomes", "refseq", "bacteria", gen, gz)
                      for gen, gz in corresp_genomes.items()]
    genomes = {genome_file:["gname", "ori_name", "path_annotate", 12, 10, 1]
               for genome_file in sorted_genomes}
    genomes["toto"] = ['totoname', 'ori toto', 'path_toto', 13, 5, 6]
    genomes_removed = {"genome": ["ref", 10]}

    # Define output directory for generated files
    outdir = os.path.join(DATA_TEST_DIR, "test_filter_write_output")
    os.makedirs(outdir)
    gspecies = "Acetobacter_fabarum"
    min_dist = 0.06

    # Check everything works without error
    assert filterg.write_outputfiles(genomes, sorted_genomes,
                                     genomes_removed, outdir, gspecies, min_dist) == 0

    # Check outfiles exist
    list_file = os.path.join(outdir, "LSTINFO-Acetobacter_fabarum-filtered-0.06.txt")
    discard_file = os.path.join(outdir, "discarded-by-minhash-Acetobacter_fabarum-0.06.txt")
    assert os.path.isfile(list_file)
    assert os.path.isfile(discard_file)

    # Check content of out files
    with open(list_file) as lf:
        # Check header
        assert "to_annotate\tgsize\tnb_conts\tL90" in lf.readline()
        # Check there are 3 genomes, with expected information
        assert "path_annotate 12\t10\t1" in lf.readline()
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert lf.readline() == "\n"

    with open(discard_file) as df:
        # Check header
        assert "genome_name\tproblem_compared_with\tdist" in df.readline()
        # Check genome line
        assert "genome" in df.readline()
        # Check no more genomes
        assert df.readline() == '\n'

    # Remove test folder
    # shutil.rmtree(outdir)