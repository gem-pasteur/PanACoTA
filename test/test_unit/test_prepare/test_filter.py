#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the download_genomes_func submodule in prepare module
"""
import os
import logging
import shutil
import pytest

import PanACoTA.prepare_module.filter_genomes as filterg

DATA_TEST_DIR = os.path.join("test", "data", "prepare")

def test_write_output():
    """
    Check that the files with kept genomes and discarded genomes are created
    """
    corresp_genomes = {"ACOR001": "ACOR001.0519.fna.gz", "ACOR002": "ACOR002.0519.fna.gz",
                       "ACOR003": "ACOR003.0519.fna.gz"}
    sorted_genomes = [os.path.join(DATA_TEST_DIR, "genomes", "refseq", "bacteria", gen, gz)
                      for gen, gz in corresp_genomes.items()]
    genomes = {genome_file:["gname", "ori_name", "path_annotate", 12, 10, 1]
               for genome_file in sorted_genomes}
    toto_path = os.path.join(DATA_TEST_DIR, "genomes", "refseq", "bacteria", "toto", "toto.fna.gz")
    genomes[toto_path] = ['totoname', 'ori toto', 'path_toto', 13, 5, 6]
    genomes_removed = {"genome": ["ref", 10]}
    sorted_genomes.append(toto_path)

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
    # File containing genomes kept, that will be used for annotation
    with open(list_file) as lf:
        # Check header
        assert "to_annotate\tgsize\tnb_conts\tL90" in lf.readline()
        # Check there are 3 genomes, with expected information
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_toto\t13\t5\t6" in lf.readline()
        assert lf.readline() == ""
        assert lf.readline() == ""

    # File with discarded genomes, and why they are discarded
    with open(discard_file) as df:
        # Check header
        assert "genome_name\tproblem_compared_with\tdist" in df.readline()
        # Check genome line
        assert "genome\tref\t10" in df.readline()
        # Check no more genomes
        assert df.readline() == ''
        assert df.readline() == ''

    # Remove test folder
    shutil.rmtree(outdir)



def test_write_output_no_discard():
    """
    Check that the files with kept genomes and discarded genomes are created, but
    file with discarded genomes is empty
    """
    corresp_genomes = {"ACOR001": "ACOR001.0519.fna.gz", "ACOR002": "ACOR002.0519.fna.gz",
                       "ACOR003": "ACOR003.0519.fna.gz"}
    sorted_genomes = [os.path.join(DATA_TEST_DIR, "genomes", "refseq", "bacteria", gen, gz)
                      for gen, gz in corresp_genomes.items()]
    genomes = {genome_file:["gname", "ori_name", "path_annotate", 12, 10, 1]
               for genome_file in sorted_genomes}
    toto_path = os.path.join(DATA_TEST_DIR, "genomes", "refseq", "bacteria", "toto", "toto.fna.gz")
    genomes[toto_path] = ['totoname', 'ori toto', 'path_toto', 13, 5, 6]
    genomes_removed = {}
    sorted_genomes.append(toto_path)

    # Define output directory for generated files
    outdir = os.path.join(DATA_TEST_DIR, "test_filter_write_output_no_discard")
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
    # File containing genomes kept, that will be used for annotation
    with open(list_file) as lf:
        # Check header
        assert "to_annotate\tgsize\tnb_conts\tL90" in lf.readline()
        # Check there are 3 genomes, with expected information
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_annotate\t12\t10\t1" in lf.readline()
        assert "path_toto\t13\t5\t6" in lf.readline()
        assert lf.readline() == ""
        assert lf.readline() == ""

    # File with discarded genomes, and why they are discarded
    with open(discard_file) as df:
        # Check header
        assert "genome_name\tproblem_compared_with\tdist" in df.readline()
        # Check no genome
        assert df.readline() == ''
        assert df.readline() == ''

    # Remove test folder
    shutil.rmtree(outdir)


def test_write_output_no_genome():
    """
    Check that the files with kept genomes and discarded genomes are created, but
    are empty because tehre is no genome
    """
    corresp_genomes = {"ACOR001": "ACOR001.0519.fna.gz", "ACOR002": "ACOR002.0519.fna.gz",
                       "ACOR003": "ACOR003.0519.fna.gz"}
    sorted_genomes = []
    genomes = {}
    genomes_removed = {}

    # Define output directory for generated files
    outdir = os.path.join(DATA_TEST_DIR, "test_filter_write_output_no_genome")
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
    # File containing genomes kept, that will be used for annotation
    with open(list_file) as lf:
        # Check header
        assert "to_annotate\tgsize\tnb_conts\tL90" in lf.readline()
        # Check there is not any genome
        assert lf.readline() == ""
        assert lf.readline() == ""

    # File with discarded genomes, and why they are discarded
    with open(discard_file) as df:
        # Check header
        assert "genome_name\tproblem_compared_with\tdist" in df.readline()
        # Check no genome
        assert df.readline() == ''
        assert df.readline() == ''

    # Remove test folder
    shutil.rmtree(outdir)


def test_write_output_no_outdir(caplog):
    """
    Check that when outdir does not exist, program ends with error message
    """
    corresp_genomes = {"ACOR001": "ACOR001.0519.fna.gz", "ACOR002": "ACOR002.0519.fna.gz",
                       "ACOR003": "ACOR003.0519.fna.gz"}
    sorted_genomes = []
    genomes_removed = {}

    # Define output directory for generated files
    outdir = os.path.join(DATA_TEST_DIR, "test_filter_write_output_no_outdir")
    gspecies = "Acetobacter_fabarum"

    # Check everything works without error
    with pytest.raises(SystemExit):
        filterg.write_outputfiles({}, sorted_genomes, genomes_removed, outdir, gspecies, 0.06)

    # Check outfiles do not exist
    list_file = os.path.join(outdir, "LSTINFO-Acetobacter_fabarum-filtered-0.06.txt")
    discard_file = os.path.join(outdir, "discarded-by-minhash-Acetobacter_fabarum-0.06.txt")
    assert not os.path.isfile(list_file)
    assert not os.path.isfile(discard_file)

    # Check output
    assert not os.path.isdir(outdir)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert ("The given output directory (test/data/prepare/test_filter_write_output_no_outdir) "
            "does not exist. We cannot create output files there") in caplog.text


def test_read_matrix():
    """
    Test that matrix from minhash is correctly read and converted to python matrix
    """
    genomes = {"genome1": ["g1_name", "g1_ori", "path_genome1", 1500, 5, 2],
               "genome2": ["g2_name", "g2_ori", "path_genome2", 20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", "path_genome3", 25003, 52, 50],
               "genome4": ["g4_name", "g4_ori", "path_genome4", 22012, 20, 10]
              }
    sorted_genomes = ["genome1", "genome2", "genome3", "genome4"]
    matrix_file = os.path.join(DATA_TEST_DIR, "test_files", "minhash_output.txt")

    mat_sp = filterg.read_matrix(genomes, sorted_genomes, matrix_file)

    # Check content of created matrix
    assert len(mat_sp) == 6
    assert mat_sp[0,1] ==0.5
    assert mat_sp[0,2] == 0.06
    assert mat_sp[0,3] == 0.005
    assert mat_sp[1,2] == 0.09
    assert mat_sp[1,3] == 0.7
    assert mat_sp[2,3] == 0.08


def test_read_matrix_no_genome():
    """
    Test that matrix from minhash is correctly read and converted to python matrix
    """
    genomes = {}
    sorted_genomes = []
    matrix_file = os.path.join(DATA_TEST_DIR, "test_files", "minhash_output_empty.txt")
    # Create empty matrix file
    open(matrix_file, "w").close()

    # Read matrix
    mat_sp = filterg.read_matrix(genomes, sorted_genomes, matrix_file)

    # Check content of created matrix
    assert len(mat_sp) == 0

    # Remove empty matrix file
    os.remove(matrix_file)


