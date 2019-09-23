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
GENOMES_DIR = os.path.join(DATA_TEST_DIR, "genomes", "refseq", "bacteria")

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
    are empty because there is no genome
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


def test_sort_genomes_minhash():
    """
    Test that a given set of genomes with their proprieties is sorted as expected
    (genomes with too bad quality removed, others sorted by l90 and nbcont)
    """
    genomes = {"genome1": ["g1_name", "g1_ori", "path_genome1", 123567, 200, 101],
               "genome2": ["g2_name", "g2_ori", "path_genome2", 20000, 3, 3],
               "genome3": ["g3_name", "g3_ori", "path_genome3", 25003, 52, 50],
               "genome4": ["g4_name", "g4_ori", "path_genome4", 123456, 1023, 11],
               "genome5": ["g5_name", "g5_ori", "path_genome5", 22012, 20, 10],
               "genome6": ["g6_name", "g6_ori", "path_genome6", 1500, 3, 2]
              }
    max_l90 = 100
    max_cont = 1000
    sorted_genomes = filterg.sort_genomes_minhash(genomes, max_l90, max_cont)
    assert sorted_genomes == ["genome6", "genome2", "genome5", "genome3"]


def test_sketch_all():
    """
    Test that all genomes are sketch, in the provided order
    """
    genomes = {"genome1": ["g1_name", "g1_ori", os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                           123567, 200, 101],
               "genome2": ["g2_name", "g2_ori", os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                           20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOR003.0519.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome3", "genome1"]
    outdir = os.path.join(DATA_TEST_DIR, "test_sketch_all")
    os.makedirs(outdir)
    list_reps = os.path.join(outdir, "test_list_reps.txt")
    out_msh = os.path.join(outdir, "out_mash.msh")
    mash_log = os.path.join(outdir, "mash_sketch.log")
    threads = 1
    filterg.sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Check that expected output files were created
    assert os.path.isfile(list_reps)
    assert os.path.isfile(mash_log)
    assert os.path.isfile(out_msh)

    # Check content of list file
    expected_lines = [os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                      os.path.join(GENOMES_DIR, "ACOR003.0519.fna"),
                      os.path.join(GENOMES_DIR, "ACOR001.0519.fna")]
    # 3 files to sketch, with expected paths
    with open(list_reps, "r") as lr:
        lines_found = lr.readlines()
        assert len(lines_found) == 3
        for line, expect in zip(lines_found, expected_lines):
            assert line.strip() == expect

    with open(mash_log, "r") as ml:
        assert ml.readline().strip() == f"Sketching {expected_lines[0]}..."
        assert ml.readline().strip() == f"Sketching {expected_lines[1]}..."
        assert ml.readline().strip() == f"Sketching {expected_lines[2]}..."
        assert ml.readline().strip() == f"Writing to {out_msh}..."

    shutil.rmtree(outdir)


def test_sketch_all_noout(caplog):
    """
    Test that, when given output directory does not exist, it ends program with error message
    """
    genomes = {"genome1": ["g1_name", "g1_ori", os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                           123567, 200, 101],
               "genome2": ["g2_name", "g2_ori", os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                           20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOR003.0519.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome3", "genome1"]
    outdir = os.path.join(DATA_TEST_DIR, "test_sketch_all_noout")
    list_reps = os.path.join(outdir, "test_list_reps.txt")
    out_msh = os.path.join(outdir, "out_mash.msh")
    mash_log = os.path.join(outdir, "mash_sketch.log")
    threads = 1
    # Check everything works without error
    with pytest.raises(SystemExit):
        filterg.sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Your output directory 'test/data/prepare/test_sketch_all_noout' "
            "does not exist") in caplog.text


def test_sketch_all_mash_exists(caplog):
    """
    Test that, when mash sketch file already exists, it does nothing: message to say
    that the file already exists, and returns 0
    """
    genomes = {"genome1": ["g1_name", "g1_ori", os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                           123567, 200, 101],
               "genome2": ["g2_name", "g2_ori", os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                           20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOR003.0519.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome3", "genome1"]
    outdir = os.path.join(DATA_TEST_DIR, "test_sketch_all_mash_exists")
    os.makedirs(outdir)
    list_reps = os.path.join(outdir, "test_list_reps.txt")
    out_msh = os.path.join(DATA_TEST_DIR, "test_files", "test_out_mash")
    mash_log = os.path.join(outdir, "mash_sketch.log")
    threads = 1
    filterg.sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Check that expected output files were created
    assert not os.path.isfile(list_reps)
    assert not os.path.isfile(mash_log)
    assert os.path.isfile(out_msh + ".msh")

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Mash sketch file test/data/prepare/test_files/test_out_mash.msh already exists. PanACoTA "
            "will use it for next step.") in caplog.text

    shutil.rmtree(outdir)


def test_sketch_all_error_mash(caplog):
    """
    Test that, when mash has a problem, PanACoTA exits with an error message
    """
    genomes = {"genome1": ["g1_name", "g1_ori", os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                           123567, 200, 101],
               "genome2": ["g2_name", "g2_ori", os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                           20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOR003.0519.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome3", "genome1"]
    outdir = os.path.join(DATA_TEST_DIR, "test_sketch_all_mash_exists")
    os.makedirs(outdir)
    list_reps = os.path.join(outdir, "test_files", "test_list_reps.txt")
    out_msh = os.path.join(DATA_TEST_DIR, "test_files", "test_mash_sketch")
    open(out_msh + ".msh", "w").close()
    mash_log = os.path.join(outdir, "mash_sketch.log")
    threads = 1
    filterg.sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Check that expected output files were created
    assert not os.path.isfile(list_reps)
    assert not os.path.isfile(mash_log)
    assert os.path.isfile(out_msh + ".msh")

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Mash sketch file test/data/prepare/test_files/out_mash.msh already exists. PanACoTA "
            "will use it for next step.") in caplog.text

    shutil.rmtree(outdir)


# def test_read_matrix():
#     """
#     Test that matrix from minhash is correctly read and converted to python matrix
#     """
#     genomes = {"genome2": ["g2_name", "g2_ori", "path_genome2", 20000, 3, 1],
#                "genome6": ["g6_name", "g6_ori", "path_genome6", 1500, 5, 2],
#                "genome5": ["g5_name", "g5_ori", "path_genome5", 22012, 20, 10],
#                "genome3": ["g3_name", "g3_ori", "path_genome3", 25003, 52, 50]
#               }
#     sorted_genomes = ["genome2", "genome6", "genome5", "genome3"]
#     matrix_file = os.path.join(DATA_TEST_DIR, "test_files", "minhash_output.txt")

#     mat_sp = filterg.read_matrix(genomes, sorted_genomes, matrix_file)

#     # Check content of created matrix
#     assert len(mat_sp) == 6
#     assert mat_sp[0,1] == 1e-6
#     assert mat_sp[0,2] == 0.07
#     assert mat_sp[0,3] == 1e-3
#     assert mat_sp[1,2] == 0.08
#     assert mat_sp[1,3] == 1e-4
#     assert mat_sp[2,3] == 5e-2


# def test_read_matrix_no_genome():
#     """
#     Test that matrix from minhash is correctly read and converted to python matrix
#     """
#     genomes = {}
#     sorted_genomes = []
#     matrix_file = os.path.join(DATA_TEST_DIR, "test_files", "minhash_output_empty.txt")
#     # Create empty matrix file
#     open(matrix_file, "w").close()

#     # Read matrix
#     mat_sp = filterg.read_matrix(genomes, sorted_genomes, matrix_file)

#     # Check content of created matrix
#     assert len(mat_sp) == 0

#     # Remove empty matrix file
#     os.remove(matrix_file)

# def test_mash_step():
#     """
#     Prepare mash step: for a given genome, keep only genomes at a correct distance
#     (between min_dist and max_dist). Add others to 'genomes_removed'
#     """
#     # CREATE MATRIX
#     nbgen = 4
#     mat_sp = dok_matrix((nbgen, nbgen), dtype=float)
#     mat_sp[0,1] = 1e-6
#     mat_sp[0,2] = 0.07
#     mat_sp[0,3] = 1e-3
#     mat_sp[1,2] = 0.08
#     mat_sp[1,3] = 1e-4
#     mat_sp[2,3] = 5e-2

#     # PREPARE MASH STEP
#     # List of genomes ordered by quality. genome1 is the best genome, genome2 the worst
#     to_try = ["genome2", "genome4", "genome3", "genome1"]
#     corresp = {"genome1":0, "genome3":1, "genome4":2, "genome2":3}
#     genomes_removed = {}
#     min_dist = 1e-4
#     max_dit = 0.06

#     # Check everything goes well
#     assert filterg.mash_step(to_try, corresp, mat_sp, genomes_removed, min_dist, max_dit) == 0
#     assert to_try == ["genome4"]
#     assert genomes_removed == {"genome3": ["genome1", 0.07], "genome2": ["genome1", 1e-6]}



# def test_mash_step_wrong_sp():
#     """
#     Prepare mash step: for a given genome, keep only genomes at a correct distance
#     (between min_dist and max_dist). Add others to 'genomes_removed'
#     """
#     # CREATE MATRIX
#     nbgen = 4
#     mat_sp = dok_matrix((nbgen, nbgen), dtype=float)
#     mat_sp[0,1] = 1e-6
#     mat_sp[0,2] = 0.07
#     mat_sp[0,3] = 1e-3
#     mat_sp[1,2] = 0.08
#     mat_sp[3,1] = 1e-4
#     mat_sp[2,3] = 5e-2

#     # PREPARE MASH STEP
#     # List of genomes ordered by quality. genome1 is the best genome, genome2 the worst
#     to_try = ["genome2", "genome4", "genome3", "genome1"]
#     corresp = {"genome1":0, "genome3":1, "genome4":2, "genome2":3}
#     genomes_removed = {}
#     min_dist = 1e-4
#     max_dit = 0.06

#     # Check everything goes well
#     assert filterg.mash_step(to_try, corresp, mat_sp, genomes_removed, min_dist, max_dit) == 0
#     assert to_try == ["genome2"]
#     assert genomes_removed == {"genome3": ["genome1", 0.07], "genome4": ["genome1", 1e-6]}


