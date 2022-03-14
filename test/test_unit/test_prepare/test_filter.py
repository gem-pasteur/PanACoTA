#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the download_genomes_func submodule in prepare module
"""
import os
import logging
import shutil
import pytest
from scipy.sparse import dok_matrix

import test.test_unit.utilities_for_tests as tutil
import PanACoTA.prepare_module.filter_genomes as filterg
import PanACoTA.utils as utils

# Define variables and functions used by several tests
DATA_TEST_DIR = os.path.join("test", "data", "prepare")
GENOMES_DIR = os.path.join(DATA_TEST_DIR, "genomes", "genomes_comparison")
GENEPATH = os.path.join(DATA_TEST_DIR, "generated_by_unit-tests")
EXP_GENOMES = {
               "ACOR001.0519.fna": ["ACOR001.0519",
                                    os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                                    os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                                    3013644, 269, 37],
               "ACOR001.0519-bis.fna": ["ACOR001.0519-bis",
                                        os.path.join(GENOMES_DIR, "ACOR001.0519-bis.fna"),
                                        os.path.join(GENOMES_DIR, "ACOR001.0519-bis.fna"),
                                        3013644, 269, 37],
               "ACOR001.0519-almost-same.fna": ["ACOR001.0519-almost-same",
                                        os.path.join(GENOMES_DIR, "ACOR001.0519-almost-same.fna"),
                                        os.path.join(GENOMES_DIR, "ACOR001.0519-almost-same.fna"),
                                        3012665, 261, 37],
               "ACOR002.0519.fna": ["ACOR002.0519",
                                    os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                                    os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                                    2997537, 78, 23],
               "ACOC.1019.fna": ["ACOC.1019", os.path.join(GENOMES_DIR, "ACOC.1019.fna"),
                                  os.path.join(GENOMES_DIR, "ACOC.1019.fna"),
                                  1587120, 1, 1]
               }

LOGFILE_BASE = os.path.join(GENEPATH, "logfile_test.txt")
LEVEL = logging.DEBUG
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]


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
    # utils.init_logger(LOGFILE_BASE, LEVEL, 'test_filter', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    shutil.rmtree(GENEPATH)
    print("teardown")


def test_write_output_no_outdir(caplog):
    """
    Check that when outdir does not exist, program ends with error message
    """
    corresp_genomes = {"ACOR001": "ACOR001.0519.fna.gz", "ACOR002": "ACOR002.0519.fna.gz",
                       "ACOR003": "ACOR003.0519.fna.gz"}
    sorted_genomes = []
    genomes_removed = {}

    # Define output directory for generated files
    outdir = os.path.join(GENEPATH, "test_filter_write_output_no_outdir")
    gspecies = "Acetobacter_fabarum"

    # Check everything works without error
    with pytest.raises(SystemExit):
        filterg.write_outputfiles({}, sorted_genomes, genomes_removed, outdir, gspecies, 1e-4, 0.06)

    # Check outfiles do not exist
    list_file = os.path.join(outdir, "LSTINFO-Acetobacter_fabarum-filtered-0.06.txt")
    discard_file = os.path.join(outdir, "discarded-by-minhash-Acetobacter_fabarum-0.06.txt")
    assert not os.path.isfile(list_file)
    assert not os.path.isfile(discard_file)

    # Check output
    assert not os.path.isdir(outdir)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert ("The given output directory (test/data/prepare/generated_by_unit-tests/test_filter_write_output_no_outdir) "
            "does not exist. We cannot create output files there") in caplog.text


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
    outdir = os.path.join(GENEPATH, "test_filter_write_output")
    os.mkdir(outdir)
    gspecies = "Acetobacter_fabarum"
    min_dist = 1e-4
    max_dist = 0.06

    # Check everything works without error
    list_file = os.path.join(outdir, "LSTINFO-Acetobacter_fabarum-filtered-0.0001_0.06.txt")
    assert filterg.write_outputfiles(genomes, sorted_genomes,
                                     genomes_removed, outdir, gspecies,
                                     min_dist, max_dist) == list_file

    # Check outfiles exist
    discard_file = os.path.join(outdir, "discarded-by-minhash-Acetobacter_fabarum-0.0001_0.06.txt")
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
    outdir = os.path.join(GENEPATH, "test_filter_write_output_no_discard")
    os.makedirs(outdir)
    gspecies = "Acetobacter_fabarum"
    min_dist = 1e-4
    max_dist = 0.06

    # Check everything works without error
    list_file = os.path.join(outdir, "LSTINFO-Acetobacter_fabarum-filtered-0.0001_0.06.txt")
    assert filterg.write_outputfiles(genomes, sorted_genomes,
                                     genomes_removed, outdir, gspecies,
                                     min_dist, max_dist) == list_file

    # Check outfiles exist
    discard_file = os.path.join(outdir, "discarded-by-minhash-Acetobacter_fabarum-0.0001_0.06.txt")
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
    outdir = os.path.join(GENEPATH, "test_filter_write_output_no_genome")
    os.makedirs(outdir)
    gspecies = "Acetobacter_fabarum"
    min_dist = 1e-4
    max_dist = 0.06

    # Check everything works without error
    list_file = os.path.join(outdir, "LSTINFO-Acetobacter_fabarum-filtered-0.0001_0.06.txt")
    assert filterg.write_outputfiles(genomes, sorted_genomes,
                                     genomes_removed, outdir, gspecies,
                                     min_dist, max_dist) == list_file

    # Check outfiles exist
    discard_file = os.path.join(outdir, "discarded-by-minhash-Acetobacter_fabarum-0.0001_0.06.txt")
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


def test_sketch_all(caplog):
    """
    Test that all genomes are sketch, in the provided order
    """
    utils.init_logger(LOGFILE_BASE, LEVEL, 'test_filter', verbose=1)
    caplog.set_level(logging.DEBUG)
    # We give 5 genomes
    genomes = {"genome1": ["g1_name", "g1_ori", os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                           123567, 200, 101],
               "genome1bis": ["g1bis_name", "g1bis_ori",
                              os.path.join(GENOMES_DIR, "ACOR001.0519-bis.fna"),
                              251500, 200, 101],
               "genome1diff": ["g1diff_name", "g1diff_ori",
                              os.path.join(GENOMES_DIR, "ACOR001.0519-almost-same.fna"),
                              1500, 3, 2],
               "genome2": ["g2_name", "g2_ori", os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                           20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOC.1019.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome1diff", "genome3", "genome1", "genome1bis"]
    outdir = os.path.join(GENEPATH, "test_sketch_all")
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
                      os.path.join(GENOMES_DIR, "ACOR001.0519-almost-same.fna"),
                      os.path.join(GENOMES_DIR, "ACOC.1019.fna"),
                      os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                      os.path.join(GENOMES_DIR, "ACOR001.0519-bis.fna")]
    assert tutil.compare_file_to_list(list_reps, expected_lines)

    with open(mash_log, "r") as ml:
        assert ml.readline().strip() == f"Sketching {expected_lines[0]}..."
        assert ml.readline().strip() == f"Sketching {expected_lines[1]}..."
        assert ml.readline().strip() == f"Sketching {expected_lines[2]}..."
        assert ml.readline().strip() == f"Sketching {expected_lines[3]}..."
        assert ml.readline().strip() == f"Sketching {expected_lines[4]}..."
        assert ml.readline().strip() == f"Writing to {out_msh}..."

    assert "Sketching all genomes..." in caplog.text
    assert ("mash sketch -o test/data/prepare/generated_by_unit-tests/test_sketch_all/out_mash.msh -p 1 -l "
            "test/data/prepare/generated_by_unit-tests/test_sketch_all/test_list_reps.txt -s 1e4") in caplog.text
    assert caplog.records[0].levelname == "INFO"
    assert caplog.records[1].levelname == "DETAIL"


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
    outdir = os.path.join(GENEPATH, "test_sketch_all_noout")
    list_reps = os.path.join(outdir, "test_list_reps.txt")
    out_msh = os.path.join(outdir, "out_mash.msh")
    mash_log = os.path.join(outdir, "mash_sketch.log")
    threads = 1
    # Check everything works without error
    with pytest.raises(SystemExit):
        filterg.sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Your output directory 'test/data/prepare/generated_by_unit-tests/test_sketch_all_noout' "
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
    outdir = os.path.join(GENEPATH, "test_sketch_all_mash_exists")
    os.makedirs(outdir)
    list_reps = os.path.join(outdir, "test_list_reps.txt")
    out_msh = os.path.join(outdir, "mash_exists")
    # Create empty msh file
    open(out_msh + ".msh", "w").close()
    mash_log = os.path.join(outdir, "mash_sketch.log")
    threads = 1
    # Check that out_msh already exists before running mash sketch
    assert os.path.isfile(out_msh + ".msh")
    filterg.sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Check that mash file still exists, but no other outfile was created
    assert not os.path.isfile(list_reps)
    assert not os.path.isfile(mash_log)
    assert os.path.isfile(out_msh + ".msh")

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Mash sketch file test/data/prepare/generated_by_unit-tests/test_sketch_all_mash_exists/mash_exists.msh "
            "already exists. PanACoTA will use it for next step.") in caplog.text


def test_sketch_all_error_mash(caplog):
    """
    Test that, when mash has a problem, PanACoTA exits with an error message
    """
    genomes = {"genome1": ["g1_name", "g1_ori", "genome1", 123567, 200, 101],
               "genome2": ["g2_name", "g2_ori", "genome2", 20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOR003.0519.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome3", "genome1"]
    outdir = os.path.join(GENEPATH, "test_sketch_all_mash_error")
    os.makedirs(outdir)
    list_reps = os.path.join(DATA_TEST_DIR, "test_files", "test_list_to_sketch.txt")
    out_msh = os.path.join(outdir, "out_mash.msh")
    mash_log = os.path.join(outdir, "mash_sketch.log")
    threads = 1

    # Test that it exists with sysExit error
    with pytest.raises(SystemExit):
        filterg.sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Check that expected output files were created
    assert os.path.isfile(list_reps)
    assert os.path.isfile(mash_log)
    assert not os.path.isfile(out_msh)

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Error while trying to sketch 3 genomes to combined archive. Maybe some genome "
            "sequences in 'tmp_files' are missing! Check logfile: test/data/prepare/generated_by_unit-tests/"
            "test_sketch_all_mash_error/mash_sketch.log") in caplog.text


def test_compare_all(caplog):
    """
    Check that comparison of all sketched sequences is as expected (output matrix is as expected)
    """
    out_msh = os.path.join(DATA_TEST_DIR, "test_files", "test_mash_output")
    matrix = os.path.join(GENEPATH, "matrix_from_test_compare_all.txt")
    mash_log = os.path.join(GENEPATH, "mashlog_from_test_compare_all.log")
    threads = 1

    # Check msh file exists
    assert os.path.isfile(out_msh + ".msh")

    filterg.compare_all(out_msh, matrix, "matrix", mash_log, threads)

    # Check output files are created
    assert os.path.isfile(matrix)
    assert os.path.isfile(mash_log)

    # Check content of matrix file
    expect_matrix = os.path.join(DATA_TEST_DIR, "test_files", "test_matrix_mash.txt")
    assert tutil.compare_file_content(matrix, expect_matrix)


def test_compare_all_matrix_exists(caplog):
    """
    Check that matrix file already exists, it returns 0 with a warning message sayng that the
    existing file will be used.
    """
    out_msh = os.path.join(DATA_TEST_DIR, "test_files", "test_mash_output")
    matrix = os.path.join(DATA_TEST_DIR, "test_files", "test_matrix_mash.txt")
    mash_log = os.path.join(GENEPATH, "mashlog_from_test_compare_all.log")
    threads = 1

    filterg.compare_all(out_msh, matrix, "matrix", mash_log, threads)

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Matrix file test/data/prepare/test_files/test_matrix_mash.txt already exists. "
            "The program will use this distance matrix to filter all genomes according to "
            "their distances") in caplog.text


def test_compare_all_error_mash(caplog):
    """
    Check that when mash has a problem, it gives an error message and closes the program
    """
    # mash file does not exist
    out_msh = os.path.join(GENEPATH, "mash.msh")
    matrix = os.path.join(GENEPATH, "matrix.txt")
    mash_log = os.path.join(GENEPATH, "mashlog_from_test_compare_all-error-mash.log")
    threads = 1

    # Test that it exists with sysExit error
    with pytest.raises(SystemExit):
        filterg.compare_all(out_msh, matrix, "matrix", mash_log, threads)

    # Check log
    caplog.set_level(logging.DEBUG)
    assert ("Error while trying to estimate pairwise distances between all genomes. "
            "See test/data/prepare/generated_by_unit-tests/mashlog_from_test_compare_all-error-mash.log") in caplog.text


def test_read_matrix():
    """
    Test that the matrix txt file is converted to a scipy matrix as expected
    """
    genomes = {"genome1": ["g1_name", "g1_ori", os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                           123567, 200, 101],
               "genome1bis": ["g1bis_name", "g1bis_ori",
                              os.path.join(GENOMES_DIR, "ACOR001.0519-bis.fna"),
                              251500, 200, 101],
               "genome1diff": ["g1diff_name", "g1diff_ori",
                              os.path.join(GENOMES_DIR, "ACOR001.0519-almost-same.fna"),
                              1500, 3, 2],
               "genome2": ["g2_name", "g2_ori", os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                           20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOC.1019.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome1diff", "genome3", "genome1", "genome1bis"]
    matrix_file = os.path.join(DATA_TEST_DIR, "test_files", "test_matrix_mash.txt")

    # Run matrix reading
    out_mat = filterg.read_matrix(genomes, sorted_genomes, matrix_file)

    exp_mat = dok_matrix((5, 5), dtype=float)
    exp_mat[0, 0] = 0  # genome2 vs genome2
    exp_mat[0, 1] = 0.000167546  # genome2 vs genome1diff
    exp_mat[0, 2] = 0.295981  # genome2 vs genome3
    exp_mat[0, 3] = 0.000143503  # genome2 vs genome1
    exp_mat[0, 4] = 0.000143503  # genome2 vs genome1bis
    exp_mat[1, 1] = 0 # genome1diff vs genome1diff
    exp_mat[1, 2] = 0.295981 # genome1diff vs genome3
    exp_mat[1, 3] = 2.38274e-05  # genome1diff vs genome1
    exp_mat[1, 4] = 2.38274e-05  # genome1diff vs genome1bis
    exp_mat[2, 3] = 0.295981  # genome3 vs genome1
    exp_mat[2, 4] = 0.295981  # genome3 vs genome1bis
    exp_mat[3, 4] = 0  # genome1 vs genome1bis

    out = out_mat.toarray()
    exp = exp_mat.toarray()
    import numpy as np
    assert  np.array_equal(out, exp)


def test_read_matrix_nofile(caplog):
    """
    Test that when the given matrix file does not exist, it exits with error message
    """
    genomes = {"genome1": ["g1_name", "g1_ori", os.path.join(GENOMES_DIR, "ACOR001.0519.fna"),
                           123567, 200, 101],
               "genome1bis": ["g1bis_name", "g1bis_ori",
                              os.path.join(GENOMES_DIR, "ACOR001.0519-bis.fna"),
                              251500, 200, 101],
               "genome1diff": ["g1diff_name", "g1diff_ori",
                              os.path.join(GENOMES_DIR, "ACOR001.0519-almost-same.fna"),
                              1500, 3, 2],
               "genome2": ["g2_name", "g2_ori", os.path.join(GENOMES_DIR, "ACOR002.0519.fna"),
                           20000, 3, 1],
               "genome3": ["g3_name", "g3_ori", os.path.join(GENOMES_DIR, "ACOC.1019.fna"), 25003, 52, 50]
               }
    sorted_genomes = ["genome2", "genome1diff", "genome3", "genome1", "genome1bis"]
    matrix_file = os.path.join(DATA_TEST_DIR, "mymatrix.txt")

    # Test that it exists with sysExit error
    with pytest.raises(SystemExit):
        out_mat = filterg.read_matrix(genomes, sorted_genomes, matrix_file)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert("Matrix file test/data/prepare/mymatrix.txt does not exist. We cannot "
           "read it and do the next steps. Program ending.") in caplog.text


def test_mash_step_2steps():
    """
    Test that when comparing a given reference genome to all others kept until now:
    - it updates 'genomes_removed' (genomes removed as they are not between min_dist
    and max_dist compared to the given reference genome.
    - it updates 'to_try' object, removing the given reference genome, as well as discarded
    genomes
    """
    to_try = ["genome1bis", "genome1", "genome3", "genome1diff", "genome2"]
    corresp = {"genome1":3, "genome1bis":4, "genome1diff":1, "genome2":0, "genome3":2}

    # Create matrix to read
    exp_mat = dok_matrix((5, 5), dtype=float)
    exp_mat[0, 0] = 0  # genome2 vs genome2
    exp_mat[0, 1] = 0.000167546  # genome2 vs genome1diff
    exp_mat[0, 2] = 0.295981  # genome2 vs genome3
    exp_mat[0, 3] = 0.000143503  # genome2 vs genome1
    exp_mat[0, 4] = 0.000143503  # genome2 vs genome1bis
    exp_mat[1, 1] = 0 # genome1diff vs genome1diff
    exp_mat[1, 2] = 0.295981 # genome1diff vs genome3
    exp_mat[1, 3] = 2.38274e-05  # genome1diff vs genome1
    exp_mat[1, 4] = 2.38274e-05  # genome1diff vs genome1bis
    exp_mat[2, 3] = 0.295981  # genome3 vs genome1
    exp_mat[2, 4] = 0.295981  # genome3 vs genome1bis
    exp_mat[3, 4] = 0  # genome1 vs genome1bis

    genomes_removed = {"toto": ["titi", 0]}
    min_dist = 1e-4
    max_dist = 0.06

    # First step: compare genome2 to all other genomes
    filterg.mash_step(to_try, corresp, exp_mat, genomes_removed, min_dist, max_dist)

    exp_removed = {"toto": ["titi", 0],
                   "genome3": ["genome2", 0.295981]}

    exp_to_try = ["genome1bis", "genome1", "genome1diff"]
    assert genomes_removed == exp_removed
    assert to_try == exp_to_try

    # Second step: compare genome1diff to all other genomes
    filterg.mash_step(exp_to_try, corresp, exp_mat, exp_removed, min_dist, max_dist)

    out_removed = {"toto": ["titi", 0],
                   "genome1": ["genome1diff", 2.38274e-05],
                   "genome1bis": ["genome1diff", 2.38274e-05],
                   "genome3": ["genome2", 0.295981]
                   }

    out_to_try = []

    assert exp_removed == out_removed
    assert exp_to_try == out_to_try


def test_mash_step_wrong_mat(caplog):
    """
    Test that when comparing a given reference genome to all others kept until now:
    - it updates 'genomes_removed' (genomes removed as they are not between min_dist
    and max_dist compared to the given reference genome.
    - it updates 'to_try' object, removing the given reference genome, as well as discarded
    genomes
    BUT, as the matrix is not triangular, prints a warning
    """
    to_try = ["genome1bis", "genome1", "genome3", "genome2", "genome1diff", "genome2"]
    corresp = {"genome1":3, "genome1bis":4, "genome1diff":1, "genome2":0, "genome3":2}

    # Create matrix to read
    exp_mat = dok_matrix((5, 5), dtype=float)
    exp_mat[0, 0] = 0  # genome2 vs genome2
    exp_mat[0, 1] = 0.000167546  # genome2 vs genome1diff
    exp_mat[0, 2] = 0.295981  # genome2 vs genome3
    exp_mat[0, 3] = 0.000143503  # genome2 vs genome1
    exp_mat[0, 4] = 0.000143503  # genome2 vs genome1bis
    exp_mat[1, 1] = 0 # genome1diff vs genome1diff
    exp_mat[1, 2] = 0.295981 # genome1diff vs genome3
    exp_mat[1, 3] = 2.38274e-05  # genome1diff vs genome1
    exp_mat[1, 4] = 2.38274e-05  # genome1diff vs genome1bis
    exp_mat[2, 3] = 0.295981  # genome3 vs genome1
    exp_mat[2, 4] = 0.295981  # genome3 vs genome1bis
    exp_mat[3, 4] = 0  # genome1 vs genome1bis

    genomes_removed = {"toto": ["titi", 0]}
    min_dist = 1e-4
    max_dist = 0.06

    # Compare genome2 to all other genomes
    filterg.mash_step(to_try, corresp, exp_mat, genomes_removed, min_dist, max_dist)

    exp_removed = {"toto": ["titi", 0],
                   "genome3": ["genome2", 0.295981],
                   "genome2": ["genome2", 0]
                   }
    exp_to_try = ["genome1bis", "genome1", "genome1diff"]
    assert genomes_removed == exp_removed
    assert to_try == exp_to_try

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert "Should never happen as mat_sp is a triangle matrix!" in caplog.text


def test_check_quality():
    """
    quality control of all genomes in the database
    """
    species_linked = "my-test-genomes"
    db_path = os.path.join(DATA_TEST_DIR, "genomes", "genomes_comparison")
    tmp_dir = os.path.join(GENEPATH, "tmp_dir_check_quality")
    os.mkdir(tmp_dir)
    max_l90 = 100
    max_cont = 100
    cutn = 0

    genomes = filterg.check_quality(species_linked, db_path, tmp_dir, max_l90, max_cont, cutn)

    assert genomes == EXP_GENOMES


def test_check_quality_no_dbdir(caplog):
    """
    quality control of all genomes in the database when given db folder does not exist:
    ends with error message
    """
    species_linked = "my-test-genomes"
    db_path = "dbpath_for_test"
    tmp_dir = "tmp_dir_check_quality"
    max_l90 = 100
    max_cont = 100
    cutn = 5

    # Test that it exists with sysExit error
    with pytest.raises(SystemExit):
        genomes = filterg.check_quality(species_linked, db_path, tmp_dir, max_l90, max_cont, cutn)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert "dbpath_for_test does not exist" in caplog.text


def test_check_quality_no_tmpdir(caplog):
    """
    quality control of all genomes in the database when given tmp folder does not exist:
    ends with error message
    """
    species_linked = "my-test-genomes"
    db_path = os.path.join(DATA_TEST_DIR, "genomes", "genomes_comparison")
    tmp_dir = "tmp_dir_check_quality"
    max_l90 = 100
    max_cont = 100
    cutn = 5

    # Test that it exists with sysExit error
    with pytest.raises(SystemExit):
        genomes = filterg.check_quality(species_linked, db_path, tmp_dir, max_l90, max_cont, cutn)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert "tmp_dir_check_quality does not exist" in caplog.text


def test_check_quality_same_contname(caplog):
    """
    quality control of all genomes in the database when one genome contains 2 contigs
    with the same name: this genome is discarded + warning message
    """
    caplog.set_level(logging.DEBUG)
    species_linked = "my-test-genomes"
    db_path = os.path.join(DATA_TEST_DIR, "genomes", "genomes_comparison")
    tmp_dir = os.path.join(GENEPATH, "tmp_dir_check_quality")
    os.makedirs(tmp_dir)
    # Copy all genomes to new database
    db_path_used = os.path.join(GENEPATH, "database")
    shutil.copytree(db_path, db_path_used)
    # Add another genome, with 2 identical contig names
    new_genome = os.path.join(db_path_used, "genome_duplicate.fst")
    with open(new_genome, "w") as ng:
        ng.write(">my_contig\nAACTATATAGGAAGACACACAATTAAGGGACAGG\n")
        ng.write(">my_contig2\nCCNNCCGATTCGAGCACACAATTAAGGGACAGG\n")
        ng.write(">my_contig\nCCANNCCATCTCTCTTATCTCTCTAANNNNCTCTANCCNNNNNNCCATTCA\n")
    max_l90 = 100
    max_cont = 100
    cutn = 0

    # Get genomes information
    genomes = filterg.check_quality(species_linked, db_path_used, tmp_dir, max_l90, max_cont, cutn)

    # Even if there was 1 more genome in the database, output list of genomes is the same,
    # as the supplementary genome contains identical contig names
    exp_genomes = {
               "ACOR001.0519.fna": ["ACOR001.0519",
                                    os.path.join(db_path_used, "ACOR001.0519.fna"),
                                    os.path.join(db_path_used, "ACOR001.0519.fna"),
                                    3013644, 269, 37],
               "ACOR001.0519-bis.fna": ["ACOR001.0519-bis",
                                        os.path.join(db_path_used, "ACOR001.0519-bis.fna"),
                                        os.path.join(db_path_used, "ACOR001.0519-bis.fna"),
                                        3013644, 269, 37],
               "ACOR001.0519-almost-same.fna": ["ACOR001.0519-almost-same",
                                        os.path.join(db_path_used, "ACOR001.0519-almost-same.fna"),
                                        os.path.join(db_path_used, "ACOR001.0519-almost-same.fna"),
                                        3012665, 261, 37],
               "ACOR002.0519.fna": ["ACOR002.0519",
                                    os.path.join(db_path_used, "ACOR002.0519.fna"),
                                    os.path.join(db_path_used, "ACOR002.0519.fna"),
                                    2997537, 78, 23],
               "ACOC.1019.fna": ["ACOC.1019", os.path.join(db_path_used, "ACOC.1019.fna"),
                                  os.path.join(db_path_used, "ACOC.1019.fna"),
                                  1587120, 1, 1]
               }

    assert genomes == exp_genomes

    # Check logs
    assert ("Total number of genomes for my-test-genomes: 6") in caplog.text
    assert ("In genome genome_duplicate.fst, '>my_contig' contig name is used for several "
            "contigs. Please put different names for each contig. This genome will be "
            "ignored.") in caplog.text


def test_check_quality_no_genome(caplog):
    """
    quality control of all genomes in the database where there is no genome to annotate:
    ends with error message
    """
    species_linked = "my-test-genomes"
    db_path = os.path.join(GENEPATH, "dbpath_for_test")
    os.mkdir(db_path)
    tmp_dir = os.path.join(GENEPATH, "tmp_dir_check_quality")
    os.mkdir(tmp_dir)
    max_l90 = 100
    max_cont = 100
    cutn = 5

    # Test that it exists with sysExit error
    with pytest.raises(SystemExit):
        genomes = filterg.check_quality(species_linked, db_path, tmp_dir, max_l90, max_cont, cutn)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert "There is no genome in test/data/prepare/generated_by_unit-tests/dbpath_for_test." in caplog.text


def test_iterative_mash(caplog):
    """
    Test that when we give all genomes, sorted, it compares them, and it saves the comparison
    matrix in a npz file, and return genomes removed
    """
    # Parameters for iterative_mash function
    sorted_genomes = ["ACOR002.0519.fna", "ACOR001.0519-almost-same.fna",
                      "ACOC.1019.fna", "ACOR001.0519.fna", "ACOR001.0519-bis.fna"]
    outdir = os.path.join(GENEPATH, "res_test_iterative_mash")
    species_linked = "my-test-species"
    min_dist = 1e-4
    max_dist = 0.06
    threads = 1
    quiet = False

    # Expected created files
    mash_dir = os.path.join(outdir, "mash_files")
    txt_matrix = os.path.join(mash_dir, "matrix-all-genomes-my-test-species.txt")
    npz_mat = os.path.join(mash_dir, "matrix-all-genomes-my-test-species.npz")

    removed = filterg.iterative_mash(sorted_genomes, EXP_GENOMES, outdir,
                                     species_linked, min_dist, max_dist, threads, quiet)

    # At least, ACOC.1019.fna should be remove by ACOR002.0519.fna, as we made them
    # completely different. For the others, depends on mash version...
    # mash1.1 get correct difference between ACOR002 and ACOR001
    # but mash2.2 considers them as too close
    assert "ACOC.1019.fna" in removed.keys()
    assert removed["ACOC.1019.fna"][0] == "ACOR002.0519.fna"

    # Check that txt and npz matrix were created
    # We cannot check their content as distances depend on mash version...
    assert os.path.isfile(npz_mat)
    assert os.path.isfile(txt_matrix)


def test_iterative_mash_npz_exists():
    """
    Test that when we give all genomes, sorted, and that the mash matrix was already calculated,
    and is now stored in an npz file, it does not re-calculate it, but returns directly removed
    genomes.
    """
    sorted_genomes = ["ACOR002.0519.fna", "ACOR001.0519-almost-same.fna",
                      "ACOC.1019.fna", "ACOR001.0519.fna", "ACOR001.0519-bis.fna"]
    # Create output dir where all mash result files will be stored
    outdir = os.path.join(GENEPATH, "res_test_iterative_mash_npz_exists")
    mash_dir = os.path.join(outdir, "mash_files")
    os.makedirs(mash_dir)
    # Copy existing npz matrix to mash output folder
    npz_matrix_model = os.path.join(DATA_TEST_DIR, "test_files",
                                    "test_npz_matrix_mash.npz")
    npz_matrix_out = os.path.join(mash_dir, "matrix-all-genomes-my-test-species.npz")
    shutil.copy(npz_matrix_model, npz_matrix_out)
    txt_matrix_out = os.path.join(mash_dir, "matrix-all-genomes-my-test-species.txt")
    species_linked = "my-test-species"
    min_dist = 1e-4
    max_dist = 0.06
    threads = 1
    quiet = False

    # Check npz matrix already exists, and txt matrix does not
    assert os.path.isfile(npz_matrix_out)
    assert not os.path.isfile(txt_matrix_out)

    removed = filterg.iterative_mash(sorted_genomes, EXP_GENOMES, outdir,
                                     species_linked, min_dist, max_dist, threads, quiet)

    # Compare output dict
    exp_removed = {"ACOC.1019.fna": ["ACOR002.0519.fna", 0.295981],
                   "ACOR001.0519-bis.fna": ["ACOR001.0519-almost-same.fna", 2.38274e-05],
                   "ACOR001.0519.fna": ["ACOR001.0519-almost-same.fna", 2.38274e-05]}

    assert removed == exp_removed

    # Check npz matrix exists (and is the same as the given file), and txt matrix was not created,
    # as we used the npz file
    assert os.path.isfile(npz_matrix_out)
    assert tutil.compare_files_bin(npz_matrix_out, npz_matrix_model)
    assert not os.path.isfile(txt_matrix_out)
