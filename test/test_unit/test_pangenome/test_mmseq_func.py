#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the mmseqs_functions submodule in pangenome module
"""

import os
import time
import shutil
import glob
import logging
import pytest

import PanACoTA.pangenome_module.mmseqs_functions as mmseqs
import PanACoTA.utils as utils
import test.test_unit.utilities_for_tests as tutil

LOGFILE_BASE = "logfile_test.txt"
LEVEL = logging.DEBUG
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]
# Define variables shared by several tests
PANDIR = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PANDIR, "test_files")
PATH_EXP_FILES = os.path.join(PANDIR, "exp_files")
GENEPATH = os.path.join(PANDIR, "generated_by_unit-tests")


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
    utils.init_logger(LOGFILE_BASE, logging.DEBUG, 'test_getseq', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    print("teardown")

# Clusters expected for the given databank (ref: [members]
EXP_CLUSTERS = {"GEN2.1017.00001.i0002_00004": ["GEN2.1017.00001.i0002_00004",
                                                "GEN4.1111.00001.i0001_00002",
                                                "GENO.1017.00001.b0002_00003",
                                                "GENO.1216.00002.i0001_00003"],
                "GEN2.1017.00001.b0003_00010": ["GEN2.1017.00001.b0003_00010"],
                "GEN2.1017.00001.b0004_00013": ["GEN2.1017.00001.b0004_00013"],
                "GEN4.1111.00001.b0001_00001": ["GEN2.1017.00001.i0002_00005",
                                                "GEN4.1111.00001.b0001_00001",
                                                "GENO.1017.00001.b0001_00002",
                                                "GENO.1216.00002.b0001_00001",
                                                "GENO.1216.00002.i0001_00002"],
                "GEN4.1111.00001.i0001_00003": ["GEN4.1111.00001.i0001_00003"],
                "GEN4.1111.00001.b0001_00009": ["GEN2.1017.00001.b0004_00011",
                                                "GEN4.1111.00001.b0001_00009",
                                                "GENO.1017.00001.b0002_00011",
                                                "GENO.1216.00002.b0002_00010"],
                "GENO.1017.00001.b0001_00001": ["GEN2.1017.00001.b0002_00006",
                                                "GENO.1017.00001.b0001_00001"],
                "GENO.1017.00001.i0002_00006": ["GEN2.1017.00001.b0003_00007",
                                                "GEN4.1111.00001.i0001_00006",
                                                "GENO.1017.00001.i0002_00006",
                                                "GENO.1017.00001.i0002_00007",
                                                "GENO.1216.00002.i0001_00007"],
                "GENO.1017.00001.i0002_00008": ["GEN2.1017.00001.i0004_00012",
                                                "GENO.1017.00001.i0002_00008"],
                "GENO.1017.00001.i0002_00009": ["GEN2.1017.00001.i0003_00008",
                                                "GEN4.1111.00001.i0001_00007",
                                                "GENO.1017.00001.i0002_00009",
                                                "GENO.1216.00002.b0001_00008"],
                "GENO.1017.00001.i0002_00010": ["GEN2.1017.00001.i0003_00009",
                                                "GEN4.1111.00001.i0001_00008",
                                                "GENO.1017.00001.i0002_00010",
                                                "GENO.1216.00002.b0002_00009"],
                "GENO.1216.00002.i0001_00004": ["GEN2.1017.00001.b0002_00003",
                                                "GENO.1216.00002.i0001_00004"],
                "GENO.1216.00002.i0001_00005": ["GEN2.1017.00001.b0001_00002",
                                                "GEN4.1111.00001.i0001_00004",
                                                "GENO.1017.00001.i0002_00004",
                                                "GENO.1216.00002.i0001_00005"],
                "GENO.1216.00002.i0001_00006": ["GEN2.1017.00001.b0001_00001",
                                                "GEN4.1111.00001.i0001_00005",
                                                "GENO.1017.00001.i0002_00005",
                                                "GENO.1216.00002.i0001_00006"],
                "GENO.1216.00002.b0003_00011": ["GENO.1216.00002.b0003_00011"],
                "GENO.1216.00002.b0003_00012": ["GENO.1216.00002.b0003_00012"]
                }

# protein families expected
FAMILIES4G = [["GEN2.1017.00001.i0002_00004", "GEN4.1111.00001.i0001_00002",
               "GENO.1017.00001.b0002_00003", "GENO.1216.00002.i0001_00003"],
              ["GEN2.1017.00001.b0003_00010"],
              ["GEN2.1017.00001.b0004_00013"],
              ["GEN2.1017.00001.i0002_00005", "GEN4.1111.00001.b0001_00001",
               "GENO.1017.00001.b0001_00002", "GENO.1216.00002.b0001_00001",
               "GENO.1216.00002.i0001_00002"],
              ["GEN4.1111.00001.i0001_00003"],
              ["GEN2.1017.00001.b0004_00011", "GEN4.1111.00001.b0001_00009",
               "GENO.1017.00001.b0002_00011", "GENO.1216.00002.b0002_00010"],
              ["GEN2.1017.00001.b0002_00006", "GENO.1017.00001.b0001_00001"],
              ["GEN2.1017.00001.b0003_00007", "GEN4.1111.00001.i0001_00006",
               "GENO.1017.00001.i0002_00006", "GENO.1017.00001.i0002_00007",
               "GENO.1216.00002.i0001_00007"],
              ["GEN2.1017.00001.i0004_00012", "GENO.1017.00001.i0002_00008"],
              ["GEN2.1017.00001.i0003_00008", "GEN4.1111.00001.i0001_00007",
               "GENO.1017.00001.i0002_00009", "GENO.1216.00002.b0001_00008"],
              ["GEN2.1017.00001.i0003_00009", "GEN4.1111.00001.i0001_00008",
               "GENO.1017.00001.i0002_00010", "GENO.1216.00002.b0002_00009"],
              ["GEN2.1017.00001.b0002_00003", "GENO.1216.00002.i0001_00004"],
              ["GEN2.1017.00001.b0001_00002", "GEN4.1111.00001.i0001_00004",
               "GENO.1017.00001.i0002_00004", "GENO.1216.00002.i0001_00005"],
              ["GEN2.1017.00001.b0001_00001", "GEN4.1111.00001.i0001_00005",
               "GENO.1017.00001.i0002_00005", "GENO.1216.00002.i0001_00006"],
              ["GENO.1216.00002.b0003_00011"],
              ["GENO.1216.00002.b0003_00012"]
              ]


# Start tests
def test_create_mmseqdb(caplog):
    """
    Test that mmseq DB is created. We do not check its content as it could change
    according to mmseq versions, and we are here testing PanACoTA, not mmseqs
    Just check that all expected outfiles are present.
    """
    caplog.set_level(logging.DEBUG)
    filename = os.path.join(GENEPATH, "test_create_mmseqsdb.msdb")
    prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    logfile = os.path.join(GENEPATH, "test_create_mmseqsdb.log")
    mmseqs.create_mmseqs_db(filename, prt_path, logfile)

    outext = ["", ".index", ".lookup", "_h", "_h.index", ".dbtype", "_h.dbtype"]
    for file in [filename + ext for ext in outext]:
        assert os.path.isfile(file)
    assert "Creating database" in caplog.text
    assert "Existing files: 0" in caplog.text
    assert "Expected extensions: 7" in caplog.text
    assert caplog.records[0].levelname == "INFO"
    assert caplog.records[1].levelname == "DEBUG"
    assert caplog.records[2].levelname == "DEBUG"
    assert os.path.isfile(logfile)


def test_create_mmseqdb_exist(caplog):
    """
    Check that, when trying to create mmseqdb while all output files already exist,
    it logs a warning message and quits without creating it
    """
    filename = os.path.join(GENEPATH, "test_create_mmseqsdb_exist.msdb")
    outext = ["", ".index", ".lookup", "_h", "_h.index", ".dbtype", "_h.dbtype"]
    for file in [filename + ext for ext in outext]:
        open(file, "w").close()
    prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    logfile = os.path.join(GENEPATH, "test_create_mmseqsdb_exist.log")
    mmseqs.create_mmseqs_db(filename, prt_path, logfile)

    # Check files created/existing
    for file in [filename + ext for ext in outext]:
        assert os.path.isfile(file)
    assert ("mmseqs database "
            "test/data/pangenome/generated_by_unit-tests/test_create_mmseqsdb_exist.msdb "
            "already exists. The program will use it.") in caplog.text
    assert caplog.records[0].levelname == "WARNING"


def test_create_mmseqdb_not_all_exist(caplog):
    """
    Check that, when trying to create mmseqdb while the output database exists but at least
    1 associated file is missing (here, ".dbtype") -> message saying that files already exist, except some associated so database will be recreated.
    """
    caplog.set_level(logging.DEBUG)
    filename = os.path.join(GENEPATH, "test_create_mmseqsdb_not-all-exists.msdb")
    outext_miss = ["", ".index", ".lookup", "_h.index", ".dbtype"]
    for file in [filename + ext for ext in outext_miss]:
        open(file, "w").close()
    prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    logfile = os.path.join(GENEPATH, "test_create_mmseqsdb_exist.log")

    # Run mmseqs creation
    mmseqs.create_mmseqs_db(filename, prt_path, logfile)

    # Check that all expected files have been created (not only the files created before running)
    # and remove them, as well as the logfile
    outext_exp = ["", ".index", ".lookup", "_h", "_h.index", "_h.dbtype"]
    for file in [filename + ext for ext in outext_exp]:
        assert os.path.isfile(file)

    # Get logs info
    found_text = caplog.text
    # Check log content + level
    assert ("mmseqs database test/data/pangenome/generated_by_unit-tests/"
            "test_create_mmseqsdb_not-all-exists.msdb already exists, but at least 1 "
            "associated file (.dbtype, .index etc). is missing. The program will "
            "remove existing files and recreate the database.") in found_text
    assert caplog.records[0].levelname == "WARNING"
    assert ("Removing 'test/data/pangenome/generated_by_unit-tests/"
            "test_create_mmseqsdb_not-all-exists.msdb'") in found_text
    assert ("Removing 'test/data/pangenome/generated_by_unit-tests/"
            "test_create_mmseqsdb_not-all-exists.msdb.index'.") in found_text
    assert ("Removing 'test/data/pangenome/generated_by_unit-tests/"
            "test_create_mmseqsdb_not-all-exists.msdb.lookup'") in found_text
    assert not ("Removing test/data/pangenome/generated_by_unit-tests/"
                "test_create_mmseqsdb_not-all-exists.msdb_h'.") in found_text
    assert ("Removing 'test/data/pangenome/generated_by_unit-tests/"
            "test_create_mmseqsdb_not-all-exists.msdb_h.index'.") in found_text
    assert ("Removing 'test/data/pangenome/generated_by_unit-tests/"
            "test_create_mmseqsdb_not-all-exists.msdb.dbtype'.") in found_text
    assert not ("Removing 'test/data/pangenome/generated_by_unit-tests/"
                "test_create_mmseqsdb_not-all-exists.msdb_h.dbtype'.") in found_text
    for num in range(1, 6):
        assert caplog.records[num].levelname == "DETAIL"
    assert ("Creating database") in found_text
    assert caplog.records[6].levelname == "INFO"
    assert ("Existing files: 0") in found_text
    assert caplog.records[7].levelname == "DEBUG"
    assert ("Expected extensions: 7") in found_text
    assert caplog.records[8].levelname == "DEBUG"


def test_cluster2file():
    """
    Check that given clusters are written as expected to output file
    """
    fileout = os.path.join(GENEPATH, "test_clusters2file.txt")
    clusters = {"ESCO.0216.00001.i001_00006":
                ["ESCO.0216.00002.b010_00115", "ESCO.0216.00001.i001_00006",
                 "ESCA.0216.00001.i001_00015", "ESCO.0216.00001.b010_00115",
                 "ESCO.0216.00002.i001_00300", "ESCA.0216.00001.b001_00003"],
                "ESCO.0216.00001.i001_00018":
                ["ESCO.0216.00002.b010_00130", "ESCO.0216.00001.i001_00018",
                 "ESCA.0216.00001.i001_00950", "ESCO.0216.00002.i001_00300",
                 "ESCO.0216.00001.b001_00003"],
                "ESCO.0216.00002.b010_01265":
                ["ESCO.0216.00002.b010_01265"]}
    fams = mmseqs.clusters_to_file(clusters, fileout)
    # num of fams are random. Just check that families are the same,
    # and that nums correspond to range(1, nb_families+1)
    exp_fams = [["ESCO.0216.00002.b010_01265"],
                ["ESCA.0216.00001.i001_00950", "ESCO.0216.00001.b001_00003",
                 "ESCO.0216.00001.i001_00018", "ESCO.0216.00002.b010_00130",
                 "ESCO.0216.00002.i001_00300"],
                ["ESCA.0216.00001.b001_00003", "ESCA.0216.00001.i001_00015",
                 "ESCO.0216.00001.i001_00006", "ESCO.0216.00001.b010_00115",
                 "ESCO.0216.00002.b010_00115", "ESCO.0216.00002.i001_00300"]
                ]
    # Check that families given are as expected
    for fam in list(fams.values()):
        found = False
        for exp_fam in exp_fams:
            if set(fam) == set(exp_fam):
                found = True
                break
        assert found
    # Check family numbers are as expected
    assert list(fams.keys()) == list(range(1, 4))
    # Check pangenome file
    with open(fileout, "r") as fo:
        for line in fo:
            num = int(line.split()[0])
            fam = line.split()[1:]
            assert num in list(range(1, 4))
            assert fam in exp_fams


def test_tsv2cluster():
    """
    Check that conversion from mmseq tsv file to clusters is as expected.
    """
    filein = os.path.join(PATH_TEST_FILES, "mmseq_clust-out.tsv")
    clusters = mmseqs.mmseq_tsv_to_clusters(filein)
    for res, clust in clusters.items():
        found = False
        for resx, clustx in EXP_CLUSTERS.items():
            if resx == res and set(clust) == set(clustx):
                found = True
                break
        assert found


def test_tsv2pangenome_outgiven():
    """
    From mmseq tsv file, generate output pangenome file with a given name for it
    """
    mmseqclust = os.path.join(PATH_TEST_FILES, "mmseq_clust-out")
    logmmseq = os.path.join(GENEPATH, "test_tsv2pan.log")
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    outfile1 = os.path.join(GENEPATH, "test_tsv2pan_outpangenome.txt")
    fams, outfile = mmseqs.mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start, outfile1)
    # Check that outfile nae was not modified by the function
    assert outfile1 == outfile
    # Check that the number of families is as expected
    assert len(fams) == 16
    # Check that all families found are expected families
    # + that family numbers are between 1 and nb_families (same number of families found, and
    # consistent family numbers)
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in list(EXP_CLUSTERS.values()):
            if fam == expfam:
                found = True
                break
        assert found
    # Check that all families written in output file are as expected output file
    # (with exception that familiy numbers are not necessarily in the same order)
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp += line_exp.split()[1:]
            lines_out += line.split()[1:]
    assert set(lines_exp) == set(lines_out)
    # Check content of logfile
    with open(logmmseq, "r") as logf:
        first_lines = [logf.readline() for _ in range(3)]
        assert first_lines == ["\n", "------------\n", "\n"]
        start_line = logf.readline().strip()
        assert start_line == f"Start: {start}"
        end_line = logf.readline().strip()
        assert "End: " in end_line


def test_tsv2pangenome_default():
    """
    From mmseq tsv file, generate output pangenome file with default name
    """
    mmseqclust = os.path.join(PATH_TEST_FILES, "mmseq_clust-out")
    mmseqsclust_test = os.path.join(GENEPATH, "mmseq_clust-out_for_test")
    shutil.copyfile(mmseqclust, mmseqsclust_test)
    mmseqclust_tsv = os.path.join(PATH_TEST_FILES, "mmseq_clust-out.tsv")
    mmseqsclust_test_tsv = os.path.join(GENEPATH, "mmseq_clust-out_for_test.tsv")
    shutil.copyfile(mmseqclust_tsv, mmseqsclust_test_tsv)
    logmmseq = os.path.join(GENEPATH, "test_tsv2pan-def.log")
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    fams, outfile = mmseqs.mmseqs_tsv_to_pangenome(mmseqsclust_test, logmmseq, start)
    exp_out = os.path.join(GENEPATH, "PanGenome-mmseq_clust-out_for_test.tsv.lst")
    # Check name of output file (because not given)
    assert outfile == exp_out
    # Check that the number of families is as expected
    assert len(fams) == 16
    # Check that families returned are as expected
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in list(EXP_CLUSTERS.values()):
            if fam == expfam:
                found = True
                break
        assert found
    # Check output file contains expected families
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp += line_exp.split()[1:]
            lines_out += line.split()[1:]
    assert set(lines_exp) == set(lines_out)
    # Check logfile content
    with open(logmmseq, "r") as logf:
        first_lines = [logf.readline() for _ in range(3)]
        assert first_lines == ["\n", "------------\n", "\n"]
        start_line = logf.readline().strip()
        assert start_line == "Start: {}".format(start)
        end_line = logf.readline().strip()
        assert "End: " in end_line


def test_mmseq2pan_givenout():
    """
    From mmseq clust output, convert to pangenome (with steps inside, already tested by the other
    functions called).+ write pangenome to ouput file
    """
    # file witch will contain pangenome
    outfile1 = os.path.join(GENEPATH, "test_mmseq2pan.lst")
    # files with output of mmseqs clusters (in mmseqs format)
    mmseqclust = os.path.join(PATH_TEST_FILES, "mmseq_clust-out")
    # mmseq db used for clustering
    mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
    # Initialyze logfile with start time
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    logmmseq = os.path.join(GENEPATH, "test_mmseq2pan-out.log")

    # Run mmseqs2pan
    fams, outf = mmseqs.mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, start, outfile1)

    # Check output filename was not changed
    assert outfile1 == outf
    # Check that the number of families is as expected
    assert len(fams) == 16
    # assert output filename was not changed
    for num, fam in fams.items():
        # Check that the number of families return by the function is as expected
        assert num in list(range(1, 17))
        found = False
        # Check that all expected families are found in fams
        for expfam in list(EXP_CLUSTERS.values()):
            if fam == expfam:
                found = True
                break
        assert found
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    #Check that families written in output file are as expected
    with open(exp_pan, "r") as ep, open(outf, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp.append(tuple(line_exp.split()[1:]))
            lines_out.append(tuple(line.split()[1:]))
    assert set(lines_exp) == set(lines_out)


def test_run_clust():
    """
    Checks that, when we run mmseq clust, it creates all files needed for after to do
    the pangenome. We do not check the content of the mmseq output files, as it could
    depend on its version, and we are here testing PanACoTA.
    """
    mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
    mmseqclust = os.path.join(GENEPATH, "test_mmseq_cluster-out")
    tmpdir = os.path.join(GENEPATH, "test_mmseq_tmp")
    os.makedirs(tmpdir)
    logmmseq = os.path.join(GENEPATH, "test_mmseq_cluster.log")
    min_id = 0.8
    threads = 1
    clust_mode = 1
    args = (mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode)
    # Check that output of mmseq does not already exist
    assert not os.path.isfile(mmseqclust)
    # Run run_mmseqs_clust on previous arguments
    mmseqs.run_mmseqs_clust(args)
    # Check that all expected files and temporary directory are created
    # and check that no more outfile is created, in order to remove all of them !! !!
    generated_outfiles = glob.glob(mmseqclust + "*")
    assert len(generated_outfiles) == 3
    assert set(generated_outfiles) == set([mmseqclust, mmseqclust + ".index",
                                           mmseqclust + ".dbtype"])
    assert os.path.isfile(logmmseq)
    assert os.path.isdir(tmpdir)


def test_get_logmmseq():
    """
    Check that the given log filename is as expected according to given information
    """
    outdir = "toto"
    prt_bank = "bank_prt"
    infoname = "GENO115"
    log = mmseqs.get_logmmseq(outdir, prt_bank, infoname)
    assert log == "toto/mmseq_bank_prt_GENO115.log"


def test_get_info():
    """
    Check that string given by get_info is as expected according to info given in input
    """
    threads = 1
    min_id = 0.8
    clust_mode = 1
    start = "STARTTIME"
    info = mmseqs.get_info(threads, min_id, clust_mode, start)
    assert info == "0.8-mode1_STARTTIME"


def test_get_info_parallel():
    """
    Check that string given by get_info is as expected according to info given in input
    """
    threads = 12
    min_id = 0.8
    clust_mode = 1
    start = "STARTTIME"
    info = mmseqs.get_info(threads, min_id, clust_mode, start)
    assert info == "0.8-mode1-th12_STARTTIME"


def test_do_pangenome_defaultname(caplog):
    """
    Check that expected output files are created,
    and compare output pangenome to the expected one.
    """
    caplog.set_level(logging.DEBUG)
    outdir = os.path.join(GENEPATH, "test_do_pangenome_outdir")
    prt_bank = "exp_EXEM.All.prt"
    mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
    min_id = 0.8
    clust_mode = 1
    threads = 1
    start = "STARTTIME"
    quiet = False
    assert not os.path.isdir(outdir)
    fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
                                        clust_mode, threads, start, quiet=quiet)
    # Check creation of output directory
    assert os.path.isdir(outdir)
    # Check creation of tmp directory
    tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
    assert os.path.isdir(tmp_dir)
    # Check presence of pangenome file
    exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME.tsv.lst")
    assert exp_out == outfile
    assert os.path.isfile(outfile)
    # Check families returned
    assert len(fams) == 16
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in list(EXP_CLUSTERS.values()):
            if fam == expfam:
                found = True
                break
        assert found
    # Check content of output pangenome file
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp.append(tuple(line_exp.split()[1:]))
            lines_out.append(tuple(line.split()[1:]))
    assert set(lines_exp) == set(lines_out)
    assert "Clustering proteins..." in caplog.text
    assert caplog.records[1].levelname == "INFO"


def test_do_pangenome_given_panfile(caplog):
    """
    Check that expected output files are created,
    and compare output pangenome to the expected one.
    No possibility to check if quiet was not applied...
    """
    caplog.set_level(logging.DEBUG)
    outdir = os.path.join(GENEPATH, "test_do_pangenome_outdir_pan-name")
    prt_bank = "exp_EXEM.All.prt"
    mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
    min_id = 0.8
    clust_mode = 1
    threads = 1
    start = "STARTTIME"
    quiet = False
    panfile = os.path.join(GENEPATH, "test_res_pangenome_given-name")
    # Check outdir does not exit before
    assert not os.path.isdir(outdir)
    fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
                                        clust_mode, threads, start, quiet=quiet, panfile=panfile)
    # Check creation of output directory
    assert os.path.isdir(outdir)
    # Check creation of tmp directory
    tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
    assert os.path.isdir(tmp_dir)
    # Check presence of pangenome file
    assert panfile == outfile
    assert os.path.isfile(outfile)
    # Check families returned
    assert len(fams) == 16
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in list(EXP_CLUSTERS.values()):
            if fam == expfam:
                found = True
                break
        assert found
    # Check content of output pangenome file
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp.append(tuple(line_exp.split()[1:]))
            lines_out.append(tuple(line.split()[1:]))
    assert set(lines_exp) == set(lines_out)
    # Check log
    assert "Clustering proteins..." in caplog.text
    assert caplog.records[1].levelname == "INFO"


def test_do_pangenome_quiet(caplog):
    """
    Check that expected output files are created,
    and compare output pangenome to the expected one.
    Check that no error appears when choosing quiet option.
    No possibility to check if quiet was applied...
    """
    caplog.set_level(logging.DEBUG)
    outdir = os.path.join(GENEPATH, "test_do_pangenome_outdir_quiet")
    prt_bank = "exp_EXEM.All.prt"
    mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
    min_id = 0.8
    clust_mode = 1
    threads = 1
    start = "STARTTIME"
    quiet = True
    assert not os.path.isdir(outdir)
    fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
                                        clust_mode, threads, start, quiet=quiet)
    # Check creation of output directory
    assert os.path.isdir(outdir)
    # Check creation of tmp directory
    tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
    assert os.path.isdir(tmp_dir)
    # Check presence of pangenome file
    exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME.tsv.lst")
    assert exp_out == outfile
    assert os.path.isfile(outfile)
    # Check families returned
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in list(EXP_CLUSTERS.values()):
            if fam == expfam:
                found = True
                break
        assert found
    # Check content of output pangenome file
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp.append(tuple(line_exp.split()[1:]))
            lines_out.append(tuple(line.split()[1:]))
    assert set(lines_exp) == set(lines_out)
    assert "Clustering proteins..." in caplog.text


def test_do_pangenome_exist(caplog):
    """
    Check that if the mmseq output file of clustering already exists, it does not
    run mmseq again, but just converts it to pangenome.
    """
    caplog.set_level(logging.DEBUG)
    outdir = os.path.join(GENEPATH, "test_do_pangenome_outdir_exist")
    prt_bank = "exp_EXEM.All.prt"
    mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
    min_id = 0.8
    clust_mode = 1
    threads = 1
    start = "STARTTIME"
    # Create clustering results in outdir
    os.makedirs(outdir)
    orig_clust = os.path.join(PATH_TEST_FILES, "mmseq_clust-out")
    out_clust = os.path.join(outdir, "exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME")
    shutil.copyfile(orig_clust, out_clust)
    shutil.copyfile(orig_clust + ".index", out_clust + ".index")
    fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
                                        clust_mode, threads, start)
    assert ("mmseqs clustering test/data/pangenome/generated_by_unit-tests/"
            "test_do_pangenome_outdir_exist/exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME "
            "already exists. The program will now convert it to a "
            "pangenome file.") in caplog.text
    # assert
    # Check creation of empty tmp directory
    tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
    assert os.path.isdir(tmp_dir)
    assert glob.glob(os.path.join(tmp_dir, "*")) == []
    # Check presence of pangenome file
    exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME.tsv.lst")
    assert exp_out == outfile
    assert os.path.isfile(outfile)
    # Check families returned
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in list(EXP_CLUSTERS.values()):
            if fam == expfam:
                found = True
                break
        assert found
    # Check content of output pangenome file
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp.append(tuple(line_exp.split()[1:]))
            lines_out.append(tuple(line.split()[1:]))
    assert set(lines_exp) == set(lines_out)
    assert "Clustering proteins..." not in caplog.text
    assert caplog.records[1].levelname == "WARNING"


def test_run_all_pangenome(caplog):
    """
    Check that, given a prt bank, it creates mmseq db, mmseq clustering, and
    outputs the expected pangenome file.
    """
    caplog.set_level(logging.DEBUG)
    min_id = 0.8
    clust_mode = 1
    outdir = os.path.join(GENEPATH, "test_run_allpangenome")
    os.makedirs(outdir)
    prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    threads = 1
    panfile = None
    quiet = False
    fams, outfile = mmseqs.run_all_pangenome(min_id, clust_mode, outdir, prt_path,
                                             threads, panfile=panfile, quiet=quiet)

    # check that tmp dir was created and not empty
    tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_*")
    assert glob.glob(os.path.join(tmp_dir, "*")) != []

    # check that pangenome file is present with expected name
    exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_*")
    found_out = glob.glob(exp_out)
    assert len(found_out) == 1
    found_out = found_out[0]
    assert outfile == found_out
    assert os.path.isfile(outfile)

    # Check content of output pangenome file
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    # Compare lines, ignoring the 1st field (family number), and order.
    # -> families can be in any order/associated to any family number
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        # get all lines in exp and out files, without family number
        for line_exp, line in zip(ep, pan):
            # line_exp.split() -> list
            # But set([list, list]) is not possible (list is not hashable).
            # -> use tuples instead of lists
            lines_exp.append(tuple(line_exp.split()[1:]))
            lines_out.append(tuple(line.split()[1:]))
    # Check that out lines are exactly the same (except order)
    assert set(lines_exp) == set(lines_out)

    # Check families returned in fams dict.
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in FAMILIES4G:
            if fam == expfam:
                found = True
                break
        assert found

    # Check logs
    assert ("Will run MMseqs2 with:\n\t- minimum sequence identity = 80.0%\n"
            "\t- cluster mode 1") in caplog.text
    assert "Creating database" in caplog.text
    assert "Existing files: 0" in caplog.text
    assert "Expected extensions: 7" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert caplog.records[0].levelname == "INFO"
    assert caplog.records[1].levelname == "INFO"
    assert caplog.records[2].levelname == "DEBUG"
    assert caplog.records[3].levelname == "DEBUG"


def test_run_all_pangenome_givenfile_parallel(caplog):
    """
    Check that, given a prt bank, it creates mmseq db, mmseq clustering, and
    outputs the expected pangenome file.
    """
    caplog.set_level(logging.DEBUG)
    min_id = 0.8
    clust_mode = 1
    outdir = os.path.join(GENEPATH, "test_run_allpangenome_paralell")
    os.makedirs(outdir)
    prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    threads = 2
    panfile = "pangenome_test_run-all-pan.lst"
    quiet = True
    fams, outfile = mmseqs.run_all_pangenome(min_id, clust_mode, outdir, prt_path,
                                             threads, panfile=panfile, quiet=quiet)

    # check that tmp dir was created and not empty
    tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1-th2*")
    assert glob.glob(os.path.join(tmp_dir, "*")) != []

    # check that pangenome file is present
    assert outfile == os.path.join(outdir, panfile)
    assert os.path.isfile(outfile)

    # Check content of output pangenome file
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    # Compare lines, ignoring the 1st field (family number), and order.
    # -> families can be in any order/associated to any family number
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp.append(tuple(line_exp.split()[1:]))
            lines_out.append(tuple(line.split()[1:]))
    # Check that out lines are exactly the same (except order)
    assert set(lines_exp) == set(lines_out)

    # Check families returned
    for num, fam in fams.items():
        assert num in list(range(1, 17))
        found = False
        for expfam in FAMILIES4G:
            if fam == expfam:
                found = True
                break
        assert found

    # Check logs
    assert ("Will run MMseqs2 with:\n\t- minimum sequence identity = 80.0%\n"
            "\t- cluster mode 1\n\t- 2 threads") in caplog.text
    assert "Creating database" in caplog.text
    assert "Existing files: 0" in caplog.text
    assert "Expected extensions: 7" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert caplog.records[0].levelname == "INFO"
    assert caplog.records[1].levelname == "INFO"
    assert caplog.records[2].levelname == "DEBUG"
    assert caplog.records[3].levelname == "DEBUG"
