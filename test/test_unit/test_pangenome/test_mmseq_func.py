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

import genomeAPCAT.pangenome_module.mmseqs_functions as mmseqs

# Define variables shared by several tests
PATH_TEST_PAN = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PATH_TEST_PAN, "test_files")
PATH_EXP_FILES = os.path.join(PATH_TEST_PAN, "exp_files")

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
    according to mmseq versions, and we are here testing genomeAPCAT, not mmseqs
    """
    caplog.set_level(logging.DEBUG)
    filename = "test_create_mmseqsdb.msdb"
    prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    logfile = "test_create_mmseqsdb.log"
    mmseqs.create_mmseqs_db(filename, prt_path, logfile)
    outext = ["", ".index", ".lookup", "_h", "_h.index"]
    for file in [filename + ext for ext in outext]:
        assert os.path.isfile(file)
        os.remove(file)
    if os.path.isfile(filename + ".dbtype"):
        os.remove(filename + ".dbtype")
    if os.path.isfile(filename + "_h.dbtype"):
        os.remove(filename + "_h.dbtype")
    assert os.path.isfile(logfile)
    os.remove(logfile)
    assert "Creating database" in caplog.text
    assert caplog.records[0].levelname == "INFO"


def test_create_mmseqdb_exist(caplog):
    """
    Check that, when trying to create mmseqdb while the output file already exists,
    it logs a warning message and quits without creating it
    """
    filename = "test_create_mmseqsdb.msdb"
    open(filename, "w").close()
    prt_path = os.path.join(PATH_TEST_FILES, "example_db", "Proteins")
    logfile = "test_create_mmseqsdb.log"
    mmseqs.create_mmseqs_db(filename, prt_path, logfile)
    outext = [".index", ".lookup", "_h", "_h.index"]
    for file in [filename + ext for ext in outext]:
        assert not os.path.isfile(file)
    assert not os.path.isfile(logfile)
    os.remove(filename)
    assert ("mmseq database test_create_mmseqsdb.msdb already exists. "
            "The program will use it.") in caplog.text
    assert caplog.records[0].levelname == "WARNING"


def test_cluster2file():
    """
    Check that given clusters are written as expected to output file
    """
    fileout = "test_clusters2file.txt"
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
    os.remove(fileout)


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
    logmmseq = "test_tsv2pan.log"
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    outfile1 = "test_tsv2pan_outpangenome.txt"
    fams, outfile = mmseqs.mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start, outfile1)
    # Check that outfile nae was not modified by the function
    assert outfile1 == outfile
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
    os.remove(outfile)
    # Check content of logfile
    with open(logmmseq, "r") as logf:
        first_lines = [logf.readline() for _ in range(3)]
        assert first_lines == ["\n", "------------\n", "\n"]
        start_line = logf.readline().strip()
        assert start_line == "Start: {}".format(start)
        end_line = logf.readline().strip()
        assert "End: " in end_line
    os.remove(logmmseq)


def test_tsv2pangenome_default():
    """
    From mmseq tsv file, generate output pangenome file with default name
    """
    mmseqclust = os.path.join(PATH_TEST_FILES, "mmseq_clust-out")
    logmmseq = "test_tsv2pan-def.log"
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    fams, outfile = mmseqs.mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start)
    exp_out = os.path.join(PATH_TEST_FILES, "PanGenome-mmseq_clust-out.tsv.lst")
    # Check name of output file (because not given)
    assert outfile == exp_out
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
    os.remove(outfile)
    # Check logfile content
    with open(logmmseq, "r") as logf:
        first_lines = [logf.readline() for _ in range(3)]
        assert first_lines == ["\n", "------------\n", "\n"]
        start_line = logf.readline().strip()
        assert start_line == "Start: {}".format(start)
        end_line = logf.readline().strip()
        assert "End: " in end_line
    os.remove(logmmseq)


# def test_mmseq2pan_givenout():
#     """
#     From mmseq clust output, convert to pangenome (with steps inside, already tested by the other
#     functions called).+ write pangenome to ouput file
#     """
#     outfile1 = "test_mmseq2pan.lst"
#     mmseqclust = os.path.join(PATH_TEST_FILES, "mmseq_clust-out")
#     mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
#     start = time.strftime('%Y-%m-%d_%H-%M-%S')
#     logmmseq = "test_mmseq2pan-out.log"
#     fams, outf = mmseqs.mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, start, outfile1)
#     assert outfile1 == outf
#     for num, fam in fams.items():
#         assert num in list(range(1, 17))
#         found = False
#         for expfam in list(EXP_CLUSTERS.values()):
#             if fam == expfam:
#                 found = True
#                 break
#         assert found
#     exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
#     with open(exp_pan, "r") as ep, open(outf, "r") as pan:
#         lines_exp = []
#         lines_out = []
#         for line_exp, line in zip(ep, pan):
#             lines_exp.append(tuple(line_exp.split()[1:]))
#             lines_out.append(tuple(line.split()[1:]))
#     assert set(lines_exp) == set(lines_out)
#     os.remove(outf)
#     os.remove(logmmseq)


# def test_run_clust():
#     """
#     Checks that, when we run mmseq clust, it creates all files needed for after to do
#     the pangenome. We do not check the content of the mmseq output files, as it could
#     depend on its version, and we are here testing genomeAPCAT.
#     """
#     mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
#     mmseqclust = "test_mmseq_cluster-out"
#     tmpdir = "test_mmseq_tmp"
#     os.makedirs(tmpdir)
#     logmmseq = "test_mmseq_cluster.log"
#     min_id = 0.8
#     threads = 1
#     clust_mode = 1
#     args = (mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode)
#     assert not os.path.isfile(mmseqclust)
#     mmseqs.run_mmseqs_clust(args)
#     assert os.path.isfile(mmseqclust)
#     assert os.path.isfile(mmseqclust + ".index")
#     assert os.path.isfile(logmmseq)
#     assert os.path.isdir(tmpdir)
#     shutil.rmtree(tmpdir)
#     os.remove(mmseqclust)
#     os.remove(mmseqclust + ".index")
#     os.remove(logmmseq)


# def test_get_logmmseq():
#     """
#     Check that the given log filename is as expected according to given information
#     """
#     outdir = "toto"
#     prt_bank = "bank_prt"
#     infoname = "GENO115"
#     log = mmseqs.get_logmmseq(outdir, prt_bank, infoname)
#     assert log == "toto/mmseq_bank_prt_GENO115.log"


# def test_get_info():
#     """
#     Check that string given by get_info is as expected according to info given in input
#     """
#     threads = 1
#     min_id = 0.8
#     clust_mode = 1
#     start = "STARTTIME"
#     info = mmseqs.get_info(threads, min_id, clust_mode, start)
#     assert info == "0.8-mode1_STARTTIME"


# def test_get_info_parallel():
#     """
#     Check that string given by get_info is as expected according to info given in input
#     """
#     threads = 12
#     min_id = 0.8
#     clust_mode = 1
#     start = "STARTTIME"
#     info = mmseqs.get_info(threads, min_id, clust_mode, start)
#     assert info == "0.8-mode1-th12_STARTTIME"


# def test_do_pangenome(caplog):
#     """
#     Check that expected output files are created,
#     and compare output pangenome to the expected one.
#     """
#     caplog.set_level(logging.DEBUG)
#     outdir = "test_do_pangenome_outdir"
#     prt_bank = "exp_EXEM.All.prt"
#     mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
#     min_id = 0.8
#     clust_mode = 1
#     threads = 1
#     start = "STARTTIME"
#     quiet = False
#     assert not os.path.isdir(outdir)
#     fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
#                                         clust_mode, threads, start, quiet=quiet)
#     # Check creation of output directory
#     assert os.path.isdir(outdir)
#     # Check creation of tmp directory
#     tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
#     assert os.path.isdir(tmp_dir)
#     # Check presence of pangenome file
#     exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME.tsv.lst")
#     assert exp_out == outfile
#     assert os.path.isfile(outfile)
#     # Check families returned
#     for num, fam in fams.items():
#         assert num in list(range(1, 17))
#         found = False
#         for expfam in list(EXP_CLUSTERS.values()):
#             if fam == expfam:
#                 found = True
#                 break
#         assert found
#     # Check content of output pangenome file
#     exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
#     with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
#         lines_exp = []
#         lines_out = []
#         for line_exp, line in zip(ep, pan):
#             lines_exp.append(tuple(line_exp.split()[1:]))
#             lines_out.append(tuple(line.split()[1:]))
#     assert set(lines_exp) == set(lines_out)
#     assert "Clustering proteins..." in caplog.text
#     shutil.rmtree(outdir)


# def test_do_pangenome_given_panfile(caplog):
#     """
#     Check that expected output files are created,
#     and compare output pangenome to the expected one.
#     """
#     caplog.set_level(logging.DEBUG)
#     outdir = "test_do_pangenome_outdir"
#     prt_bank = "exp_EXEM.All.prt"
#     mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
#     min_id = 0.8
#     clust_mode = 1
#     threads = 1
#     start = "STARTTIME"
#     quiet = False
#     panfile = "test_res_pangenome"
#     assert not os.path.isdir(outdir)
#     fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
#                                         clust_mode, threads, start, quiet=quiet, panfile=panfile)
#     print(outfile)
#     # Check creation of output directory
#     assert os.path.isdir(outdir)
#     # Check creation of tmp directory
#     tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
#     assert os.path.isdir(tmp_dir)
#     # Check presence of pangenome file
#     assert panfile == outfile
#     assert os.path.isfile(outfile)
#     # Check families returned
#     for num, fam in fams.items():
#         assert num in list(range(1, 17))
#         found = False
#         for expfam in list(EXP_CLUSTERS.values()):
#             if fam == expfam:
#                 found = True
#                 break
#         assert found
#     # Check content of output pangenome file
#     exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
#     with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
#         lines_exp = []
#         lines_out = []
#         for line_exp, line in zip(ep, pan):
#             lines_exp.append(tuple(line_exp.split()[1:]))
#             lines_out.append(tuple(line.split()[1:]))
#     assert set(lines_exp) == set(lines_out)
#     assert "Clustering proteins..." in caplog.text
#     shutil.rmtree(outdir)
#     os.remove(panfile)


# def test_do_pangenome_quiet(caplog):
#     """
#     Check that expected output files are created,
#     and compare output pangenome to the expected one.
#     Check that no error appears when choosing quiet option.
#     """
#     caplog.set_level(logging.DEBUG)
#     outdir = "test_do_pangenome_outdir"
#     prt_bank = "exp_EXEM.All.prt"
#     mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
#     min_id = 0.8
#     clust_mode = 1
#     threads = 1
#     start = "STARTTIME"
#     quiet = True
#     assert not os.path.isdir(outdir)
#     fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
#                                         clust_mode, threads, start, quiet=quiet)
#     # Check creation of output directory
#     assert os.path.isdir(outdir)
#     # Check creation of tmp directory
#     tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
#     assert os.path.isdir(tmp_dir)
#     # Check presence of pangenome file
#     exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME.tsv.lst")
#     assert exp_out == outfile
#     assert os.path.isfile(outfile)
#     # Check families returned
#     for num, fam in fams.items():
#         assert num in list(range(1, 17))
#         found = False
#         for expfam in list(EXP_CLUSTERS.values()):
#             if fam == expfam:
#                 found = True
#                 break
#         assert found
#     # Check content of output pangenome file
#     exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
#     with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
#         lines_exp = []
#         lines_out = []
#         for line_exp, line in zip(ep, pan):
#             lines_exp.append(tuple(line_exp.split()[1:]))
#             lines_out.append(tuple(line.split()[1:]))
#     assert set(lines_exp) == set(lines_out)
#     assert "Clustering proteins..." in caplog.text
#     shutil.rmtree(outdir)


# def test_do_pangenome_exist(caplog):
#     """
#     Check that if the mmseq output file of clustering already exists, it does not
#     run mmseq again, but just converts it to pangenome.
#     """
#     caplog.set_level(logging.DEBUG)
#     outdir = "test_do_pangenome_outdir_exist"
#     prt_bank = "exp_EXEM.All.prt"
#     mmseqdb = os.path.join(PATH_TEST_FILES, "mmseq_db")
#     min_id = 0.8
#     clust_mode = 1
#     threads = 1
#     start = "STARTTIME"
#     # Create clustering results in outdir
#     os.makedirs(outdir)
#     orig_clust = os.path.join(PATH_TEST_FILES, "mmseq_clust-out")
#     out_clust = os.path.join(outdir, "exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME")
#     shutil.copyfile(orig_clust, out_clust)
#     shutil.copyfile(orig_clust + ".index", out_clust + ".index")
#     fams, outfile = mmseqs.do_pangenome(outdir, prt_bank, mmseqdb, min_id,
#                                         clust_mode, threads, start)
#     assert ("mmseqs clustering test_do_pangenome_outdir_exist/exp_EXEM.All.prt-clust-0.8-"
#             "mode1_STARTTIME already exists. The program will now convert it to a "
#             "pangenome file.") in caplog.text
#     # Check creation of empty tmp directory
#     tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_STARTTIME")
#     assert os.path.isdir(tmp_dir)
#     assert glob.glob(os.path.join(tmp_dir, "*")) == []
#     # Check presence of pangenome file
#     exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_STARTTIME.tsv.lst")
#     assert exp_out == outfile
#     assert os.path.isfile(outfile)
#     # Check families returned
#     for num, fam in fams.items():
#         assert num in list(range(1, 17))
#         found = False
#         for expfam in list(EXP_CLUSTERS.values()):
#             if fam == expfam:
#                 found = True
#                 break
#         assert found
#     # Check content of output pangenome file
#     exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
#     with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
#         lines_exp = []
#         lines_out = []
#         for line_exp, line in zip(ep, pan):
#             lines_exp.append(tuple(line_exp.split()[1:]))
#             lines_out.append(tuple(line.split()[1:]))
#     assert set(lines_exp) == set(lines_out)
#     assert "Clustering proteins..." not in caplog.text
#     shutil.rmtree(outdir)


# def test_run_all_pangenome(caplog):
#     """
#     Check that, given a prt bank, it creates mmseq db, mmseq clustering, and
#     outputs the expected pangenome file.
#     """
#     caplog.set_level(logging.DEBUG)
#     min_id = 0.8
#     clust_mode = 1
#     outdir = "test_run_allpangenome"
#     os.makedirs(outdir)
#     prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
#     threads = 1
#     panfile = None
#     quiet = False
#     fams, outfile = mmseqs.run_all_pangenome(min_id, clust_mode, outdir, prt_path,
#                                              threads, panfile=panfile, quiet=quiet)
#     # check that tmp dir was created and not empty
#     tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1_*")
#     assert glob.glob(os.path.join(tmp_dir, "*")) != []
#     # check that pangenome file is present
#     exp_out = os.path.join(outdir, "PanGenome-exp_EXEM.All.prt-clust-0.8-mode1_*")
#     found_out = glob.glob(exp_out)
#     assert len(found_out) == 1
#     found_out = found_out[0]
#     assert outfile == found_out
#     assert os.path.isfile(outfile)
#     # Check content of output pangenome file
#     exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
#     with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
#         lines_exp = []
#         lines_out = []
#         for line_exp, line in zip(ep, pan):
#             lines_exp.append(tuple(line_exp.split()[1:]))
#             lines_out.append(tuple(line.split()[1:]))
#     assert set(lines_exp) == set(lines_out)
#     # Check families returned
#     for num, fam in fams.items():
#         assert num in list(range(1, 17))
#         found = False
#         for expfam in FAMILIES4G:
#             if fam == expfam:
#                 found = True
#                 break
#         assert found
#     assert ("Will run MMseqs2 with:\n\t- minimum sequence identity = 0.8\n"
#             "\t- cluster mode 1") in caplog.text
#     shutil.rmtree(outdir)


# def test_run_all_pangenome_givenfile_parallel(caplog):
#     """
#     Check that, given a prt bank, it creates mmseq db, mmseq clustering, and
#     outputs the expected pangenome file.
#     """
#     caplog.set_level(logging.DEBUG)
#     min_id = 0.8
#     clust_mode = 1
#     outdir = "test_run_allpangenome"
#     os.makedirs(outdir)
#     prt_path = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
#     threads = 2
#     panfile = "pangenome_test_run-all-pan.lst"
#     quiet = True
#     fams, outfile = mmseqs.run_all_pangenome(min_id, clust_mode, outdir, prt_path,
#                                              threads, panfile=panfile, quiet=quiet)
#     # check that tmp dir was created and not empty
#     tmp_dir = os.path.join(outdir, "tmp_exp_EXEM.All.prt_0.8-mode1-th2*")
#     assert glob.glob(os.path.join(tmp_dir, "*")) != []
#     # check that pangenome file is present
#     assert outfile == os.path.join(outdir, panfile)
#     assert os.path.isfile(outfile)
#     # Check content of output pangenome file
#     exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
#     with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
#         lines_exp = []
#         lines_out = []
#         for line_exp, line in zip(ep, pan):
#             lines_exp.append(tuple(line_exp.split()[1:]))
#             lines_out.append(tuple(line.split()[1:]))
#     assert set(lines_exp) == set(lines_out)
#     # Check families returned
#     for num, fam in fams.items():
#         assert num in list(range(1, 17))
#         found = False
#         for expfam in FAMILIES4G:
#             if fam == expfam:
#                 found = True
#                 break
#         assert found
#     assert ("Will run MMseqs2 with:\n\t- minimum sequence identity = 0.8\n"
#             "\t- cluster mode 1\n\t- 2 threads") in caplog.text
#     shutil.rmtree(outdir)
