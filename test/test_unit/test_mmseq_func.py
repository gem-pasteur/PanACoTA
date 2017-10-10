#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the mmseqs_functions submodule in pangenome module
"""

import pytest
import os
import time
import shutil

import genomeAPCAT.pangenome_module.mmseqs_functions as mmseqs


@pytest.fixture(scope="function")
def path_test_pan():
    return os.path.join("test", "data", "pangenome")


@pytest.fixture(scope="function")
def path_test_files():
    return os.path.join(path_test_pan(), "test_files")


@pytest.fixture(scope="function")
def path_exp_files():
    return os.path.join(path_test_pan(), "exp_files")


@pytest.fixture(scope="function")
def exp_clusters():
    """
    Expected clusters generated my mmseq
    """
    clusters = {"GEN4.1111.00001.b0001_00001":
                ["GEN4.1111.00001.b0001_00001", "GENO.0817.00001.b0001_00002",
                 "GENO.1216.00002.b0001_00001", "GENO.1216.00002.i0001_00002"],
                "GEN4.1111.00001.i0001_00003": ["GEN4.1111.00001.i0001_00003"],
                "GEN4.1111.00001.b0001_00009":
                 ["GEN4.1111.00001.b0001_00009", "GENO.0817.00001.b0002_00011",
                  "GENO.1216.00002.b0002_00010"],
                "GENO.0817.00001.b0001_00001": ["GENO.0817.00001.b0001_00001"],
                "GENO.0817.00001.b0002_00003":
                 ["GENO.0817.00001.b0002_00003", "GEN4.1111.00001.i0001_00002",
                  "GENO.1216.00002.i0001_00003"],
                "GENO.0817.00001.i0002_00006":
                 ["GENO.0817.00001.i0002_00006", "GENO.0817.00001.i0002_00007",
                  "GENO.1216.00002.i0001_00007", "GEN4.1111.00001.i0001_00006"],
                "GENO.0817.00001.i0002_00008": ["GENO.0817.00001.i0002_00008"],
                "GENO.0817.00001.i0002_00009":
                 ["GENO.0817.00001.i0002_00009", "GENO.1216.00002.b0001_00008",
                  "GEN4.1111.00001.i0001_00007"],
                "GENO.1216.00002.i0001_00004": ["GENO.1216.00002.i0001_00004"],
                "GENO.1216.00002.i0001_00005":
                 ["GENO.1216.00002.i0001_00005", "GENO.0817.00001.i0002_00004",
                  "GEN4.1111.00001.i0001_00004"],
                "GENO.1216.00002.i0001_00006":
                 ["GENO.1216.00002.i0001_00006", "GENO.0817.00001.i0002_00005",
                  "GEN4.1111.00001.i0001_00005"],
                "GENO.1216.00002.b0002_00009":
                 ["GENO.1216.00002.b0002_00009", "GEN4.1111.00001.i0001_00008",
                  "GENO.0817.00001.i0002_00010"],
                "GENO.1216.00002.b0003_00011": ["GENO.1216.00002.b0003_00011"],
                "GENO.1216.00002.b0003_00012": ["GENO.1216.00002.b0003_00012"]
                }
    return clusters


def test_create_mmseqdb(path_test_files, caplog):
    """
    Test that mmseq DB is created. We do not check its content as it could change
    according to mmseq versions, and we are here testing genomeAPCAT, not mmseqs
    """
    filename = "test_create_mmseqsdb.msdb"
    prt_path = os.path.join(path_test_files, "example_db", "Proteins")
    logfile = "test_create_mmseqsdb.log"
    mmseqs.create_mmseqs_db(filename, prt_path, logfile)
    outext = ["", ".index", ".lookup", "_h", "_h.index"]
    for file in [filename + ext for ext in outext]:
        assert os.path.isfile(file)
        os.remove(file)
    assert os.path.isfile(logfile)
    os.remove(logfile)
    assert "Creating database" in caplog.text
    assert caplog.records[0].levelname == "INFO"


def test_create_mmseqdb_exist(path_test_files, caplog):
    """
    Check that, when trying to create mmseqdb while the output file already exists,
    it logs a warning message and quits without creating it
    """
    filename = "test_create_mmseqsdb.msdb"
    open(filename, "w").close()
    prt_path = os.path.join(path_test_files, "example_db", "Proteins")
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


def test_tsv2cluster(path_test_files, exp_clusters):
    """
    Check that conversion from mmseq tsv file to clusters is as expected.
    """
    filein = os.path.join(path_test_files, "mmseq_tsvfile.tsv")
    clusters = mmseqs.mmseq_tsv_to_clusters(filein)
    assert clusters == exp_clusters


def test_tsv2pangenome_outgiven(path_test_files, path_exp_files, exp_clusters):
    """
    From mmseq tsv file, generate output pangenome file
    """
    mmseqclust = os.path.join(path_test_files, "mmseq_tsvfile")
    logmmseq = "test_tsv2pan.log"
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    outfile1 = "test_tsv2pan_outpangenome.txt"
    fams, outfile = mmseqs.mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start, outfile1)
    assert outfile1 == outfile
    for num, fam in fams.items():
        assert num in list(range(1, 15))
        found = False
        for expfam in list(exp_clusters.values()):
            if fam == expfam:
                found = True
                break
        assert found
    exp_pan = os.path.join(path_exp_files, "exp_pangenome.txt")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp += line_exp.split()[1:]
            lines_out += line.split()[1:]
    assert set(lines_exp) == set(lines_out)
    os.remove(outfile)
    with open(logmmseq, "r") as logf:
        first_lines = [logf.readline() for _ in range(3)]
        assert first_lines == ["\n", "------------\n", "\n"]
        start_line = logf.readline().strip()
        assert start_line == "Start: {}".format(start)
        end_line = logf.readline().strip()
        assert "End: " in end_line
    os.remove(logmmseq)


def test_tsv2pangenome_default(path_test_files, path_exp_files, exp_clusters):
    """
    From mmseq tsv file, generate output pangenome file
    """
    mmseqclust = os.path.join(path_test_files, "mmseq_tsvfile")
    logmmseq = "test_tsv2pan-def.log"
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    fams, outfile = mmseqs.mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start)
    exp_out = os.path.join(path_test_files, "PanGenome-mmseq_tsvfile.tsv.lst")
    assert outfile == exp_out
    for num, fam in fams.items():
        assert num in list(range(1, 15))
        found = False
        for expfam in list(exp_clusters.values()):
            if fam == expfam:
                found = True
                break
        assert found
    exp_pan = os.path.join(path_exp_files, "exp_pangenome.txt")
    with open(exp_pan, "r") as ep, open(outfile, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp += line_exp.split()[1:]
            lines_out += line.split()[1:]
    assert set(lines_exp) == set(lines_out)
    os.remove(outfile)
    with open(logmmseq, "r") as logf:
        first_lines = [logf.readline() for _ in range(3)]
        assert first_lines == ["\n", "------------\n", "\n"]
        start_line = logf.readline().strip()
        assert start_line == "Start: {}".format(start)
        end_line = logf.readline().strip()
        assert "End: " in end_line
    os.remove(logmmseq)


def test_mmseq2pan_givenout(path_test_files, path_exp_files, exp_clusters):
    """
    From mmseq output, convert to pangenome (with steps inside, already tested by the other
    functions called)
    """
    outfile1 = "test_mmseq2pan.lst"
    mmseqclust = os.path.join(path_test_files, "mmseq_clust-out")
    mmseqdb = os.path.join(path_test_files, "mmseq_db")
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    logmmseq = "test_mmseq2pan-out.log"
    fams, outf = mmseqs.mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, start, outfile1)
    assert outfile1 == outf
    for num, fam in fams.items():
        assert num in list(range(1, 15))
        found = False
        for expfam in list(exp_clusters.values()):
            if fam == expfam:
                found = True
                break
        assert found
    exp_pan = os.path.join(path_exp_files, "exp_pangenome.txt")
    with open(exp_pan, "r") as ep, open(outf, "r") as pan:
        lines_exp = []
        lines_out = []
        for line_exp, line in zip(ep, pan):
            lines_exp += line_exp.split()[1:]
            lines_out += line.split()[1:]
    assert set(lines_exp) == set(lines_out)
    os.remove(outf)
    os.remove(logmmseq)


def test_run_clust(path_test_files):
    """
    Checks that, when we run mmseq clust, it creates all files needed for after to do
    the pangenome. We do not check the content of the mmseq output files, as it could
    depend on its version, and we are here testing genomeAPCAT.
    """
    mmseqdb = os.path.join(path_test_files, "mmseq_db")
    mmseqclust = "test_mmseq_cluster-out"
    tmpdir = "test_mmseq_tmp"
    os.makedirs(tmpdir)
    logmmseq = "test_mmseq_cluster.log"
    min_id = 0.8
    threads = 1
    clust_mode = 1
    args = (mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode)
    assert not os.path.isfile(mmseqclust)
    mmseqs.run_mmseqs_clust(args)
    assert os.path.isfile(mmseqclust)
    assert os.path.isfile(mmseqclust + ".index")
    assert os.path.isfile(logmmseq)
    assert os.path.isdir(tmpdir)
    shutil.rmtree(tmpdir)
    os.remove(mmseqclust)
    os.remove(mmseqclust + ".index")
    os.remove(logmmseq)


def test_get_logmmseq():
    """
    Check that the given log filename is as expected according to given information
    """
    outdir = "toto"
    prt_bank = "bank_prt"
    infoname = "GENO115"
    log = mmseqs.get_logmmseq(outdir, prt_bank, infoname)
    assert log == "toto/mmseq_bank_prt_GENO115.log"
