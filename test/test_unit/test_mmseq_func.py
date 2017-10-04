#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the mmseqs_functions submodule in pangenome module
"""

import pytest
import os

import genomeAPCAT.pangenome_module.mmseqs_functions as mmseqs


@pytest.fixture(scope="function")
def path_test_pan():
    return os.path.join("test", "data", "pangenome")


@pytest.fixture(scope="function")
def path_test_files():
    return os.path.join(path_test_pan(), "test_files")


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

