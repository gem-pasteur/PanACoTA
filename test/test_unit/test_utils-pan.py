#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the utils_pangenome submodule of genomeAPCAT
"""
import logging
import os

import shutil

from genomeAPCAT import utils_pangenome as upan, utils


# Define functions and variables shared by several tests
def my_logger():
    """
    logger given to function called by a subprocess
    """

    def make_logger(name="test_post_mmseq"):
        """
        Create logger according to name given
        """
        logfile_base = "log_" + name
        level = logging.DEBUG
        utils.init_logger(logfile_base, level, name, verbose=0, quiet=False)
        return logfile_base

    return make_logger


# Define common variables

PAN_FILE = os.path.join("test", "data", "pangenome", "exp_files", "exp_pangenome.txt")

FAMS_BY_STRAIN = \
    {'1': {'GEN4.1111.00001': ['GEN4.1111.00001.b0001_00001'],
           'GENO.0817.00001': ['GENO.0817.00001.b0001_00002'],
           'GENO.1216.00002': ['GENO.1216.00002.b0001_00001', 'GENO.1216.00002.i0001_00002']
           },
     '10': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00004'],
            'GENO.0817.00001': ['GENO.0817.00001.i0002_00004'],
            'GENO.1216.00002': ['GENO.1216.00002.i0001_00005']
            },
     '11': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00005'],
            'GENO.0817.00001': ['GENO.0817.00001.i0002_00005'],
            'GENO.1216.00002': ['GENO.1216.00002.i0001_00006']
            },
     '12': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00008'],
            'GENO.0817.00001': ['GENO.0817.00001.i0002_00010'],
            'GENO.1216.00002': ['GENO.1216.00002.b0002_00009']
            },
     '13': {'GENO.1216.00002': ['GENO.1216.00002.b0003_00011']},
     '14': {'GENO.1216.00002': ['GENO.1216.00002.b0003_00012']},
     '2': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00003']},
     '3': {'GEN4.1111.00001': ['GEN4.1111.00001.b0001_00009'],
           'GENO.0817.00001': ['GENO.0817.00001.b0002_00011'],
           'GENO.1216.00002': ['GENO.1216.00002.b0002_00010']
           },
     '4': {'GENO.0817.00001': ['GENO.0817.00001.b0001_00001']},
     '5': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00002'],
           'GENO.0817.00001': ['GENO.0817.00001.b0002_00003'],
           'GENO.1216.00002': ['GENO.1216.00002.i0001_00003']
           },
     '6': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00006'],
           'GENO.0817.00001': ['GENO.0817.00001.i0002_00006', 'GENO.0817.00001.i0002_00007'],
           'GENO.1216.00002': ['GENO.1216.00002.i0001_00007']
           },
     '7': {'GENO.0817.00001': ['GENO.0817.00001.i0002_00008']},
     '8': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00007'],
           'GENO.0817.00001': ['GENO.0817.00001.i0002_00009'],
           'GENO.1216.00002': ['GENO.1216.00002.b0001_00008']
           },
     '9': {'GENO.1216.00002': ['GENO.1216.00002.i0001_00004']}
     }

FAMILIES = {'1': ['GEN4.1111.00001.b0001_00001', 'GENO.0817.00001.b0001_00002',
                  'GENO.1216.00002.b0001_00001', 'GENO.1216.00002.i0001_00002'],
            '10': ['GEN4.1111.00001.i0001_00004', 'GENO.0817.00001.i0002_00004',
                   'GENO.1216.00002.i0001_00005'],
            '11': ['GEN4.1111.00001.i0001_00005', 'GENO.0817.00001.i0002_00005',
                   'GENO.1216.00002.i0001_00006'],
            '12': ['GEN4.1111.00001.i0001_00008', 'GENO.0817.00001.i0002_00010',
                   'GENO.1216.00002.b0002_00009'],
            '13': ['GENO.1216.00002.b0003_00011'],
            '14': ['GENO.1216.00002.b0003_00012'],
            '2': ['GEN4.1111.00001.i0001_00003'],
            '3': ['GEN4.1111.00001.b0001_00009', 'GENO.0817.00001.b0002_00011',
                  'GENO.1216.00002.b0002_00010'],
            '4': ['GENO.0817.00001.b0001_00001'],
            '5': ['GEN4.1111.00001.i0001_00002', 'GENO.0817.00001.b0002_00003',
                  'GENO.1216.00002.i0001_00003'],
            '6': ['GEN4.1111.00001.i0001_00006', 'GENO.0817.00001.i0002_00006',
                  'GENO.0817.00001.i0002_00007', 'GENO.1216.00002.i0001_00007'],
            '7': ['GENO.0817.00001.i0002_00008'],
            '8': ['GEN4.1111.00001.i0001_00007', 'GENO.0817.00001.i0002_00009',
                  'GENO.1216.00002.b0001_00008'],
            '9': ['GENO.1216.00002.i0001_00004']
            }

ALL_STRAINS = ['GEN4.1111.00001', 'GENO.0817.00001', 'GENO.1216.00002']


def test_read_gene():
    """
    Check that when reading a given gene name, it extracts expected information
    """
    gene = "ESCO.1016.00012.i012_00015"
    num = "1"
    fams_by_strain = {"1": {}}
    all_strains = []
    upan.read_gene(gene, num, fams_by_strain, all_strains)
    assert all_strains == ["ESCO.1016.00012"]
    assert fams_by_strain == {"1": {"ESCO.1016.00012": [gene]}}


def test_read_gene_other_format():
    """
    Check that when reading a gene name which does not have the gembase
    format, it still returns the right information
    """
    gene = "my_gene-name_other_0001"
    num = "1"
    fams_by_strain = {"1": {}}
    all_strains = []
    upan.read_gene(gene, num, fams_by_strain, all_strains)
    assert all_strains == ["my_gene-name_other"]
    assert fams_by_strain == {"1": {"my_gene-name_other": [gene]}}


def test_read_gene_strain_known():
    """
    Check that when reading a gene name, and its corresponding
    strain already exists in fams_by_strain and all_strains, it just adds
    this new gene to fams_by_strain, and does nothing to all_strains
    """
    gene = "ESCO.1016.00012.i012_00015"
    num = "1"
    fams_by_strain = {"1": {"ESCO.1016.00012": ["ESCO.1016.00012.i001_01"]}}
    all_strains = ["ESCO.1016.00012"]
    upan.read_gene(gene, num, fams_by_strain, all_strains)
    assert all_strains == ["ESCO.1016.00012"]
    assert fams_by_strain == {"1": {"ESCO.1016.00012": ["ESCO.1016.00012.i001_01",
                                                        gene]}}


def test_read_panfile(caplog):
    """
    check that it reads the pangenome file and returns the expected objects
    """
    fbs, fams, sas = upan.read_pan_file(PAN_FILE)
    assert fbs == FAMS_BY_STRAIN
    assert fams == FAMILIES
    assert sas == ALL_STRAINS
    assert "Reading and getting information from pangenome file" in caplog.text


def test_get_fams(caplog):
    """
    Test that when giving all members for each family, it returns all strains involved in
    those families, and all families with members sorted by strain
    """
    bystrain, sortedf = upan.get_fams_info(FAMILIES)
    assert bystrain == FAMS_BY_STRAIN
    assert sortedf == ALL_STRAINS
    assert "Retrieving information from pan families" in caplog.text


def test_read_pan_filetxt(caplog):
    """
    Test that when giving a pangenome file, it returns all families as expected.
    """
    fbs, fams, ass = upan.read_pangenome(PAN_FILE)
    assert fbs == FAMS_BY_STRAIN
    assert fams == FAMILIES
    assert ass == ALL_STRAINS
    assert "Reading and getting information from pangenome file" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert os.path.isfile(PAN_FILE + ".bin")
    # os.remove(PAN_FILE + ".bin")


def test_read_pan_filebin(caplog):
    """
    Test that when giving only a pangenome filename, and the corresponding bin file exists,
    it reads the binary file, and returns expected objects.
    """
    panbin = os.path.join("test", "data", "pangenome", "test_files", "pangenome.bin")
    panbinres = PAN_FILE + ".bin"
    shutil.copyfile(panbin, panbinres)
    fbs, fams, ass = upan.read_pangenome(PAN_FILE)
    assert fbs == FAMS_BY_STRAIN
    assert fams == FAMILIES
    assert ass == ALL_STRAINS
    assert "Retrieving info from binary file" in caplog.text
    os.remove(panbinres)


def test_read_pan_fams(caplog):
    """
    Test that when giving a pangenome file, and families, it directly extracts strain information
    from the families: pangenome file does not need to exist, and a binary file is created
    """
    logger = my_logger()
    name = "test_read_panfams"
    logger(name)
    panfile = "test_readpanbin.txt"
    fbs, fams, ass = upan.read_pangenome(panfile, FAMILIES)
    assert fbs == FAMS_BY_STRAIN
    assert fams == FAMILIES
    assert ass == ALL_STRAINS
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert os.path.isfile(panfile + ".bin")
    os.remove(panfile + ".bin")
    os.remove("log_" + name + ".log")
    os.remove("log_" + name + ".log.details")
    os.remove("log_" + name + ".log.err")


def test_read_pan_fams_binok(caplog):
    """
    Test that when giving a pangenome file, and families, it directly extracts strain information
    from the families: pangenome file does not need to exist. However, the pangenome.bin file
    already exists (whatever its content), and is then not recreated.
    """
    logger = my_logger()
    name = "test_read_panfams_bin"
    logger(name)
    panfile = "test_readpanbinok.txt"
    # Create bn pangenome file (which is empty
    open(panfile + ".bin", "w").close()
    fbs, fams, ass = upan.read_pangenome(panfile, FAMILIES)
    assert fbs == FAMS_BY_STRAIN
    assert fams == FAMILIES
    assert ass == ALL_STRAINS
    with open(panfile + ".bin", "r") as panf:
        all_lines = panf.readlines()
    assert all_lines == []
    assert "Retrieving information from pan families" in caplog.text
    assert os.path.isfile(panfile + ".bin")
    os.remove(panfile + ".bin")
    os.remove("log_" + name + ".log")
    os.remove("log_" + name + ".log.details")
    os.remove("log_" + name + ".log.err")