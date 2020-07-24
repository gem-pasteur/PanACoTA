#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the mmseqs_functions submodule in pangenome module
"""

import os
import io
import logging
import pytest
import shutil

import PanACoTA.pangenome_module.post_treatment as post
import PanACoTA.utils as utils
import test.test_unit.utilities_for_tests as tutil

# Define variables shared by several tests
PANDIR = os.path.join("test", "data", "pangenome")
PATH_TEST_FILES = os.path.join(PANDIR, "test_files")
PATH_EXP_FILES = os.path.join(PANDIR, "exp_files")
GENEPATH = os.path.join(PANDIR, "generated_by_unit-tests")
# log files
LOGFILE_BASE = "logfile_test.txt"
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
    utils.init_logger(LOGFILE_BASE, logging.DEBUG, 'test_post_mmseq', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    print("teardown")

# Define common variables
FAMS_BY_STRAIN = \
    {'1': {'GEN2.1017.00001': ['GEN2.1017.00001.i0002_00004'],
           'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00002'],
           'GENO.1017.00001': ['GENO.1017.00001.b0002_00003'],
           'GENO.1216.00002': ["GENO.1216.00002.i0001_00003"]
           },
     '2': {'GEN2.1017.00001': ['GEN2.1017.00001.b0003_00010']},
     '3': {'GEN2.1017.00001': ['GEN2.1017.00001.b0004_00013']},
     '4': {'GEN2.1017.00001': ['GEN2.1017.00001.i0002_00005'],
           'GEN4.1111.00001': ["GEN4.1111.00001.b0001_00001"],
           'GENO.1017.00001': ["GENO.1017.00001.b0001_00002"],
           'GENO.1216.00002': ["GENO.1216.00002.b0001_00001", "GENO.1216.00002.i0001_00002"]},
     '5': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00003']},
     '6': {'GEN2.1017.00001': ['GEN2.1017.00001.b0004_00011'],
           'GEN4.1111.00001': ['GEN4.1111.00001.b0001_00009'],
           'GENO.1017.00001': ['GENO.1017.00001.b0002_00011'],
           'GENO.1216.00002': ['GENO.1216.00002.b0002_00010']
           },
     '7': {'GEN2.1017.00001': ['GEN2.1017.00001.b0002_00006'],
           'GENO.1017.00001': ["GENO.1017.00001.b0001_00001"]},
     '8': {'GEN2.1017.00001': ['GEN2.1017.00001.b0003_00007'],
           'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00006'],
           'GENO.1017.00001': ['GENO.1017.00001.i0002_00006', 'GENO.1017.00001.i0002_00007'],
           'GENO.1216.00002': ['GENO.1216.00002.i0001_00007']
           },
     '9': {'GEN2.1017.00001': ['GEN2.1017.00001.i0004_00012'],
           'GENO.1017.00001': ['GENO.1017.00001.i0002_00008']},
     '10': {'GEN2.1017.00001': ['GEN2.1017.00001.i0003_00008'],
            'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00007'],
            'GENO.1017.00001': ['GENO.1017.00001.i0002_00009'],
            'GENO.1216.00002': ['GENO.1216.00002.b0001_00008']
            },
     '11': {'GEN2.1017.00001': ['GEN2.1017.00001.i0003_00009'],
            'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00008'],
            'GENO.1017.00001': ['GENO.1017.00001.i0002_00010'],
            'GENO.1216.00002': ['GENO.1216.00002.b0002_00009']
            },
     '12': {'GEN2.1017.00001': ['GEN2.1017.00001.b0002_00003'],
            'GENO.1216.00002': ['GENO.1216.00002.i0001_00004']
            },
     '13': {'GEN2.1017.00001': ['GEN2.1017.00001.b0001_00002'],
            'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00004'],
            'GENO.1017.00001': ['GENO.1017.00001.i0002_00004'],
            'GENO.1216.00002': ['GENO.1216.00002.i0001_00005']},
     '14': {'GEN2.1017.00001': ['GEN2.1017.00001.b0001_00001'],
            'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00005'],
            'GENO.1017.00001': ['GENO.1017.00001.i0002_00005'],
            'GENO.1216.00002': ['GENO.1216.00002.i0001_00006']},
     '15': {'GENO.1216.00002': ['GENO.1216.00002.b0003_00011']},
     '16': {'GENO.1216.00002': ['GENO.1216.00002.b0003_00012']}
     }

FAMILIES = {'1': ['GEN2.1017.00001.i0002_00004', 'GEN4.1111.00001.i0001_00002',
                  'GENO.1017.00001.b0002_00003', 'GENO.1216.00002.i0001_00003'],
            '2': ['GEN2.1017.00001.b0003_00010'],
            '3': ['GEN2.1017.00001.b0004_00013'],
            '4': ['GEN2.1017.00001.i0002_00005', 'GEN4.1111.00001.b0001_00001',
                  'GENO.1017.00001.b0001_00002', 'GENO.1216.00002.b0001_00001',
                  'GENO.1216.00002.i0001_00002'],
            '5': ['GEN4.1111.00001.i0001_00003'],
            '6': ['GEN2.1017.00001.b0004_00011', 'GEN4.1111.00001.b0001_00009',
                  'GENO.1017.00001.b0002_00011', 'GENO.1216.00002.b0002_00010'],
            '7': ['GEN2.1017.00001.b0002_00006', 'GENO.1017.00001.b0001_00001'],
            '8': ['GEN2.1017.00001.b0003_00007', 'GEN4.1111.00001.i0001_00006',
                  'GENO.1017.00001.i0002_00006', 'GENO.1017.00001.i0002_00007',
                  'GENO.1216.00002.i0001_00007'],
            '9': ['GEN2.1017.00001.i0004_00012', 'GENO.1017.00001.i0002_00008'],
            '10': ['GEN2.1017.00001.i0003_00008', 'GEN4.1111.00001.i0001_00007',
                   'GENO.1017.00001.i0002_00009', 'GENO.1216.00002.b0001_00008'],
            '11': ['GEN2.1017.00001.i0003_00009', 'GEN4.1111.00001.i0001_00008',
                   'GENO.1017.00001.i0002_00010', 'GENO.1216.00002.b0002_00009'],
            '12': ['GEN2.1017.00001.b0002_00003', 'GENO.1216.00002.i0001_00004'],
            '13': ['GEN2.1017.00001.b0001_00002', 'GEN4.1111.00001.i0001_00004',
                   'GENO.1017.00001.i0002_00004', 'GENO.1216.00002.i0001_00005'],
            '14': ['GEN2.1017.00001.b0001_00001', 'GEN4.1111.00001.i0001_00005',
                   'GENO.1017.00001.i0002_00005', 'GENO.1216.00002.i0001_00006'],
            '15': ['GENO.1216.00002.b0003_00011'],
            '16': ['GENO.1216.00002.b0003_00012']
            }

ALL_STRAINS = ['GEN2.1017.00001', 'GEN4.1111.00001', 'GENO.1017.00001', 'GENO.1216.00002']

EXP_QUALIS = {'1': [1, 1, 1, 1], '2': [1, 0, 0, 0], '3': [1, 0, 0, 0], '4': [1, 1, 1, 1],
              '5': [0, 1, 0, 0], '6': [1, 1, 1, 1], '7': [1, 0, 1, 0], '8': [1, 1, 1, 1],
              '9': [1, 0, 1, 0], '10': [1, 1, 1, 1], '11': [1, 1, 1, 1], '12': [1, 0, 0, 1],
              '13': [1, 1, 1, 1], '14': [1, 1, 1, 1], '15': [0, 0, 0, 1], '16': [0, 0, 0, 1]}

EXP_QUANTIS = {'1': [1, 1, 1, 1], '2': [1, 0, 0, 0], '3': [1, 0, 0, 0], '4': [1, 1, 1, 2],
               '5': [0, 1, 0, 0], '6': [1, 1, 1, 1], '7': [1, 0, 1, 0], '8': [1, 1, 2, 1],
               '9': [1, 0, 1, 0], '10': [1, 1, 1, 1], '11': [1, 1, 1, 1], '12': [1, 0, 0, 1],
               '13': [1, 1, 1, 1], '14': [1, 1, 1, 1], '15': [0, 0, 0, 1], '16': [0, 0, 0, 1]}
                # nb sumquant sumqual nb0 nbmono nbmulti sum maxmulti
EXP_SUMS = {'1': [4, 4, 4, 0, 4, 0, 4, 1], '2': [1, 1, 1, 3, 1, 0, 4, 1],
            '3': [1, 1, 1, 3, 1, 0, 4, 1], '4': [5, 5, 4, 0, 3, 1, 4, 2],
            '5': [1, 1, 1, 3, 1, 0, 4, 1], '6': [4, 4, 4, 0, 4, 0, 4, 1],
            '7': [2, 2, 2, 2, 2, 0, 4, 1], '8': [5, 5, 4, 0, 3, 1, 4, 2],
            '9': [2, 2, 2, 2, 2, 0, 4, 1], '10': [4, 4, 4, 0, 4, 0, 4, 1],
            '11': [4, 4, 4, 0, 4, 0, 4, 1], '12': [2, 2, 2, 2, 2, 0, 4, 1],
            '13': [4, 4, 4, 0, 4, 0, 4, 1], '14': [4, 4, 4, 0, 4, 0, 4, 1],
            '15': [1, 1, 1, 3, 1, 0, 4, 1], '16': [1, 1, 1, 3, 1, 0, 4, 1]}

EXP_QUALIF = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst.quali_transpose.txt")
EXP_QUANTIF = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst.quanti_transpose.txt")
EXP_SUMF = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst.summary.txt")


def test_write_outputs(caplog):
    """
    Check that given some families, the qualitative and quantitative matrices,
    as well as the summary file are as expected.
    """
    caplog.set_level(logging.DEBUG)
    base = "test_write_out"
    pqlf = io.StringIO(base + ".quali_transpose.txt")
    pqtf = io.StringIO(base + ".quanti_transpose.txt")
    psf = io.StringIO(base + ".sum.txt")
    # run cmd
    res = post.generate_and_write_outputs(FAMS_BY_STRAIN, FAMILIES, 
                                          ALL_STRAINS, pqlf, pqtf, psf)
    # Check returned outputs
    (qualis, quantis, sums) = res
    assert qualis == EXP_QUALIS
    assert quantis == EXP_QUANTIS
    assert sums == EXP_SUMS
    # Check generated files
    # check content of matrix quali file
    with open(EXP_QUALIF, "r") as eq:
        for line_out, line_exp in zip(pqlf.getvalue().split("\n"), eq):
            assert line_out == line_exp.strip()
    # Check content of matrix quanti file
    with open(EXP_QUANTIF, "r") as eq:
        for line_out, line_exp in zip(pqtf.getvalue().split("\n"), eq):
            assert line_out == line_exp.strip()
    # Check content of summary file
    with open(EXP_SUMF, "r") as eq:
        next(eq)  # skip header
        for line_out, line_exp in zip(psf.getvalue().split("\n"), eq):
            assert line_out == line_exp.strip()

    # Check logs
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text
    assert caplog.records[0].levelname == "INFO"
    # Close io objects and discard memory buffers
    pqlf.close()
    pqtf.close()
    psf.close()


def test_open_out():
    """
    Check that given some families and a pagenome file, it creates 3 output files,
    with the expected content (quanti, quali, summary)
    """
    pangenome = os.path.join(GENEPATH, "test_open_out_pangenome.txt")
    res = post.open_outputs_to_write(FAMS_BY_STRAIN, FAMILIES, ALL_STRAINS, pangenome)

    # Check function output
    qualis, quantis, sums = res
    assert qualis == EXP_QUALIS
    assert quantis == EXP_QUANTIS
    assert sums == EXP_SUMS

    # Check presence and content of quali matrix file
    assert os.path.isfile(pangenome + ".quali.txt")
    assert tutil.compare_order_content(pangenome + ".quali.txt", EXP_QUALIF)   
    # Check presence and content of quanti matrix file
    assert os.path.isfile(pangenome + ".quanti.txt")
    assert tutil.compare_order_content(pangenome + ".quanti.txt", EXP_QUANTIF) 
    # Check presence and content of summary file
    assert os.path.isfile(pangenome + ".summary.txt")
    assert tutil.compare_order_content(pangenome + ".summary.txt", EXP_SUMF) 


def test_all_post():
    """
    Check that when running main method of post-treatment, it creates the 3 output files
    expected, with the expected content.
    """
    pangenome = os.path.join(GENEPATH, "test_all_post")
    post.post_treat(FAMILIES, pangenome)

    # Check presence and content of quali matrix file
    assert os.path.isfile(pangenome + ".quali.txt")
    assert tutil.compare_order_content(pangenome + ".quali.txt", EXP_QUALIF)   

    # Check presence and content of quanti matrix file
    assert os.path.isfile(pangenome + ".quanti.txt")
    assert tutil.compare_order_content(pangenome + ".quanti.txt", EXP_QUANTIF) 

    # Check presence and content of summary file
    assert os.path.isfile(pangenome + ".summary.txt")
    assert tutil.compare_order_content(pangenome + ".summary.txt", EXP_SUMF) 

    # Check that bin pangenome file was created (as it did not exist before)
    assert os.path.isfile(pangenome + ".bin")
    