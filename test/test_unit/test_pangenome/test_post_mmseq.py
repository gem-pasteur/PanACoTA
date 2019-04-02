#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the mmseqs_functions submodule in pangenome module
"""

import os
import io
import logging

import genomeAPCAT.pangenome_module.post_treatment as post
import genomeAPCAT.utils as utils


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

EXP_QUALIF = os.path.join("test", "data", "pangenome", "exp_files",
                          "exp_pangenome-4genomes.lst.quali.txt")
EXP_QUANTIF = os.path.join("test", "data", "pangenome", "exp_files",
                           "exp_pangenome-4genomes.lst.quanti.txt")
EXP_SUMF = os.path.join("test", "data", "pangenome", "exp_files",
                        "exp_pangenome-4genomes.lst.summary.txt")


def test_write_outputs():
    """
    Check that given some families, the qualitative and quantitative matrices,
    as well as the summary file are as expected.
    """
    base = "test_write_out"
    pqlf = io.StringIO(base + ".quali.txt")
    pqtf = io.StringIO(base + ".quanti.txt")
    psf = io.StringIO(base + ".sum.txt")
    res = post.generate_and_write_outputs(FAMS_BY_STRAIN, FAMILIES, ALL_STRAINS, pqlf, pqtf, psf)
    (qualis, quantis, sums) = res
    assert qualis == EXP_QUALIS
    assert quantis == EXP_QUANTIS
    assert sums == EXP_SUMS
    # check content of matrix quali file
    with open(EXP_QUALIF, "r") as eq:
        eq.readline()  # skip header
        for line_out, line_exp in zip(pqlf.getvalue().split("\n"), eq):
            assert line_out == line_exp.strip()
    # Check content of matrix quanti file
    with open(EXP_QUANTIF, "r") as eq:
        eq.readline()  # skip header
        for line_out, line_exp in zip(pqtf.getvalue().split("\n"), eq):
            assert line_out == line_exp.strip()
    # Check content of summary file
    with open(EXP_SUMF, "r") as eq:
        eq.readline()  # skip header
        for line_out, line_exp in zip(psf.getvalue().split("\n"), eq):
            assert line_out == line_exp.strip()
    # Close io objects and discard memory buffers
    pqlf.close()
    pqtf.close()
    psf.close()


def test_open_out():
    """
    Check that given some families and a pagenome file, it creates 3 output files,
    with the expected content (quanti, quali, summary)
    """
    pangenome = "test_open_out_pangenome.txt"
    res = post.open_outputs_to_write(FAMS_BY_STRAIN, FAMILIES, ALL_STRAINS, pangenome)
    qualis, quantis, sums = res
    assert qualis == EXP_QUALIS
    assert quantis == EXP_QUANTIS
    assert sums == EXP_SUMS

    # Check presence and content of quali matrix file
    assert os.path.isfile(pangenome + ".quali.txt")
    with open(pangenome + ".quali.txt", "r") as panf, open(EXP_QUALIF, "r") as eq:
        for line_out, line_exp in zip(panf, eq):
            assert line_out == line_exp
    os.remove(pangenome + ".quali.txt")

    # Check presence and content of quanti matrix file
    assert os.path.isfile(pangenome + ".quanti.txt")
    with open(pangenome + ".quanti.txt", "r") as panf, open(EXP_QUANTIF, "r") as eq:
        for line_out, line_exp in zip(panf, eq):
            assert line_out == line_exp
    os.remove(pangenome + ".quanti.txt")

    # Check presence and content of summary file
    assert os.path.isfile(pangenome + ".summary.txt")
    with open(pangenome + ".summary.txt", "r") as panf, open(EXP_SUMF, "r") as eq:
        for line_out, line_exp in zip(panf, eq):
            assert line_out == line_exp
    os.remove(pangenome + ".summary.txt")


def test_all_post():
    """
    Check that when running main method of post-treatment, it creates the 3 output files
    expected, with the expected content.
    """
    logger = my_logger()
    logname = logger(name="test_all_post")
    pangenome = "test_all_post"
    post.post_treat(FAMILIES, pangenome)

    # Check presence and content of quali matrix file
    assert os.path.isfile(pangenome + ".quali.txt")
    with open(pangenome + ".quali.txt", "r") as panf, open(EXP_QUALIF, "r") as eq:
        for line_out, line_exp in zip(panf, eq):
            assert line_out == line_exp
    os.remove(pangenome + ".quali.txt")

    # Check presence and content of quanti matrix file
    assert os.path.isfile(pangenome + ".quanti.txt")
    with open(pangenome + ".quanti.txt", "r") as panf, open(EXP_QUANTIF, "r") as eq:
        for line_out, line_exp in zip(panf, eq):
            assert line_out == line_exp
    os.remove(pangenome + ".quanti.txt")

    # Check presence and content of summary file
    assert os.path.isfile(pangenome + ".summary.txt")
    with open(pangenome + ".summary.txt", "r") as panf, open(EXP_SUMF, "r") as eq:
        for line_out, line_exp in zip(panf, eq):
            assert line_out == line_exp
    os.remove(pangenome + ".summary.txt")

    # Check that bin pangenome file was created (as it did not exist before)
    assert os.path.isfile(pangenome + ".bin")
    os.remove(pangenome + ".bin")

    os.remove(logname + ".log")
    os.remove(logname + ".log.err")
    os.remove(logname + ".log.details")
