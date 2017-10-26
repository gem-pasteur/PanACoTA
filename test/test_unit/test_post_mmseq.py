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

EXP_QUALIS = {'1': [1, 1, 1], '2': [1, 0, 0], '3': [1, 1, 1], '4': [0, 1, 0], '5': [1, 1, 1],
              '6': [1, 1, 1], '7': [0, 1, 0], '8': [1, 1, 1], '9': [0, 0, 1], '10': [1, 1, 1],
              '11': [1, 1, 1], '12': [1, 1, 1], '13': [0, 0, 1], '14': [0, 0, 1]}

EXP_QUANTIS = {'1': [1, 1, 2], '2': [1, 0, 0], '3': [1, 1, 1], '4': [0, 1, 0], '5': [1, 1, 1],
               '6': [1, 2, 1], '7': [0, 1, 0], '8': [1, 1, 1], '9': [0, 0, 1], '10': [1, 1, 1],
               '11': [1, 1, 1], '12': [1, 1, 1], '13': [0, 0, 1], '14': [0, 0, 1]}

EXP_SUMS = {'1': [4, 4, 3, 0, 2, 1, 3, 2], '2': [1, 1, 1, 2, 1, 0, 3, 1],
            '3': [3, 3, 3, 0, 3, 0, 3, 1], '4': [1, 1, 1, 2, 1, 0, 3, 1],
            '5': [3, 3, 3, 0, 3, 0, 3, 1], '6': [4, 4, 3, 0, 2, 1, 3, 2],
            '7': [1, 1, 1, 2, 1, 0, 3, 1], '8': [3, 3, 3, 0, 3, 0, 3, 1],
            '9': [1, 1, 1, 2, 1, 0, 3, 1], '10': [3, 3, 3, 0, 3, 0, 3, 1],
            '11': [3, 3, 3, 0, 3, 0, 3, 1], '12': [3, 3, 3, 0, 3, 0, 3, 1],
            '13': [1, 1, 1, 2, 1, 0, 3, 1], '14': [1, 1, 1, 2, 1, 0, 3, 1]}

EXP_QUALIF = os.path.join("test", "data", "pangenome", "exp_files", "exp_pangenome.txt.quali.txt")
EXP_QUANTIF = os.path.join("test", "data", "pangenome", "exp_files",
                           "exp_pangenome.txt.quanti.txt")
EXP_SUMF = os.path.join("test", "data", "pangenome", "exp_files", "exp_pangenome.txt.summary.txt")


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

# missing tests for utils_pangenome
# and functional tests for subcommands/pangenome
