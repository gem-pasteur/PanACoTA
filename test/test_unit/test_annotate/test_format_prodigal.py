#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for general_format_functions.py
"""

import os
import logging
import shutil
from io import StringIO
import pytest

import PanACoTA.annotate_module.format_prodigal as prodigalfunc
import PanACoTA.utils as utils
import test.test_unit.utilities_for_tests as tutil

ANNOTEDIR = os.path.join("test", "data", "annotate")
GENOMES_DIR = os.path.join(ANNOTEDIR, "genomes")
EXP_ANNOTE = os.path.join(ANNOTEDIR, "exp_files")
TEST_ANNOTE = os.path.join(ANNOTEDIR, "test_files")
GENEPATH = os.path.join(ANNOTEDIR, "generated_by_unit-tests")

LOGFILE_BASE = os.path.join(GENEPATH, "logs", "logfile")
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
    os.mkdir(GENEPATH)
    # utils.init_logger(LOGFILE_BASE, 0, 'test_fastme', verbose=1)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    # for f in LOGFILES:
    #     if os.path.exists(f):
    #         os.remove(f)
    print("teardown")


def test_create_gen_lst(caplog):
    """
    Check that generated gen and lst files are as expected.
    In the test file, all genomes have names different from gembase name
    This test file contains the following aspects:
    - gene in D strand (start < end)
    - gene in C strand (start > end)
    - CDS features
    - contigs with more than 2 genes
    - contig with only 2 genes (both 'b' loc)
    - contig with 1 gene ('b' loc)
    - contig without gene (should be skipped)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    genfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.ffn")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007"
               }
    name = "test.0417.00002"
    res_gen_file = os.path.join(GENEPATH, "prodigal_res.gen")
    res_lst_file = os.path.join(GENEPATH, "prodigal_res.lst")
    gpath = "original_genome_name"
    assert prodigalfunc.create_gene_lst(contigs, genfile, res_gen_file, res_lst_file, gpath,
                                        name, logger)
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    assert tutil.compare_order_content(exp_lst, res_lst_file)
    exp_gen = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.gen")
    assert tutil.compare_order_content(exp_gen, res_gen_file)


def test_create_gen_lst_cont_unknown(caplog):
    """
    A contig name in the gen file does not exist -> error message, and all result files
    must be removed for this genome
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    genfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.faa")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007"
               }
    name = "test.0417.00002"
    res_gen_file = os.path.join(GENEPATH, "prodigal_res.gen")
    res_lst_file = os.path.join(GENEPATH, "prodigal_res.lst")
    gpath = "original_genome_name"
    assert not prodigalfunc.create_gene_lst(contigs, genfile, res_gen_file, res_lst_file,
                                            gpath, name, logger)
    assert ("my contig found in test/data/annotate/test_files/original_name.fna-prodigalRes/"
            "prodigal_out_for_test.faa does not exist in original_genome_name") in caplog.text


def test_create_gff(caplog):
    """
    Check generated gff file. Must have all sequences in header (even replicons without gene),
    and then 1 line per gene
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.gff")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007"
               }
    sizes = {"test.0417.00002.0001": 14000,
             "test.0417.00002.0002": 5000,
             "test.0417.00002.0003": 4600,
             "test.0417.00002.0004": 8000,
             "test.0417.00002.0005": 1,
             "test.0417.00002.0006": 10,
             "test.0417.00002.0007": 15000,
            }
    name = "test.0417.00002"
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    gpath = "original_genome_name"
    assert prodigalfunc.create_gff(gpath, gfffile, res_gff_file, exp_lst, contigs, sizes)
    exp_gff = os.path.join(EXP_ANNOTE, "res_create_gff_prodigal.gff")
    assert tutil.compare_order_content(exp_gff, res_gff_file)


def test_create_gff_wrong_start(caplog):
    """
    Check that when trying to generate gff, if problem, exits with expected message
    Here, the start position in the gff file generated by prodigal is not the same as
    the start position in the lstinfo file (which was taken from the ffn
    file generated by prodigal)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test-wrong-start.gff")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007"
               }
    sizes = {"test.0417.00002.0001": 14000,
             "test.0417.00002.0002": 5000,
             "test.0417.00002.0003": 4600,
             "test.0417.00002.0004": 8000,
             "test.0417.00002.0005": 1,
             "test.0417.00002.0006": 10,
             "test.0417.00002.0007": 15000,
            }
    name = "test.0417.00002"
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    gpath = "original_genome_name"
    assert not prodigalfunc.create_gff(gpath, gfffile, res_gff_file, exp_lst, contigs, sizes)
    assert ("Files prodigal_out_for_test-wrong-start.ffn and "
            "prodigal_out_for_test-wrong-start.gff "
            "(in prodigal tmp_files: original_genome_name-prodigalRes) "
            "do not have the same start value for gene EPKOMDHM_00008 "
            "(78 in gff, 77 in ffn") in caplog.text


def test_create_prt(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    assert prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    exp_prt = os.path.join(EXP_ANNOTE, "res_create_prt_prodigal.faa")
    assert tutil.compare_order_content(exp_prt, res_prt_file)


def test_create_prt_wrong_lst(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-wrongformat.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Problem in format of lstline (1279\t2346\tCDS\ttest.0417.00002.0002i_00005\tNA\t"
            "| NA | NA | fefer | NA") in caplog.text


def test_create_prt_short_lst(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-shortlst.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("No more protein in lst file. We cannot get information on this "
            "protein (>toto_00011 # 2419 # 3000 # 1 # a)! Check that you do not have "
            "more proteins than genes in prodigal results") in caplog.text


def test_create_prt_end_not_int_lst(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-notint.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Start and/or end of protein test.0417.00002.0001i_00002 position is not a number "
            "(start = 4416; end = a6068)") in caplog.text


def test_create_prt_not_divisible3_lst(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-not-divisible3.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Gene test.0417.00002.0001b_00003 has a number of nucleotides (3001) "
            "that is not divisible by 3.") in caplog.text


def test_create_prt_not_moreprots_lst(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal_out_for_test.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-more-proteins.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Protein test.0417.00002.0007b_00015 is in .lst file but its sequence is not "
            "in the protein file generated by prodigal.") in caplog.text


def test_format_1genome(caplog):
    """
    Test that when prodigal results are ok, all files are
    generated as expected.
    """
# # def test_tbl_to_lst_not_changed_names(caplog):
#     """
#     Check that generated lstinfo file is as expected, when the genome name is the same as
#     it already was in the genome given to prokka.
#     The test tblfile contains the following aspects:
#     - gene in D strand (start < end)
#     - gene in C strand (start > end)
#     - CDS features (some with all info = ECnumber, gene name, product etc. ;
#     some with missing info)
#     - tRNA type
#     - repeat_region type (*2)
#     - contigs with more than 2 genes
#     - contig with only 2 genes (both 'b' loc)
#     - contig with 1 gene ('b' loc)
#     - contig without gene (should be skipped)
#     """
#     caplog.set_level(logging.DEBUG)
#     logger = logging.getLogger("test_prokka")
#     tblfile = os.path.join(TEST_ANNOTE, "original_name.fna-prokkaRes", "prokka_out_for_test.tbl")
#     lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
#     contigs = ["test.0417.00002.0001\t50", "test.0417.00002.0002\t50", "test.0417.00002.0003\t50",
#                "test.0417.00002.0004\t50", "test.0417.00002.0005\t50", "test.0417.00002.0006\t50",
#                "test.0417.00002.0007\t50"]
#     name = "test.0417.00002"
#     assert prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, logger, changed_name=False)
#     exp_lst = os.path.join(EXP_ANNOTE, "res_tbl2lst.lst")
#     assert tutil.compare_order_content(exp_lst, lstfile)

