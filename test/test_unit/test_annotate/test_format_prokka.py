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

from PanACoTA.annotate_module import format_prokka as prokkafunc
import PanACoTA.utils as utils
import test.test_unit.utilities_for_tests as tutil

ANNOTEDIR = os.path.join("test", "data", "annotate")
EXP_ANNOTE = os.path.join(ANNOTEDIR, "exp_files")
TEST_ANNOTE = os.path.join(ANNOTEDIR, "test_files")
GENEPATH = os.path.join(ANNOTEDIR, "generated_by_unit-tests")

LOGFILE_BASE = os.path.join(GENEPATH, "logfile")
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
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    print("teardown")


def test_tbl_to_lst_not_changed_names(caplog):
    """
    Check that generated lstinfo file is as expected, when the genome name is the same as
    it already was in the genome given to prokka.
    The test tblfile contains the following aspects:
    - gene in D strand (start < end)
    - gene in C strand (start > end)
    - CDS features (some with all info = ECnumber, gene name, product etc. ;
    some with missing info)
    - tRNA type
    - repeat_region type (*2) -> should be ignored in .lst
            * 1 in prokka1 version (start end repeat_region
                                        rpt_family CRISPR
                                        score 7)
            * 1 in prokka2 version (start end CRISPR
                                        note CRISPR with x repeat units
                                        rpt_family CRISPR)
    - contigs with more than 2 genes
    - contig with only 2 genes (both 'b' loc)
    - contig with 1 gene ('b' loc)
    - contig without gene (should be skipped)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prokka")
    tblfile = os.path.join(TEST_ANNOTE, "original_name.fna-prokkaRes", "prokka_out_for_test.tbl")
    lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007",
              }
    name = "test.0417.00002"
    gpath = "genome_init"
    assert prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, gpath)
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    assert tutil.compare_order_content(exp_lst, lstfile)


def test_tbl_to_lst_changed_names(caplog):
    """
    Check that generated lstinfo file is as expected, when the genome name is not the same as
    it already was in the genome given to prokka.
    The test tblfile contains the following aspects:
    - gene in D strand (start < end)
    - gene in C strand (start > end)
    - CDS features (some with all info = ECnumber, gene name, product etc. ;
    some with missing info)
    - tRNA type
    - repeat_region type (*2)
    - contigs with more than 2 genes
    - contig with only 2 genes (both 'b' loc)
    - contig with 1 gene ('b' loc)
    - contig without gene (should be skipped)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prokka")
    tblfile = os.path.join(TEST_ANNOTE, "prokka_out_tbl_changed-contnames.tbl")
    lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
    contigs = {"toto_1": "test.0417.00002.0001",
               "toto_2": "test.0417.00002.0002",
               "toto_3": "test.0417.00002.0003",
               "toto_4": "test.0417.00002.0004",
               "toto_5": "test.0417.00002.0005",
               "toto_6": "test.0417.00002.0006",
               "toto_7": "test.0417.00002.0007",
              }
    name = "test.0417.00002"
    gpath = "path_to_genome"
    assert prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, gpath)
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    assert tutil.compare_order_content(exp_lst, lstfile)


def test_tbl_to_lst_unknown(caplog):
    """
    A contig name in the fna file does not exist -> error message, and all result files
    must be removed for this genome
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prokka")
    tblfile = os.path.join(TEST_ANNOTE, "prokka_out_tbl_changed-contnames.tbl")
    lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
    contigs = {"toto_1": "test.0417.00002.0001",
               "toto_2": "test.0417.00002.0002",
               "toto_3": "test.0417.00002.0003",
               "toto_5": "test.0417.00002.0004",
               "toto_6": "test.0417.00002.0006",
               "toto_7": "test.0417.00002.0007",
              }
    name = "test.0417.00002"
    gpath = "path_to_genome"
    assert not prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, gpath)
    assert ("'toto_4' found in test/data/annotate/test_files/prokka_out_tbl_changed-contnames.tbl "
            "does not exist in path_to_genome") in caplog.text


def test_tbl_to_lst_empty_tbl(caplog):
    """
    The tbl file is empty. Error message saying that the file is empty, and returns False
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prokka")
    tblfile = os.path.join(GENEPATH, "prokka_out_empty_tbl.tbl")
    open(tblfile, "w").close()
    lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
    contigs = {"toto_1": "test.0417.00002.0001",
               "toto_2": "test.0417.00002.0002",
               "toto_3": "test.0417.00002.0003",
               "toto_5": "test.0417.00002.0004",
               "toto_6": "test.0417.00002.0006",
               "toto_7": "test.0417.00002.0007",
              }
    name = "test.0417.00002"
    gpath = "path_to_genome"
    assert not prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, gpath)
    assert ("test/data/annotate/generated_by_unit-tests/"
            "prokka_out_empty_tbl.tbl is empty.") in caplog.text


def test_tbl_to_lst_wrong_tbl(caplog):
    """
    The tbl file is in wrong format.
    Error message saying that it is not the expected format, and returns False
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prokka")
    tblfile = os.path.join(GENEPATH, "prokka_out_empty_tbl.tbl")
    with open(tblfile, "w") as tblf:
        tblf.write("wrong format")
    lstfile = os.path.join(GENEPATH, "res_test_tbl2lst.lst")
    contigs = {"toto_1": "test.0417.00002.0001",
               "toto_2": "test.0417.00002.0002",
               "toto_3": "test.0417.00002.0003",
               "toto_5": "test.0417.00002.0004",
               "toto_6": "test.0417.00002.0006",
               "toto_7": "test.0417.00002.0007",
              }
    name = "test.0417.00002"
    gpath = "path_to_genome"
    assert not prokkafunc.tbl2lst(tblfile, lstfile, contigs, name, gpath)
    assert ("Wrong format for test/data/annotate/generated_by_unit-tests/"
            "prokka_out_empty_tbl.") in caplog.text


def test_create_gff(caplog):
    """
    Check generated gff file. Must have all sequences in header (even replicons without gene),
    and then 1 line per gene
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "original_name.fna-prokkaRes",
                           "prokka_out_for_test.gff")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007",
              }
    sizes = {"test.0417.00002.0001": 84,
             "test.0417.00002.0002": 103,
             "test.0417.00002.0003": 122,
             "test.0417.00002.0004": 35,
             "test.0417.00002.0005": 198,
             "test.0417.00002.0006": 128,
             "test.0417.00002.0007": 85
            }
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    res_lst = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    gpath = "original_genome_name"
    assert prokkafunc.generate_gff(gpath, gfffile, res_gff_file, res_lst, sizes, contigs)
    exp_gff = os.path.join(EXP_ANNOTE, "res_create_gff-prokka.gff")
    assert tutil.compare_order_content(exp_gff, res_gff_file)


def test_create_gff_wrong_start(caplog):
    """
    Check generated gff file. A sequence has not the same start in gff and in tbl
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "prokka_out_gff-wrong-start.gff")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007",
              }
    sizes = {"test.0417.00002.0001": 84,
             "test.0417.00002.0002": 103,
             "test.0417.00002.0003": 122,
             "test.0417.00002.0004": 35,
             "test.0417.00002.0005": 198,
             "test.0417.00002.0006": 128,
             "test.0417.00002.0007": 85
            }
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    res_lst = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    gpath = "original_genome_name"
    assert not prokkafunc.generate_gff(gpath, gfffile, res_gff_file, res_lst, sizes, contigs)
    assert ("Files prokka_out_gff-wrong-start.tbl and prokka_out_gff-wrong-start.gff "
            "(in prokka tmp_files: original_genome_name-prokkaRes) do not have the same "
            "start value for gene EPKOMDHM_00011 (3000 in gff, 3500 in tbl)") in caplog.text


def test_create_gff_error_gff(caplog):
    """
    Check generated gff file. The prokka output gff file has a problem (locus_tag != ID)
    -> returns False with error message
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "prokka_out_gff-error.gff")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007",
              }
    sizes = {"test.0417.00002.0001": 84,
             "test.0417.00002.0002": 103,
             "test.0417.00002.0003": 122,
             "test.0417.00002.0004": 35,
             "test.0417.00002.0005": 198,
             "test.0417.00002.0006": 128,
             "test.0417.00002.0007": 85
            }
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    res_lst = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    gpath = "original_genome_name"
    assert not prokkafunc.generate_gff(gpath, gfffile, res_gff_file, res_lst, sizes, contigs)
    assert ("Problem in prokka_out_gff-error.gff: ID=EPKOMDHM_00005 whereas "
            "locus_tag=toto") in caplog.text


def test_create_gff_wrong_format(caplog):
    """
    Check generated gff file. The prokka output gff file has a problem (locus_tag != ID)
    -> returns False with error message
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(GENEPATH, "prokka_out_wrong_format.gff")
    with open(gfffile, "w") as gfff:
        gfff.write("##gff-version3\n")
        gfff.write("##sequence-region bis 1 600\n")
        gfff.write("JGIKIPgffgIJ    Prodigal:2.6    CDS 287 787\n")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007",
              }
    sizes = {"test.0417.00002.0001": 84,
             "test.0417.00002.0002": 103,
             "test.0417.00002.0003": 122,
             "test.0417.00002.0004": 35,
             "test.0417.00002.0005": 198,
             "test.0417.00002.0006": 128,
             "test.0417.00002.0007": 85
            }
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    res_lst = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    gpath = "original_genome_name"
    assert not prokkafunc.generate_gff(gpath, gfffile, res_gff_file, res_lst, sizes, contigs)
    assert ("Wrong format for test/data/annotate/generated_by_unit-tests/"
            "prokka_out_wrong_format.gff") in caplog.text


def test_create_gen(caplog):
    """
    Check create gen file.
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    ffnfile = os.path.join(TEST_ANNOTE, "original_name.fna-prokkaRes",
                           "prokka_out_for_test.ffn")
    lstfile = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    res_gen_file = os.path.join(GENEPATH, "prodigal_res.gen")
    assert prokkafunc.create_gen(ffnfile, lstfile, res_gen_file)
    exp_gen = os.path.join(EXP_ANNOTE, "res_create_gene_prokka.gen")
    assert tutil.compare_order_content(exp_gen, res_gen_file)


def test_create_gen_supgen(caplog):
    """
    Check create gen file. But there is a gene in ffn that is not in lst -> error
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    ffnfile = os.path.join(TEST_ANNOTE, "prokka_out_for_test-supGene.ffn")
    lstfile = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    res_gen_file = os.path.join(GENEPATH, "prodigal_res.gen")
    assert not prokkafunc.create_gen(ffnfile, lstfile, res_gen_file)
    assert ("Missing info for gene >JGIKIPIJ_03050 (from test/data/annotate/test_files/prokka_out_for_test-supGene.ffn) in test/data/annotate/exp_files/res_create_lst-prokka.lst. If it is actually present in the lst file, check that genes are ordered by increasing number in both lst and ffn files.") in caplog.text


def test_create_gen_missingSeq(caplog):
    """
    Check create gen file. A gene in lst does not have a sequence in ffn.
    Just skip it, and go to next sequence for gen file.
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    ffnfile = os.path.join(TEST_ANNOTE, "prokka_out_for_test-noSeqFor1gene.ffn")
    lstfile = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    res_gen_file = os.path.join(GENEPATH, "prodigal_res.gen")
    assert prokkafunc.create_gen(ffnfile, lstfile, res_gen_file)
    exp_gen = os.path.join(EXP_ANNOTE, "res_create_gene_prokka-missGene.gen")
    assert tutil.compare_order_content(exp_gen, res_gen_file)


def test_create_prt(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prokkaRes",
                           "prokka_out_for_test.faa")
    res_prt_file = os.path.join(GENEPATH, "prokka_res.prt")
    lstfile = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    assert prokkafunc.create_prt(protfile, lstfile, res_prt_file)
    exp_prt = os.path.join(EXP_ANNOTE, "res_create_prt_prokka.faa")
    assert tutil.compare_order_content(exp_prt, res_prt_file)


def test_create_prt_wrong_header_int(caplog):
    """
    Test creating prt file, but the faa file has a header with wrong format (>JGIKIPIJ_d0008)
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "prokka_out_for_test-wrongHeaderInt.faa")
    res_prt_file = os.path.join(GENEPATH, "prokka_res.prt")
    lstfile = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    assert not prokkafunc.create_prt(protfile, lstfile, res_prt_file)
    assert ("Unknown header format >JGIKIPIJ_d0008 in test/data/annotate/test_files/"
            "prokka_out_for_test-wrongHeaderInt.faa. Gene ID is not a number.") in caplog.text


def test_create_prt_wrong_header_sep(caplog):
    """
    Test creating prt file, but the faa file has a header with wrong format,
    no '_' between name and gene ID (>JGIKIPIJ00008)
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "prokka_out_for_test-wrongHeaderSep.faa")
    res_prt_file = os.path.join(GENEPATH, "prokka_res.prt")
    lstfile = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    assert not prokkafunc.create_prt(protfile, lstfile, res_prt_file)
    assert ("Unknown header format >JGIKIPIJ00008 in test/data/annotate/test_files/"
            "prokka_out_for_test-wrongHeaderSep.faa. Gene ID is not a number.") in caplog.text


def test_create_prt_wrong_unknown_prot(caplog):
    """
    Test creating prt file, but the faa file has a protein (>sup-prot_00012)
    which is not in the lst file
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "prokka_out_for_test-supHeader.faa")
    res_prt_file = os.path.join(GENEPATH, "prokka_res.prt")
    lstfile = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    assert not prokkafunc.create_prt(protfile, lstfile, res_prt_file)
    assert ("Missing info for protein >sup-prot_00012 (from test/data/annotate/test_files/"
            "prokka_out_for_test-supHeader.faa) in test/data/annotate/exp_files/"
            "res_create_lst-prokka.lst. If it is actually present in the lst file, check that "
            "proteins are ordered by increasing number in both lst and faa files.") in caplog.text


def test_format_1genome(caplog):
    """
    Test that when prokka results are ok, all files are generated as expected.
    """
    caplog.set_level(logging.DEBUG)
    name = "test.0417.00002"
    # path to original genome, given to prodigal for annotation
    gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    prok_path = TEST_ANNOTE
    # Create result directories
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gene_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff")
    os.makedirs(prot_dir)
    os.makedirs(lst_dir)
    os.makedirs(rep_dir)
    os.makedirs(gene_dir)
    os.makedirs(gff_dir)

    assert prokkafunc.format_one_genome(gpath, name, prok_path, lst_dir, prot_dir, gene_dir,
                                        rep_dir, gff_dir)

    # Check output files content
    # Replicons
    exp_rep = os.path.join(EXP_ANNOTE, "res_created_rep-prokka.fna")
    res_rep_file = os.path.join(rep_dir, "test.0417.00002.fna")
    assert tutil.compare_order_content(exp_rep, res_rep_file)
    # Proteins
    exp_prt = os.path.join(EXP_ANNOTE, "res_create_prt_prokka.faa")
    res_prt_file = os.path.join(prot_dir, "test.0417.00002.prt")
    assert tutil.compare_order_content(exp_prt, res_prt_file)
    # Genes
    exp_gen = os.path.join(EXP_ANNOTE, "res_create_gene_prokka.gen")
    res_gen_file = os.path.join(gene_dir, "test.0417.00002.gen")
    assert tutil.compare_order_content(exp_gen, res_gen_file)
    # LSTINFO
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_lst-prokka.lst")
    res_lst_file = os.path.join(lst_dir, "test.0417.00002.lst")
    assert tutil.compare_order_content(exp_lst, res_lst_file)
    # gff
    exp_gff = os.path.join(EXP_ANNOTE, "res_create_gff-prokka.gff")
    res_gff_file = os.path.join(gff_dir, "test.0417.00002.gff")
    assert tutil.compare_order_content(exp_gff, res_gff_file)


def test_format_1genome_emptygpath(caplog):
    """
    Test on formatting prokka results, when original ffn file is empty -> error message,
    and no file generated
    """
    caplog.set_level(logging.DEBUG)
    name = "prokka_out_for_test"
    # Create empty file, that we give to prodigal for formatting step
    gpath =  os.path.join(GENEPATH, "original_name-empty.fna")
    open(gpath, "w").close()
    # Create prokka result files (empty files, will not be read)
    gpath_prokres =  gpath + "-prokkaRes"
    os.makedirs(gpath_prokres)
    fna_prokres = os.path.join(gpath_prokres, "prokka_out_for_test.fna")
    open(fna_prokres, "w").close()
    tbl_prokres = os.path.join(gpath_prokres, "prokka_out_for_test.tbl")
    open(tbl_prokres, "w").close()
    gff_prokres = os.path.join(gpath_prokres, "prokka_out_for_test.gff")
    open(gff_prokres, "w").close()
    ffn_prokres = os.path.join(gpath_prokres, "prokka_out_for_test.ffn")
    open(ffn_prokres, "w").close()
    faa_prokres = os.path.join(gpath_prokres, "prokka_out_for_test.faa")
    open(faa_prokres, "w").close()
    # Create result directories
    prok_path = GENEPATH
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gen_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff")
    os.makedirs(rep_dir)
    os.makedirs(gff_dir)
    os.makedirs(lst_dir)
    os.makedirs(gen_dir)
    # Add empty res lst, gff and gen files, to check that it is removed at the end
    res_gff_file = os.path.join(gff_dir, "prokka_out_for_test.gff")
    open(res_gff_file, "w").close()
    assert len(os.listdir(gff_dir) ) == 1
    res_lst_file = os.path.join(lst_dir, "prokka_out_for_test.lst")
    open(res_lst_file, "w").close()
    assert len(os.listdir(lst_dir) ) == 1
    res_gen_file = os.path.join(gen_dir, "prokka_out_for_test.gen")
    open(res_gen_file, "w").close()
    assert len(os.listdir(gen_dir) ) == 1
    # res_gen_file = os.path.join(gen_dir, "prokka_out_for_test.gen")
    # open(res_gen_file, "w").close()
    # assert len(os.listdir(gen_dir) ) == 1

    assert not prokkafunc.format_one_genome(gpath, name, prok_path, lst_dir, prot_dir, gen_dir,
                                              rep_dir, gff_dir)
    # Check that all files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(gen_dir) ) == 0
    # Check log
    assert ("Problems while generating Replicon file for prokka_out_for_test") in caplog.text


def test_format_1genome_pb_tbl(caplog):
    """
    Test on formatting prokka results, when prokka output tbl file does not have
    the expected format -> error message, and no file generated
    """
    caplog.set_level(logging.DEBUG)
    name = "test.0417.00002"
    # path to original genome, given to prodigal for annotation
    orig_gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    # In generated_by_tests folder, create the original genome given to prokka
    # (copy from test_file)
    used_gpath = os.path.join(GENEPATH, "original_name.fna")
    used_respath = used_gpath + "-prokkaRes"
    os.makedirs(used_respath)
    shutil.copyfile(orig_gpath, used_gpath)

    # Create tbl_file with a wrong format
    with open(os.path.join(used_respath, "prokka_out_for_test.tbl"), "w") as ori:
        ori.write(">wrongheader # 1 # 2 # 1 # toto")
    # Add empty prokka res gff ffn and faa files (they won't be read, as it will stop
    # at tbl2lst)
    orig_fna =  os.path.join(TEST_ANNOTE, "original_name.fna-prokkaRes", "prokka_out_for_test.fna")
    fna_prokres = os.path.join(used_respath, "prokka_out_for_test.fna")
    shutil.copyfile(orig_fna, fna_prokres)
    res_gff_file = os.path.join(used_respath, "prokka_out_for_test.gff")
    open(res_gff_file, "w").close()
    res_ffn_file = os.path.join(used_respath, "prokka_out_for_test.ffn")
    open(res_ffn_file, "w").close()
    res_faa_file = os.path.join(used_respath, "prokka_out_for_test.faa")
    open(res_faa_file, "w").close()

    # Create output directories
    prok_path = GENEPATH
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gen_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff")
    os.makedirs(rep_dir)
    os.makedirs(gff_dir)
    os.makedirs(lst_dir)
    os.makedirs(gen_dir)
    # Add empty res lst, gff and gen files, to check that it is removed at the end
    res_gff_file = os.path.join(gff_dir, "test.0417.00002.gff")
    open(res_gff_file, "w").close()
    assert len(os.listdir(gff_dir) ) == 1
    res_lst_file = os.path.join(lst_dir, "test.0417.00002.lst")
    open(res_lst_file, "w").close()
    assert len(os.listdir(lst_dir) ) == 1
    res_gen_file = os.path.join(gen_dir, "test.0417.00002.gen")
    open(res_gen_file, "w").close()
    assert len(os.listdir(gen_dir) ) == 1

    # Run formatting
    assert not prokkafunc.format_one_genome(used_gpath, name, prok_path, lst_dir, prot_dir,
                                            gen_dir, rep_dir, gff_dir)

    # Check that all files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(gen_dir) ) == 0
    assert("Wrong format for test/data/annotate/generated_by_unit-tests/"
           "original_name.fna-prokkaRes/prokka_out_for_test.tbl.") in caplog.text
    assert ("Problems while generating LSTINFO file for test.0417.00002") in caplog.text


def test_format_1genome_pb_gff(caplog):
    """
    Test on formatting prokka results, when prokka output gff file does not have
    the expected format -> error message, and no file generated
    """
    caplog.set_level(logging.DEBUG)
    name = "test.0417.00002"
    # path to original genome, given to prodigal for annotation
    orig_gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    # In generated_by_tests folder, create the original genome given to prokka
    # (copy from test_file)
    used_gpath = os.path.join(GENEPATH, "original_name.fna")
    used_respath = used_gpath + "-prokkaRes"
    os.makedirs(used_respath)
    shutil.copyfile(orig_gpath, used_gpath)
    # Copy tbl file, which is as expected (tbl2lst must succeed)
    orig_tbl = os.path.join(orig_gpath + "-prokkaRes",
                            "prokka_out_for_test.tbl")
    used_tbl = os.path.join(used_respath, "prokka_out_for_test.tbl")
    shutil.copyfile(orig_tbl, used_tbl)
    orig_fna =  os.path.join(orig_gpath + "-prokkaRes", "prokka_out_for_test.fna")
    fna_prokres = os.path.join(used_respath, "prokka_out_for_test.fna")
    shutil.copyfile(orig_fna, fna_prokres)

    # Create gff_file with a wrong format
    with open(os.path.join(used_respath, "prokka_out_for_test.gff"), "w") as ori:
        ori.write("wrongheader # 1 # 2 # 1 # toto")
    # Add empty prokka res ffn and faa files
    res_ffn_file = os.path.join(used_respath, "prokka_out_for_test.ffn")
    open(res_ffn_file, "w").close()
    res_faa_file = os.path.join(used_respath, "prokka_out_for_test.faa")
    open(res_faa_file, "w").close()

    # Create output directories
    prok_path = GENEPATH
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gen_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff")
    os.makedirs(rep_dir)
    os.makedirs(gff_dir)
    os.makedirs(lst_dir)
    os.makedirs(gen_dir)
    # Add empty res lst, gff and gen files, to check that it is removed at the end
    res_gff_file = os.path.join(gff_dir, "test.0417.00002.gff")
    open(res_gff_file, "w").close()
    assert len(os.listdir(gff_dir) ) == 1
    res_lst_file = os.path.join(lst_dir, "test.0417.00002.lst")
    open(res_lst_file, "w").close()
    assert len(os.listdir(lst_dir) ) == 1
    res_gen_file = os.path.join(gen_dir, "test.0417.00002.gen")
    open(res_gen_file, "w").close()
    assert len(os.listdir(gen_dir) ) == 1

    # Run formatting
    assert not prokkafunc.format_one_genome(used_gpath, name, prok_path, lst_dir, prot_dir,
                                            gen_dir, rep_dir, gff_dir)

    # Check that all files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(gen_dir) ) == 0
    assert("Wrong format for test/data/annotate/generated_by_unit-tests/"
           "original_name.fna-prokkaRes/prokka_out_for_test.gff.") in caplog.text
    assert ("Problems while generating .gff file for test.0417.00002") in caplog.text


def test_format_1genome_pb_ffn(caplog):
    """
    Test on formatting prokka results, when prokka output ffn file does not have
    the expected format -> error message, and no file generated
    """
    caplog.set_level(logging.DEBUG)
    name = "test.0417.00002"
    # path to original genome, given to prodigal for annotation
    orig_gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    # In generated_by_tests folder, create the original genome given to prokka
    # (copy from test_file)
    used_gpath = os.path.join(GENEPATH, "original_name.fna")
    used_respath = used_gpath + "-prokkaRes"
    os.makedirs(used_respath)
    shutil.copyfile(orig_gpath, used_gpath)
    # Copy tbl and gff files, which is as expected (tbl2lst and generate_gff must succeed)
    orig_tbl = os.path.join(orig_gpath + "-prokkaRes",
                            "prokka_out_for_test.tbl")
    used_tbl = os.path.join(used_respath, "prokka_out_for_test.tbl")
    shutil.copyfile(orig_tbl, used_tbl)
    orig_gff = os.path.join(orig_gpath + "-prokkaRes",
                            "prokka_out_for_test.gff")
    used_gff = os.path.join(used_respath, "prokka_out_for_test.gff")
    shutil.copyfile(orig_gff, used_gff)
    orig_fna =  os.path.join(orig_gpath + "-prokkaRes", "prokka_out_for_test.fna")
    fna_prokres = os.path.join(used_respath, "prokka_out_for_test.fna")
    shutil.copyfile(orig_fna, fna_prokres)

    # Create ffn_file with a wrong format
    orig_ffn = os.path.join(TEST_ANNOTE, "prokka_out_for_test-supGene.ffn")
    used_ffn = os.path.join(used_respath, "prokka_out_for_test.ffn")
    shutil.copyfile(orig_ffn, used_ffn)
    # Add empty prokka res faa file
    res_faa_file = os.path.join(used_respath, "prokka_out_for_test.faa")
    open(res_faa_file, "w").close()

    # Create output directories
    prok_path = GENEPATH
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gen_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff")
    os.makedirs(rep_dir)
    os.makedirs(gff_dir)
    os.makedirs(lst_dir)
    os.makedirs(gen_dir)
    # Add empty res lst, gff and gen files, to check that it is removed at the end
    res_gff_file = os.path.join(gff_dir, "test.0417.00002.gff")
    open(res_gff_file, "w").close()
    assert len(os.listdir(gff_dir) ) == 1
    res_lst_file = os.path.join(lst_dir, "test.0417.00002.lst")
    open(res_lst_file, "w").close()
    assert len(os.listdir(lst_dir) ) == 1
    res_gen_file = os.path.join(gen_dir, "test.0417.00002.gen")
    open(res_gen_file, "w").close()
    assert len(os.listdir(gen_dir) ) == 1

    # Run formatting
    assert not prokkafunc.format_one_genome(used_gpath, name, prok_path, lst_dir, prot_dir,
                                            gen_dir, rep_dir, gff_dir)

    # Check that all files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(gen_dir) ) == 0
    assert("Missing info for gene >JGIKIPIJ_03050 (from test/data/annotate/"
           "generated_by_unit-tests/original_name.fna-prokkaRes/prokka_out_for_test.ffn) "
           "in test/data/annotate/generated_by_unit-tests/LSTINFO/test.0417.00002.lst. "
           "If it is actually present in the lst file, "
           "check that genes are ordered by increasing number in both lst and "
           "ffn files.") in caplog.text
    assert ("Problems while generating .gen file for test.0417.00002") in caplog.text


def test_format_1genome_pb_faa(caplog):
    """
    Test on formatting prokka results, when prokka output faa file does not have
    the expected format -> error message, and no file generated
    """
    caplog.set_level(logging.DEBUG)
    name = "test.0417.00002"
    # path to original genome, given to prodigal for annotation
    orig_gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    # In generated_by_tests folder, create the original genome given to prokka
    # (copy from test_file)
    used_gpath = os.path.join(GENEPATH, "original_name.fna")
    used_respath = used_gpath + "-prokkaRes"
    os.makedirs(used_respath)
    shutil.copyfile(orig_gpath, used_gpath)
    # Copy tbl and gff files, which is as expected (tbl2lst and generate_gff must succeed)
    orig_fna =  os.path.join(orig_gpath + "-prokkaRes", "prokka_out_for_test.fna")
    fna_prokres = os.path.join(used_respath, "prokka_out_for_test.fna")
    shutil.copyfile(orig_fna, fna_prokres)
    orig_tbl = os.path.join(orig_gpath + "-prokkaRes",
                            "prokka_out_for_test.tbl")
    used_tbl = os.path.join(used_respath, "prokka_out_for_test.tbl")
    shutil.copyfile(orig_tbl, used_tbl)
    orig_gff = os.path.join(orig_gpath + "-prokkaRes",
                            "prokka_out_for_test.gff")
    used_gff = os.path.join(used_respath, "prokka_out_for_test.gff")
    shutil.copyfile(orig_gff, used_gff)
    orig_ffn = os.path.join(orig_gpath + "-prokkaRes",
                            "prokka_out_for_test.ffn")
    used_ffn = os.path.join(used_respath, "prokka_out_for_test.ffn")
    shutil.copyfile(orig_ffn, used_ffn)
    # Create faa_file with a wrong format
    orig_faa = os.path.join(TEST_ANNOTE, "prokka_out_for_test-wrongHeaderInt.faa")
    used_faa = os.path.join(used_respath, "prokka_out_for_test.faa")
    shutil.copyfile(orig_faa, used_faa)

    # Create output directories
    prok_path = GENEPATH
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gen_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff")
    os.makedirs(rep_dir)
    os.makedirs(gff_dir)
    os.makedirs(lst_dir)
    os.makedirs(gen_dir)
    os.makedirs(prot_dir)
    # Add empty res lst, gff and gen files, to check that it is removed at the end
    res_gff_file = os.path.join(gff_dir, "test.0417.00002.gff")
    open(res_gff_file, "w").close()
    assert len(os.listdir(gff_dir) ) == 1
    res_lst_file = os.path.join(lst_dir, "test.0417.00002.lst")
    open(res_lst_file, "w").close()
    assert len(os.listdir(lst_dir) ) == 1
    res_gen_file = os.path.join(gen_dir, "test.0417.00002.gen")
    open(res_gen_file, "w").close()
    assert len(os.listdir(gen_dir) ) == 1

    # Run formatting
    assert not prokkafunc.format_one_genome(used_gpath, name, prok_path, lst_dir, prot_dir,
                                            gen_dir, rep_dir, gff_dir)

    # Check that all files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(gen_dir) ) == 0
    assert("Unknown header format >JGIKIPIJ_d0008 in test/data/annotate/generated_by_unit-tests/"
           "original_name.fna-prokkaRes/prokka_out_for_test.faa. "
           "Gene ID is not a number.") in caplog.text
    assert ("Problems while generating .prt file for test.0417.00002") in caplog.text
