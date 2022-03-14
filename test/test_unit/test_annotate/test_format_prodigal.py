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

import test.test_unit.utilities_for_tests as tutil
import PanACoTA.utils as utils
from PanACoTA.annotate_module import format_prodigal as prodigalfunc

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
    genfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.ffn")
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
    assert prodigalfunc.create_gene_lst(contigs, genfile, res_gen_file, res_lst_file, gpath,
                                        name)
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
    genfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.faa")
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
    assert not prodigalfunc.create_gene_lst(contigs, genfile, res_gen_file, res_lst_file,
                                            gpath, name)
    assert ("'my_contig' found in test/data/annotate/test_files/original_name.fna-prodigalRes/"
            "prodigal.outtest.ok.faa does not exist in original_genome_name") in caplog.text


def test_create_gff(caplog):
    """
    Check generated gff file. Must have all sequences in header (even replicons without gene),
    and then 1 line per gene
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.gff")
    contigs = {"JGIKIPgffgIJ": "test.0417.00002.0001",
               "toto": "test.0417.00002.0002",
               "other_header": "test.0417.00002.0003",
               "my_contig": "test.0417.00002.0004",
               "bis": "test.0417.00002.0005",
               "ter": "test.0417.00002.0006",
               "contname": "test.0417.00002.0007"
               }
    sizes = {"test.0417.00002.0001": 84,
             "test.0417.00002.0002": 103,
             "test.0417.00002.0003": 122,
             "test.0417.00002.0004": 35,
             "test.0417.00002.0005": 198,
             "test.0417.00002.0006": 128,
             "test.0417.00002.0007": 85,
            }
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
    gfffile = os.path.join(TEST_ANNOTE, "prodigal.outtest.wrong-start.gff")
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
    assert ("Files prodigal.outtest.wrong-start.ffn and "
            "prodigal.outtest.wrong-start.gff "
            "(in prodigal tmp_files: original_genome_name-prodigalRes) "
            "do not have the same start value for gene EPKOMDHM_00008 "
            "(78 in gff, 77 in ffn") in caplog.text


def test_create_gff_wrong_end(caplog):
    """
    Check that when trying to generate gff, if problem, exits with expected message
    Here, the start position in the gff file generated by prodigal is not the same as
    the start position in the lstinfo file (which was taken from the ffn
    file generated by prodigal)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "prodigal.outtest.wrong-end.gff")
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
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    gpath = "original_genome_name"
    assert not prodigalfunc.create_gff(gpath, gfffile, res_gff_file, exp_lst, contigs, sizes)
    assert ("Files prodigal.outtest.wrong-end.ffn and "
            "prodigal.outtest.wrong-end.gff "
            "(in prodigal tmp_files: original_genome_name-prodigalRes) "
            "do not have the same end value for gene EPKOMDHM_00009 "
            "(2347 in gff, 2346 in ffn") in caplog.text


def test_create_gff_wrong_type(caplog):
    """
    Check that when trying to generate gff, if problem, exits with expected message
    Here, the start position in the gff file generated by prodigal is not the same as
    the start position in the lstinfo file (which was taken from the ffn
    file generated by prodigal)
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    gfffile = os.path.join(TEST_ANNOTE, "prodigal.outtest.wrong-type.gff")
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
    res_gff_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    gpath = "original_genome_name"
    assert not prodigalfunc.create_gff(gpath, gfffile, res_gff_file, exp_lst, contigs, sizes)
    assert ("Files prodigal.outtest.wrong-type.ffn and "
            "prodigal.outtest.wrong-type.gff "
            "(in prodigal tmp_files: original_genome_name-prodigalRes) "
            "do not have the same type value for gene EPKOMDHM_00008 "
            "(tRNA in gff, CDS in ffn") in caplog.text


def test_create_prt(caplog):
    """
    Check that prt file is generated as expected
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.prt")
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    assert prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    exp_prt = os.path.join(EXP_ANNOTE, "res_create_prt_prodigal.faa")
    assert tutil.compare_order_content(exp_prt, res_prt_file)


def test_create_prt_wrong_lst(caplog):
    """
    Check that prt file is not generated if there is a problem in the lst file
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-wrongformat.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Problem in format of lstline (1279\t2346\tCDS\ttest.0417.00002.0002i_00005\tNA\t"
            "| NA | NA | fefer | NA") in caplog.text


def test_create_prt_short_lst(caplog):
    """
    Check that prt file is not generated when there are more proteins in faa than in lst
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-shortlst.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("No more protein in lst file. We cannot get information on this "
            "protein (>toto_00011 # 2419 # 3000 # 1 # a)! Check that you do not have "
            "more proteins than genes in prodigal results") in caplog.text


def test_create_prt_end_not_int_lst(caplog):
    """
    Check that prt file is not generated when an end value of a protein in lst file is not an int
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-notint.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Start and/or end of protein test.0417.00002.0001i_00002 position is not a number "
            "(start = 4416; end = a6068)") in caplog.text


def test_create_prt_not_divisible3_lst(caplog):
    """
    Check that prt file is not generated when there is a gene which length is not
    divisible by 3
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                           "prodigal.outtest.ok.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-not-divisible3.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Gene test.0417.00002.0001b_00003 has a number of nucleotides (3001) "
            "that is not divisible by 3.") in caplog.text


def test_create_prt_not_moreprots_lst(caplog):
    """
    Check that prt file is not generated when there are more proteins in lst than in faa
    """
    caplog.set_level(logging.DEBUG)
    protfile = os.path.join(TEST_ANNOTE, "original_name.fna-prodigalRes",
                        "prodigal.outtest.ok.faa")
    res_prt_file = os.path.join(GENEPATH, "prodigal_res.gff")
    exp_lst = os.path.join(TEST_ANNOTE, "test_create_prt_prodigal-more-proteins.lst")
    assert not prodigalfunc.create_prt(protfile, res_prt_file, exp_lst)
    assert ("Protein test.0417.00002.0007b_00016 is in .lst file but its sequence is not "
            "in the protein file generated by prodigal.") in caplog.text


def test_format_1genome(caplog):
    """
    Test that when prodigal results are ok, all files are generated as expected.
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    name = "prodigal.outtest.ok"
    gpath =  os.path.join(TEST_ANNOTE, "original_name.fna") # path to original genome, given to prodigal for annotation
    prod_path = TEST_ANNOTE
    # Generate result folders
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

    assert prodigalfunc.format_one_genome(gpath, name, prod_path, lst_dir, prot_dir, gene_dir,
                                          rep_dir, gff_dir)
    assert os.path.isfile(os.path.join(prot_dir, "prodigal.outtest.ok.prt"))
    assert os.path.isfile(os.path.join(lst_dir, "prodigal.outtest.ok.lst"))
    assert os.path.isfile(os.path.join(rep_dir, "prodigal.outtest.ok.fna"))
    assert os.path.isfile(os.path.join(gene_dir, "prodigal.outtest.ok.gen"))
    assert os.path.isfile(os.path.join(gff_dir, "prodigal.outtest.ok.gff"))


def test_format_1genome_emptygpath(caplog):
    """
    Test on formatting prodigal results, when ffn file is empty -> error message,
    and no file generated
    """
    caplog.set_level(logging.DEBUG)
    name = "prodigal.outtest.ok"
    # Create empty file, that we give to prodigal for formatting step
    gpath =  os.path.join(GENEPATH, "original_name-empty.fna")
    open(gpath, "w").close()
    prod_path = GENEPATH
    # Create prodigal result files (empty, then won't be read)
    prodigal_dir = gpath + "-prodigalRes"
    os.makedirs(prodigal_dir)
    prodigal_faa = os.path.join(gpath + "-prodigalRes", "notread.faa")
    prodigal_ffn = os.path.join(gpath + "-prodigalRes", "notread.ffn")
    prodigal_gff = os.path.join(gpath + "-prodigalRes", "notread.gff")
    for file in [prodigal_faa, prodigal_gff, prodigal_ffn]:
        open(file, "w").close()
    # Generate result folders
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
    res_gff_file = os.path.join(gff_dir, "prodigal.outtest.ok.gff")
    open(res_gff_file, "w").close()
    assert len(os.listdir(gff_dir) ) == 1
    res_lst_file = os.path.join(lst_dir, "prodigal.outtest.ok.lst")
    open(res_lst_file, "w").close()
    assert len(os.listdir(lst_dir) ) == 1
    res_gen_file = os.path.join(gen_dir, "prodigal.outtest.ok.gen")
    open(res_gen_file, "w").close()
    assert len(os.listdir(gen_dir) ) == 1

    assert not prodigalfunc.format_one_genome(gpath, name, prod_path, lst_dir, prot_dir, gen_dir,
                                              rep_dir, gff_dir)
    # Check that all files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(gen_dir) ) == 0
    assert ("Your genome test/data/annotate/generated_by_unit-tests/original_name-empty.fna "
            "does not contain any sequence, or is not in fasta format.") in caplog.text
    assert ("Problems while generating Replicon file for prodigal.outtest.ok") in caplog.text


def test_format_1genome_wrongffn(caplog):
    """
    Test on formatting prodigal results, when the ffn file generated by prodigal does
    not have the same contig name as in original fna file
    -> error message,
    -> no file generated
    """
    caplog.set_level(logging.DEBUG)
    name = "prodigal.outtest.ok"
    # path to original genome, given to prodigal for annotation
    orig_gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    orig_prodpath = orig_gpath + "-prodigalRes"
    used_gpath = os.path.join(GENEPATH, "original_name.fna")
    used_respath = used_gpath + "-prodigalRes"
    # Add original genome, and prodigal results to result folder
    shutil.copyfile(orig_gpath, used_gpath)
    shutil.copytree(orig_prodpath, used_respath)
    # In GENEPATH folder, create the original genome given to prodigal
    # (copy from test_file)
    # Create gen_file with a header not existing
    with open(os.path.join(used_respath, "prodigal.outtest.ok.ffn"), "w") as ori:
        ori.write(">wrongheader # 1 # 2 # 1 # toto")
    prod_path = GENEPATH
    # Generate result folders
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gene_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff")
    os.makedirs(rep_dir)
    os.makedirs(gene_dir)
    os.makedirs(lst_dir)
    os.makedirs(gff_dir)
    os.makedirs(prot_dir)
    # Add empty res gff file, to check that it is removed at the end
    res_gff_file = os.path.join(gff_dir, "prodigal.outtest.ok.gff")
    open(res_gff_file, "w").close()
    assert len(os.listdir(gff_dir) ) == 1
    # Run formatting
    assert not prodigalfunc.format_one_genome(used_gpath, name, prod_path, lst_dir, prot_dir,
                                              gene_dir, rep_dir, gff_dir)
    # Check that all files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(gene_dir) ) == 0
    assert len(os.listdir(prot_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert ("'wrongheader' found in test/data/annotate/generated_by_unit-tests/"
            "original_name.fna-prodigalRes/prodigal.outtest.ok.ffn does not exist in "
            "test/data/annotate/generated_by_unit-tests/original_name.fna") in caplog.text
    assert ("Problems while generating .gen and .lst files for prodigal.outtest.ok") in caplog.text


def test_format_1genome_wronggff(caplog):
    """
    Test on formatting prodigal results, when the gff file generated by prodigal
    does not have same info as ffn file generated by prodigal (= lst file generated by panacota)
    -> error message,
    -> no file generated
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    name = "prodigal.outtest.ok"

    # path to original genome, given to prodigal for annotation
    orig_gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    orig_prodpath = orig_gpath + "-prodigalRes"
    # In generated_by_tests folder, create the original genome given to prodigal
    # (copy from test_file)
    used_gpath = os.path.join(GENEPATH, "original_name.fna")
    used_respath = used_gpath + "-prodigalRes"
    # Add original genome, and prodigal results to result folder
    shutil.copyfile(orig_gpath, used_gpath)
    shutil.copytree(orig_prodpath, used_respath)
    # Copy gff file, but modify to get wrong start position
    orig_gff = os.path.join(orig_gpath + "-prodigalRes", "prodigal.outtest.ok.gff")
    used_gff = os.path.join(used_respath, "prodigal.outtest.ok.gff")
    with open(orig_gff, "r") as orig, open(used_gff, "w") as gffw:
        for line in orig:
            if line.startswith("#"):
                gffw.write(line)
            else:
                fields_g = line.split("\t")
                (contig_name, source, type_g, start_g, end_g,
                 score, strand_g, phase, attributes) = fields_g
                gffw.write("\t".join([contig_name, source, type_g, "12", end_g, score, strand_g,
                                      phase, attributes]))
                break
    prod_path = GENEPATH
    # Generate result folders
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gene_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff3")
    os.makedirs(rep_dir)
    os.makedirs(gene_dir)
    os.makedirs(prot_dir)
    os.makedirs(lst_dir)
    os.makedirs(gff_dir)
    assert not prodigalfunc.format_one_genome(used_gpath, name, prod_path, lst_dir, prot_dir,
                                              gene_dir, rep_dir, gff_dir)
    # Check that replicon and gff files were removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(gene_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(prot_dir) ) == 0
    assert ("Files prodigal.outtest.ok.ffn and prodigal.outtest.ok.gff (in prodigal tmp_files: "
            "test/data/annotate/generated_by_unit-tests/original_name.fna-prodigalRes) do "
            "not have the same start value for gene EPKOMDHM_00001 "
            "(12 in gff, 287 in ffn)") in caplog.text
    assert ("Problems while generating .gff (gff3 folder) file for "
            "prodigal.outtest.ok") in caplog.text



def test_format_1genome_wrongprt(caplog):
    """
    Test on formatting prodigal results, when the faa file generated by prodigal does
    not have as many proteins as the ffn file generated by prodigal (=lst file generated
    by panacota)
    -> error message,
    -> no file generated
    """
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger("test_prodigal")
    name = "prodigal.outtest.ok"

    # path to original genome, given to prodigal for annotation
    orig_gpath =  os.path.join(TEST_ANNOTE, "original_name.fna")
    orig_prodpath = orig_gpath + "-prodigalRes"
    # In generated_by_tests folder, create the original genome given to prodigal
    # (copy from test_file)
    used_gpath = os.path.join(GENEPATH, "original_name.fna")
    used_respath = used_gpath + "-prodigalRes"
    # Add original genome, and prodigal results to result folder
    shutil.copyfile(orig_gpath, used_gpath)
    shutil.copytree(orig_prodpath, used_respath)
    orig_faa = os.path.join(orig_gpath + "-prodigalRes", "prodigal.outtest.ok.faa")
    used_faa = os.path.join(used_respath, "prodigal.outtest.ok.faa")
    with open(orig_faa, "r") as faa, open(used_faa, "w") as faar:
        for i in range(7):
            faa.readline()
        for line in faa:
            faar.write(line)
    # Generate result folders
    prod_path = GENEPATH
    prot_dir = os.path.join(GENEPATH, "Proteins")
    lst_dir = os.path.join(GENEPATH, "LSTINFO")
    rep_dir = os.path.join(GENEPATH, "Replicons")
    gene_dir = os.path.join(GENEPATH, "Genes")
    gff_dir = os.path.join(GENEPATH, "gff3")
    os.makedirs(rep_dir)
    os.makedirs(gene_dir)
    os.makedirs(lst_dir)
    os.makedirs(gff_dir)
    os.makedirs(prot_dir)
    # Copy generated lstfile, but modify first line to have a difference between gff and lst starts

    assert not prodigalfunc.format_one_genome(used_gpath, name, prod_path, lst_dir, prot_dir,
                                              gene_dir, rep_dir, gff_dir)
    # Check that replicon file was removed
    assert len(os.listdir(rep_dir) ) == 0
    assert len(os.listdir(gene_dir) ) == 0
    assert len(os.listdir(lst_dir) ) == 0
    assert len(os.listdir(gff_dir) ) == 0
    assert len(os.listdir(prot_dir) ) == 0
    assert ("Protein prodigal.outtest.ok.0007b_00013 is in .lst file but its sequence is not in "
            "the protein file generated by prodigal.") in caplog.text
    assert ("Problems while generating .prt file (Proteins folder) for "
            "prodigal.outtest.ok") in caplog.text
