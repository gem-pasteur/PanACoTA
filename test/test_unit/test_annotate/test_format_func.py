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

import PanACoTA.utils as utils
import PanACoTA.annotate_module.general_format_functions as ffunc
import test.test_unit.utilities_for_tests as tutil


ANNOTEDIR = os.path.join("test", "data", "annotate")
EXP_ANNOTE = os.path.join(ANNOTEDIR, "exp_files")
GENEPATH = os.path.join(ANNOTEDIR, "generated_by_unit-tests")


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
    # utils.init_logger(LOGFILE_BASE, 0, 'test_fastme', verbose=1)
    if os.path.isdir(GENEPATH):
        content = os.listdir(GENEPATH)
        for f in content:
            assert f.startswith(".fuse")
    else:
        os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    print("teardown")


# Define variables and functions used by several tests
def my_logger():
    """
    logger given to function called by a subprocess
    """
    import multiprocessing
    m = multiprocessing.Manager()
    q = m.Queue()
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger('test_gene-format')
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    return q, logging.getLogger('process')


# Start tests
def test_write_gene():
    """
    Test that lstinfo line is written as expected when writing info for
    a gene (CDS).
    """
    gtype = "CDS"
    locus_num = "5621221"
    gene_name = "abc"
    product = "new product"
    cont_loc = "i"
    genome = "ESCO.0216.00005"
    cont_num = 15
    ecnum = "454.12.5"
    inf2 = "more information... dfd | with | pipe|characters..."
    db_xref = "mydb|pipe"
    strand = "C"
    start = str(154)
    end = str(656)
    lstfile = os.path.join(GENEPATH, "toto.lst")
    lstopenfile = open(lstfile, "w")
    lst_line = ffunc.write_gene(gtype, locus_num, gene_name, product,
                                cont_loc, genome, cont_num, ecnum, inf2, db_xref, strand,
                                start, end, lstopenfile)
    lstopenfile.close()
    assert lst_line == ("154\t656\tC\tCDS\tESCO.0216.00005.0015i_5621221\tabc\t| new product "
                        "| 454.12.5 | more information... dfd _ with _ pipe_characters... | "
                        "mydb_pipe")
    exp_file = os.path.join(EXP_ANNOTE, "res_test_write_geneCDS.lst")
    assert tutil.compare_order_content(exp_file, lstfile)


def test_write_header_gene():
    """
    From a given line of lstinfo file, giving info for a gene (start, end, gene name,
    product, EC number, more information), check that the header line of the protein and
    gene files are generated as expected.
    """
    outfile = StringIO()
    lstline = ("4416\t6068\tD\tCDS\ttest.0417.00002.0001i_00005\tyiaD\t| "
               "putative lipoprotein YiaD | 6.3.2.- | similar to AA sequence:UniProtKB:P37665")
    ffunc.write_header(lstline, outfile)
    res = outfile.getvalue()
    outfile.close()
    exp = (">test.0417.00002.0001i_00005 1653 yiaD | putative lipoprotein YiaD | 6.3.2.- "
           "| similar to AA sequence:UniProtKB:P37665\n")
    assert res == exp


def test_write_header_gene_no_name():
    """
    From a given line of lstinfo file, giving info for a gene with many unknown parts (gene
    name, product, EC number and more information are NAs), check that the header line of the
    protein and gene files are generated as expected.
    """
    outfile = StringIO()
    lstline = ("4632\t5000\tC\tCDS\ttest.0417.00002.0002b_00011\tNA\t| hypothetical protein "
               "| NA | NA")
    ffunc.write_header(lstline, outfile)
    res = outfile.getvalue()
    exp = ">test.0417.00002.0002b_00011 369 NA | hypothetical protein | NA | NA\n"
    assert res == exp
    outfile.close()


def test_handle_genome_badprok():
    """
    Test that when we try to format a genome which was annotated by prokka, but original genome
    is empty -> cannot format genome -> returns False associated with genome name
    """
    # Create prokka output dir and pur expected files (empty, we want to generate an error)
    gpath = os.path.join(GENEPATH, "toto.fasta")
    open(gpath, "w").close()
    prok_path = gpath + "-prokkaRes"
    os.makedirs(prok_path)
    fna_res = os.path.join(prok_path, " toto.fna")
    tbl_res = os.path.join(prok_path, " toto.tbl")
    gff_res = os.path.join(prok_path, "toto.gff")
    ffn_res = os.path.join(prok_path, "toto.ffn")
    faa_res = os.path.join(prok_path, "toto.faa")
    for file in [fna_res, tbl_res, gff_res, ffn_res, faa_res]:
        open(file, "w").close()
    # Create output directory for .fna file
    rep_dir = os.path.join(GENEPATH, "Replicons")
    os.makedirs(rep_dir)
    # Get args for function
    args = ("toto.fasta", "name", gpath, GENEPATH, "lst/dir", "prot/dir",
            "gene/dir", rep_dir, "gff/dir", False, my_logger()[0])
    ok_format, genome = ffunc.handle_genome(args)
    assert ok_format == False
    assert genome == "toto.fasta"


def test_handle_genome_badprodigal():
    """
    Test that when we try to format a genome which was annotated by prokka, but original genome
    is empty -> cannot format genome -> returns False associated with genome name
    """
    # Create prokka output dir and pur expected files (empty, we want to generate an error)
    gpath = os.path.join(GENEPATH, "wrong.fasta")
    open(gpath, "w").close()
    prodi_path = gpath + "-prodigalRes"
    os.makedirs(prodi_path)
    gff_res = os.path.join(prodi_path, "toto.gff")
    ffn_res = os.path.join(prodi_path, "toto.ffn")
    faa_res = os.path.join(prodi_path, "toto.faa")
    for file in [gff_res, ffn_res, faa_res]:
        open(file, "w").close()
    # Create output directory for .fna file
    rep_dir = os.path.join(GENEPATH, "Replicons")
    os.makedirs(rep_dir)
    # Get args for function
    args = ("wrong.fasta", "name", gpath, GENEPATH, "lst/dir", "prot/dir",
            "gene/dir", rep_dir, "gff/dir", True, my_logger()[0])
    ok_format, genome = ffunc.handle_genome(args)
    assert ok_format == False
    assert genome == "wrong.fasta"


def test_handle_genome_formatok(caplog):
    """
    Test that when we try to format a genome which was annotated by prokka without any problem
    It returns True associated with the genome name
    """
    caplog.set_level(logging.DEBUG)
    name = "test.0417.00002"
    # path to original genome, given to prodigal for annotation
    gpath =  os.path.join(ANNOTEDIR, "test_files", "original_name.fna")
    prok_path = os.path.join(ANNOTEDIR, "test_files")
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
    # Get args for function
    args = ("original_name", name, gpath, prok_path, lst_dir, prot_dir,
            gene_dir, rep_dir, gff_dir, False, my_logger()[0])
    ok_format, genome = ffunc.handle_genome(args)
    assert ok_format == True
    assert genome == "original_name"
    # Check generated files
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


def test_handle_genome_formatok_prodigal(caplog):
    """
    Test that when we try to format a genome which was annotated by prodigal without any problem
    It returns True associated with the genome name
    """
    caplog.set_level(logging.DEBUG)
    name_orig = "prodigal.outtest.ok"
    name = "test.0417.00002"
    # path to original genome, given to prodigal for annotation
    gpath =  os.path.join(ANNOTEDIR, "test_files", "original_name.fna")
    prodi_path = os.path.join(ANNOTEDIR, "test_files")
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
    # Get args for function
    args = (name_orig, name, gpath, prodi_path, lst_dir, prot_dir,
            gene_dir, rep_dir, gff_dir, True, my_logger()[0])
    ok_format, genome = ffunc.handle_genome(args)
    assert ok_format == True
    assert genome == name_orig
    # Check generated files
    exp_rep = os.path.join(EXP_ANNOTE, "res_created_rep-prokka.fna")
    res_rep_file = os.path.join(rep_dir, "test.0417.00002.fna")
    assert tutil.compare_order_content(exp_rep, res_rep_file)
    # Proteins
    exp_prt = os.path.join(EXP_ANNOTE, "res_create_prt_prodigal.faa")
    res_prt_file = os.path.join(prot_dir, "test.0417.00002.prt")
    assert tutil.compare_order_content(exp_prt, res_prt_file)
    # Genes
    exp_gen = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.gen")
    res_gen_file = os.path.join(gene_dir, "test.0417.00002.gen")
    assert tutil.compare_order_content(exp_gen, res_gen_file)
    # LSTINFO
    exp_lst = os.path.join(EXP_ANNOTE, "res_create_gene_lst_prodigal.lst")
    res_lst_file = os.path.join(lst_dir, "test.0417.00002.lst")
    assert tutil.compare_order_content(exp_lst, res_lst_file)
    # gff
    exp_gff = os.path.join(EXP_ANNOTE, "res_create_gff_prodigal.gff")
    res_gff_file = os.path.join(gff_dir, "test.0417.00002.gff")
    assert tutil.compare_order_content(exp_gff, res_gff_file)


def test_format_all_prokka(caplog):
    """
    Test that when giving a list of genomes, for which prokka ran without problem,
    they are formatted, with all expected files created.
    """
    caplog.set_level(logging.DEBUG)
    # genomes = {genome: [name, gpath, to_annot, size, nbcont, l90]}
    # Get genome names we want to format (with their path)
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
    gpaths = [os.path.join(ANNOTEDIR, "genomes", name) for name in gnames]
    onames = ["test_runprokka_H299", "test.0417.00002"]
    genomes = {gnames[0]: [onames[0], gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: [onames[1], gpaths[1], gpaths[1], 456464645, 5, 1]
               }
    res_path = GENEPATH
    annotated_path = os.path.join(ANNOTEDIR, "exp_files")
    # Format both genomes
    skipped_format = ffunc.format_genomes(genomes, res_path, annotated_path, False, threads=2)
    assert skipped_format == []
    # Get all names of expected output files
    exp_dir = os.path.join(EXP_ANNOTE, "res_formatAll", "prokka")
    exp_folders = ["LSTINFO", "Proteins", "Genes", "Replicons", "gff3"]
    exp_extensions = [".lst", ".prt", ".gen", ".fna", ".gff"]
    # Check that output files are created, and contain what is expected
    for fol, ext in zip(exp_folders, exp_extensions):
        exp_files = [os.path.join(exp_dir, fol, name + ext) for name in onames]
        res_files = [os.path.join(res_path, fol, name + ext) for name in onames]
        for res, exp in zip(res_files, exp_files):
            assert os.path.isfile(res)
            assert tutil.compare_order_content(res, exp)
    # Check log
    assert "Formatting all genomes" in caplog.text


def test_format_all_prodigal(caplog):
    """
    Test that when giving a list of genomes, for which prokka ran without problem,
    they are formatted, with all expected files created.
    """
    caplog.set_level(logging.DEBUG)
    # genomes = {genome: [name, gpath, to_annot, size, nbcont, l90]}
    # Get genome names we want to format (with their path)
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
    gpaths = [os.path.join(ANNOTEDIR, "genomes", name) for name in gnames]
    onames = ["test_runprokka_H299", "test.0417.00002"]
    genomes = {gnames[0]: [onames[0], gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: [onames[1], gpaths[1], gpaths[1], 456464645, 5, 1]
               }
    res_path = GENEPATH
    annotated_path = os.path.join(ANNOTEDIR, "exp_files")
    # Format both genomes
    skipped_format = ffunc.format_genomes(genomes, res_path, annotated_path, True, threads=2)
    assert skipped_format == []
    # Get all names of expected output files
    exp_dir = os.path.join(EXP_ANNOTE, "res_formatAll", "prodigal")
    exp_folders = ["LSTINFO", "Proteins", "Genes", "Replicons", "gff3"]
    exp_extensions = [".lst", ".prt", ".gen", ".fna", ".gff"]
    # Check that output files are created, and contain what is expected
    for fol, ext in zip(exp_folders, exp_extensions):
        exp_files = [os.path.join(exp_dir, fol, name + ext) for name in onames]
        res_files = [os.path.join(res_path, fol, name + ext) for name in onames]
        for res, exp in zip(res_files, exp_files):
            assert os.path.isfile(res)
            assert tutil.compare_order_content(res, exp)
    # Check log
    assert "Formatting all genomes" in caplog.text


def test_format_allpb_prokka(caplog):
    """
    Test that when giving a list of genomes, 1 that is correctly formatted, and 1 has a pb,
    it returns the last one in skipped_format
    """
    caplog.set_level(logging.DEBUG)
    # Create empty original sequence files
    genomes = ["wrong.fasta", "error.fasta"]
    gpaths = [os.path.join(GENEPATH, name) for name in genomes]
    for file in gpaths:
        open(file, "w").close()
    # Add prokka (empty) result files to prokkaRes directory
    prok_paths = [gpath + "-prokkaRes" for gpath in gpaths]
    for prok_path in prok_paths:
        os.makedirs(prok_path)
        tbl_res = os.path.join(prok_path, "toto.tbl")
        gff_res = os.path.join(prok_path, "toto.gff")
        ffn_res = os.path.join(prok_path, "toto.ffn")
        faa_res = os.path.join(prok_path, "toto.faa")
        fna_res = os.path.join(prok_path, "toto.fna")
        for file in [fna_res, tbl_res, gff_res, ffn_res, faa_res]:
            open(file, "w").close()
    # Create output directory for .fna files
    rep_dir = os.path.join(GENEPATH, "Replicons")
    os.makedirs(rep_dir)
    # genomes = {genome: [name, gpath, to_annot, size, nbcont, l90]}
    genomes = {genomes[0]: ["test_wrong-fasta", gpaths[0], gpaths[0], 12656, 3, 1],
               genomes[1]: ["test_error-fasta", gpaths[1], gpaths[1], 456464645, 5, 1]
              }
    res_path = GENEPATH
    annotated_path = GENEPATH
    # Try to format both genomes
    skipped_format = ffunc.format_genomes(genomes, res_path, annotated_path, False, threads=1)
    assert skipped_format == ["wrong.fasta", "error.fasta"]
    # Get all names of expected output files
    exp_folders = ["LSTINFO", "Proteins", "Genes", "Replicons", "gff3"]
    for res_folder in [os.path.join(res_path, folder) for folder in exp_folders]:
        assert len(os.listdir(res_folder)) == 0
    # Check log
    assert "Formatting all genomes" in caplog.text
    assert ("Your genome test/data/annotate/generated_by_unit-tests/"
            "wrong.fasta-prokkaRes/toto.fna does not "
            "contain any sequence, or is not in fasta format.") in caplog.text
    assert ("Your genome test/data/annotate/generated_by_unit-tests/"
            "error.fasta-prokkaRes/toto.fna does not "
            "contain any sequence, or is not in fasta format.") in caplog.text
    assert "Problems while generating Replicon file for test_wrong-fasta" in caplog.text
    assert "Problems while generating Replicon file for test_error-fasta" in caplog.text


def test_format_1pb_prodigal(caplog):
    """
    Test that when giving a list of genomes, 1 that is correctly formatted, and 1 has a pb,
    it returns the last one in skipped_format
    """
    caplog.set_level(logging.DEBUG)
    # GENOME 2: Create empty original genome file
    genome1 = "wrong.fasta"
    gpath1 = os.path.join(GENEPATH, "wrong.fasta")
    open(gpath1, "w").close()
    # Add prodigal (empty) result files to prodigalRes directory
    prodi_path = gpath1 + "-prodigalRes"
    os.makedirs(prodi_path)
    gff_res = os.path.join(prodi_path, "toto.gff")
    ffn_res = os.path.join(prodi_path, "toto.ffn")
    faa_res = os.path.join(prodi_path, "toto.faa")
    for file in [gff_res, ffn_res, faa_res]:
        open(file, "w").close()
    # Create output directory for .fna file
    rep_dir = os.path.join(GENEPATH, "Replicons")
    os.makedirs(rep_dir)
    # GENOME 2
    genome2 = "H299_H561.fasta"
    gpath2 = os.path.join(ANNOTEDIR, "genomes", genome2)
    # Copy results of prodigal for this genome to output dir (GENEPATH)
    orig_res_files = os.path.join(EXP_ANNOTE, genome2 + '-prodigalRes')
    used_res_path = os.path.join(GENEPATH, genome2 + "-prodigalRes")
    shutil.copytree(orig_res_files, used_res_path)
    # genomes = {genome: [name, gpath, to_annot, size, nbcont, l90]}
    genomes = {genome1: ["test_genome1", gpath1, gpath1, 12656, 3, 1],
               genome2: ["test_runprokka_H299", gpath2, gpath2, 456464645, 5, 1]
               }
    res_path = GENEPATH
    annotated_path = GENEPATH
    # Format both genomes
    skipped_format = ffunc.format_genomes(genomes, res_path, annotated_path, True, threads=2)
    assert skipped_format == ["wrong.fasta"]
    # Get all names of expected output files
    exp_dir = os.path.join(EXP_ANNOTE, "res_formatAll", "prodigal")
    exp_folders = ["LSTINFO", "Proteins", "Genes", "Replicons", "gff3"]
    exp_extensions = [".lst", ".prt", ".gen", ".fna", ".gff"]
    # Check that output files are created, and contain what is expected
    for fol, ext in zip(exp_folders, exp_extensions):
        exp_file = os.path.join(exp_dir, fol, "test_runprokka_H299" + ext)
        res_file = os.path.join(res_path, fol, "test_runprokka_H299" + ext)
        assert os.path.isfile(res_file)
        assert tutil.compare_order_content(res_file, exp_file)
    # Check log
    assert "Formatting all genomes" in caplog.text
    assert ("Your genome test/data/annotate/generated_by_unit-tests/wrong.fasta does not "
            "contain any sequence, or is not in fasta format.") in caplog.text
    assert "Problems while generating Replicon file for test_genome1" in caplog.text
