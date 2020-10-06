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


def test_contig_name():
    """
    test that when we give a genome name and a contig number, it returns the expected fasta header
    for Replicon files (no gene number)
    """
    genome = "ESCO.1218.00005"
    cont_num = 30
    head_line = ffunc.get_contig_name(genome, cont_num)
    assert head_line == ">ESCO.1218.00005.0030"


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
    tbl_res = os.path.join(prok_path, " toto.tbl")
    gff_res = os.path.join(prok_path, "toto.gff")
    ffn_res = os.path.join(prok_path, "toto.ffn")
    faa_res = os.path.join(prok_path, "toto.faa")
    for file in [tbl_res, gff_res, ffn_res, faa_res]:
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


# def test_format_all_prokka():
#     """
#     Test that when giving a list of genomes, for which prokka ran without problem,
#     they are formatted, with all expected files created.
#     """
#     # genomes = {genome: [name, gpath, to_annot, size, nbcont, l90]}
#     initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
#     initpaths = [os.path.join(ANNOTEDIR, "genomes", name) for name in initnames]
#     gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-changeName.fna-short-contig.fna"]
#     onames = ["test_runprokka_H299", "test.0417.00002"]
#     gpaths = [os.path.join(ANNOTEDIR, "genomes", name) for name in gnames]
#     for f1, f2 in zip(initpaths, gpaths):
#         shutil.copyfile(f1, f2)
#     genomes = {gnames[0]: [onames[0], gpaths[0], gpaths[0], 12656, 3, 1],
#                gnames[1]: [onames[1], gpaths[1], gpaths[1], 456464645, 5, 1]
#                }
#     prok_path = os.path.join(ANNOTEDIR, "exp_files")
#     res_path = GENEPATH
#     skipped_format = ffunc.format_genomes(genomes, res_path,
#                                                    prok_path, False, threads=4)
#     assert skipped_format == []
#     lstfiles = [os.path.join(res_path, "LSTINFO", name + ".lst") for name in onames]
#     prtfiles = [os.path.join(res_path, "Proteins", name + ".prt") for name in onames]
#     genfiles = [os.path.join(res_path, "Genes", name + ".gen") for name in onames]
#     repfiles = [os.path.join(res_path, "Replicons", name + ".fna") for name in onames]
#     gfffiles = [os.path.join(res_path, "gff3", name + ".gff") for name in onames]
#     for f in lstfiles + prtfiles + genfiles + repfiles + gfffiles:
#         assert os.path.isfile(f)
#     shutil.rmtree(os.path.join(res_path, "LSTINFO"))
#     shutil.rmtree(os.path.join(res_path, "Proteins"))
#     shutil.rmtree(os.path.join(res_path, "Genes"))
#     shutil.rmtree(os.path.join(res_path, "Replicons"))
#     shutil.rmtree(os.path.join(res_path, "gff3"))


# def test_format_all_result_false():
#     """
#     Test that when giving a list of 2 genomes, 1 for which prokka ran without problem,
#     1 for which prokka had problems (given with False in results),
#     the correct genome is formatted, with all
#     expected files created, and the genome with problems is not formatted.
#     """
#     # genomes = {genome: [name, gpath, size, nbcont, l90]}
#     initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
#     initpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in initnames]
#     gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-changeName.fna-short-contig.fna"]
#     onames = ["test_runprokka_H299", "test.0417.00002"]
#     gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
#     for f1, f2 in zip(initpaths, gpaths):
#         shutil.copyfile(f1, f2)
#     genomes = {gnames[0]: [onames[0], gpaths[0], 12656, 3, 1],
#                gnames[1]: [onames[1], gpaths[1], 456464645, 5, 1]
#                }
#     prok_path = os.path.join("test", "data", "annotate", "exp_files")
#     res_path = os.path.join("test", "data", "annotate")
#     results = {gnames[0]: True, gnames[1]: False}
#     skipped, skipped_format = ffunc.format_genomes(genomes, results, res_path, prok_path)
#     assert skipped == ["B2_A3_5.fasta-changeName.fna-short-contig.fna"]
#     assert skipped_format == []
#     lstfiles = os.path.join(res_path, "LSTINFO")
#     prtfiles = os.path.join(res_path, "Proteins")
#     genfiles = os.path.join(res_path, "Genes")
#     repfiles = os.path.join(res_path, "Replicons")
#     gfffiles = os.path.join(res_path, "gff3")
#     assert os.path.isfile(os.path.join(lstfiles, onames[0] + ".lst"))
#     assert not os.path.isfile(os.path.join(lstfiles, onames[1] + ".lst"))
#     assert os.path.isfile(os.path.join(prtfiles, onames[0] + ".prt"))
#     assert not os.path.isfile(os.path.join(prtfiles, onames[1] + ".prt"))
#     assert os.path.isfile(os.path.join(genfiles, onames[0] + ".gen"))
#     assert not os.path.isfile(os.path.join(genfiles, onames[1] + ".gen"))
#     assert os.path.isfile(os.path.join(repfiles, onames[0] + ".fna"))
#     assert not os.path.isfile(os.path.join(repfiles, onames[1] + ".fna"))
#     assert os.path.isfile(os.path.join(gfffiles, onames[0] + ".gff"))
#     assert not os.path.isfile(os.path.join(gfffiles, onames[1] + ".gff"))
#     shutil.rmtree(os.path.join(res_path, "LSTINFO"))
#     shutil.rmtree(os.path.join(res_path, "Proteins"))
#     shutil.rmtree(os.path.join(res_path, "Genes"))
#     shutil.rmtree(os.path.join(res_path, "Replicons"))
#     shutil.rmtree(os.path.join(res_path, "gff3"))
#     for f in gpaths:
#         os.remove(f)


# def test_format_all_not_result():
#     """
#     Test that when giving a list of 2 genomes, but only 1 is in the results list (and prokka ran
#     without problems for it), the correct genome is formatted, with all
#     expected files created, and the other genome is not formatted, and does not appear in
#     skipped list (as it was removed from the study before annotation step, probably by QC).
#     """
#     # genomes = {genome: [name, gpath, size, nbcont, l90]}
#     initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
#     initpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in initnames]
#     gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-changeName.fna-short-contig.fna"]
#     onames = ["test_runprokka_H299", "test.0417.00002"]
#     gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
#     for f1, f2 in zip(initpaths, gpaths):
#         shutil.copyfile(f1, f2)
#     genomes = {gnames[0]: [onames[0], gpaths[0], 12656, 3, 1],
#                gnames[1]: [onames[1], gpaths[1], 456464645, 5, 1]
#                }
#     prok_path = os.path.join("test", "data", "annotate", "exp_files")
#     res_path = os.path.join("test", "data", "annotate")
#     results = {gnames[0]: True}
#     skipped, skipped_format = ffunc.format_genomes(genomes, results, res_path, prok_path)
#     assert skipped == []
#     assert skipped_format == []
#     lstfiles = os.path.join(res_path, "LSTINFO")
#     prtfiles = os.path.join(res_path, "Proteins")
#     genfiles = os.path.join(res_path, "Genes")
#     repfiles = os.path.join(res_path, "Replicons")
#     gfffiles = os.path.join(res_path, "gff3")
#     assert os.path.isfile(os.path.join(lstfiles, onames[0] + ".lst"))
#     assert not os.path.isfile(os.path.join(lstfiles, onames[1] + ".lst"))
#     assert os.path.isfile(os.path.join(prtfiles, onames[0] + ".prt"))
#     assert not os.path.isfile(os.path.join(prtfiles, onames[1] + ".prt"))
#     assert os.path.isfile(os.path.join(genfiles, onames[0] + ".gen"))
#     assert not os.path.isfile(os.path.join(genfiles, onames[1] + ".gen"))
#     assert os.path.isfile(os.path.join(repfiles, onames[0] + ".fna"))
#     assert not os.path.isfile(os.path.join(repfiles, onames[1] + ".fna"))
#     assert os.path.isfile(os.path.join(gfffiles, onames[0] + ".gff"))
#     assert not os.path.isfile(os.path.join(gfffiles, onames[1] + ".gff"))
#     shutil.rmtree(os.path.join(res_path, "LSTINFO"))
#     shutil.rmtree(os.path.join(res_path, "Proteins"))
#     shutil.rmtree(os.path.join(res_path, "Genes"))
#     shutil.rmtree(os.path.join(res_path, "Replicons"))
#     shutil.rmtree(os.path.join(res_path, "gff3"))
#     for f in gpaths:
#         os.remove(f)

#         # probleme avec .fna de onames[0] qui n'est pas créé...


# def test_format_all_error():
#     """
#     Test that when giving a list of 2 genomes, prokka ran without problem for both.
#     But a problem appears while formatting the 2nd one. So, the 2nd one is not formatted,
#     and appears in skipped_format. The first one is formated, and check that all
#     output files are created.
#     """
#     # genomes = {genome: [name, gpath, size, nbcont, l90]}
#     name = "test.0417.00002"
#     initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
#     initpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in initnames]
#     gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-problems.fna-short-contig.fna"]
#     onames = ["test_runprokka_H299", "test.0417.00002"]
#     gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
#     for f1, f2 in zip(initpaths, gpaths):
#         shutil.copyfile(f1, f2)
#     genomes = {gnames[0]: [onames[0], gpaths[0], 12656, 3, 1],
#                gnames[1]: [onames[1], gpaths[1], 456464645, 5, 1]
#                }
#     prok_path = os.path.join("test", "data", "annotate", "exp_files")
#     res_path = os.path.join("test", "data", "annotate")
#     tbl_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
#                             name + ".tbl")
#     tblout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
#                           name + ".tbl")
#     shutil.copyfile(tbl_init, tblout)
#     gff_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
#                             name + ".gff")
#     gffout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
#                           name + ".gff")
#     shutil.copyfile(gff_init, gffout)
#     results = {gnames[0]: True, gnames[1]: True}
#     skipped, skipped_format = ffunc.format_genomes(genomes, results, res_path, prok_path)
#     assert skipped == []
#     assert skipped_format == ["B2_A3_5.fasta-problems.fna-short-contig.fna"]
#     lstfiles = os.path.join(res_path, "LSTINFO")
#     prtfiles = os.path.join(res_path, "Proteins")
#     genfiles = os.path.join(res_path, "Genes")
#     repfiles = os.path.join(res_path, "Replicons")
#     gfffiles = os.path.join(res_path, "gff3")
#     assert os.path.isfile(os.path.join(lstfiles, onames[0] + ".lst"))
#     assert not os.path.isfile(os.path.join(lstfiles, onames[1] + ".lst"))
#     assert os.path.isfile(os.path.join(prtfiles, onames[0] + ".prt"))
#     assert not os.path.isfile(os.path.join(prtfiles, onames[1] + ".prt"))
#     assert os.path.isfile(os.path.join(genfiles, onames[0] + ".gen"))
#     assert not os.path.isfile(os.path.join(genfiles, onames[1] + ".gen"))
#     assert os.path.isfile(os.path.join(repfiles, onames[0] + ".fna"))
#     assert not os.path.isfile(os.path.join(repfiles, onames[1] + ".fna"))
#     assert os.path.isfile(os.path.join(gfffiles, onames[0] + ".gff"))
#     assert not os.path.isfile(os.path.join(gfffiles, onames[1] + ".gff"))
#     shutil.rmtree(os.path.join(res_path, "LSTINFO"))
#     shutil.rmtree(os.path.join(res_path, "Proteins"))
#     shutil.rmtree(os.path.join(res_path, "Genes"))
#     shutil.rmtree(os.path.join(res_path, "Replicons"))
#     shutil.rmtree(os.path.join(res_path, "gff3"))
#     os.remove(tblout)
#     os.remove(gffout)
#     for f in gpaths:
#         os.remove(f)
