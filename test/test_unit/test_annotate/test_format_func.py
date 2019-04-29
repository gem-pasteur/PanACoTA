#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for format_functions.py
"""

import os
import logging
import shutil
from io import StringIO
import PanACoTA.annote_module.format_functions as ffunc
import PanACoTA.utils as utils


# Define variables and functions used by several tests
def my_logger():
    """
    logger given to function called by a subprocess
    """
    import multiprocessing
    m = multiprocessing.Manager()
    q = m.Queue()
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    return q, logging.getLogger('process')


# Start tests
def test_write_gene():
    """
    Test that lstinfo line is written as expected when writing info for
    a gene (CDS). Also check that crispr number is not changed
    """
    gtype = "CDS"
    locus_num = "5621221"
    gene_name = "abc"
    product = "new product"
    crispr_num = 1
    cont_loc = "i"
    genome = "ESCO.0216.00005"
    cont_num = 15
    ecnum = "454.12.5"
    inf2 = "more information... dfd | with | pipe|characters..."
    strand = "C"
    start = 154
    end = 656
    lstfile = "toto.lst"
    lstopenfile = open(lstfile, "w")
    res = ffunc.write_gene(gtype, locus_num, gene_name, product, crispr_num,
                           cont_loc, genome, cont_num, ecnum, inf2, strand,
                           start, end, lstopenfile)
    lstopenfile.close()
    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_test_write_geneCDS.lst")
    assert res == crispr_num
    with open(exp_file, "r") as expf, open(lstfile, "r") as lstf:
        for line_exp, line_out in zip(expf, lstf):
            assert line_exp == line_out
    os.remove(lstfile)


def test_write_crispr():
    """
    Test that lstinfo line is written as expected when writing info for CRISPR,
    and that crispr num increased by 1
    """
    gtype = "repeat_region"
    locus_num = "465"
    gene_name = "NA"
    product = "NA"
    crispr_num = 1
    cont_loc = "b"
    genome = "ESCO.0216.00005"
    cont_num = 15
    ecnum = "NA"
    inf2 = "more information... dfd | with | pipe|characters..."
    strand = "D"
    start = 154
    end = 656
    lstfile = "toto.lst"
    lstopenfile = open(lstfile, "w")
    res = ffunc.write_gene(gtype, locus_num, gene_name, product, crispr_num,
                           cont_loc, genome, cont_num, ecnum, inf2, strand,
                           start, end, lstopenfile)
    lstopenfile.close()
    exp_file = os.path.join("test", "data", "annotate", "exp_files",
                            "res_test_write_geneCRISPR.lst")
    assert res == 2
    with open(exp_file, "r") as expf, open(lstfile, "r") as lstf:
        for line_exp, line_out in zip(expf, lstf):
            assert line_exp == line_out
    os.remove(lstfile)


def test_tbl_to_lst():
    """
    Check that generated lstinfo file is as expected, when the genome name is the same as
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
    tblfile = os.path.join("test", "data", "annotate", "test_files", "original_name.fna-prokkaRes",
                           "prokka_out_for_test.tbl")
    lstfile = os.path.join("test", "data", "annotate", "test_tbl2lst.lst")
    exp_lst = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    name = "test.0417.00002"
    assert ffunc.tbl2lst(tblfile, lstfile, name)
    with open(exp_lst, "r") as expf, open(lstfile, "r") as lstf:
        for line_exp, line_out in zip(expf, lstf):
            assert line_exp == line_out
    os.remove(lstfile)


def test_tbl_to_lst_new_name():
    """
    Check that generated lstinfo file is as expected, when the genome name has changed between
    the one given to prokka, and the name given now.
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
    tblfile = os.path.join("test", "data", "annotate", "test_files", "original_name.fna-prokkaRes",
                           "prokka_out_for_test.tbl")
    lstfile = os.path.join("test", "data", "annotate", "test_tbl2lstNewName.lst")
    exp_lst = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst-newName.lst")
    name = "test.0417.00010"
    assert not ffunc.tbl2lst(tblfile, lstfile, name)
    with open(exp_lst, "r") as expf, open(lstfile, "r") as lstf:
        for line_exp, line_out in zip(expf, lstf):
            assert line_exp == line_out
    os.remove(lstfile)


def test_write_header_gene():
    """
    From a given line of lstinfo file, giving info for a gene (start, end, gene name,
    product, EC number, more information), check that the header line of the protein and
    gene files are generated as expected.
    """
    outfile = StringIO()
    lstline = ("4416\t6068\tD\tCDS\ttest.0417.00002.i0001_00005\tyiaD\t| "
               "putative lipoprotein YiaD | 6.3.2.- | similar to AA sequence:UniProtKB:P37665")
    ffunc.write_header(lstline, outfile)
    res = outfile.getvalue()
    exp = (">test.0417.00002.i0001_00005 1653 yiaD | putative lipoprotein YiaD | 6.3.2.- "
           "| similar to AA sequence:UniProtKB:P37665\n")
    assert res == exp
    outfile.close()


def test_write_header_gene_no_name():
    """
    From a given line of lstinfo file, giving info for a gene with many unknown parts (gene
    name, product, EC number and more information are NAs), check that the header line of the
    protein and gene files are generated as expected.
    """
    outfile = StringIO()
    lstline = ("4632\t5000\tC\tCDS\ttest.0417.00002.b0002_00011\tNA\t| hypothetical protein "
               "| NA | NA")
    ffunc.write_header(lstline, outfile)
    res = outfile.getvalue()
    exp = ">test.0417.00002.b0002_00011 369 NA | hypothetical protein | NA | NA\n"
    assert res == exp
    outfile.close()


def test_write_header_crispr():
    """
    From a given line of lstinfo file, giving info for a CRISPR check that the header
    line of the protein and gene files are generated as expected.
    """
    outfile = StringIO()
    lstline = ("296902\t2968265\tC\tCRISPR\ttest.0417.00002.b0003_CRISPR1\tcrispr\t| "
               "crispr-array | NA | NA")
    ffunc.write_header(lstline, outfile)
    res = outfile.getvalue()
    exp = ">test.0417.00002.b0003_CRISPR1 2671364 crispr | crispr-array | NA | NA\n"
    assert res == exp
    outfile.close()


def test_create_prt_wrong_header_sep():
    """
    Test that, when creating prt file from faa and lst, if a header of faa file is
    not in the right format (protein name and number are not separated by '_'),
    it writes an error, erases prt file, and returns False.
    """
    logger = my_logger()
    faaseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-wrongHeaderSep.faa")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    prtseq = os.path.join("test", "data", "annotate", "test_create_prt-wrongHeadSep.prt")
    assert not ffunc.create_prt(faaseq, lstfile, prtseq, logger[1])
    assert not os.path.isfile(prtseq)
    msg = ("Unknown header format >JGIKIPIJ00008 in test/data/annotate/test_files/"
           "prokka_out_for_test-wrongHeaderSep.faa. Error: invalid literal for int() "
           "with base 10: '>JGIKIPIJ00008'\nPrt file not created from "
           "test/data/annotate/test_files/prokka_out_for_test-wrongHeaderSep.faa.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_prt_wrong_header_int():
    """
    Test that, when creating prt file from faa and lst, if a header of faa file is
    not in the right format (protein name and number are separated by '_', but protein num
    contains a letter), it writes an error, erases prt file, and returns False.
    """
    logger = my_logger()
    faaseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-wrongHeaderInt.faa")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    prtseq = os.path.join("test", "data", "annotate", "test_create_prt-wrongHeadInt.prt")
    assert not ffunc.create_prt(faaseq, lstfile, prtseq, logger[1])
    assert not os.path.isfile(prtseq)
    msg = ("Unknown header format >JGIKIPIJ_d0008 in test/data/annotate/test_files/"
           "prokka_out_for_test-wrongHeaderInt.faa. Error: invalid literal for int() "
           "with base 10: 'd0008'\nPrt file not created from "
           "test/data/annotate/test_files/prokka_out_for_test-wrongHeaderInt.faa.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_prt_miss_lst():
    """
    Test that, when creating prt file from faa and lst, if a protein of faa file is not present in
    the lst file, it writes an error, removes the prt file, and returns False.
    """
    logger = my_logger()
    faaseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-supHeader.faa")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    prtseq = os.path.join("test", "data", "annotate", "test_create_prt-missLst.prt")
    assert not ffunc.create_prt(faaseq, lstfile, prtseq, logger[1])
    assert not os.path.isfile(prtseq)
    msg = ("Missing info for protein >sup-prot_00012 in "
           "test/data/annotate/exp_files/res_tbl2lst.lst. If it is "
           "actually present in the lst file, check that proteins are ordered by "
           "increasing number in both lst and faa files.\n"
           "Prt file not created from test/data/annotate/test_files/"
           "prokka_out_for_test-supHeader.faa.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_prt_wrong_order():
    """
    Test that, when creating prt file from faa and lst, if a protein of faa file is not in
    increasing protein number, so that it does not correspond to the protein in the lstinfo file,
    it writes an error, removes the prt file, and returns False.
    """
    logger = my_logger()
    faaseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-wrongOrder.faa")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    prtseq = os.path.join("test", "data", "annotate", "test_create_prt-wrongOrder.prt")
    assert not ffunc.create_prt(faaseq, lstfile, prtseq, logger[1])
    assert not os.path.isfile(prtseq)
    msg = ("Missing info for protein >appears_after_13_00011 in "
           "test/data/annotate/exp_files/res_tbl2lst.lst. If it is "
           "actually present in the lst file, check that proteins are ordered by "
           "increasing number in both lst and faa files.\n"
           "Prt file not created from test/data/annotate/test_files/"
           "prokka_out_for_test-wrongOrder.faa.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_prt_ok():
    """
    Test that when everything is ok in both faa and lst files, the prt file is
    created as expected.
    """
    logger = my_logger()
    faaseq = os.path.join("test", "data", "annotate", "test_files", "original_name.fna-prokkaRes",
                          "prokka_out_for_test.faa")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    prtseq = os.path.join("test", "data", "annotate", "test_create_prt.prt")
    assert ffunc.create_prt(faaseq, lstfile, prtseq, logger[1])
    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_create_prt.faa")
    with open(exp_file, "r") as expf, open(prtseq, "r") as prtf:
        for line_exp, line_out in zip(expf, prtf):
            assert line_exp == line_out
    os.remove(prtseq)


def test_create_gen_sup_crispr():
    """
    Test that when there is a CRISPR in the ffn file, but not in lstinfo,
    it generates an error, because the CRISPR ID does not correspond to the gene ID in lstinfo.
    It should return False, write an error message, and remove the .gen file.
    Moreover, the CRISPR ID is not in the same format as a gene ID, so the error should
    be on the format.
    """
    logger = my_logger()
    ffnseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-supCRISPR.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    genseq = os.path.join("test", "data", "annotate", "test_create_gen_supCRISPR.gen")
    assert not ffunc.create_gen(ffnseq, lstfile, genseq, logger[1])
    assert not os.path.isfile(genseq)
    msg = ("Unknown header format >prokka_out_for_test in test/data/annotate/test_files/"
           "prokka_out_for_test-supCRISPR.ffn.\nGen file will not be created.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_gen_sup_gene():
    """
    Test that, when creating gen file from ffn and lst, if a gene of ffn file is not present in
    the lst file, it writes an error, removes the gen file, and returns False.
    """
    logger = my_logger()
    ffnseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-supGene.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    prtseq = os.path.join("test", "data", "annotate", "test_create_gen-supgene.prt")
    assert not ffunc.create_gen(ffnseq, lstfile, prtseq, logger[1])
    assert not os.path.isfile(prtseq)
    msg = ("Missing info for gene >sup_gene_00012 in test/data/annotate/exp_files/"
           "res_tbl2lst.lst. If it is actually present "
           "in the lst file, check that genes are ordered by increasing "
           "number in both lst and ffn files.\nGen file not created"
           " from test/data/annotate/test_files/prokka_out_for_test-supGene.ffn.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_gen_miss_crispr():
    """
    Test for situation where there are 2 CRISPRs in the lstinfo file, the first one
    is not in the ffn file, while the second one is in the ffn file.
    It should return an error message about CRISPR number not corresponding, as
    the first CRISPR found in ffn (CRISPR 1) corresponds to CRISPR2 in lstinfo file.
    Gene file should be removed, and the function should return False
    """
    logger = my_logger()
    ffnseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-missCRISPR.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    genseq = os.path.join("test", "data", "annotate", "test_create_gen_missCRISPR.gen")
    assert not ffunc.create_gen(ffnseq, lstfile, genseq, logger[1])
    assert not os.path.isfile(genseq)
    msg = ("Problem with CRISPR numbers in test/data/annotate/exp_files/res_tbl2lst.lst. "
           "CRISPR >prokka_out_for_test in ffn is CRISPR num 1, whereas it is annotated "
           "as CRISPR num 2 in lst file.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_gen_no_crispr_ffn():
    """
    Test that when the there is a CRISPR in the lst file, but not in the ffn file,
    everything goes well, as in some versions of prokka (1.12)n CRISPRs are not in ffn while
    they are specified in lst. Function should return True, and gene file created.
    """
    ffnseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-noCRISPRffn.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    genseq = os.path.join("test", "data", "annotate", "test_create_gen_noCRISPRffn.gen")
    assert ffunc.create_gen(ffnseq, lstfile, genseq, my_logger()[1])
    assert os.path.isfile(genseq)
    exp_file = os.path.join("test", "data", "annotate", "exp_files",
                            "res_create_gen_noCRISPRffn.gen")
    with open(exp_file, "r") as expf, open(genseq, "r") as prtf:
        for line_exp, line_out in zip(expf, prtf):
            assert line_exp == line_out
    os.remove(genseq)


def test_create_gen_wrong_header_sep():
    """
    Test that, when creating gen file from ffn and lst, if a header of ffn file is
    not in the right format (gene name and number are not separated by '_'),
    it writes an error, erases gen file, and returns False.
    """
    logger = my_logger()
    ffnseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-wrongFormat.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    genseq = os.path.join("test", "data", "annotate", "test_create_gen_wrongHeadSep.gen")
    assert not ffunc.create_gen(ffnseq, lstfile, genseq, logger[1])
    assert not os.path.isfile(genseq)
    msg = ("Unknown header format >JGIKIPIJ-00005 in test/data/annotate/test_files/"
           "prokka_out_for_test-wrongFormat.ffn.\n"
           "Gen file will not be created.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_gen_wrong_header_int():
    """
    Test that, when creating gen file from ffn and lst, if a header of ffn file is
    not in the right format (gene name and number are separated by '_', but gene num
    contains a letter), it writes an error, erases gen file, and returns False.
    """
    logger = my_logger()
    ffnseq = os.path.join("test", "data", "annotate", "test_files",
                          "prokka_out_for_test-wrongInt.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    genseq = os.path.join("test", "data", "annotate", "test_create_gen_wrongHeadInt.gen")
    assert not ffunc.create_gen(ffnseq, lstfile, genseq, logger[1])
    assert not os.path.isfile(genseq)
    msg = ("Unknown header format >JGIKIPIJ_a00005 in test/data/annotate/test_files/"
           "prokka_out_for_test-wrongInt.ffn.\n"
           "Gen file will not be created.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_gen_wrong_lst_int():
    """
    Test that, when creating gen file from ffn and lst, if a gene name in lst file is
    not in the right format (gene name and number are separated by '_', but gene num
    contains a letter), it writes an error, erases gen file, and returns False.
    Because the gene name in ffn won't be found in lst (a it contains an error in lst).
    """
    logger = my_logger()
    ffnseq = os.path.join("test", "data", "annotate", "test_files", "original_name.fna-prokkaRes",
                          "prokka_out_for_test.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst-wrongGeneName.lst")
    genseq = os.path.join("test", "data", "annotate", "test_create_gen_wrongLstHeadInt.gen")
    assert not ffunc.create_gen(ffnseq, lstfile, genseq, logger[1])
    assert not os.path.isfile(genseq)
    msg = ("Missing info for gene >JGIKIPIJ_00009 in test/data/annotate/exp_files/"
           "res_tbl2lst-wrongGeneName.lst. If it is actually present "
           "in the lst file, check that genes are ordered by increasing "
           "number in both lst and ffn files.\nGen file not created "
           "from test/data/annotate/test_files/original_name.fna-prokkaRes/"
           "prokka_out_for_test.ffn.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_create_gen_ok():
    """
    Test that when everything is ok in both ffn and lst files, the gen file is
    created as expected.
    """
    faaseq = os.path.join("test", "data", "annotate", "test_files", "original_name.fna-prokkaRes",
                          "prokka_out_for_test.ffn")
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "res_tbl2lst.lst")
    genseq = os.path.join("test", "data", "annotate", "test_create_gen.gen")
    assert ffunc.create_gen(faaseq, lstfile, genseq, my_logger()[1])
    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_create_gen.gen")
    with open(exp_file, "r") as expf, open(genseq, "r") as prtf:
        for line_exp, line_out in zip(expf, prtf):
            assert line_exp == line_out
    os.remove(genseq)


def test_handle_line_gff():
    """
    Check that given a line in lstinfo with no information (gene name, ecnumber, product,
    inference), it returns the good gff line
    """
    line_lst = ("201\t743\tC\tCDS\tESCO.1015.00001.b0001_00001\tNA\t"
                "| NA | NA | NA")
    line_gff = ("H561_S27_L001__1\tProdigal:2.6\tCDS\t201\t743\t.\t-\t0\t"
                "ID=ONKACNIE_00001;inference=ab initio "
                "prediction:Prodigal:2.6;locus_tag=ONKACNIE_00001;product=hypothetical protein")
    outgff = "test_handle_line_gff.gff"
    gfff = open(outgff, "w")
    ffunc.handle_line_gff(line_lst, line_gff, gfff)
    gfff.close()
    with open(outgff, "r") as gf:
        lines = gf.readlines()
        assert len(lines) == 1
        exp_line = ("ESCO.1015.00001.0001\tProdigal:2.6\tCDS\t201\t743\t.\t-\t.\t"
                    "ID=ESCO.1015.00001.b0001_00001;locus_tag=ESCO.1015.00001.b0001_00001\n")
        assert lines[0] == exp_line
    os.remove(outgff)


def test_handle_line_gff_ecnum():
    """
    Check that given a line in lstinfo with only ec_number information,
    it returns the good gff line
    """
    line_lst = ("201\t743\tC\tCDS\tESCO.1015.00001.b0001_00001\tNA\t"
                "| NA | 123.995.4546 | NA")
    line_gff = ("H561_S27_L001__1\tProdigal:2.6\tCDS\t201\t743\t.\t-\t0\t"
                "ID=ONKACNIE_00001;inference=ab initio "
                "prediction:Prodigal:2.6;locus_tag=ONKACNIE_00001;product=hypothetical protein")
    outgff = "test_handle_line_gff-ecnum.gff"
    gfff = open(outgff, "w")
    ffunc.handle_line_gff(line_lst, line_gff, gfff)
    gfff.close()
    with open(outgff, "r") as gf:
        lines = gf.readlines()
        assert len(lines) == 1
        exp_line = ("ESCO.1015.00001.0001\tProdigal:2.6\tCDS\t201\t743\t.\t-\t.\t"
                    "ID=ESCO.1015.00001.b0001_00001;eC_number=123.995.4546;"
                    "locus_tag=ESCO.1015.00001.b0001_00001\n")
        assert lines[0] == exp_line
    os.remove(outgff)


def test_handle_line_gff_gene():
    """
    Check that given a line in lstinfo with only gene name information,
    it returns the good gff line
    """
    line_lst = ("201\t743\tC\tCDS\tESCO.1015.00001.b0001_00001\tge4a\t"
                "| NA | NA | NA")
    line_gff = ("H561_S27_L001__1\tProdigal:2.6\tCDS\t201\t743\t.\t-\t0\t"
                "ID=ONKACNIE_00001;inference=ab initio "
                "prediction:Prodigal:2.6;locus_tag=ONKACNIE_00001;product=hypothetical protein")
    outgff = "test_handle_line_gff-gene.gff"
    gfff = open(outgff, "w")
    ffunc.handle_line_gff(line_lst, line_gff, gfff)
    gfff.close()
    with open(outgff, "r") as gf:
        lines = gf.readlines()
        assert len(lines) == 1
        exp_line = ("ESCO.1015.00001.0001\tProdigal:2.6\tCDS\t201\t743\t.\t-\t.\t"
                    "ID=ESCO.1015.00001.b0001_00001;Name=ge4a;gene=ge4a;"
                    "locus_tag=ESCO.1015.00001.b0001_00001\n")
        assert lines[0] == exp_line
    os.remove(outgff)


def test_handle_line_gff_inf():
    """
    Check that given a line in lstinfo with only inference information,
    it returns the good gff line
    """
    line_lst = ("201\t743\tC\tCDS\tESCO.1015.00001.b0001_00001\tNA\t"
                "| NA | NA | ab initio prediction:Prodigal:2.6")
    line_gff = ("H561_S27_L001__1\tProdigal:2.6\tCDS\t201\t743\t.\t-\t0\t"
                "ID=ONKACNIE_00001;inference=ab initio "
                "prediction:Prodigal:2.6;locus_tag=ONKACNIE_00001;product=hypothetical protein")
    outgff = "test_handle_line_gff-inf.gff"
    gfff = open(outgff, "w")
    ffunc.handle_line_gff(line_lst, line_gff, gfff)
    gfff.close()
    with open(outgff, "r") as gf:
        lines = gf.readlines()
        assert len(lines) == 1
        exp_line = ("ESCO.1015.00001.0001\tProdigal:2.6\tCDS\t201\t743\t.\t-\t.\t"
                    "ID=ESCO.1015.00001.b0001_00001;inference=ab initio "
                    "prediction:Prodigal:2.6;locus_tag=ESCO.1015.00001.b0001_00001\n")
        assert lines[0] == exp_line
    os.remove(outgff)


def test_handle_line_gff_prod():
    """
    Check that given a line in lstinfo with only product information,
    it returns the good gff line
    """
    line_lst = ("201\t743\tC\tCDS\tESCO.1015.00001.b0001_00001\tNA\t"
                "| hypothetical protein | NA | NA")
    line_gff = ("H561_S27_L001__1\tProdigal:2.6\tCDS\t201\t743\t.\t-\t0\t"
                "ID=ONKACNIE_00001;inference=ab initio "
                "prediction:Prodigal:2.6;locus_tag=ONKACNIE_00001;product=hypothetical protein")
    outgff = "test_handle_line_gff-prod.gff"
    gfff = open(outgff, "w")
    ffunc.handle_line_gff(line_lst, line_gff, gfff)
    gfff.close()
    with open(outgff, "r") as gf:
        lines = gf.readlines()
        assert len(lines) == 1
        exp_line = ("ESCO.1015.00001.0001\tProdigal:2.6\tCDS\t201\t743\t.\t-\t.\t"
                    "ID=ESCO.1015.00001.b0001_00001;locus_tag=ESCO.1015.00001.b0001_00001;"
                    "product=hypothetical protein\n")
        assert lines[0] == exp_line
    os.remove(outgff)


def test_generate_gff():
    """
    Test creating gff file.
    """
    prokgff = os.path.join("test", "data", "annotate", "test_files", "prokka_out_gff.gff")
    gffout = os.path.join("test", "data", "annotate", "test_creategff.gff")
    lstgenome = os.path.join("test", "data", "annotate", "test_files", "lstinfo_for_gff.lst")
    assert ffunc.generate_gff(prokgff, gffout, lstgenome, my_logger()[1])
    exp_gff = os.path.join("test", "data", "annotate", "exp_files", "res_create_gff.gff")
    with open(gffout, "r") as gffo, open(exp_gff, "r") as expf:
        for line_exp, line_out in zip(expf, gffo):
            assert line_exp == line_out
    os.remove(gffout)


def test_generate_gff_error():
    """
    Test creating gff file.
    """
    logger = my_logger()
    prokgff = os.path.join("test", "data", "annotate", "test_files", "prokka_out_gff-error.gff")
    lstgenome = os.path.join("test", "data", "annotate", "test_files", "lstinfo_for_gff.lst")
    gffout = os.path.join("test", "data", "annotate", "test_creategff.gff")
    assert not ffunc.generate_gff(prokgff, gffout, lstgenome, logger[1])
    os.remove(gffout)
    q = logger[0]
    assert q.qsize() == 1
    logfound = q.get()
    msg = ("Problem with your gff file. '##FASTA' is not a gff entry line, whereas it should "
           "correspond to '863\t1795\tD\tCDS\tESCO.1015.00001.b0003_00016\tNA\t"
           "| hypothetical protein | NA | NA'")
    assert msg in logfound.message
    assert logfound.levelname == "ERROR"


def test_handle_genome_nores():
    """
    Test that when we try to format a genome which is not in results,
    it returns a tuple with "no_res" and the genome name.
    """
    results = {"abcd.fasta": True}
    args = ("toto.fasta", "name", "genome/path", "prokka/path", "lst/dir", "prot/dir",
            "gene/dir", "rep/dir", "gff/dir", results, my_logger()[0])
    res = ffunc.handle_genome(args)
    assert res == ("no_res", "toto.fasta")


def test_handle_genome_badprok():
    """
    Test that when we try to format a genome which is in results, but with False,
    it returns a tuple with "bad_prokka" and the genome name.
    """
    results = {"abcd.fasta": True, "toto.fasta": False}
    args = ("toto.fasta", "name", "genome/path", "prokka/path", "lst/dir", "prot/dir",
            "gene/dir", "rep/dir", "gff/dir", results, my_logger()[0])
    res = ffunc.handle_genome(args)
    assert res == ("bad_prokka", "toto.fasta")


def test_handle_genome_formatok():
    """
    Test that when we try to format a genome which is in results, with True,
    it returns a tuple with "True" and the genome name.
    """
    gpath = os.path.join("test", "data", "annotate", "genomes",
                         "B2_A3_5.fasta-split5N.fna-short-contig.fna")
    name = "test.0417.00002"
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    lst_dir = os.path.join("test", "data", "annotate")
    prot_dir = lst_dir
    gene_dir = lst_dir
    rep_dir = lst_dir
    gff_dir = lst_dir
    results = {"B2_A3_5.fasta-split5N.fna-short-contig.fna": True, "toto.fasta": False}
    args = ("B2_A3_5.fasta-split5N.fna-short-contig.fna", name, gpath, prok_path,
            lst_dir, prot_dir,
            gene_dir, rep_dir, gff_dir, results, my_logger()[0])
    res = ffunc.handle_genome(args)
    assert res == (True, "B2_A3_5.fasta-split5N.fna-short-contig.fna")
    os.remove(os.path.join(lst_dir, name + ".prt"))
    os.remove(os.path.join(lst_dir, name + ".fna"))
    os.remove(os.path.join(lst_dir, name + ".gen"))
    os.remove(os.path.join(lst_dir, name + ".lst"))
    os.remove(os.path.join(lst_dir, name + ".gff"))


def test_handle_genome_formaterror():
    """
    Test that when we try to format a genome which is in results, but with False,
    it returns a tuple with "bad_prokka" and the genome name.
    """
    logger = my_logger()
    gpath = os.path.join("test", "data", "annotate", "genomes",
                         "B2_A3_5.fasta-problems.fna-short-contig.fna")
    name = "test.0417.00002"
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    tbl_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".tbl")
    tblout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".tbl")
    shutil.copyfile(tbl_init, tblout)
    gff_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".gff")
    gffout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".gff")
    shutil.copyfile(gff_init, gffout)
    lst_dir = os.path.join("test", "data", "annotate")
    prot_dir = lst_dir
    gene_dir = lst_dir
    rep_dir = lst_dir
    gff_dir = lst_dir
    results = {"B2_A3_5.fasta-problems.fna-short-contig.fna": True, "toto.fasta": False}
    args = ("B2_A3_5.fasta-problems.fna-short-contig.fna", name, gpath,
            prok_path, lst_dir, prot_dir, gene_dir, rep_dir, gff_dir, results, logger[0])
    res = ffunc.handle_genome(args)
    assert res == (False, "B2_A3_5.fasta-problems.fna-short-contig.fna")
    msg = ("Unknown header format >EPKOMDHM_i00002 hypothetical protein in "
           "test/data/annotate/exp_files/B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes/"
           "test.0417.00002.ffn.\n"
           "Gen file will not be created.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    # remove tblout which was copied for this test
    assert not os.path.isfile(os.path.join(lst_dir, name + ".prt"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gen"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".fna"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".lst"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gff"))
    os.remove(tblout)
    os.remove(gffout)


def test_format1genome():
    """
    Test that formatting a genome (making .prt, .gen, .fna, .lst) works, with a genome
    which did not change name between prokka run and format step.
    """
    gpath = os.path.join("test", "data", "annotate", "genomes",
                         "B2_A3_5.fasta-split5N.fna-short-contig.fna")
    name = "test.0417.00002"
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    lst_dir = os.path.join("test", "data", "annotate")
    prot_dir = lst_dir
    gene_dir = lst_dir
    rep_dir = lst_dir
    gff_dir = lst_dir
    assert ffunc.format_one_genome(gpath, name, prok_path, lst_dir, prot_dir,
                                   gene_dir, rep_dir, gff_dir, my_logger()[1])
    # Check that all files were created
    assert os.path.isfile(os.path.join(lst_dir, name + ".lst"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".fna"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".prt"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".gen"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".gff"))
    # Check the contents of the files
    explst = os.path.join(prok_path, "res_format-B2.lst")
    expprt = os.path.join(prok_path, "res_format-B2.prt")
    expgen = os.path.join(prok_path, "res_format-B2.gen")
    expgff = os.path.join(prok_path, "res_format-B2.gff")
    with open(explst, "r") as lstf, open(os.path.join(lst_dir, name + ".lst"), "r") as lsto:
        for line_exp, line_out in zip(lstf, lsto):
            assert line_exp == line_out
    with open(expprt, "r") as expf, open(os.path.join(lst_dir, name + ".prt"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    with open(expgen, "r") as expf, open(os.path.join(lst_dir, name + ".gen"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    with open(gpath, "r") as expf, open(os.path.join(lst_dir, name + ".fna"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    with open(expgff, "r") as expf, open(os.path.join(lst_dir, name + ".gff"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    os.remove(os.path.join(lst_dir, name + ".lst"))
    os.remove(os.path.join(lst_dir, name + ".prt"))
    os.remove(os.path.join(lst_dir, name + ".fna"))
    os.remove(os.path.join(lst_dir, name + ".gen"))
    os.remove(os.path.join(lst_dir, name + ".gff"))


def test_format1genome_change_head():
    """
    Test that formatting a genome (making .prt, .gen, .fna, .lst) works, with a genome
    which changed its name between prokka and format step.
    """
    ginit = os.path.join("test", "data", "annotate", "genomes", "B2_A3_5.fasta-changeName.fna")
    gpath = os.path.join("test", "data", "annotate", "genomes",
                         "B2_A3_5.fasta-changeName.fna-short-contig.fna")
    shutil.copyfile(ginit, gpath)
    name = "test.0417.00002"
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    lst_dir = os.path.join("test", "data", "annotate")
    prot_dir = lst_dir
    gene_dir = lst_dir
    rep_dir = lst_dir
    gff_dir = lst_dir
    assert ffunc.format_one_genome(gpath, name, prok_path, lst_dir,
                                   prot_dir, gene_dir, rep_dir, gff_dir, my_logger()[1])
    # Check that all files were created
    assert os.path.isfile(os.path.join(lst_dir, name + ".lst"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".fna"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".prt"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".gen"))
    assert os.path.isfile(os.path.join(lst_dir, name + ".gff"))
    # Check the contents of the files
    explst = os.path.join(prok_path, "res_format-B2.lst")
    expprt = os.path.join(prok_path, "res_format-B2.prt")
    expgen = os.path.join(prok_path, "res_format-B2.gen")
    expgff = os.path.join(prok_path, "res_format-B2.gff")
    expreplicons = os.path.join("test", "data", "annotate", "genomes",
                                "B2_A3_5.fasta-split5N.fna-short-contig.fna")
    with open(explst, "r") as lstf, open(os.path.join(lst_dir, name + ".lst"), "r") as lsto:
        for line_exp, line_out in zip(lstf, lsto):
            assert line_exp == line_out
    with open(expprt, "r") as expf, open(os.path.join(lst_dir, name + ".prt"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    with open(expgen, "r") as expf, open(os.path.join(lst_dir, name + ".gen"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    with open(expgff, "r") as expf, open(os.path.join(lst_dir, name + ".gff"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    with open(expreplicons, "r") as expf, open(os.path.join(lst_dir, name + ".fna"), "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    os.remove(os.path.join(lst_dir, name + ".lst"))
    os.remove(os.path.join(lst_dir, name + ".prt"))
    os.remove(os.path.join(lst_dir, name + ".fna"))
    os.remove(os.path.join(lst_dir, name + ".gen"))
    os.remove(os.path.join(lst_dir, name + ".gff"))
    os.remove(gpath)


def test_format1genome_problemgen():
    """
    Test that formatting a genome (making .prt, .gen, .fna, .lst) returns an error message
    and does not create any output file if there is a problem while converting the
    .ffn to .gen
    """
    logger = my_logger()
    gpath = os.path.join("test", "data", "annotate", "genomes",
                         "B2_A3_5.fasta-problems.fna-short-contig.fna")
    name = "test.0417.00002"
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    tbl_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".tbl")
    tblout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".tbl")
    shutil.copyfile(tbl_init, tblout)
    gff_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".gff")
    gffout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".gff")
    shutil.copyfile(gff_init, gffout)
    lst_dir = os.path.join("test", "data", "annotate")
    prot_dir = lst_dir
    gene_dir = lst_dir
    rep_dir = lst_dir
    gff_dir = lst_dir
    assert not ffunc.format_one_genome(gpath, name, prok_path, lst_dir, prot_dir,
                                       gene_dir, rep_dir, gff_dir, logger[1])
    # Check that all files were not created
    assert not os.path.isfile(os.path.join(lst_dir, name + ".lst"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".fna"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".prt"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gen"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gff"))
    msg = ("Unknown header format >EPKOMDHM_i00002 hypothetical protein in "
           "test/data/annotate/exp_files/B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes/"
           "test.0417.00002.ffn.\n"
           "Gen file will not be created.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    # remove tblout which was copied for this test
    os.remove(tblout)
    os.remove(gffout)


def test_format1genome_problemprt():
    """
    Test that formatting a genome (making .prt, .gen, .fna, .lst) works, with a genome
    which did not change name between prokka run and format step.
    """
    logger = my_logger()
    gpath = os.path.join("test", "data", "annotate", "genomes",
                         "B2_A3_5.fasta-problems.fna-short-contig.fna")
    name = "test.0417.00002"
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    # copy tbl without errors to error prokka dir
    tbl_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".tbl")
    tblout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".tbl")
    shutil.copyfile(tbl_init, tblout)
    # copy gff without errors to error prokka dir
    gff_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".gff")
    gffout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".gff")
    shutil.copyfile(gff_init, gffout)
    # copy ffn without error to error prokka dir
    ffn_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".ffn")
    ffn_ok = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".ffn")
    ffn_error = ffn_ok + "-error"
    # change name of ffn file with error to keep it for later (used for tests)
    shutil.copyfile(ffn_ok, ffn_error)
    # copy ffn without error to prokka res (erasing ffn with error)
    shutil.copyfile(ffn_init, ffn_ok)
    lst_dir = os.path.join("test", "data", "annotate")
    prot_dir = lst_dir
    gene_dir = lst_dir
    rep_dir = lst_dir
    gff_dir = lst_dir
    assert not ffunc.format_one_genome(gpath, name, prok_path, lst_dir, prot_dir,
                                       gene_dir, rep_dir, gff_dir, logger[1])
    # Check that all files were not created
    assert not os.path.isfile(os.path.join(lst_dir, name + ".lst"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".fna"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".prt"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gen"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gff"))
    msg = ("Unknown header format >EPKOMDHM00003 hypothetical protein in "
           "test/data/annotate/exp_files/B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes/"
           "test.0417.00002.faa. "
           "Error: invalid literal for int() with base 10: '>EPKOMDHM00003'\n"
           "Prt file not created from test/data/annotate/exp_files/"
           "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes/test.0417.00002.faa.")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    # remove files which were copied for this test (tblout). And rename ffn with errors
    # to its original name.
    os.rename(ffn_error, ffn_ok)
    os.remove(tblout)
    os.remove(gffout)


def test_format1genome_problemgff():
    """
    Test that formatting a genome (making .prt, .gen, .fna, .lst) returns an error message
    and does not create any output file if there is a problem while converting the
    .ffn to .gen
    """
    logger = my_logger()
    gpath = os.path.join("test", "data", "annotate", "genomes",
                         "B2_A3_5.fasta-problems.fna-short-contig.fna")
    name = "test.0417.00002"
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    tbl_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".tbl")
    tblout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".tbl")
    shutil.copyfile(tbl_init, tblout)
    gff_init = os.path.join("test", "data", "annotate", "test_files", "prokka_out_gff-error.gff")
    gffout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".gff")
    shutil.copyfile(gff_init, gffout)
    lst_dir = os.path.join("test", "data", "annotate")
    prot_dir = lst_dir
    gene_dir = lst_dir
    rep_dir = lst_dir
    gff_dir = lst_dir
    assert not ffunc.format_one_genome(gpath, name, prok_path, lst_dir, prot_dir,
                                       gene_dir, rep_dir, gff_dir, logger[1])
    # Check that all files were not created
    assert not os.path.isfile(os.path.join(lst_dir, name + ".lst"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".fna"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".prt"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gen"))
    assert not os.path.isfile(os.path.join(lst_dir, name + ".gff"))
    msg = ("Problem with your gff file. '##FASTA' is not a gff entry line, whereas it should "
           "correspond to '11249\t12328\tD\tCDS\ttest.0417.00002.i0002_00016\tNA\t"
           "| hypothetical protein | NA | NA'")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    # remove tblout which was copied for this test
    os.remove(tblout)
    os.remove(gffout)


def test_format_all():
    """
    Test that when giving a list of genomes, for which prokka ran without problem,
    they are formatted, with all expected files created.
    """
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
    initpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in initnames]
    gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-changeName.fna-short-contig.fna"]
    onames = ["test_runprokka_H299", "test.0417.00002"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    for f1, f2 in zip(initpaths, gpaths):
        shutil.copyfile(f1, f2)
    genomes = {gnames[0]: [onames[0], gpaths[0], 12656, 3, 1],
               gnames[1]: [onames[1], gpaths[1], 456464645, 5, 1]
               }
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    res_path = os.path.join("test", "data", "annotate")
    results = {gname: True for gname in gnames}
    skipped, skipped_format = ffunc.format_genomes(genomes, results, res_path,
                                                   prok_path, threads=4)
    assert skipped == []
    assert skipped_format == []
    lstfiles = [os.path.join(res_path, "LSTINFO", name + ".lst") for name in onames]
    prtfiles = [os.path.join(res_path, "Proteins", name + ".prt") for name in onames]
    genfiles = [os.path.join(res_path, "Genes", name + ".gen") for name in onames]
    repfiles = [os.path.join(res_path, "Replicons", name + ".fna") for name in onames]
    gfffiles = [os.path.join(res_path, "gff3", name + ".gff") for name in onames]
    for f in lstfiles + prtfiles + genfiles + repfiles + gfffiles:
        assert os.path.isfile(f)
    shutil.rmtree(os.path.join(res_path, "LSTINFO"))
    shutil.rmtree(os.path.join(res_path, "Proteins"))
    shutil.rmtree(os.path.join(res_path, "Genes"))
    shutil.rmtree(os.path.join(res_path, "Replicons"))
    shutil.rmtree(os.path.join(res_path, "gff3"))
    for f in gpaths:
        os.remove(f)


def test_format_all_result_false():
    """
    Test that when giving a list of 2 genomes, 1 for which prokka ran without problem,
    1 for which prokka had problems (given with False in results),
    the correct genome is formatted, with all
    expected files created, and the genome with problems is not formatted.
    """
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
    initpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in initnames]
    gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-changeName.fna-short-contig.fna"]
    onames = ["test_runprokka_H299", "test.0417.00002"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    for f1, f2 in zip(initpaths, gpaths):
        shutil.copyfile(f1, f2)
    genomes = {gnames[0]: [onames[0], gpaths[0], 12656, 3, 1],
               gnames[1]: [onames[1], gpaths[1], 456464645, 5, 1]
               }
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    res_path = os.path.join("test", "data", "annotate")
    results = {gnames[0]: True, gnames[1]: False}
    skipped, skipped_format = ffunc.format_genomes(genomes, results, res_path, prok_path)
    assert skipped == ["B2_A3_5.fasta-changeName.fna-short-contig.fna"]
    assert skipped_format == []
    lstfiles = os.path.join(res_path, "LSTINFO")
    prtfiles = os.path.join(res_path, "Proteins")
    genfiles = os.path.join(res_path, "Genes")
    repfiles = os.path.join(res_path, "Replicons")
    gfffiles = os.path.join(res_path, "gff3")
    assert os.path.isfile(os.path.join(lstfiles, onames[0] + ".lst"))
    assert not os.path.isfile(os.path.join(lstfiles, onames[1] + ".lst"))
    assert os.path.isfile(os.path.join(prtfiles, onames[0] + ".prt"))
    assert not os.path.isfile(os.path.join(prtfiles, onames[1] + ".prt"))
    assert os.path.isfile(os.path.join(genfiles, onames[0] + ".gen"))
    assert not os.path.isfile(os.path.join(genfiles, onames[1] + ".gen"))
    assert os.path.isfile(os.path.join(repfiles, onames[0] + ".fna"))
    assert not os.path.isfile(os.path.join(repfiles, onames[1] + ".fna"))
    assert os.path.isfile(os.path.join(gfffiles, onames[0] + ".gff"))
    assert not os.path.isfile(os.path.join(gfffiles, onames[1] + ".gff"))
    shutil.rmtree(os.path.join(res_path, "LSTINFO"))
    shutil.rmtree(os.path.join(res_path, "Proteins"))
    shutil.rmtree(os.path.join(res_path, "Genes"))
    shutil.rmtree(os.path.join(res_path, "Replicons"))
    shutil.rmtree(os.path.join(res_path, "gff3"))
    for f in gpaths:
        os.remove(f)


def test_format_all_not_result():
    """
    Test that when giving a list of 2 genomes, but only 1 is in the results list (and prokka ran
    without problems for it), the correct genome is formatted, with all
    expected files created, and the other genome is not formatted, and does not appear in
    skipped list (as it was removed from the study before annotation step, probably by QC).
    """
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
    initpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in initnames]
    gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-changeName.fna-short-contig.fna"]
    onames = ["test_runprokka_H299", "test.0417.00002"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    for f1, f2 in zip(initpaths, gpaths):
        shutil.copyfile(f1, f2)
    genomes = {gnames[0]: [onames[0], gpaths[0], 12656, 3, 1],
               gnames[1]: [onames[1], gpaths[1], 456464645, 5, 1]
               }
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    res_path = os.path.join("test", "data", "annotate")
    results = {gnames[0]: True}
    skipped, skipped_format = ffunc.format_genomes(genomes, results, res_path, prok_path)
    assert skipped == []
    assert skipped_format == []
    lstfiles = os.path.join(res_path, "LSTINFO")
    prtfiles = os.path.join(res_path, "Proteins")
    genfiles = os.path.join(res_path, "Genes")
    repfiles = os.path.join(res_path, "Replicons")
    gfffiles = os.path.join(res_path, "gff3")
    assert os.path.isfile(os.path.join(lstfiles, onames[0] + ".lst"))
    assert not os.path.isfile(os.path.join(lstfiles, onames[1] + ".lst"))
    assert os.path.isfile(os.path.join(prtfiles, onames[0] + ".prt"))
    assert not os.path.isfile(os.path.join(prtfiles, onames[1] + ".prt"))
    assert os.path.isfile(os.path.join(genfiles, onames[0] + ".gen"))
    assert not os.path.isfile(os.path.join(genfiles, onames[1] + ".gen"))
    assert os.path.isfile(os.path.join(repfiles, onames[0] + ".fna"))
    assert not os.path.isfile(os.path.join(repfiles, onames[1] + ".fna"))
    assert os.path.isfile(os.path.join(gfffiles, onames[0] + ".gff"))
    assert not os.path.isfile(os.path.join(gfffiles, onames[1] + ".gff"))
    shutil.rmtree(os.path.join(res_path, "LSTINFO"))
    shutil.rmtree(os.path.join(res_path, "Proteins"))
    shutil.rmtree(os.path.join(res_path, "Genes"))
    shutil.rmtree(os.path.join(res_path, "Replicons"))
    shutil.rmtree(os.path.join(res_path, "gff3"))
    for f in gpaths:
        os.remove(f)

        # probleme avec .fna de onames[0] qui n'est pas créé...


def test_format_all_error():
    """
    Test that when giving a list of 2 genomes, prokka ran without problem for both.
    But a problem appears while formatting the 2nd one. So, the 2nd one is not formatted,
    and appears in skipped_format. The first one is formated, and check that all
    output files are created.
    """
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    name = "test.0417.00002"
    initnames = ["H299_H561.fasta", "B2_A3_5.fasta-changeName.fna"]
    initpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in initnames]
    gnames = ["H299_H561.fasta-short-contig.fna", "B2_A3_5.fasta-problems.fna-short-contig.fna"]
    onames = ["test_runprokka_H299", "test.0417.00002"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    for f1, f2 in zip(initpaths, gpaths):
        shutil.copyfile(f1, f2)
    genomes = {gnames[0]: [onames[0], gpaths[0], 12656, 3, 1],
               gnames[1]: [onames[1], gpaths[1], 456464645, 5, 1]
               }
    prok_path = os.path.join("test", "data", "annotate", "exp_files")
    res_path = os.path.join("test", "data", "annotate")
    tbl_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".tbl")
    tblout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".tbl")
    shutil.copyfile(tbl_init, tblout)
    gff_init = os.path.join(prok_path, "B2_A3_5.fasta-split5N.fna-short-contig.fna-prokkaRes",
                            name + ".gff")
    gffout = os.path.join(prok_path, "B2_A3_5.fasta-problems.fna-short-contig.fna-prokkaRes",
                          name + ".gff")
    shutil.copyfile(gff_init, gffout)
    results = {gnames[0]: True, gnames[1]: True}
    skipped, skipped_format = ffunc.format_genomes(genomes, results, res_path, prok_path)
    assert skipped == []
    assert skipped_format == ["B2_A3_5.fasta-problems.fna-short-contig.fna"]
    lstfiles = os.path.join(res_path, "LSTINFO")
    prtfiles = os.path.join(res_path, "Proteins")
    genfiles = os.path.join(res_path, "Genes")
    repfiles = os.path.join(res_path, "Replicons")
    gfffiles = os.path.join(res_path, "gff3")
    assert os.path.isfile(os.path.join(lstfiles, onames[0] + ".lst"))
    assert not os.path.isfile(os.path.join(lstfiles, onames[1] + ".lst"))
    assert os.path.isfile(os.path.join(prtfiles, onames[0] + ".prt"))
    assert not os.path.isfile(os.path.join(prtfiles, onames[1] + ".prt"))
    assert os.path.isfile(os.path.join(genfiles, onames[0] + ".gen"))
    assert not os.path.isfile(os.path.join(genfiles, onames[1] + ".gen"))
    assert os.path.isfile(os.path.join(repfiles, onames[0] + ".fna"))
    assert not os.path.isfile(os.path.join(repfiles, onames[1] + ".fna"))
    assert os.path.isfile(os.path.join(gfffiles, onames[0] + ".gff"))
    assert not os.path.isfile(os.path.join(gfffiles, onames[1] + ".gff"))
    shutil.rmtree(os.path.join(res_path, "LSTINFO"))
    shutil.rmtree(os.path.join(res_path, "Proteins"))
    shutil.rmtree(os.path.join(res_path, "Genes"))
    shutil.rmtree(os.path.join(res_path, "Replicons"))
    shutil.rmtree(os.path.join(res_path, "gff3"))
    os.remove(tblout)
    os.remove(gffout)
    for f in gpaths:
        os.remove(f)
