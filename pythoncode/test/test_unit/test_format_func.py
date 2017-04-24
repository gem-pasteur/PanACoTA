#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for format_functions.py
"""

import pytest
import os
import pipelinepackage.format_functions as ffunc
from io import StringIO


def test_write_gene():
    """
    Test that lstinfo line is written as expected when writting info for
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
    exp_file = os.path.join("test", "data", "exp_files", "res_test_write_geneCDS.lst")
    assert res == crispr_num
    with open(exp_file, "r") as expf, open(lstfile, "r") as lstf:
        for line_exp, line_out in zip(expf, lstf):
            assert line_exp == line_out
    os.remove(lstfile)


def test_write_CRISPR():
    """
    Test that lstinfo line is written as expected when writting info for CRISPR,
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
    exp_file = os.path.join("test", "data", "exp_files", "res_test_write_geneCRISPR.lst")
    assert res == 2
    with open(exp_file, "r") as expf, open(lstfile, "r") as lstf:
        for line_exp, line_out in zip(expf, lstf):
            assert line_exp == line_out
    os.remove(lstfile)


def test_tbl_to_lst():
    """
    Check that generated lstinfo file is as expected.
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
    tblfile = os.path.join("test", "data", "test_files", "prokka_out_for_test.tbl")
    lstfile = os.path.join("test", "data", "test_tbl2lst.lst")
    exp_lst = os.path.join("test", "data", "exp_files", "res_tbl2lst.lst")
    ffunc.tbl2lst(tblfile, lstfile)
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
           "| similar to AA sequence:UniProtKB:P37665")
    assert res == exp
    outfile.close()


def test_write_header_geneNoName():
    """
    From a given line of lstinfo file, giving info for a gene with many unknowk parts (gene
    name, product, EC number and more information are NAs), check that the header line of the
    protein and gene files are generated as expected.
    """
    outfile = StringIO()
    lstline = ("4632\t5000\tC\tCDS\ttest.0417.00002.b0002_00011\tNA\t| hypothetical protein "
               "| NA | NA")
    ffunc.write_header(lstline, outfile)
    res = outfile.getvalue()
    exp = (">test.0417.00002.b0002_00011 369 NA | hypothetical protein | NA | NA")
    assert res == exp
    outfile.close()


def test_write_header_geneNoName():
    """
    From a given line of lstinfo file, giving info for a CRISPR check that the header
    line of the protein and gene files are generated as expected.
    """
    outfile = StringIO()
    lstline = ("296902\t2968265\tC\tCRISPR\ttest.0417.00002.b0003_CRISPR1\tcrispr\t| "
               "crispr-array | NA | NA")
    ffunc.write_header(lstline, outfile)
    res = outfile.getvalue()
    exp = (">test.0417.00002.b0003_CRISPR1 2671364 crispr | crispr-array | NA | NA")
    assert res == exp
    outfile.close()


