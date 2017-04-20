#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pipelinepackage.genome_seq_functions as gfunc
import pytest
import os

def test_sort_genomes():
    """
    Test the function sorting genomes by L90 and nb contigs.
    genome = name, path, gsize, nbcont, L90]
    """
    genome1 = ["SAEN.1015.", "path/to/genome1", 10000, 11, 2]
    genome2 = ["SAEN.1015.", "path/to/genome2", 10000, 12, 2]
    genome3 = ["SAEN.1015.", "path/to/genome3", 10000, 12, 1]
    genome4 = ["ESCO.0216.", "path/to/genome4", 10000, 12, 1]

    genomes = {1:genome1, 2:genome2, 3:genome3, 4:genome4}
    sorted_g = sorted(genomes.items(), key=gfunc.sort_genomes)
    exp = [(4, genome4), (3, genome3), (1, genome1), (2, genome2)]
    assert sorted_g == exp


def test_save_contig_5N():
    """
    Test that the given contig is split at each stretch of at least 5 'N', and not at
    stretches of less than 5 'N'.
    Check that the contigs in the output file are named as expected, and that
    the contig sizes are well reported.
    """
    pat = "NNNNN+"  # at least 5 'N'
    cur_cont = ("AACCGTGTCTCTCGGAGCNNNNCCGTTCGGCTCNCGGTCNNNNNCCGTTATNNCGGTTCGCNNNCTGGTC"
                "GGCTTATNNNNNNNNNNNNCCTGGTATTCGGCGCTTCNC")
    cur_cont_name = ">ESCO.0216.00001_cont2"
    # one contig saved before running, check that it is not erased
    contig_sizes = {">ESCO.0216.00001_cont1": 1623}
    seq_file = os.path.join("test", "data", "test_save_contig5N.faa")
    resf = open(seq_file, "w")
    gfunc.save_contig(pat, cur_cont, cur_cont_name, contig_sizes, resf)
    resf.close()
    exp = {">ESCO.0216.00001_cont1": 1623, ">ESCO.0216.00001_cont2_0": 39,
           ">ESCO.0216.00001_cont2_1": 33, ">ESCO.0216.00001_cont2_2": 20}
    assert contig_sizes == exp

    exp_file = os.path.join("test", "data", "exp_files", "res_save_contig5N.faa")
    with open(exp_file, "r") as expf, open(seq_file, "r") as seqf:
        for line_exp, line_seq in zip(expf, seqf):
            assert line_exp == line_seq
    os.remove(seq_file)


def test_save_contig_ATCG():
    """
    Test that the given contig is split at each pattern 'ATCG' (just for test)
    Check that the contigs in the output file are named as expected, and that
    the contig sizes are well reported.
    """
    pat = "ATCG"  # split each time those 4 letters are found in the sequence
    cur_cont = ("AAATGGTCTCGATGATCGATCGAGGGATTCGGAATCGGGCTCTGAATTCGATCGGTAGCTCTCGGGA"
                "GCTCTAGGCTCGTACGCCGTGATCGCATCGGTTCGTATCGATCGATCGATCGGGGGG")
    cur_cont_name = ">ESCO.0216.00001_cont2"
    # one contig saved before running, check that it is not erased
    contig_sizes = {">ESCO.0216.00001_cont1": 1623}
    seq_file = os.path.join("test", "data", "test_save_contigATCG.faa")
    resf = open(seq_file, "w")
    gfunc.save_contig(pat, cur_cont, cur_cont_name, contig_sizes, resf)
    resf.close()
    exp = {">ESCO.0216.00001_cont1": 1623, ">ESCO.0216.00001_cont2_0": 14,
           ">ESCO.0216.00001_cont2_1": 11, ">ESCO.0216.00001_cont2_2": 13,
           ">ESCO.0216.00001_cont2_3": 34, ">ESCO.0216.00001_cont2_4": 1,
           ">ESCO.0216.00001_cont2_5": 6, ">ESCO.0216.00001_cont2_6": 5}
    assert contig_sizes == exp

    exp_file = os.path.join("test", "data", "exp_files", "res_save_contigATCG.faa")
    with open(exp_file, "r") as expf, open(seq_file, "r") as seqf:
        for line_exp, line_seq in zip(expf, seqf):
            assert line_exp == line_seq
    os.remove(seq_file)


def test_calc_l90_exact():
    """
    Calculate L90 according to the given genome size and contig sizes
    2 contigs get exactly 90% of the genome
    """
    cont_size = {1: 3, 2:800, 3:100, 4:90, 5:7}
    l90 = gfunc.calc_l90(cont_size)
    assert l90 == 2


def test_calc_l90_more():
    """
    Calculate L90 according to the given genome size and contig sizes
    3 contigs get exactly more than 90%, but 2 contigs get less -> l90 = 3
    """
    cont_size = {1: 3, 2: 800, 3: 90, 4: 90, 5: 17}
    l90 = gfunc.calc_l90(cont_size)
    assert l90 == 3

