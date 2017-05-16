#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pytest
import os
import genomeAPCAT.qc_annote_module.genome_seq_functions as gfunc
import test.test_unit.util_tests as util_tests


def test_sort_genomes():
    """
    Test the function sorting genomes by L90 and nb contigs.
    genome = name, path, gsize, nbcont, L90]
    """
    genome1 = ["SAEN.1116.", "path/to/genome1", 10000, 11, 2]
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


def test_rename_contigs():
    """
    From a given sequence, rename all its contigs with the given gembase name + a number,
    and save the output sequence to the given res_path.
    Check that the output file is as expected.
    """
    gpath = os.path.join("test", "data", "genomes", "H299_H561.fasta")
    gembase_name = "ESCO.0216.00005"
    res_path = os.path.join("test", "data")
    outfile = os.path.join(res_path, "H299_H561.fasta-gembase.fna")
    exp_file = os.path.join("test", "data", "exp_files", "res_H299_H561-ESCO00005.fna")
    gfunc.rename_genome_contigs(gembase_name, gpath, outfile)
    with open(exp_file, "r") as expf, open(outfile, "r") as of:
        for line_exp, line_seq in zip(expf, of):
            assert line_exp == line_seq
    os.remove(outfile)


def test_rename_genomes():
    """
    From a list of genomes ({genome: [name.date, path, gsize, nbcont, L90]}),
    order them by species, and by decreasing quality (L90, nb_cont), and rename them,
    as well as their contigs.
    """
    genomes_dir = os.path.join("test", "data", "genomes")
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
          "genome5.fasta", "genome6.fasta", "genome7.fasta"]

    genomes = {gs[0]: ["SAEN.1113", os.path.join(genomes_dir, gs[0]), 51, 4, 2],
               gs[1]: ["SAEN.1114", os.path.join(genomes_dir, gs[1]), 67, 3, 3],
               gs[2]: ["ESCO.0416", os.path.join(genomes_dir, gs[2]), 70, 4, 1],
               gs[3]: ["ESCO.0216", os.path.join(genomes_dir, gs[3]), 114, 5, 2],
               gs[4]: ["SAEN.1115", os.path.join(genomes_dir, gs[4]), 106, 3, 1],
               gs[5]: ["ESCO.0216", os.path.join(genomes_dir, gs[5]), 116, 4, 2],
               gs[6]: ["SAEN.1115", os.path.join(genomes_dir, gs[6]), 137, 3, 2]}
    res_path = os.path.join("test", "data")
    out_f = [os.path.join(res_path, gname + "-gembase.fna") for gname in gs]
    gfunc.rename_all_genomes(genomes, res_path)
    # SAEN genomes 1 and 2 have same characteristics. Their place will be chosen randomly,
    # so take into account both choices
    exp_genomes =  {gs[0]: ["SAEN.1113.00003", out_f[0], 51, 4, 2],
                    gs[1]: ["SAEN.1114.00004", out_f[1], 67, 3, 3],
                    gs[2]: ["ESCO.0416.00001", out_f[2], 70, 4, 1],
                    gs[3]: ["ESCO.0216.00003", out_f[3], 114, 5, 2],
                    gs[4]: ["SAEN.1115.00001", out_f[4], 106, 3, 1],
                    gs[5]: ["ESCO.0216.00002", out_f[5], 116, 4, 2],
                    gs[6]: ["SAEN.1115.00002", out_f[6], 137, 3, 2]}
    exp_f = [os.path.join("test", "data", "exp_files", "res_" + gname + "-gembase.fna")
             for gname in gs]
    assert genomes == exp_genomes
    for exp, out in zip(exp_f, out_f):
        with open(exp, "r") as expf, open(out, "r") as outf:
            for line_exp, line_out in zip(expf, outf):
                assert line_exp == line_out
        os.remove(out)


def test_analyse1genome_nocut():
    """
    Analyse the given genome: without cutting at stretches of N, calculate
    its genome size, nb contigs and L90, and add it to the genomes dict, as well as
    the path to the genome file.
    """
    dbpath = os.path.join("test", "data", "genomes")
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = gs[1]
    tmp_path = os.path.join("plop")
    cut = False
    pat = "NNNNN+"
    gfunc.analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes)
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114", os.path.join(dbpath, gs[1]), 67, 3, 3],
                   gs[2]: ["ESCO.0416"]}
    assert genomes == exp_genomes


# def test_analyse1genome_nocut_empty():
#     """
#     Analyse the given genome: without cutting at stretches of N. The genome is an empty
#     file, so it is not possible to calculate L90
#     """
#     dbpath = os.path.join("test", "data", "genomes")
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
#     open(os.path.join(dbpath, gs[3]), "w").close()
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"]}
#     genome = gs[1]
#     tmp_path = os.path.join("plop")
#     cut = False
#     pat = "NNNNN+"
#     gfunc.analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes)
#     exp_genomes = {gs[0]: ["SAEN.1113"],
#                    gs[1]: ["SAEN.1114", os.path.join(dbpath, gs[1]), 67, 3, 3],
#                    gs[2]: ["ESCO.0416"]}
#     assert genomes == exp_genomes


def test_analyse1genome_cut():
    """
    Analyse the given genome: cut at each stretch of 5 N, put it to a new file,
    and then calculate its genome size, nb contigs and L90. Add this information
    to the genomes dict, as well as the path to the genome file (cut).
    """
    dbpath = os.path.join("test", "data", "genomes")
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = gs[1]
    tmp_path = os.path.join("test", "data")
    cut = True
    pat = "NNNNN+"
    gfunc.analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes)
    out_f = os.path.join(tmp_path, gs[1] + "-split5N.fna")
    exp_f = os.path.join(tmp_path, "exp_files", "res_genome2.fasta-split5N.fna")
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114", out_f, 55, 5, 4],
                   gs[2]: ["ESCO.0416"]}
    assert genomes == exp_genomes
    with open(out_f, "r") as outf, open(exp_f, "r") as expf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    os.remove(out_f)


def test_analyseAllGenomes_nocut():
    """
    Analyze all given genomes: don't cut at stretches of N, but look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    dbpath = os.path.join("test", "data", "genomes")
    gpaths = [os.path.join(dbpath, gname) for gname in gs]
    tmp_path = os.path.join("test", "data")
    nbn = 0
    # Run analysis
    gfunc.analyse_all_genomes(genomes, dbpath, tmp_path, nbn)
    # construct expected results
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], 67, 3, 3],
                   gs[2]: ["ESCO.0416", gpaths[2], 70, 4, 1]}
    assert exp_genomes == genomes


def test_analyseAllGenomes_cut():
    """
    Analyze all given genomes: cut at each stretch of 5 'N', look at the output sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    dbpath = os.path.join("test", "data", "genomes")
    tmp_path = os.path.join("test", "data")
    gpaths = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
    nbn = 5
    # Run analysis
    gfunc.analyse_all_genomes(genomes, dbpath, tmp_path, nbn)
    # construct expected results
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], 55, 5, 4],
                   gs[2]: ["ESCO.0416", gpaths[2], 70, 4, 1]}
    out_f = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
    exp_f = [os.path.join(tmp_path, "exp_files", "res_" + gname + "-split5N.fna") for gname in gs]
    assert exp_genomes == genomes
    for out, exp in zip(out_f, exp_f):
        with open(out, "r") as outf, open(exp, "r") as expf:
            for line_exp, line_out in zip(expf, outf):
                assert line_exp == line_out
        os.remove(out)


def test_plot_dist():
    """
    For all genomes, plot the distribution of their L90 values, and their number of contigs.
    Add a vertical line at the given threshold.
    genomes: {genome: [name, path, size, nbcont, l90]}
    """
    genomes_dir = os.path.join("test", "data", "genomes")
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
          "genome5.fasta", "genome6.fasta", "genome7.fasta"]
    genomes = {gs[0]: ["SAEN.1113", os.path.join(genomes_dir, gs[0]), 51, 2, 2],
               gs[1]: ["SAEN.1114", os.path.join(genomes_dir, gs[1]), 67, 15, 13],
               gs[2]: ["ESCO.0416", os.path.join(genomes_dir, gs[2]), 70, 15, 11],
               gs[3]: ["ESCO.0216", os.path.join(genomes_dir, gs[3]), 114, 17, 11],
               gs[4]: ["SAEN.1115", os.path.join(genomes_dir, gs[4]), 106, 17, 12],
               gs[5]: ["ESCO.0216", os.path.join(genomes_dir, gs[5]), 116, 60, 50],
               gs[6]: ["SAEN.1115", os.path.join(genomes_dir, gs[6]), 137, 20, 12]}
    res_path = os.path.join("test", "data")
    exp_path = os.path.join(res_path, "exp_files")
    listfile_base = "test_plot_dist"
    l90 = 13
    nbconts = 19
    gfunc.plot_distributions(genomes, res_path, listfile_base, l90, nbconts)
    outfiles = [os.path.join(res_path, "QC_L90-" + listfile_base + ".png"),
                os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")]
    expfiles = [os.path.join(exp_path, "res_QC_L90-" + listfile_base + ".png"),
                os.path.join(exp_path, "res_QC_nb-contigs-" + listfile_base + ".png")]
    for out, exp in zip(outfiles, expfiles):
        assert util_tests.compare_files(out, exp)
        os.remove(out)

