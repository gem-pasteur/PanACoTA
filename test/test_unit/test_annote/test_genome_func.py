#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pytest
import os
import PanACoTA.annote_module.genome_seq_functions as gfunc
import matplotlib
matplotlib.use('AGG')


# Define variables used by several tests
DBPATH = os.path.join("test", "data", "annotate", "genomes")
BASELINE_DIR = os.path.join("..", "..", "data", "annotate", "exp_files", "baseline")


# Start tests
def test_sort_genomes():
    """
    Test the function sorting genomes by L90 and nb contigs.
    genome = name, path, gsize, nbcont, L90]
    """
    genome1 = ["SAEN.1116.", "path/to/genome1", 10000, 11, 2]
    genome2 = ["SAEN.1015.", "path/to/genome2", 10000, 12, 2]
    genome3 = ["SAEN.1015.", "path/to/genome3", 10000, 12, 1]
    genome4 = ["ESCO.0216.", "path/to/genome4", 10000, 12, 1]

    genomes = {1: genome1, 2: genome2, 3: genome3, 4: genome4}
    sorted_g = sorted(genomes.items(), key=gfunc.sort_genomes)
    exp = [(4, genome4), (3, genome3), (1, genome1), (2, genome2)]
    assert sorted_g == exp


def test_save_contig_5n():
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
    seq_file = os.path.join("test", "data", "annotate", "test_save_contig5N.faa")
    resf = open(seq_file, "w")
    gfunc.save_contig(pat, cur_cont, cur_cont_name, contig_sizes, resf, -1)
    resf.close()
    exp = {">ESCO.0216.00001_cont1": 1623, ">ESCO.0216.0000_0": 39,
           ">ESCO.0216.0000_1": 33, ">ESCO.0216.0000_2": 20}
    assert contig_sizes == exp

    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_save_contig5N.faa")
    with open(exp_file, "r") as expf, open(seq_file, "r") as seqf:
        for line_exp, line_seq in zip(expf, seqf):
            assert line_exp == line_seq
    os.remove(seq_file)


def test_save_contig_atcg():
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
    seq_file = os.path.join("test", "data", "annotate", "test_save_contigATCG.faa")
    resf = open(seq_file, "w")
    gfunc.save_contig(pat, cur_cont, cur_cont_name, contig_sizes, resf, -1)
    resf.close()
    exp = {">ESCO.0216.00001_cont1": 1623, ">ESCO.0216.0000_0": 14,
           ">ESCO.0216.0000_1": 11, ">ESCO.0216.0000_2": 13,
           ">ESCO.0216.0000_3": 34, ">ESCO.0216.0000_4": 1,
           ">ESCO.0216.0000_5": 6, ">ESCO.0216.0000_6": 5}
    assert contig_sizes == exp

    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_save_contigATCG.faa")
    with open(exp_file, "r") as expf, open(seq_file, "r") as seqf:
        for line_exp, line_seq in zip(expf, seqf):
            assert line_exp == line_seq
    os.remove(seq_file)


def test_calc_l90_exact():
    """
    Calculate L90 according to the given genome size and contig sizes
    2 contigs get exactly 90% of the genome
    """
    cont_size = {1: 3, 2: 800, 3: 100, 4: 90, 5: 7}
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


def test_rename_genomes():
    """
    From a list of genomes ({genome: [name.date, path, gsize, nbcont, L90]}),
    order them by species, and by decreasing quality (L90, nb_cont), and rename them,
    as well as their contigs.
    """
    genomes_dir = os.path.join("test", "data", "annotate", "genomes")
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
          "genome5.fasta", "genome6.fasta", "genome7.fasta"]

    genomes = {gs[0]: ["SAEN.1113", os.path.join(genomes_dir, gs[0]), 51, 4, 2],
               gs[1]: ["SAEN.1114", os.path.join(genomes_dir, gs[1]), 67, 3, 3],
               gs[2]: ["ESCO.0416", os.path.join(genomes_dir, gs[2]), 70, 4, 1],
               gs[3]: ["ESCO.0216", os.path.join(genomes_dir, gs[3]), 114, 5, 2],
               gs[4]: ["SAEN.1115", os.path.join(genomes_dir, gs[4]), 106, 3, 1],
               gs[5]: ["ESCO.0216", os.path.join(genomes_dir, gs[5]), 116, 4, 2],
               gs[6]: ["SAEN.1115", os.path.join(genomes_dir, gs[6]), 137, 3, 2]}
    gfunc.rename_all_genomes(genomes)
    # SAEN genomes 1 and 2 have same characteristics. Their place will be chosen randomly,
    # so take into account both choices
    exp_genomes = {gs[0]: ["SAEN.1113.00003", os.path.join(genomes_dir, gs[0]), 51, 4, 2],
                   gs[1]: ["SAEN.1114.00004", os.path.join(genomes_dir, gs[1]), 67, 3, 3],
                   gs[2]: ["ESCO.0416.00001", os.path.join(genomes_dir, gs[2]), 70, 4, 1],
                   gs[3]: ["ESCO.0216.00003", os.path.join(genomes_dir, gs[3]), 114, 5, 2],
                   gs[4]: ["SAEN.1115.00001", os.path.join(genomes_dir, gs[4]), 106, 3, 1],
                   gs[5]: ["ESCO.0216.00002", os.path.join(genomes_dir, gs[5]), 116, 4, 2],
                   gs[6]: ["SAEN.1115.00002", os.path.join(genomes_dir, gs[6]), 137, 3, 2]}
    assert genomes == exp_genomes


def test_analyse1genome_nocut():
    """
    Analyse the given genome: without cutting at stretches of N, calculate
    its genome size, nb contigs and L90, and add it to the genomes dict, as well as
    the path to the genome file.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = gs[1]
    tmp_path = os.path.join("test", "data", "annotate")
    # Put genome file in tmppath instead of dbpath, as if it was
    # the result of concatenation of several files for the same genome,
    # done in the first step.
    orig_file = os.path.join(DBPATH, genome)
    out_file = os.path.join(tmp_path, genome)
    os.rename(orig_file, out_file)
    cut = False
    pat = "NNNNN+"
    assert gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
    outf = os.path.join(tmp_path, gs[1] + "-short-contig.fna")
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114", outf, 67, 3, 3],
                   gs[2]: ["ESCO.0416"]}
    assert genomes == exp_genomes
    exp_file = os.path.join(tmp_path, "exp_files", "res_test_analyse-genome2.fna")
    with open(exp_file, "r") as expf, open(outf, "r") as of:
        for linee, lineo in zip(expf, of):
            assert linee == lineo
    os.remove(outf)
    os.rename(out_file, orig_file)


def test_analyse1genome_nocut_empty():
    """
    Analyse the given genome: without cutting at stretches of N. The genome is an empty
    file, so it is not possible to calculate L90
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
    open(os.path.join(DBPATH, gs[3]), "w").close()
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["ESCO.0415"]}
    genome = gs[3]
    tmp_path = os.path.join("test", "data", "annotate")
    cut = False
    pat = "NNNNN+"
    assert not gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114"],
                   gs[2]: ["ESCO.0416"],
                   gs[3]: ["ESCO.0415"]}
    assert genomes == exp_genomes
    os.remove(os.path.join(DBPATH, gs[3]))
    os.remove(os.path.join(tmp_path, gs[3] + "-short-contig.fna"))


def test_analyse1genome_cut():
    """
    Analyse the given genome: cut at each stretch of 5 N, put it to a new file,
    and then calculate its genome size, nb contigs and L90. Add this information
    to the genomes dict, as well as the path to the genome file (cut).
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = gs[1]
    tmp_path = os.path.join("test", "data", "annotate")
    cut = True
    pat = "NNNNN+"
    assert gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
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


def test_analyse1genome_cut_same_names():
    """
    Analyse a genome. Its contig names all have the same first 20 characters. There is no
    stretch of at least 5N, so contigs are not split.
    New contig names should be uniq, and not all ending with _0!
    """
    genome = "genome_long_header.fst"
    genomes = {genome: ["SAEN.1015.0117"]}
    tmp_path = os.path.join("test", "data", "annotate")
    cut = True
    pat = "NNNNN+"
    assert gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
    out_f = os.path.join(tmp_path, genome + "-split5N.fna")
    exp_f = os.path.join("test", "data", "annotate", "exp_files",
                         "res_genome_short-long_header.fst")
    exp_genomes = {genome: ["SAEN.1015.0117", out_f, 151, 3, 3]}
    assert genomes == exp_genomes
    with open(out_f, "r") as outf, open(exp_f, "r") as expf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    os.remove(out_f)


def test_analyse1genome_cut_empty():
    """
    Analyse the given genome: cut at each stretch of 5 N, but the file is empty.
    Check that it returns False
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
    open(os.path.join(DBPATH, gs[3]), "w").close()
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["ESCO.0415"]}
    genome = gs[3]
    tmp_path = os.path.join("test", "data", "annotate")
    cut = True
    pat = "NNNNN+"
    assert not gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114"],
                   gs[2]: ["ESCO.0416"],
                   gs[3]: ["ESCO.0415"]}
    assert genomes == exp_genomes
    out_f = os.path.join(tmp_path, gs[3] + "-split5N.fna")
    with open(out_f, "r") as outf:
        assert outf.readlines() == []
    os.remove(os.path.join(DBPATH, gs[3]))
    os.remove(out_f)


def test_analyse_all_genomes_nocut():
    """
    Analyze all given genomes: don't cut at stretches of N, but look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    tmp_path = os.path.join("test", "data", "annotate")
    opaths = [os.path.join(tmp_path, gname + "-short-contig.fna") for gname in gs]
    nbn = 0
    # Run analysis
    gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
    # construct expected results
    exp_genomes = {gs[0]: ["SAEN.1113", opaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", opaths[1], 67, 3, 3],
                   gs[2]: ["ESCO.0416", opaths[2], 70, 4, 1]}
    assert exp_genomes == genomes
    for f in opaths:
        os.remove(f)


def test_analyse_all_genomes_nocut_empty():
    """
    Analyze all given genomes: don't cut at stretches of N, but look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
    open(os.path.join(DBPATH, gs[3]), "w").close()
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["ESCO.0123"]}
    tmp_path = os.path.join("test", "data", "annotate")
    opaths = [os.path.join(tmp_path, gname + "-short-contig.fna") for gname in gs]
    nbn = 0
    # Run analysis
    gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
    # construct expected results
    exp_genomes = {gs[0]: ["SAEN.1113", opaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", opaths[1], 67, 3, 3],
                   gs[2]: ["ESCO.0416", opaths[2], 70, 4, 1]}
    assert exp_genomes == genomes
    os.remove(os.path.join(DBPATH, gs[3]))
    for f in opaths:
        os.remove(f)


def test_analyse_all_genomes_cut():
    """
    Analyze all given genomes: cut at each stretch of 5 'N', look at the output sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    tmp_path = os.path.join("test", "data", "annotate")
    gpaths = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
    nbn = 5
    # Run analysis
    gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
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


def test_analyse_all_genomes_cut_empty():
    """
    Analyze all given genomes: cut at each stretch of 5 'N', look at the output sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
    open(os.path.join(DBPATH, gs[3]), "w").close()
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["ESCO.0123"]}
    tmp_path = os.path.join("test", "data", "annotate")
    gpaths = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
    nbn = 5
    # Run analysis
    gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
    # construct expected results
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], 55, 5, 4],
                   gs[2]: ["ESCO.0416", gpaths[2], 70, 4, 1]}
    out_f = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
    exp_f = [os.path.join(tmp_path, "exp_files", "res_" + gname + "-split5N.fna") for gname in gs]
    assert exp_genomes == genomes
    for out, exp in zip(out_f[:-1], exp_f[:-1]):
        with open(out, "r") as outf, open(exp, "r") as expf:
            for line_exp, line_out in zip(expf, outf):
                assert line_exp == line_out
        os.remove(out)
    with open(out_f[-1], "r") as outf:
        assert outf.readlines() == []
    os.remove(out_f[-1])
    os.remove(os.path.join(DBPATH, gs[3]))


def get_plot_distribs():
    """
    For all genomes, plot the distribution of their L90 values, and their number of contigs.
    Add a vertical line at the given threshold.
    genomes: {genome: [name, path, size, nbcont, l90]}
    output of plot_distributions is L90_vals, nbcont_vals, l90_dist, nbcont_dist
    these outputs will be compared to expected results in tests
    """
    genomes_dir = os.path.join("test", "data", "annotate", "genomes")
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
          "genome5.fasta", "genome6.fasta", "genome7.fasta"]
    genomes = {gs[0]: ["SAEN.1113", os.path.join(genomes_dir, gs[0]), 51, 2, 2],
               gs[1]: ["SAEN.1114", os.path.join(genomes_dir, gs[1]), 67, 15, 13],
               gs[2]: ["ESCO.0416", os.path.join(genomes_dir, gs[2]), 70, 15, 11],
               gs[3]: ["ESCO.0216", os.path.join(genomes_dir, gs[3]), 114, 17, 11],
               gs[4]: ["SAEN.1115", os.path.join(genomes_dir, gs[4]), 106, 17, 12],
               gs[5]: ["ESCO.0216", os.path.join(genomes_dir, gs[5]), 116, 60, 50],
               gs[6]: ["SAEN.1115", os.path.join(genomes_dir, gs[6]), 137, 20, 12]}
    res_path = os.path.join("test", "data", "annotate")
    listfile_base = "test_plot_dist"
    l90 = 13
    nbconts = 19
    outdist = gfunc.plot_distributions(genomes, res_path, listfile_base, l90, nbconts)
    return outdist


@pytest.mark.mpl_image_compare(baseline_dir=BASELINE_DIR, tolerance=6, backend="agg")
def test_dist_l90():
    """
    For created L90 graph, check that calculated L90 values are as expected,
    and graph is also as expected
    """
    res_path = os.path.join("test", "data", "annotate")
    listfile_base = "test_plot_dist"
    outfile1 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
    outfile2 = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
    l90, _, dist, _ = get_plot_distribs()
    # Check that png file was created
    assert os.path.isfile(outfile1)
    assert os.path.isfile(outfile2)
    os.remove(outfile1)
    os.remove(outfile2)
    # Check values calculated for l90
    assert set(l90) == {2, 13, 11, 11, 12, 50, 12}
    # Check that output plot is as expected
    return dist


@pytest.mark.mpl_image_compare(baseline_dir=BASELINE_DIR, tolerance=6, backend="agg")
def test_dist_nbcont():
    """
    For created L90 graph, check that calculated L90 values are as expected,
    and graph is also as expected
    """
    res_path = os.path.join("test", "data", "annotate")
    listfile_base = "test_plot_dist"
    outfile1 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
    outfile2 = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
    _, nbcont, _, dist = get_plot_distribs()
    # Check that png file was created
    assert os.path.isfile(outfile1)
    assert os.path.isfile(outfile2)
    os.remove(outfile1)
    os.remove(outfile2)
    # Check values calculated for l90
    assert set(nbcont) == {2, 15, 15, 17, 17, 60, 20}
    # Check that output plot is as expected
    return dist
