#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pytest
import os
import logging

import test.test_unit.utilities_for_tests as util
import PanACoTA.annotate_module.genome_seq_functions as gfunc

import matplotlib
matplotlib.use('AGG')

# Define variables used by several tests
DBPATH = os.path.join("test", "data", "annotate", "genomes")
TMP_PATH = os.path.join('test', 'data', 'annotate', "tmp_files")
EXP_DIR = os.path.join('test', 'data', 'annotate', 'exp_files')
logger = logging.getLogger('test_genome_func')
# BASELINE_DIR = os.path.join("..", "..", "data", "annotate", "exp_files", "baseline")


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

    genomes = {gs[0]: ["SAEN.1113", os.path.join(genomes_dir, gs[0]), "pathtoseq1", 51, 4, 2],
               gs[1]: ["SAEN.1114", os.path.join(genomes_dir, gs[1]), "pathToSeq2", 67, 3, 3],
               gs[2]: ["ESCO.0416", os.path.join(genomes_dir, gs[2]), "pathToSeq3", 70, 4, 1],
               gs[3]: ["ESCO.0216", os.path.join(genomes_dir, gs[3]), "pathToSeq4", 114, 5, 2],
               gs[4]: ["SAEN.1115", os.path.join(genomes_dir, gs[4]), "path_to_seq5", 106, 3, 1],
               gs[5]: ["ESCO.0216", os.path.join(genomes_dir, gs[5]), "pathtoseq6", 116, 4, 2],
               gs[6]: ["SAEN.1115", os.path.join(genomes_dir, gs[6]), "pathtoseq7", 137, 3, 2]}
    gfunc.rename_all_genomes(genomes)
    # SAEN genomes 1 and 2 have same characteristics. Their place will be chosen randomly,
    # so take into account both choices
    exp_genomes = {gs[0]: ["SAEN.1113.00003",
                           os.path.join(genomes_dir, gs[0]), "pathtoseq1", 51, 4, 2],
                   gs[1]: ["SAEN.1114.00004",
                           os.path.join(genomes_dir, gs[1]), "pathToSeq2", 67, 3, 3],
                   gs[2]: ["ESCO.0416.00001",
                           os.path.join(genomes_dir, gs[2]), "pathToSeq3", 70, 4, 1],
                   gs[3]: ["ESCO.0216.00003",
                           os.path.join(genomes_dir, gs[3]), "pathToSeq4", 114, 5, 2],
                   gs[4]: ["SAEN.1115.00001",
                           os.path.join(genomes_dir, gs[4]), "path_to_seq5", 106, 3, 1],
                   gs[5]: ["ESCO.0216.00002",
                           os.path.join(genomes_dir, gs[5]), "pathtoseq6", 116, 4, 2],
                   gs[6]: ["SAEN.1115.00002",
                           os.path.join(genomes_dir, gs[6]), "pathtoseq7", 137, 3, 2]}
    assert genomes == exp_genomes


def test_get_outdir_prodigal_nocut():
  """
  When we use prodigal, and do not cut at each stretch of 5N, no need to create a modified
  sequence. So, check that get_output_dir returns an empty res_path (as we will not create
  any new sequence)
  """
  soft = "prodigal"
  genome = "genome1.fasta"
  cut = False
  pat = None
  gpath, grespath = gfunc.get_output_dir(soft, DBPATH, TMP_PATH, genome, cut, pat)
  assert not grespath
  assert gpath == os.path.join(DBPATH, "genome1.fasta")


def test_get_outdir_prodigal_cut():
  """
  When using prodigal, and cut at each stretch of 4 'N', check that the output
  dir where modified sequence must be store is as expected (in tmp_dir, and adding "-split4N").
  """
  soft = "prodigal"
  genome = "genome1.fasta"
  cut = True
  pat = "NNNN+"
  gpath, grespath = gfunc.get_output_dir(soft, DBPATH, TMP_PATH, genome, cut, pat)
  assert gpath == os.path.join(DBPATH, "genome1.fasta")
  assert grespath == os.path.join(TMP_PATH, "genome1.fasta_prodigal-split4N.fna")


def test_get_outdir_prokka_nocut():
  """
  When using prokka, and don't cut sequence, check that the output
  dir where modified sequence must be store is as expected (in tmp dir, indicating
  "shorter_contigs" because we need to shorten the contig headers for prokka)
  """
  soft = "prokka"
  genome = "genome1.fasta"
  cut = False
  pat = None
  gpath, grespath = gfunc.get_output_dir(soft, DBPATH, TMP_PATH, genome, cut, pat)
  assert gpath == os.path.join(DBPATH, "genome1.fasta")
  assert grespath == os.path.join(TMP_PATH, "genome1.fasta_prokka-shorter-contigs.fna")


def test_get_outdir_prokka_cut():
  """
  When using prokka, and cut at each stretch of 3 'N', check that the output
  dir where modified sequence must be store is as expected (in tmp_di, indicating "-split3N"
  because sequence will be cut)
  """
  soft = "prokka"
  genome = "genome1.fasta"
  cut = True
  pat = "NNN+"
  gpath, grespath = gfunc.get_output_dir(soft, DBPATH, TMP_PATH, genome, cut, pat)
  assert gpath == os.path.join(DBPATH, "genome1.fasta")
  assert grespath == os.path.join(TMP_PATH, "genome1.fasta_prokka-split3N.fna")


def test_get_outdir_no_input_seq():
    """
    Test that when the given genome name does not exist in the dbpath, it exists in tmp path.
    It means that it is a concatenation of 2 files of dbpath, which was saved in tmppath
    """
    soft = "prodigal"
    genome = "prokka_out_for_test-supHeader.faa"
    tmp_path = os.path.join("test", "data", "annotate", "test_files")
    cut = True
    pat = "NNNNNNN+"
    gpath, grespath = gfunc.get_output_dir(soft, DBPATH, tmp_path, genome, cut, pat)
    assert gpath == os.path.join(tmp_path, "prokka_out_for_test-supHeader.faa")
    assert grespath == os.path.join(tmp_path,
                                  "prokka_out_for_test-supHeader.faa_prodigal-split7N.fna")


def test_split_contig_nocut():
    """
    Test that when a contig must not be cut, it returns the current number of contigs + 1
    (we just add this contig, which is not split), and writes
    this new contig to the new file
    """
    pat = None
    whole_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name_for_my_sequence"
    contig_sizes = {"contig_1": 10}
    resfile = os.path.join("test", "data", "annotate", "test_split_contig_nocut.fna")
    gresf = open(resfile, "w")
    num = 2

    num = gfunc.split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num)
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_nocut.fna")
    assert num == 3
    assert os.path.exists(resfile)
    assert util.compare_order_content(resfile, exp_file)

    # Remove created file
    os.remove(resfile)


def test_split_contig_cut():
    """
    Test that when a contig must not be cut at each stretch of 3 'N', and the sequence
    contains 1 stretch of 5 'N', it cuts the sequence into 2 contigs, ignoring all 'N' of the
    stretch (even if there are more than 3 'N'). Then, it returns the current number of contigs + 2
    and writes those new contigs to the new file
    """
    pat = "NNN+"
    whole_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name_for_my_sequence"
    contig_sizes = {">contig_1": 10}
    resfile = os.path.join("test", "data", "annotate", "test_split_contig_nocut.fna")
    gresf = open(resfile, "w")
    num = 2

    num = gfunc.split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num)
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_cut3N.fna")
    assert num == 4
    assert os.path.exists(resfile)
    assert util.compare_order_content(resfile, exp_file)

    # Remove created file
    os.remove(resfile)


def test_split_empty_contig():
    """
    Test that when we want to split an empty contig (for exemple, if the sequence starts with
    'NNNN', it returns 2 contigs: 1 empty and 1 with the rest of the sequence), this contig is
    just ignored (not written nor taken into account in contig num).
    """
    pat = "NNN+"
    whole_seq = "NNNNNAACTGCTTTTTAAGCGCGCTCCTGCGNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name_for_my_sequence"
    contig_sizes = {"contig_1": 10}
    resfile = os.path.join("test", "data", "annotate", "test_split_contig_nocut.fna")
    gresf = open(resfile, "w")
    num = 2

    num = gfunc.split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num)
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_empty_contig.fna")
    assert num == 3
    assert os.path.exists(resfile)
    assert util.compare_order_content(resfile, exp_file)

    # Remove created file
    os.remove(resfile)


def test_format_contig_cut():
    """
    For a given contig, if we want to annotate it with prodigal, and cut at each stretch of 5 'N'
    check that it writes this contig, split, in the expected file
    """
    cut = True
    pat = 'NNNNN+'
    cur_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name_for_my_sequence"
    contig_sizes = {}
    resfile = os.path.join("test", "data", "annotate", "test_format_cont_cut5N.fna")
    gresf = open(resfile, "w")
    num = 2

    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, contig_sizes, gresf,
                               num, logger=None) == 4
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_cut3N.fna")
    assert os.path.exists(resfile)
    assert util.compare_order_content(resfile, exp_file)
    assert contig_sizes == {">my_contig_name_for_my_sequence_2\n": 26,
                            ">my_contig_name_for_my_sequence_3\n": 25}

    # Remove created file
    os.remove(resfile)


def test_format_contig_nocut_prokka():
    """
    For a given contig, if we want to annotate it with prodigal, and cut at each stretch of 5 'N'
    check that it writes this contig, split, in the expected file
    """
    cut = False
    pat = None
    cur_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name_for_my_sequence"
    contig_sizes = {}
    resfile = os.path.join("test", "data", "annotate", "test_format_cont_nocut_prokka.fna")
    gresf = open(resfile, "w")
    num = 2

    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, contig_sizes, gresf,
                               num, logger=None) == 3
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_nocut.fna")
    assert os.path.exists(resfile)
    assert util.compare_order_content(resfile, exp_file)
    assert contig_sizes == {">my_contig_name_for_my_sequence_2\n": 56}

    # Remove created file
    os.remove(resfile)


def test_format_contig_nocut_prodigal_notSameName():
    """
    For a given contig, if we want to annotate it with prodigal, and do not cut, then we keep the same file. However, we must check that contig names are all different.
    Add 2 contigs, to be sure the 'num' parameter is not increased.
    """
    cut = False
    pat = None
    cur_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_seq2 = 'AACGTGGTCAGAGCGTG'
    cur_contig_name = ">my_contig_name_for_my_sequence"
    cur_contig_name2 = ">mycontigname"
    contig_sizes = {">mycontig": 155}
    num = 1

    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, contig_sizes, None,
                               num, logger=None) == 1
    assert gfunc.format_contig(cut, pat, cur_seq2, cur_contig_name2, contig_sizes, None,
                               num, logger=None) == 1

    assert contig_sizes == {">my_contig_name_for_my_sequence": 56,
                            ">mycontigname": 17,
                            ">mycontig": 155}


def test_format_contig_nocut_prodigal_SameName(caplog):
    """
    For a given contig, if we want to annotate it with prodigal, and do not cut, then we keep the same file. However, we must check that contig names are all different.
    Try to add a contig which name is already used, check that it prints the expected error,
    and returns -1
    """
    cut = False
    pat = None
    cur_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_seq2 = 'AACGTGGTCAGAGCGTG'
    cur_contig_name = ">my_contig_name_for_my_sequence"
    cur_contig_name2 = ">my_contig2"
    contig_sizes = {">mycontig": 155, ">my_contig_name_for_my_sequence":45}
    num = 1

    # Try to add a contig already existing -> error log
    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, contig_sizes, None,
                               num, logger) == -1
    assert contig_sizes == {">my_contig_name_for_my_sequence": 45,
                            ">mycontig": 155}
    # Check logs
    caplog.set_level(logging.DEBUG)
    assert (">my_contig_name_for_my_sequence contig name is used for several contigs. "
            "Please put different names for each contig. This genome will be "
            "ignored.") in caplog.text

    # Add a contig with new name. contig_sizes is completed
    assert gfunc.format_contig(cut, pat, cur_seq2, cur_contig_name2, contig_sizes, None,
                               num, logger) == 1
    assert contig_sizes == {">my_contig_name_for_my_sequence": 45,
                            ">my_contig2": 17,
                            ">mycontig": 155}


def test_analyse1genome_nocut_prodigal():
    """
    Analyse the given genome, without cutting at stretches of N, will be annotated by prodigal:
    -> no new file created.
    Calculate its genome size, nb contigs and L90, and add it to the genomes dict, as well as
    the path to the genome file.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = gs[1]
    cut = False
    pat = None
    soft = "prodigal"
    assert gfunc.analyse_genome(genome, DBPATH, TMP_PATH, cut, pat, genomes, soft)

    # Check that information on analyzed genome are correct. And path to 'genome to annotate'
    # is the same as the path to the genome itself
    outf = os.path.join(DBPATH, "genome2.fasta")
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114", outf, outf, 67, 3, 3],
                   gs[2]: ["ESCO.0416"]}
    assert genomes == exp_genomes


# def test_analyse1genome_cut_prodigal():
#     # New file created in tmp path, with same contig names, but contigs seqs cut at each stretch
#     # -> check 'genomes' + output file
#     assert False


# def test_analyse1genome_nocut_prokka():
#     # New file created in tmp with shortened contig names, but same content (not split)
#     # -> check 'genomes' + output file
#     assert False


# def test_analyse1genome_cut_prokka():
#     # New file created in tmp with shortened contig names, and contigs cut at each stretch
#     # -> check 'genomes' + output file
#     assert False

# def test_analyse1genome_empty():
#     assert False


# tests
# -> analyse genome
# -> analyse all genomes
# -> plot_distributions

# def test_analyse1genome_nocut_empty():
#     """
#     Analyse the given genome: without cutting at stretches of N. The genome is an empty
#     file, so it is not possible to calculate L90
#     """
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
#     open(os.path.join(DBPATH, gs[3]), "w").close()
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"],
#                gs[3]: ["ESCO.0415"]}
#     genome = gs[3]
#     tmp_path = os.path.join("test", "data", "annotate")
#     cut = False
#     pat = "NNNNN+"
#     assert not gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
#     exp_genomes = {gs[0]: ["SAEN.1113"],
#                    gs[1]: ["SAEN.1114"],
#                    gs[2]: ["ESCO.0416"],
#                    gs[3]: ["ESCO.0415"]}
#     assert genomes == exp_genomes
#     os.remove(os.path.join(DBPATH, gs[3]))
#     os.remove(os.path.join(tmp_path, gs[3] + "-short-contig.fna"))


# def test_analyse1genome_cut():
#     """
#     Analyse the given genome: cut at each stretch of 5 N, put it to a new file,
#     and then calculate its genome size, nb contigs and L90. Add this information
#     to the genomes dict, as well as the path to the genome file (cut).
#     """
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"]}
#     genome = gs[1]
#     tmp_path = os.path.join("test", "data", "annotate")
#     cut = True
#     pat = "NNNNN+"
#     assert gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
#     out_f = os.path.join(tmp_path, gs[1] + "-split5N.fna")
#     exp_f = os.path.join(tmp_path, "exp_files", "res_genome2.fasta-split5N.fna")
#     exp_genomes = {gs[0]: ["SAEN.1113"],
#                    gs[1]: ["SAEN.1114", out_f, 55, 5, 4],
#                    gs[2]: ["ESCO.0416"]}
#     assert genomes == exp_genomes
#     with open(out_f, "r") as outf, open(exp_f, "r") as expf:
#         for line_exp, line_out in zip(expf, outf):
#             assert line_exp == line_out
#     os.remove(out_f)


# def test_analyse1genome_cut_same_names():
#     """
#     Analyse a genome. Its contig names all have the same first 20 characters. There is no
#     stretch of at least 5N, so contigs are not split.
#     New contig names should be uniq, and not all ending with _0!
#     """
#     genome = "genome_long_header.fst"
#     genomes = {genome: ["SAEN.1015.0117"]}
#     tmp_path = os.path.join("test", "data", "annotate")
#     cut = True
#     pat = "NNNNN+"
#     assert gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
#     out_f = os.path.join(tmp_path, genome + "-split5N.fna")
#     exp_f = os.path.join("test", "data", "annotate", "exp_files",
#                          "res_genome_short-long_header.fst")
#     exp_genomes = {genome: ["SAEN.1015.0117", out_f, 151, 3, 3]}
#     assert genomes == exp_genomes
#     with open(out_f, "r") as outf, open(exp_f, "r") as expf:
#         for line_exp, line_out in zip(expf, outf):
#             assert line_exp == line_out
#     os.remove(out_f)


# def test_analyse1genome_cut_empty():
#     """
#     Analyse the given genome: cut at each stretch of 5 N, but the file is empty.
#     Check that it returns False
#     """
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
#     open(os.path.join(DBPATH, gs[3]), "w").close()
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"],
#                gs[3]: ["ESCO.0415"]}
#     genome = gs[3]
#     tmp_path = os.path.join("test", "data", "annotate")
#     cut = True
#     pat = "NNNNN+"
#     assert not gfunc.analyse_genome(genome, DBPATH, tmp_path, cut, pat, genomes)
#     exp_genomes = {gs[0]: ["SAEN.1113"],
#                    gs[1]: ["SAEN.1114"],
#                    gs[2]: ["ESCO.0416"],
#                    gs[3]: ["ESCO.0415"]}
#     assert genomes == exp_genomes
#     out_f = os.path.join(tmp_path, gs[3] + "-split5N.fna")
#     with open(out_f, "r") as outf:
#         assert outf.readlines() == []
#     os.remove(os.path.join(DBPATH, gs[3]))
#     os.remove(out_f)


# def test_analyse_all_genomes_nocut():
#     """
#     Analyze all given genomes: don't cut at stretches of N, but look at their sequence
#     file, to calculate L90, genome size and nb contigs. Add this information, as well as the
#     path to the genomic sequence, to the genomes dict.
#     """
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"]}
#     tmp_path = os.path.join("test", "data", "annotate")
#     opaths = [os.path.join(tmp_path, gname + "-short-contig.fna") for gname in gs]
#     nbn = 0
#     # Run analysis
#     gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
#     # construct expected results
#     exp_genomes = {gs[0]: ["SAEN.1113", opaths[0], 51, 4, 2],
#                    gs[1]: ["SAEN.1114", opaths[1], 67, 3, 3],
#                    gs[2]: ["ESCO.0416", opaths[2], 70, 4, 1]}
#     assert exp_genomes == genomes
#     for f in opaths:
#         os.remove(f)


# def test_analyse_all_genomes_nocut_empty():
#     """
#     Analyze all given genomes: don't cut at stretches of N, but look at their sequence
#     file, to calculate L90, genome size and nb contigs. Add this information, as well as the
#     path to the genomic sequence, to the genomes dict.
#     """
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
#     open(os.path.join(DBPATH, gs[3]), "w").close()
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"],
#                gs[3]: ["ESCO.0123"]}
#     tmp_path = os.path.join("test", "data", "annotate")
#     opaths = [os.path.join(tmp_path, gname + "-short-contig.fna") for gname in gs]
#     nbn = 0
#     # Run analysis
#     gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
#     # construct expected results
#     exp_genomes = {gs[0]: ["SAEN.1113", opaths[0], 51, 4, 2],
#                    gs[1]: ["SAEN.1114", opaths[1], 67, 3, 3],
#                    gs[2]: ["ESCO.0416", opaths[2], 70, 4, 1]}
#     assert exp_genomes == genomes
#     os.remove(os.path.join(DBPATH, gs[3]))
#     for f in opaths:
#         os.remove(f)


# def test_analyse_all_genomes_cut():
#     """
#     Analyze all given genomes: cut at each stretch of 5 'N', look at the output sequence
#     file, to calculate L90, genome size and nb contigs. Add this information, as well as the
#     path to the genomic sequence, to the genomes dict.
#     """
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"]}
#     tmp_path = os.path.join("test", "data", "annotate")
#     gpaths = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
#     nbn = 5
#     # Run analysis
#     gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
#     # construct expected results
#     exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], 51, 4, 2],
#                    gs[1]: ["SAEN.1114", gpaths[1], 55, 5, 4],
#                    gs[2]: ["ESCO.0416", gpaths[2], 70, 4, 1]}
#     out_f = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
#     exp_f = [os.path.join(tmp_path, "exp_files", "res_" + gname + "-split5N.fna") for gname in gs]
#     assert exp_genomes == genomes
#     for out, exp in zip(out_f, exp_f):
#         with open(out, "r") as outf, open(exp, "r") as expf:
#             for line_exp, line_out in zip(expf, outf):
#                 assert line_exp == line_out
#         os.remove(out)


# def test_analyse_all_genomes_cut_empty():
#     """
#     Analyze all given genomes: cut at each stretch of 5 'N', look at the output sequence
#     file, to calculate L90, genome size and nb contigs. Add this information, as well as the
#     path to the genomic sequence, to the genomes dict.
#     """
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "empty.fasta"]
#     open(os.path.join(DBPATH, gs[3]), "w").close()
#     genomes = {gs[0]: ["SAEN.1113"],
#                gs[1]: ["SAEN.1114"],
#                gs[2]: ["ESCO.0416"],
#                gs[3]: ["ESCO.0123"]}
#     tmp_path = os.path.join("test", "data", "annotate")
#     gpaths = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
#     nbn = 5
#     # Run analysis
#     gfunc.analyse_all_genomes(genomes, DBPATH, tmp_path, nbn)
#     # construct expected results
#     exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], 51, 4, 2],
#                    gs[1]: ["SAEN.1114", gpaths[1], 55, 5, 4],
#                    gs[2]: ["ESCO.0416", gpaths[2], 70, 4, 1]}
#     out_f = [os.path.join(tmp_path, gname + "-split5N.fna") for gname in gs]
#     exp_f = [os.path.join(tmp_path, "exp_files", "res_" + gname + "-split5N.fna") for gname in gs]
#     assert exp_genomes == genomes
#     for out, exp in zip(out_f[:-1], exp_f[:-1]):
#         with open(out, "r") as outf, open(exp, "r") as expf:
#             for line_exp, line_out in zip(expf, outf):
#                 assert line_exp == line_out
#         os.remove(out)
#     with open(out_f[-1], "r") as outf:
#         assert outf.readlines() == []
#     os.remove(out_f[-1])
#     os.remove(os.path.join(DBPATH, gs[3]))


# def get_plot_distribs():
#     """
#     For all genomes, plot the distribution of their L90 values, and their number of contigs.
#     Add a vertical line at the given threshold.
#     genomes: {genome: [name, path, size, nbcont, l90]}
#     output of plot_distributions is L90_vals, nbcont_vals, l90_dist, nbcont_dist
#     these outputs will be compared to expected results in tests
#     """
#     genomes_dir = os.path.join("test", "data", "annotate", "genomes")
#     gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
#           "genome5.fasta", "genome6.fasta", "genome7.fasta"]
#     genomes = {gs[0]: ["SAEN.1113", os.path.join(genomes_dir, gs[0]), 51, 2, 2],
#                gs[1]: ["SAEN.1114", os.path.join(genomes_dir, gs[1]), 67, 15, 13],
#                gs[2]: ["ESCO.0416", os.path.join(genomes_dir, gs[2]), 70, 15, 11],
#                gs[3]: ["ESCO.0216", os.path.join(genomes_dir, gs[3]), 114, 17, 11],
#                gs[4]: ["SAEN.1115", os.path.join(genomes_dir, gs[4]), 106, 17, 12],
#                gs[5]: ["ESCO.0216", os.path.join(genomes_dir, gs[5]), 116, 60, 50],
#                gs[6]: ["SAEN.1115", os.path.join(genomes_dir, gs[6]), 137, 20, 12]}
#     res_path = os.path.join("test", "data", "annotate")
#     listfile_base = "test_plot_dist"
#     l90 = 13
#     nbconts = 19
#     outdist = gfunc.plot_distributions(genomes, res_path, listfile_base, l90, nbconts)
#     return outdist


# @pytest.mark.mpl_image_compare(baseline_dir=BASELINE_DIR, tolerance=6, backend="agg")
# def test_dist_l90():
#     """
#     For created L90 graph, check that calculated L90 values are as expected,
#     and graph is also as expected
#     """
#     res_path = os.path.join("test", "data", "annotate")
#     listfile_base = "test_plot_dist"
#     outfile1 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
#     outfile2 = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
#     l90, _, dist, _ = get_plot_distribs()
#     # Check that png file was created
#     assert os.path.isfile(outfile1)
#     assert os.path.isfile(outfile2)
#     os.remove(outfile1)
#     os.remove(outfile2)
#     # Check values calculated for l90
#     assert set(l90) == {2, 13, 11, 11, 12, 50, 12}
#     # Check that output plot is as expected
#     return dist


# @pytest.mark.mpl_image_compare(baseline_dir=BASELINE_DIR, tolerance=6, backend="agg")
# def test_dist_nbcont():
#     """
#     For created L90 graph, check that calculated L90 values are as expected,
#     and graph is also as expected
#     """
#     res_path = os.path.join("test", "data", "annotate")
#     listfile_base = "test_plot_dist"
#     outfile1 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
#     outfile2 = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
#     _, nbcont, _, dist = get_plot_distribs()
#     # Check that png file was created
#     assert os.path.isfile(outfile1)
#     assert os.path.isfile(outfile2)
#     os.remove(outfile1)
#     os.remove(outfile2)
#     # Check values calculated for l90
#     assert set(nbcont) == {2, 15, 15, 17, 17, 60, 20}
#     # Check that output plot is as expected
#     return dist
