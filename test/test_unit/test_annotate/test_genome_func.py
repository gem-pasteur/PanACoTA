#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pytest
import os
import logging
import shutil

import test.test_unit.utilities_for_tests as tutil
import PanACoTA.annotate_module.genome_seq_functions as gfunc

import matplotlib
matplotlib.use('AGG')

# Define variables used by several tests
DBDIR = os.path.join("test", "data", "annotate")
GEN_PATH = os.path.join(DBDIR, "genomes")
TMP_PATH = os.path.join(DBDIR, "tmp_files")
EXP_DIR = os.path.join(DBDIR, 'exp_files')
BASELINE_DIR = os.path.abspath(os.path.join(EXP_DIR, "baseline"))
GENEPATH = os.path.join(DBDIR, "generated_by_unit-tests")
logger = logging.getLogger('test_genome_func')


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
    # utils.init_logger(LOGFILE_BASE, 0, 'test_postalign', verbose=1)
    if os.path.isdir(GENEPATH):
        content = os.listdir(GENEPATH)
        for f in content:
            assert f.startswith(".fuse")
    else:
        os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")



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


def test_calc_l90_error():
    """
    Calculate L90 according to the given genome size and contig sizes
    3 contigs get exactly more than 90%, but 2 contigs get less -> l90 = 3
    """
    cont_size = {}
    assert not gfunc.calc_l90(cont_size)


def test_rename_genomes():
    """
    From a list of genomes ({genome: [name.date, path, path_to_seq, gsize, nbcont, L90]}),
    order them by species, and by decreasing quality (L90, nb_cont), and rename them,
    as well as their contigs.
    """
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
          "genome5.fasta", "genome6.fasta", "genome7.fasta"]

    genomes = {gs[0]: ["SAEN.1113", os.path.join(GEN_PATH, gs[0]), "pathtoseq1", 51, 4, 2],
               gs[1]: ["SAEN.1114", os.path.join(GEN_PATH, gs[1]), "pathToSeq2", 67, 3, 3],
               gs[2]: ["ESCO.0416", os.path.join(GEN_PATH, gs[2]), "pathToSeq3", 70, 4, 1],
               gs[3]: ["ESCO.0216", os.path.join(GEN_PATH, gs[3]), "pathToSeq4", 114, 5, 2],
               gs[4]: ["SAEN.1115", os.path.join(GEN_PATH, gs[4]), "path_to_seq5", 106, 3, 1],
               gs[5]: ["ESCO.0216", os.path.join(GEN_PATH, gs[5]), "pathtoseq6", 116, 4, 2],
               gs[6]: ["SAEN.1115", os.path.join(GEN_PATH, gs[6]), "pathtoseq7", 137, 3, 2]}
    gfunc.rename_all_genomes(genomes)
    # SAEN genomes 1 and 2 have same characteristics. Their place will be chosen randomly,
    # so take into account both choices
    exp_genomes = {gs[0]: ["SAEN.1113.00003",
                           os.path.join(GEN_PATH, gs[0]), "pathtoseq1", 51, 4, 2],
                   gs[1]: ["SAEN.1114.00004",
                           os.path.join(GEN_PATH, gs[1]), "pathToSeq2", 67, 3, 3],
                   gs[2]: ["ESCO.0416.00001",
                           os.path.join(GEN_PATH, gs[2]), "pathToSeq3", 70, 4, 1],
                   gs[3]: ["ESCO.0216.00003",
                           os.path.join(GEN_PATH, gs[3]), "pathToSeq4", 114, 5, 2],
                   gs[4]: ["SAEN.1115.00001",
                           os.path.join(GEN_PATH, gs[4]), "path_to_seq5", 106, 3, 1],
                   gs[5]: ["ESCO.0216.00002",
                           os.path.join(GEN_PATH, gs[5]), "pathtoseq6", 116, 4, 2],
                   gs[6]: ["SAEN.1115.00002",
                           os.path.join(GEN_PATH, gs[6]), "pathtoseq7", 137, 3, 2]}
    assert genomes == exp_genomes


def test_split_contig_nocut():
    """
    Test that when a contig must not be cut, it returns the current number of contigs + 1
    (we just add this contig, which is not split), and writes
    this new contig to the new file
    """
    pat = None
    whole_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig name for my_sequence"
    contig_sizes = {"contig_1": 10}
    resfile = os.path.join(GENEPATH, "test_split_contig_nocut.fna")
    gresf = open(resfile, "w")
    num = 2

    num = gfunc.split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num)
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_nocut.fna")
    assert num == 3
    assert os.path.exists(resfile)
    assert tutil.compare_order_content(resfile, exp_file)


def test_split_contig_cut():
    """
    Test that when a contig must not be cut at each stretch of 3 'N', and the sequence
    contains 1 stretch of 5 'N', it cuts the sequence into 2 contigs, ignoring all 'N' of the
    stretch (even if there are more than 3 'N'). Then, it returns the current number of contigs + 2
    and writes those new contigs to the new file
    """
    pat = "NNN+"
    whole_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name for_my_sequence"
    contig_sizes = {">contig_1": 10}
    resfile = os.path.join(GENEPATH, "test_split_contig_nocut.fna")
    gresf = open(resfile, "w")
    num = 2

    num = gfunc.split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num)
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_cut3N.fna")
    assert num == 4
    assert os.path.exists(resfile)
    assert tutil.compare_order_content(resfile, exp_file)


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
    resfile = os.path.join(GENEPATH, "test_split_contig_nocut.fna")
    gresf = open(resfile, "w")
    num = 2

    num = gfunc.split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num)
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_empty_contig.fna")
    assert num == 3
    assert os.path.exists(resfile)
    assert tutil.compare_order_content(resfile, exp_file)


def test_format_contig_cut():
    """
    For a given contig, if we want to annotate it, and cut at each stretch of 5 'N'
    check that it writes this contig, split, in the expected file
    """
    cut = True
    pat = 'NNNNN+'
    cur_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name for_my_sequence"
    contig_sizes = {}
    resfile = os.path.join(GENEPATH, "test_format_cont_cut5N.fna")
    gresf = open(resfile, "w")
    num = 2

    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, "genome", contig_sizes, gresf,
                               num, logger=None) == 4
    gresf.close()

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_cut3N.fna")
    assert os.path.exists(resfile)
    assert tutil.compare_order_content(resfile, exp_file)
    assert contig_sizes == {">2_my_contig_name for_my_sequence\n": 26,
                            ">3_my_contig_name for_my_sequence\n": 25}


def test_format_contig_nocut():
    """
    For a given contig, if we want to annotate it with prokka, and do not cut at each stretch of
    5 'N'check that it writes this contig as given
    """
    cut = False
    pat = None
    cur_seq = "AACTGCTTTTTAAGCGCGCTCCTGCGNNNNNGGTTGTGTGGGCCCAGAGCGAGNCG"
    cur_contig_name = ">my_contig_name_for_my_sequence\n"
    contig_sizes = {}
    resfile = os.path.join(GENEPATH, "test_format_cont_nocut_prokka.fna")
    gresf = None
    num = 2

    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, "genome", contig_sizes, gresf,
                               num, logger=None) == 2

    exp_file = os.path.join(EXP_DIR, "exp_split_contig_nocut.fna")
    assert not os.path.exists(resfile)
    assert contig_sizes == {">my_contig_name_for_my_sequence\n": 56}


def test_format_contig_nocut_notDuplicateName():
    """
    For a given contig, if we want to annotate it with prodigal, and do not cut,
    then we keep the same file (no need to split at 20 characters)
    However, we must check that contig names are all different.
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

    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, "genome", contig_sizes, None,
                               num, logger=None) == 1
    assert gfunc.format_contig(cut, pat, cur_seq2, cur_contig_name2, "genome", contig_sizes, None,
                               num, logger=None) == 1

    assert contig_sizes == {">my_contig_name_for_my_sequence": 56,
                            ">mycontigname": 17,
                            ">mycontig": 155}


def test_format_contig_nocut_DuplicateName(caplog):
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
    assert gfunc.format_contig(cut, pat, cur_seq, cur_contig_name, "genome", contig_sizes, None,
                               num, logger) == -1
    assert contig_sizes == {">my_contig_name_for_my_sequence": 45,
                            ">mycontig": 155}
    # Check logs
    caplog.set_level(logging.DEBUG)
    assert ("In genome genome, '>my_contig_name_for_my_sequence' contig name is used for "
            "several contigs. "
            "Please put different names for each contig. This genome will be "
            "ignored.") in caplog.text

    # Add a contig with new name. contig_sizes is completed
    assert gfunc.format_contig(cut, pat, cur_seq2, cur_contig_name2, "genome", contig_sizes, None,
                               num, logger) == 1
    assert contig_sizes == {">my_contig_name_for_my_sequence": 45,
                            ">my_contig2": 17,
                            ">mycontig": 155}


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
  gpath, grespath = gfunc.get_output_dir(soft, GEN_PATH, GENEPATH, genome, cut, pat)
  assert not grespath
  assert gpath == os.path.join(GEN_PATH, "genome1.fasta")


def test_get_outdir_prodigal_cut():
  """
  When using prodigal, and cut at each stretch of 4 'N', check that the output
  dir where modified sequence must be store is as expected (in tmp_dir, and adding "-split4N").
  """
  soft = "prodigal"
  genome = "genome1.fasta"
  cut = True
  pat = "NNNN+"
  gpath, grespath = gfunc.get_output_dir(soft, GEN_PATH, GENEPATH, genome, cut, pat)
  assert gpath == os.path.join(GEN_PATH, "genome1.fasta")
  assert grespath == os.path.join(GENEPATH, "genome1.fasta_prodigal-split4N.fna")


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
  gpath, grespath = gfunc.get_output_dir(soft, GEN_PATH, GENEPATH, genome, cut, pat)
  assert gpath == os.path.join(GEN_PATH, "genome1.fasta")
  assert not grespath


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
  gpath, grespath = gfunc.get_output_dir(soft, GEN_PATH, GENEPATH, genome, cut, pat)
  assert gpath == os.path.join(GEN_PATH, "genome1.fasta")
  assert grespath == os.path.join(GENEPATH, "genome1.fasta_prokka-split3N.fna")


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
    gpath, grespath = gfunc.get_output_dir(soft, GEN_PATH, tmp_path, genome, cut, pat)
    assert gpath == os.path.join(tmp_path, "prokka_out_for_test-supHeader.faa")
    assert grespath == os.path.join(tmp_path,
                                  "prokka_out_for_test-supHeader.faa_prodigal-split7N.fna")


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
    assert gfunc.analyse_genome(genome, GEN_PATH, GENEPATH, cut, pat, genomes, soft, logger)

    # Check that information on analyzed genome are correct. And path to 'genome to annotate'
    # is the same as the path to the genome itself
    outf = os.path.join(GEN_PATH, "genome2.fasta")
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114", outf, outf, 67, 3, 3],
                   gs[2]: ["ESCO.0416"]}
    assert genomes == exp_genomes


def test_analyse1genome_cut_prodigal():
    '''
    Analyse the given genome, cutting at stretches of 5N, in order to annotate it
    Create new genome file in outdir, calculate genome size, nb contigs and L90, and add it
    to the genomes dict, as well as the path to the genome file.
    '''
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = gs[1]
    cut = True
    pat = "NNNNN+"
    soft = "prodigal"
    assert gfunc.analyse_genome(genome, GEN_PATH, GENEPATH, cut, pat, genomes, soft, logger)

    # Check that information on analyzed genome are correct. And path to 'genome to annotate'
    # is the same as the path to the genome itself
    initf = os.path.join(GEN_PATH, "genome2.fasta")  # initial genome path
    outf = os.path.join(GENEPATH, "genome2.fasta_prodigal-split5N.fna")  # path to geerated genome
    exp_out = os.path.join(EXP_DIR, "genome2-split5N.fna") # expected generated genome
    assert os.path.isfile(outf)
    assert tutil.compare_order_content(outf, exp_out)
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114", initf, outf, 55, 5, 4],
                   gs[2]: ["ESCO.0416"]}
    assert genomes == exp_genomes


def test_analyse1genome_cut_prokka():
    '''
    Analyse the given genome, cutting at stretches of 5N, in order to annotate it with prokka
    Create new genome file in outdir, with shortened contig names, calculate genome size,
    nb contigs and L90, and add it
    to the genomes dict, as well as the path to the genome file.
    '''
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = gs[1]
    cut = True
    pat = "NNNNN+"
    soft = "prokka"
    assert gfunc.analyse_genome(genome, GEN_PATH, GENEPATH, cut, pat, genomes, soft, logger)

    # Check that information on analyzed genome are correct. And path to 'genome to annotate'
    # is the same as the path to the genome itself
    initf = os.path.join(GEN_PATH, "genome2.fasta")  # initial genome path
    outf = os.path.join(GENEPATH, "genome2.fasta_prokka-split5N.fna")  # path to geerated genome
    exp_out = os.path.join(EXP_DIR, "genome2-split5N.fna") # expected generated genome
    assert os.path.isfile(outf)
    assert tutil.compare_order_content(outf, exp_out)
    exp_genomes = {gs[0]: ["SAEN.1113"],
                   gs[1]: ["SAEN.1114", initf, outf, 55, 5, 4],
                   gs[2]: ["ESCO.0416"]}
    assert genomes == exp_genomes


def test_analyse1genome_cut_same_names():
    """
    Analyse a genome. Its contig names all have the same first 20 characters. There is no
    stretch of at least 5N, so contigs are not split.
    New contig names should be uniq, and not all starting with '1_'!
    """
    genome = "genome_long_header.fst"
    genomes = {genome: ["SAEN.1015.0117"]}
    cut = True
    pat = "NNNNN+"
    assert gfunc.analyse_genome(genome, GEN_PATH, GENEPATH, cut, pat, genomes, "prokka", logger)
    in_f = os.path.join(GEN_PATH, genome)  # genome given
    out_f = os.path.join(GENEPATH, genome + "_prokka-split5N.fna")  # genome generated
    exp_f = os.path.join(EXP_DIR, "res_genome_short-long_header.fst") # expected genome generated
    exp_genomes = {genome: ["SAEN.1015.0117", in_f, out_f, 151, 3, 3]}
    assert genomes == exp_genomes
    assert tutil.compare_order_content(out_f, exp_f)


def test_analyse1genome_same_names_nocut(caplog):
    """
    Analyse a genome. 2 contigs have the same name, and we do not generate a new file (no cut)
    should return false with corresponding error message
    """
    caplog.set_level(logging.DEBUG)
    genome = "genome_2_identical_headers.fst"
    genomes = {genome: ["SAEN.1015.0117"]}
    cut = False
    pat = None
    assert not gfunc.analyse_genome(genome, GEN_PATH, GENEPATH, cut, pat, genomes,
                                    "prodigal", logger)
    assert ("In genome genome_2_identical_headers.fst, '>myheader' contig name is used for "
            "several contigs. Please put different names for "
            "each contig. This genome will be ignored") in caplog.text


def test_analyse1genome_same_last_name_nocut(caplog):
    """
    Analyse a genome. 2 contigs have the same name, and we do not generate a new file (no cut)
    should return false with corresponding error message
    """
    caplog.set_level(logging.DEBUG)
    genome = "genome_2_identical_headers-lastContig.fst"
    genomes = {genome: ["SAEN.1015.0117"]}
    cut = False
    pat = None
    assert not gfunc.analyse_genome(genome, GEN_PATH, GENEPATH, cut, pat, genomes,
                                    "prodigal", logger)
    assert ("In genome genome_2_identical_headers-lastContig.fst, '>myheader' contig name "
            "is used for several contigs. Please put different names for "
            "each contig. This genome will be ignored") in caplog.text


def test_analyse1genome_nofile(caplog):
    '''
    Test that when we ask to analyse a genome whose sequence file does not exist, it returns false
    with corresponding error message.
    '''
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = "toto"
    cut = True
    pat = "NNNNN+"
    soft = "prodigal"
    assert not gfunc.analyse_genome(genome, GEN_PATH, GENEPATH, cut, pat, genomes, soft, logger)
    assert "The file toto does not exist"


def test_analyse1genome_empty(caplog):
    '''
    Test that when we ask to analyse a genome whose sequence file does not exist, it returns false
    with corresponding error message.
    '''
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    genome = "genome1.fasta"
    # create empty genome file
    open(os.path.join(GENEPATH, genome), "w").close()
    cut = True
    pat = "NNNNN+"
    soft = "prodigal"
    assert not gfunc.analyse_genome(genome, GENEPATH, GENEPATH, cut, pat, genomes, soft, logger)
    assert ("Your file test/data/annotate/generated_by_unit-tests/genome1.fasta does not contain "
            "any gene. Please check that you really gave a fasta sequence file") in caplog.text
    gres = os.path.join(GENEPATH, "toto.fna_prodigal-split5N.fna")
    assert not os.path.isfile(gres)


def test_analyse_all_genomes_nocut(caplog):
    """
    Analyze all given genomes: don't cut at stretches of N, but look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    nbn = 0
    # Run analysis
    gfunc.analyse_all_genomes(genomes, GEN_PATH, GENEPATH, nbn, "prokka", logger, quiet=False)
    # construct expected results
    gpaths = [os.path.join(GEN_PATH, gname) for gname in gs]
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], gpaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], gpaths[1], 67, 3, 3],
                   gs[2]: ["ESCO.0416", gpaths[2], gpaths[2], 70, 4, 1]}
    assert exp_genomes == genomes
    assert ("Calculating genome size, number of contigs, L90") in caplog.text


def test_analyse_all_genomes_binary(caplog):
    """
    Analyze all given genomes: don't cut at stretches of N, but look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    1 file is a binary file: write warning message and remove it from analysis.
    """
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome.fna.bin"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["BIN.1234"]}
    nbn = 0
    # Run analysis
    gfunc.analyse_all_genomes(genomes, GEN_PATH, GENEPATH, nbn, "prokka", logger, quiet=False)
    # construct expected results
    gpaths = [os.path.join(GEN_PATH, gname) for gname in gs]
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], gpaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], gpaths[1], 67, 3, 3],
                   gs[2]: ["ESCO.0416", gpaths[2], gpaths[2], 70, 4, 1]}
    assert exp_genomes == genomes
    assert ("Calculating genome size, number of contigs, L90") in caplog.text
    assert ("'genome.fna.bin' does not seem to be a fasta file. It "
            "will be ignored.") in caplog.text


def test_analyse_all_genomes_cut(caplog):
    """
    Analyze all given genomes: cut at stretches of 3N, and look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict.
    """
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"]}
    nbn = 3
    # Run analysis
    gfunc.analyse_all_genomes(genomes, GEN_PATH, GENEPATH, nbn, "prokka", logger, quiet=False)
    # construct expected results
    gpaths = [os.path.join(GEN_PATH, gname) for gname in gs]
    opaths = [os.path.join(GENEPATH, gname + "_prokka-split3N.fna") for gname in gs]
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], opaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], opaths[1], 51, 6, 5],
                   gs[2]: ["ESCO.0416", gpaths[2], opaths[2], 70, 4, 1]}
    assert exp_genomes == genomes
    assert ("Cutting genomes at each time there are at least 3 'N' in a row, "
            "and then, calculating genome size, number of contigs and L90.") in caplog.text


def test_analyse_all_genomes_nocut_empty(caplog):
    """
    Analyze all given genomes: don't cut at stretches of N, but look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict. 1 genome is empty -> should be removed
    """
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "empty.fasta", "genome3.fasta"]
    empty_genome = os.path.join(GEN_PATH, gs[2])
    # Add an empty genome to the original database
    open(empty_genome, "w").close()
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["ESCO.0123"]}
    nbn = 0
    # Run analysis
    gfunc.analyse_all_genomes(genomes, GEN_PATH, GENEPATH, nbn, "prokka", logger, quiet=True)
    # construct expected results
    gpaths = [os.path.join(GEN_PATH, gname) for gname in gs]
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], gpaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], gpaths[1], 67, 3, 3],
                   gs[3]: ["ESCO.0123", gpaths[3], gpaths[3], 70, 4, 1]}
    assert exp_genomes == genomes
    assert ("Calculating genome size, number of contigs, L90") in caplog.text
    assert ("Your file test/data/annotate/genomes/empty.fasta "
            "does not contain any gene. Please check that you really gave a "
            "fasta sequence file") in caplog.text

    # remove the empty genome
    os.remove(empty_genome)


def test_analyse_all_genomes_cut_empty(caplog):
    """
    Analyze all given genomes: cut at stretches of 3N, and look at their sequence
    file, to calculate L90, genome size and nb contigs. Add this information, as well as the
    path to the genomic sequence, to the genomes dict. 1 genome is empty -> should be removed
    """
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "empty.fasta", "genome3.fasta"]
    empty_genome = os.path.join(GEN_PATH, gs[2])
    # Add an empty genome to the original database
    open(empty_genome, "w").close()
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["ESCO.0123"]}
    nbn = 3
    # Run analysis
    gfunc.analyse_all_genomes(genomes, GEN_PATH, GENEPATH, nbn, "prokka", logger, quiet=True)
    # construct expected results
    gpaths = [os.path.join(GEN_PATH, gname) for gname in gs]
    opaths = [os.path.join(GENEPATH, gname + "_prokka-split3N.fna") for gname in gs]
    exp_genomes = {gs[0]: ["SAEN.1113", gpaths[0], opaths[0], 51, 4, 2],
                   gs[1]: ["SAEN.1114", gpaths[1], opaths[1], 51, 6, 5],
                   gs[3]: ["ESCO.0123", gpaths[3], opaths[3], 70, 4, 1]}
    assert exp_genomes == genomes
    assert ("Cutting genomes at each time there are at least 3 'N' in a row, "
            "and then, calculating genome size, number of contigs and L90.") in caplog.text
    assert ("Your file test/data/annotate/genomes/empty.fasta "
            "does not contain any gene. Please check that you really gave a "
            "fasta sequence file") in caplog.text

    # remove the empty genome
    os.remove(empty_genome)


def test_analyse_all_genomes_noseq(caplog):
    """
    Analyze all given genomes: no given sequence file exists
    -> Exits with error message
    """
    caplog.set_level(logging.DEBUG)
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta"]
    genomes = {gs[0]: ["SAEN.1113"],
               gs[1]: ["SAEN.1114"],
               gs[2]: ["ESCO.0416"],
               gs[3]: ["ESCO.0123"]}
    nbn = 3
    # Run analysis
    with pytest.raises(SystemExit):
        gfunc.analyse_all_genomes(genomes, "toto", GENEPATH, nbn, "prokka", logger, quiet=True)
    assert ("No genome was found in the database folder toto. See logfile "
            "for more information.") in caplog.text


@pytest.mark.mpl_image_compare(baseline_dir=BASELINE_DIR, tolerance=6, backend="agg")
def test_dist_l90():
    """
    For created L90 graph, check that calculated L90 values are as expected,
    and graph is also as expected
    """
    listfile_base = "test_plot_dist"
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
          "genome5.fasta", "genome6.fasta", "genome7.fasta"]
    genomes = {gs[0]: ["SAEN.1113", 'orig_path',  'annot_path', 51, 2, 2],
               gs[1]: ["SAEN.1114", 'orig_path', 'annot_path', 67, 15, 13],
               gs[2]: ["ESCO.0416", 'orig_path', 'annot_path', 70, 15, 11],
               gs[3]: ["ESCO.0216", 'orig_path', 'annot_path', 114, 17, 11],
               gs[4]: ["SAEN.1115", 'orig_path', 'annot_path', 106, 17, 12],
               gs[5]: ["ESCO.0216", 'orig_path', 'annot_path', 116, 60, 50],
               gs[6]: ["SAEN.1115", 'orig_path', 'annot_path', 137, 20, 12]}
    l90 = 13
    nbconts = 19
    l90, _, dist, _ = gfunc.plot_distributions(genomes, GENEPATH, listfile_base, l90, nbconts)
    outfile1 = os.path.join(GENEPATH, "QC_L90-" + listfile_base + ".png")
    outfile2 = os.path.join(GENEPATH, "QC_nb-contigs-" + listfile_base + ".png")
    # Check that png file was created
    assert os.path.isfile(outfile1)
    assert os.path.isfile(outfile2)
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
    listfile_base = "test_plot_dist"
    gs = ["genome1.fasta", "genome2.fasta", "genome3.fasta", "genome4.fasta",
          "genome5.fasta", "genome6.fasta", "genome7.fasta"]
    genomes = {gs[0]: ["SAEN.1113", 'orig_path',  'annot_path', 51, 2, 2],
               gs[1]: ["SAEN.1114", 'orig_path', 'annot_path', 67, 15, 13],
               gs[2]: ["ESCO.0416", 'orig_path', 'annot_path', 70, 15, 11],
               gs[3]: ["ESCO.0216", 'orig_path', 'annot_path', 114, 17, 11],
               gs[4]: ["SAEN.1115", 'orig_path', 'annot_path', 106, 17, 12],
               gs[5]: ["ESCO.0216", 'orig_path', 'annot_path', 116, 60, 50],
               gs[6]: ["SAEN.1115", 'orig_path', 'annot_path', 137, 20, 12]}
    l90 = 13
    nbconts = 19
    _, nbcont, _, dist = gfunc.plot_distributions(genomes, GENEPATH, listfile_base, l90, nbconts)
    outfile1 = os.path.join(GENEPATH, "QC_L90-" + listfile_base + ".png")
    outfile2 = os.path.join(GENEPATH, "QC_nb-contigs-" + listfile_base + ".png")
    # Check that png file was created
    assert os.path.isfile(outfile1)
    assert os.path.isfile(outfile2)
    # Check values calculated for nbconts
    assert set(nbcont) == {2, 15, 15, 17, 17, 60, 20}
    # Check that output plot is as expected
    return dist
