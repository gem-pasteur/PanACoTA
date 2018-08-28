#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import genomeAPCAT.utils as utils
import pytest
import os
import logging
import shutil
import matplotlib

matplotlib.use('AGG')

# Define variables used by several tests
BASELINE_DIR = os.path.join("..", "data", "annotate", "exp_files", "baseline")


# Start tests
def test_check_install():
    """
    Try to run prokka, which is installed, and check that there is no problem
    """
    assert utils.check_installed("prokka")


def test_check_install_error():
    """
    Try to run a command which does not exist, and check that it closes the program
    with exit code 1
    """
    assert not utils.check_installed("plop false command...")


@pytest.mark.mpl_image_compare(baseline_dir=BASELINE_DIR, tolerance=5, backend="agg")
def test_plot_dist():
    """
    Plot a given distribution, and check that output is as expected
    """
    values = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 10]
    limit = 3
    res_dir = os.path.join("test", "data", "annotate")
    os.makedirs(res_dir, exist_ok=True)
    # reffile = os.path.join("test", "data", "annotate", "exp_files", "res_plot_distr.png")
    title = "Distribution test"
    text = "Max L90 ="
    myfig = utils.plot_distr(values, limit, title, text)
    return myfig


def test_skipped_prokka(capsys):
    """
    Test that when the list of skipped genomes (because of prokka run) is not empty,
    it writes the right message.
    """
    logfile_base = "test_prokka"
    utils.init_logger(logfile_base, 0, '', verbose=1)
    skipped = ["toto", "genome", "genome2"]
    utils.write_warning_skipped(skipped)
    out, err = capsys.readouterr()
    assert ("Prokka had problems while annotating some genomes, or did not "
            "find any gene. Hence, they are not "
            "formatted, and absent from your output database. Please look at their "
            "Prokka logs (<output_directory>/tmp_files/<genome_name>-prokka.log) and "
            "to the current error log (<output_directory>/<input_filename>.log.err)"
            " to get more information, and run again to annotate and format them. "
            "Here are the genomes (problem with prokka or no "
            "gene found):") in err
    assert ("\\n\\t- toto\\n\\t- genome\\n\\t- genome2" in err or
            "\n\t- toto\n\t- genome\n\t- genome2" in err)
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_skipped_format(capsys):
    """
    Test that when the list of skipped genomes (format step could not run) is not empty,
    it writes the right message.
    """
    logfile_base = "test_prokka"
    utils.init_logger(logfile_base, 0, '', verbose=1)
    skipped_format = ["toto", "genome", "genome2"]
    utils.write_warning_skipped(skipped_format, do_format=True)
    out, err = capsys.readouterr()
    assert ("Some genomes were annotated by prokka, but could not be formatted, "
            "and are hence absent from your output database. Please look at log "
            "files to get more information about why they could not be ") in err
    assert ("formatted.\n\t- toto\n\t- genome\n\t- genome2\n" in err or
            "formatted.\\n\\t- toto\\n\\t- genome\\n\\t- genome2" in err)
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_write_discarded():
    """
    Test that the list of discarded genomes is written as expected.
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], 9876546, 6, 2]
               }
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3"]
    list_file = os.path.join("test", "data", "annotate", "list_genomes.txt")
    res_path = os.path.join("test", "data", "annotate")
    utils.write_discarded(genomes, kept_genomes, list_file, res_path)
    outfile = os.path.join("test", "data", "annotate", "discarded-list_genomes.lst")
    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_test_write_discard.lst")
    # There is no order in the discarded file. So, just check that the lines
    # written are as expected.
    with open(outfile, "r") as outf, open(exp_file, "r") as expf:
        exp_lines = expf.readlines()
        out_lines = outf.readlines()
    assert set(out_lines) == set(exp_lines)
    os.remove(outfile)


def test_write_discarded_qc():
    """
    Test that the list with information on all genomes when we run with QC only is
    written as expected
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], 9876546, 6, 2]
               }
    kept_genomes = []
    list_file = os.path.join("test", "data", "annotate", "list_genomes.txt")
    res_path = os.path.join("test", "data", "annotate")
    utils.write_discarded(genomes, kept_genomes, list_file, res_path, qc=True)
    outfile = os.path.join("test", "data", "annotate", "info-genomes-list_genomes.lst")
    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_test_write_info_qc.lst")
    # There is no order in the discarded file. So, just check that the lines
    # written are as expected.
    with open(outfile, "r") as outf, open(exp_file, "r") as expf:
        exp_lines = expf.readlines()
        out_lines = outf.readlines()
    assert set(out_lines) == set(exp_lines)
    os.remove(outfile)


def test_write_discarded_empty():
    """
    Test that when the list of genomes is empty, but the list of kept-genomes is not
    (should never happen...), it writes only the header of discarded lst file.
    """
    genomes = {}
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3"]
    list_file = os.path.join("test", "data", "annotate", "list_genomes.txt")
    res_path = os.path.join("test", "data", "annotate")
    utils.write_discarded(genomes, kept_genomes, list_file, res_path)
    outfile = os.path.join("test", "data", "annotate", "discarded-list_genomes.lst")
    with open(outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "orig_name\tgsize\tnb_conts\tL90\n"
    os.remove(outfile)


def test_write_discarded_all_kept():
    """
    Test that when all genomes are kept, the discarded lst file only contains the
    header line.
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], 9876546, 6, 2]
               }
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3", "genome2", "genome1"]
    list_file = os.path.join("test", "data", "annotate", "list_genomes.txt")
    res_path = os.path.join("test", "data", "annotate")
    utils.write_discarded(genomes, kept_genomes, list_file, res_path)
    outfile = os.path.join("test", "data", "annotate", "discarded-list_genomes.lst")
    with open(outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "orig_name\tgsize\tnb_conts\tL90\n"
    os.remove(outfile)


def test_write_lstinfo():
    """
    Test that lstinfo file is written as expected.
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417.00010", gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417.00006", gpaths[1], 456464645, 5, 1],
               gnames[2]: ["genome.0417.00008", gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417.00008", gpaths[3], 6549, 16, 8],
               gnames[4]: ["genome.0417.00001", gpaths[4], 9876546, 6, 2]
               }
    list_file = os.path.join("test", "data", "annotate", "list_genomes.txt")
    outdir = os.path.join("test", "data", "annotate")
    utils.write_lstinfo(list_file, genomes, outdir)
    outfile = os.path.join(outdir, "LSTINFO-list_genomes.lst")
    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_test_write_lstinfo.lst")
    with open(outfile, "r") as outf, open(exp_file, "r") as expf:
        for line_out, line_exp in zip(outf, expf):
            assert line_out == line_exp
    os.remove(outfile)


def test_write_lstinfo_nogenome():
    """
    Test that when there is no genome fully annotated, lstinfo contains
    only header.
    """
    genomes = {}
    list_file = os.path.join("test", "data", "annotate", "list_genomes.txt")
    outdir = os.path.join("test", "data", "annotate")
    utils.write_lstinfo(list_file, genomes, outdir)
    outfile = os.path.join(outdir, "LSTINFO-list_genomes.lst")
    with open(outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "gembase_name\torig_name\tgsize\tnb_conts\tL90\n"
    os.remove(outfile)


def test_sort_gene():
    """
    Test that genomes are sorted by species first, and then by strain number.
    """
    # genomes = {genome_orig, [gembase, path, gsize, nbcont, L90]}
    genomes = ["genome.0417.00010", "toto.0417.00010", "genome1.0417.00002",
               "genome.0417.00015", "totn.0417.00010", "genome.0417.00009",
               "genome.0517.00001", "toto.0417.00011"]
    sorted_genomes = sorted(genomes, key=utils.sort_genomes)
    exp = ["genome.0517.00001", "genome.0417.00009", "genome.0417.00010",
           "genome.0417.00015", "genome1.0417.00002", "totn.0417.00010",
           "toto.0417.00010", "toto.0417.00011"]
    assert sorted_genomes == exp


def test_sort_gene_tuple():
    """
    Test that genomes are sorted by species first, and then by strain number.
    """
    # genomes = {genome_orig, [gembase, path, gsize, nbcont, L90]}
    genomes = {"name1": ["genome.0417.00010", "path/to/genome", 123456, 50, 3],
               "name2": ["toto.0417.00010", "path/to/genome", 123456, 50, 3],
               "name3": ["genome1.0417.00002", "path/to/genome", 123456, 50, 3],
               "name4": ["genome.0417.00015", "path/to/genome", 123456, 50, 3],
               "name5": ["totn.0417.00010", "path/to/genome", 123456, 50, 3],
               "name6": ["genome.0417.00009", "path/to/genome", 123456, 50, 3],
               "name7": ["genome.0517.00001", "path/to/genome", 123456, 50, 3],
               "name8": ["toto.0417.00011", "path/to/genome", 123456, 50, 3], }
    sorted_genomes = sorted(genomes.items(), key=utils.sort_genomes)
    exp = [("name7", genomes["name7"]), ("name6", genomes["name6"]),
           ("name1", genomes["name1"]), ("name4", genomes["name4"]),
           ("name3", genomes["name3"]), ("name5", genomes["name5"]),
           ("name2", genomes["name2"]), ("name8", genomes["name8"])]
    assert sorted_genomes == exp


def test_sort_gene_noformat():
    """
    Test that when genomes are not in the gembase format, it returns
    genomes in alphabetical order
    """
    # genomes = {genome_orig, [gembase, path, gsize, nbcont, L90]}
    genomes = {"name1": ["genome.0417.00010", "path/to/genome", 123456, 50, 3],
               "name2": ["toto.0417.00010", "path/to/genome", 123456, 50, 3],
               "name3": ["genome1.0417.00002", "path/to/genome", 123456, 50, 3],
               "name4": ["genome.0417.00015", "path/to/genome", 123456, 50, 3],
               "name5": ["mygenome.0416", "path/to/genome", 123456, 50, 3],
               "name6": ["genome.0417", "path/to/genome", 123456, 50, 3],
               "name7": ["genome.0517.00001", "path/to/genome", 123456, 50, 3],
               "name8": ["toto.0417.00011", "path/to/genome", 123456, 50, 3], }
    sorted_genomes = sorted(genomes.items(), key=utils.sort_genomes)
    exp = [("name7", genomes["name7"]), ("name1", genomes["name1"]),
           ("name4", genomes["name4"]), ("name6", genomes["name6"]),
           ("name3", genomes["name3"]), ("name5", genomes["name5"]),
           ("name2", genomes["name2"]), ("name8", genomes["name8"])]
    assert sorted_genomes == exp


def test_sort_proteins():
    """
    Test that when a list of proteins is given, it returns the list sorted by :
    - species
    - strain number
    - protein number
    """
    proteins = ["ESCO.0417.00010.i0001_12354", "ESCA.0617.00001.i0001_005",
                "ESCO.0517.00001.i0001_12354", "ESCO.0417.00010.i0001_1",
                "TOTO.0417.00001.i0001_12", "TOTO.0717.00001.i0008_9"]
    sorted_prot = sorted(proteins, key=utils.sort_proteins)
    exp = ["ESCA.0617.00001.i0001_005", "ESCO.0517.00001.i0001_12354",
           "ESCO.0417.00010.i0001_1", "ESCO.0417.00010.i0001_12354",
           "TOTO.0717.00001.i0008_9", "TOTO.0417.00001.i0001_12"]
    assert sorted_prot == exp


def test_sort_proteins_other_format():
    """
    Test that when a list of proteins is given, it returns the list sorted by :
    - species
    - strain number
    - protein number
    for all proteins where this information is available,
    and then, the other proteins, sorted by :
    - species (string before "_")
    - protein number (num after "_")
    """
    proteins = ["ESCO.0417.00010.i0001_12354", "ESCO.0617.00001.i0001_005", "my_protein_0002",
                "ESCO.0517.00001.i0001_12354", "my_protein_0001", "ESCO.0417.00010.i0001_1",
                "a_my_prot_15", "a_my_prot_2",
                "TOTO.0417.00001.i0001_12", "TOTO.0717.00001.i0008_9"]
    sorted_prot = sorted(proteins, key=utils.sort_proteins)
    exp = ["ESCO.0617.00001.i0001_005", "ESCO.0517.00001.i0001_12354",
           "ESCO.0417.00010.i0001_1", "ESCO.0417.00010.i0001_12354",
           "TOTO.0717.00001.i0008_9", "TOTO.0417.00001.i0001_12",
           "a_my_prot_2", "a_my_prot_15", "my_protein_0001", "my_protein_0002"]
    assert sorted_prot == exp


def test_sort_proteins_error_format1(capsys):
    """
    Test that when a protein name does not follow the format <alpha_num>_<num>,
    it gives an error.
    """
    proteins = ["ESCO.0417.00010.i0001_12354", "ESCO.0617.00001.i0001_005",
                "ESCO.0517.00001.i0001_12354", "error-protein"]
    with pytest.raises(SystemExit):
        sorted(proteins, key=utils.sort_proteins)
    out, err = capsys.readouterr()
    assert ("ERROR: Protein error-protein does not have the required format. It must contain, "
            "at least <alpha-num>_<num_only>, and at best "
            "<name>.<date>.<strain_num>.<contig_info>_<prot_num>. "
            "Please change its name.") in err


def test_read_genomes_nofile(capsys):
    """
    Test that when the genome list file provided does not exist, it
    ends the program with an error message
    """
    with pytest.raises(SystemExit):
        utils.read_genomes("toto.txt", "TOTO", "0417", "db/path", "tmppath")
    out, err = capsys.readouterr()
    assert "ERROR: Your list file " in err
    assert "toto.txt" in err
    assert "does not exist. Please provide a list file." in err
    assert "Ending program." in err


def test_read_genomes_wrongname():
    """
    Test that when the list file contains only genome names which do not exist,
    it returns an empty list of genomes to annotate/format.
    """
    name = "ESCO"
    date = "0417"
    dbpath = os.path.join("test", "data", "annotate", "genomes")
    tmppath = "tmppath"
    list_file = os.path.join("test", "data", "annotate", "test_files",
                             "list_genomes-wrongNames.txt")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath)
    assert genomes == {}


def test_read_genomes_ok(capsys):
    """
    Test that when the list file contains genomes existing, it returns the expected list
    of genomes
    """
    logfile_base = "test_utils"
    utils.init_logger(logfile_base, 0, '', verbose=1)
    name = "ESCO"
    date = "0417"
    dbpath = os.path.join("test", "data", "annotate", "genomes")
    tmppath = "tmppath"
    list_file = os.path.join("test", "data", "annotate", "test_files", "list_genomes.lst")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath)
    exp = {"A_H738.fasta": ["ESCO.0417"],
           "B2_A3_5.fasta-split5N.fna-short-contig.fna": ["ESCO.0417"],
           "H299_H561.fasta": ["ABCD.0417"], "genome2.fasta": ["TOTO.0417"],
           "genome3.fasta": ["ESCO.0512"], "genome4.fasta": ["TOTO.0417"],
           "genome5.fasta": ["TOTO.0114"]}
    assert exp == genomes
    _, err = capsys.readouterr()
    assert "genome.fst genome file does not exist. It will be ignored." in err
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_read_genomes_errors(capsys):
    """
    Test that when the list file contains errors in name and date provided,
    it returns the expected errors, and the expected genome list.
    """
    logfile_base = "test_utils"
    utils.init_logger(logfile_base, 0, '', verbose=1)
    name = "ESCO"
    date = "0417"
    dbpath = os.path.join("test", "data", "annotate", "genomes")
    tmppath = "tmppath"
    list_file = os.path.join("test", "data", "annotate", "test_files", "list_genomes-errors.txt")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath)
    exp = {"A_H738.fasta": ["ESCO.0417"],
           "B2_A3_5.fasta-split5N.fna-short-contig.fna": ["ESCO.0417"],
           "H299_H561.fasta": ["ABCD.0417"], "genome2.fasta": ["TOTO.0417"],
           "genome3.fasta": ["ESCO.0512"], "genome4.fasta": ["ESCO.0417"],
           "genome5.fasta": ["ESCO.0417"]}
    assert genomes == exp
    _, err = capsys.readouterr()
    assert ("Invalid name/date given for genome A_H738.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default name (ESCO) and date (0417) will "
            "be used.") in err
    assert (
           "Invalid name abc given for genome B2_A3_5.fasta-split5N.fna-short-contig.fna. Only put "
           "4 alphanumeric characters in your date and name. For "
           "this genome, the default name (ESCO) will "
           "be used.") in err
    assert ("Invalid date 152 given for genome H299_H561.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default date (0417) will "
            "be used.") in err
    assert ("Invalid date 1-03 given for genome genome2.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default date (0417) will "
            "be used.") in err
    assert ("genome.fst genome file does not exist. "
            "It will be ignored.") in err
    assert ("Invalid name a/b2 given for genome genome3.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default name (ESCO) will "
            "be used.") in err
    assert ("Invalid name #esc given for genome genome5.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default name (ESCO) will "
            "be used.") in err
    assert ("Invalid date 1_14 given for genome genome5.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default date (0417) will "
            "be used.") in err
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_read_genomes_multi_files(capsys):
    """
    Test that when the list file contains several filenames for 1 same genome,
    it returns the expected genome list, the expected errors (when some genome
    files do not exist) and the expected concatenated files.
    """
    logfile_base = "test_utils"
    utils.init_logger(logfile_base, 0, '', verbose=1)
    name = "ESCO"
    date = "0417"
    dbpath = os.path.join("test", "data", "annotate", "genomes")
    tmppath = os.path.join("test", "data", "annotate")
    list_file = os.path.join("test", "data", "annotate", "test_files",
                             "list_genomes-multi-files.txt")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath)
    exp = {"A_H738.fasta-all.fna": ["ESCO.0417"],
           "H299_H561.fasta-all.fna": ["ABCD.0417"], "genome2.fasta": ["TOTO.0417"],
           "genome3.fasta": ["ESCO.0512"], "genome4.fasta": ["TOTO.0417"],
           "genome5.fasta": ["TOTO.0114"]}
    assert exp == genomes
    # Check error messages
    _, err = capsys.readouterr()
    assert ("genome.fna genome file does not exist. Its file will be ignored "
            "when concatenating ['H299_H561.fasta', 'genome6.fasta', 'genome.fna']") in err
    assert "genome.fst genome file does not exist. It will be ignored." in err
    assert ("toto.fst genome file does not exist. Its file will be ignored "
            "when concatenating ['toto.fst', 'toto.fasta', 'genome.fst']") in err
    assert ("toto.fasta genome file does not exist. Its file will be ignored "
            "when concatenating ['toto.fst', 'toto.fasta', 'genome.fst']") in err
    assert ("genome.fst genome file does not exist. Its file will be ignored "
            "when concatenating ['toto.fst', 'toto.fasta', 'genome.fst']") in err
    assert ("None of the genome files in ['toto.fst', 'toto.fasta', 'genome.fst'] exist. "
            "This genome will be ignored.") in err
    # Check that files were concatenated as expected
    concat1 = os.path.join(tmppath, "A_H738.fasta-all.fna")
    exp_concat1 = os.path.join(dbpath, "A_H738-and-B2_A3_5.fna")
    concat2 = os.path.join(tmppath, "H299_H561.fasta-all.fna")
    exp_concat2 = os.path.join(dbpath, "H299_H561-and-genome6.fna")
    with open(concat1, "r") as outf, open(exp_concat1, "r") as expf:
        for line_out, line_exp in zip(outf, expf):
            assert line_out == line_exp
    with open(concat2, "r") as outf, open(exp_concat2, "r") as expf:
        for line_out, line_exp in zip(outf, expf):
            assert line_out == line_exp
    os.remove(concat1)
    os.remove(concat2)
    os.remove(logfile_base + ".log")
    os.remove(logfile_base + ".log.details")
    os.remove(logfile_base + ".log.err")


def test_check_resdirlst(capsys):
    """
    Test that when the result directory already contains .lst files in LSTINFO,
    program ends with an error message.
    """
    resdir = os.path.join("test", "data", "annotate", "test_check_resdir")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "LSTINFO"))
    open(os.path.join(resdir, "LSTINFO", "toto.lst"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(resdir)
    out, err = capsys.readouterr()
    assert ("ERROR: Your output directory already has .lst files in the "
            "LSTINFO folder. Provide another result directory, or remove the "
            "files in this one.") in err
    shutil.rmtree(resdir)


def test_check_resdirprt(capsys):
    """
    Test that when the result directory already contains .prt files in Proteins,
    program ends with an error message.
    """
    resdir = os.path.join("test", "data", "annotate", "test_check_resdir")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "Proteins"))
    open(os.path.join(resdir, "Proteins", "toto.prt"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(resdir)
    out, err = capsys.readouterr()
    assert ("ERROR: Your output directory already has .prt files in the "
            "Proteins folder. Provide another result directory, or remove the "
            "files in this one.") in err
    shutil.rmtree(resdir)


def test_check_resdirgen(capsys):
    """
    Test that when the result directory already contains .gen files in Genes,
    program ends with an error message.
    """
    resdir = os.path.join("test", "data", "annotate", "test_check_resdir")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "Genes"))
    open(os.path.join(resdir, "Genes", "toto.gen"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(resdir)
    out, err = capsys.readouterr()
    assert ("ERROR: Your output directory already has .gen files in the "
            "Genes folder. Provide another result directory, or remove the "
            "files in this one.") in err
    shutil.rmtree(resdir)


def test_check_resdirrep(capsys):
    """
    Test that when the result directory already contains .fna files in Replicons,
    program ends with an error message.
    """
    resdir = os.path.join("test", "data", "annotate", "test_check_resdir")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "Replicons"))
    open(os.path.join(resdir, "Replicons", "toto.fna"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(resdir)
    out, err = capsys.readouterr()
    assert ("ERROR: Your output directory already has .fna files in the "
            "Replicons folder. Provide another result directory, or remove the "
            "files in this one.") in err
    shutil.rmtree(resdir)


def test_check_resdirgff(capsys):
    """
    Test that when the result directory already contains .fna files in Replicons,
    program ends with an error message.
    """
    resdir = os.path.join("test", "data", "annotate", "test_check_resdir_gff")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "gff3"))
    open(os.path.join(resdir, "gff3", "toto.gff"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(resdir)
    out, err = capsys.readouterr()
    assert ("ERROR: Your output directory already has .gff files in the "
            "gff3 folder. Provide another result directory, or remove the "
            "files in this one.") in err
    shutil.rmtree(resdir)


def test_check_resdirotherext():
    """
    Test that when the result directory contains txt files in Replicons dir,
    there is no problem.
    """
    resdir = os.path.join("test", "data", "annotate", "test_check_resdir")
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(resdir, "Replicons"))
    open(os.path.join(resdir, "Replicons", "toto.txt"), "w").close()
    utils.check_out_dirs(resdir)
    shutil.rmtree(resdir)


def test_check_resdirnodir():
    """
    Test that when the result directory does not already exist, there is no problem.
    """
    utils.check_out_dirs("totoresdir")


def test_run_cmd_error_noquit(capsys):
    """
    Test that when we try to run a command which does not exist, it returns an error message,
    but does not exit the program (eof=False).
    """
    cmd = "toto"
    error = "error trying to run toto"
    assert utils.run_cmd(cmd, error) == -1
    out, err = capsys.readouterr()
    assert "error trying to run toto: toto does not exist" in err


def test_run_cmd_error_noquit_logger(capsys):
    """
    Test that when we try to run a command which does not exist, it returns an error message,
    but does not exit the program (eof=False). With a given logger where error is written.
    """
    cmd = "toto"
    logger = logging.getLogger("default")
    error = "error trying to run toto"
    assert utils.run_cmd(cmd, error, logger=logger) == -1
    out, err = capsys.readouterr()
    assert "error trying to run toto: toto does not exist" in err


def test_run_cmd_error_quit(capsys):
    """
    Test that when we try to run a command which does not exist, it returns an error message,
    and exits the program (eof=True)
    """
    cmd = "toto"
    error = "error trying to run toto"
    with pytest.raises(SystemExit):
        utils.run_cmd(cmd, error, eof=True)
    out, err = capsys.readouterr()
    assert "error trying to run toto: toto does not exist" in err


def test_run_cmd_retcode_non0(caplog):
    """
    Test that when the command fails, it returns a non-zero int, writes an error message,
    but does not quit the program (eof=False by default).
    """
    cmd = "prokka"
    error = "error trying to run prodigal"
    assert utils.run_cmd(cmd, error).returncode != 0
    assert error in caplog.text


def test_run_cmd_retcode_non0_quit(caplog):
    """
    Test that when the command fails, it returns a non-zero int, writes an error message,
    and exits the program (eof=True).
    """
    cmd = "prokka"
    error = "error trying to run prodigal"
    with pytest.raises(SystemExit):
        utils.run_cmd(cmd, error, eof=True)
    assert error in caplog.text


def test_run_cmd_error_stderrfile():
    """
    Test that when we try to run a command which does not exist, and direct its output to
    a file instead of stderr, we have the expected error written in the given error file.
    """
    cmd = "toto"
    error = "error trying to run toto"
    outfile = open(os.path.join("test", "data", "annotate", "stderr_run_cmd.txt"), "w")
    assert utils.run_cmd(cmd, error, stderr=outfile) == -1
    outfile.close()
    os.remove(os.path.join("test", "data", "annotate", "stderr_run_cmd.txt"))


def test_run_cmd_error_stdoutfile():
    """
    Test that when we try to run a command which does not exist, and direct its output to
    a file instead of stderr, we have the expected error written in the given error file.
    """
    cmd = "toto"
    error = "error trying to run toto"
    outfile = open(os.path.join("test", "data", "annotate", "stdout_run_cmd.txt"), "w")
    assert utils.run_cmd(cmd, error, stdout=outfile) == -1
    outfile.close()
    os.remove(os.path.join("test", "data", "annotate", "stdout_run_cmd.txt"))


def test_rename_contigs():
    """
    From a given sequence, rename all its contigs with the given gembase name + a number,
    and save the output sequence to the given res_path.
    Check that the output file is as expected.
    """
    gpath = os.path.join("test", "data", "annotate", "genomes", "H299_H561.fasta")
    gembase_name = "ESCO.0216.00005"
    res_path = os.path.join("test", "data", "annotate")
    outfile = os.path.join(res_path, "H299_H561.fasta-short-contig.fna")
    exp_file = os.path.join("test", "data", "annotate", "exp_files", "res_H299_H561-ESCO00005.fna")
    utils.rename_genome_contigs(gembase_name, gpath, outfile)
    with open(exp_file, "r") as expf, open(outfile, "r") as of:
        for line_exp, line_seq in zip(expf, of):
            assert line_exp == line_seq
    os.remove(outfile)


def test_cat_nobar(capsys):
    """
    Check that when cat is called on a list of several files, the output file
    contains what is expected (concatenation of content of all input files)
    """
    import glob
    list_files = glob.glob(os.path.join("test", "data", "annotate", "genomes", "*.fasta"))
    outfile = "test_catfile.txt"
    utils.cat(list_files, outfile)
    exp_file = os.path.join("test", "data", "annotate", "exp_files",
                            "res_test_cat_genomes_fasta.fst")
    with open(exp_file, 'r') as expf, open(outfile, 'r') as outf:
        lines_exp = expf.readlines()
        lines_out = outf.readlines()
    assert set(lines_exp) == set(lines_out)
    _, err = capsys.readouterr()
    assert "/{} (".format(len(list_files)) not in err
    os.remove(outfile)


def test_cat_bar():
    """
    Check that when cat is called on a list of several files, the output file
    contains what is expected (concatenation of content of all input files)
    """
    import glob
    list_files = glob.glob(os.path.join("test", "data", "annotate", "genomes", "*.fasta"))
    outfile = "test_catfile.txt"
    title = "test cat progressbar"
    utils.cat(list_files, outfile, title=title)
    exp_file = os.path.join("test", "data", "annotate", "exp_files",
                            "res_test_cat_genomes_fasta.fst")
    with open(exp_file, 'r') as expf, open(outfile, 'r') as outf:
        lines_exp = expf.readlines()
        lines_out = outf.readlines()
    assert set(lines_exp) == set(lines_out)
    os.remove(outfile)


def test_detail():
    """
    Check value returned for detail level
    """
    assert utils.detail_lvl() == 15


def test_grep():
    """
    Check that it returns all lines containing the given pattern
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = "filein.txt"
    with open(filein, "w") as ff:
        ff.write(lines)
    pattern = "[0-9]+toto\.txt [a-z]{6}"
    lines_grep = utils.grep(filein, pattern)
    exp = ["1toto.txt otherletters",
           "123toto.txt otherletters"]
    assert exp == lines_grep
    os.remove(filein)


def test_grep_count():
    """
    Check that it returns the number of lines containing the given pattern
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = "filein.txt"
    with open(filein, "w") as ff:
        ff.write(lines)
    pattern = "[0-9]+toto\.txt [a-z]{6}"
    lines_grep = utils.grep(filein, pattern, counts=True)
    assert lines_grep == 2
    os.remove(filein)


def test_count_lines():
    """
    Test that it returns the number of lines in given file
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = "filein.txt"
    with open(filein, "w") as ff:
        ff.write(lines)
    nbline = utils.count(filein)
    assert nbline == 5
    os.remove(filein)


def test_count_words():
    """
    Test that it returns the number of words in given file
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = "filein.txt"
    with open(filein, "w") as ff:
        ff.write(lines)
    nbword = utils.count(filein, get="words")
    assert nbword == 10
    os.remove(filein)


def test_count_error(capsys):
    """
    test that when we want to count something else than 'lines' or 'words', it returns an error
    """
    with pytest.raises(SystemExit):
        utils.count("filein", get="letters")
    _, err = capsys.readouterr()
    assert "Choose what you want to count among ['lines', 'words']" in err


def test_save_bin():
    """
    Test that python objects are correctly saved into binary file
    """
    obj1 = {1: "toto", "key": 455, "key2": "val"}
    obj2 = [1, 2, 5, "toto", "plop"]
    obj3 = "a string"
    objects = [obj1, obj2, obj3]
    fileout = "filout.bin"
    utils.save_bin(objects, fileout)
    import _pickle as pickle
    with open(fileout, "rb") as binf:
        obj = pickle.load(binf)
    assert objects == obj
    os.remove(fileout)


def test_load_bin():
    """
    Test that python objects are correctly recovered from binary file
    """
    fileout = os.path.join("test", "data", "utils", "objects_saved.bin")
    found = utils.load_bin(fileout)
    obj1 = {1: "toto", "key": 455, "key2": "val"}
    obj2 = [1, 2, 5, "toto", "plop"]
    obj3 = "a string"
    objects = [obj1, obj2, obj3]
    assert found == objects


def test_write_list():
    """
    Test that the given list is returned as expected string
    """
    inlist = [1, 2, "toto", {1: 2, "1": 5}, 1e-6]
    tostr = utils.write_list(inlist)
    exp1 = "1 2 toto {1: 2, '1': 5} 1e-06\n"
    exp2 = "1 2 toto {'1': 5, 1: 2} 1e-06\n"
    assert tostr == exp1 or tostr == exp2


def test_remove_exits():
    """
    Test that a given file is removed if it exists
    """
    infile = "toto.txt"
    open(infile, "w").close()
    assert os.path.isfile(infile)
    utils.remove(infile)
    assert not os.path.isfile(infile)


def test_remove_not_exist():
    """
    Test that removing a file which does not exist brings no error
    """
    infile = "toto.txt"
    assert not os.path.isfile(infile)
    utils.remove(infile)
    assert not os.path.isfile(infile)
