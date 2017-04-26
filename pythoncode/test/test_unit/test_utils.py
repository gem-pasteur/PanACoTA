#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pipelinepackage.utils as utils
import pytest
import os
import logging
import test.test_unit.util_tests as util_tests


def test_check_install():
    """
    Try to run prokka, which is installed, and check that there is no problem
    """
    ret = utils.check_installed(["prokka"])
    assert ret == 1


def test_check_install_error(capsys):
    """
    Try to run a command which does not exist, and check that it closes the program
    with exit code 1
    """
    with pytest.raises(SystemExit):
        utils.check_installed(["plop false command..."])
    out, err = capsys.readouterr()
    assert "plop false command... failed:" in err


def test_class_filter():
    """
    Check that for a class LessThanFilter(warning), info and debug are allowed,
    but warning, error and critical are not.
    """
    a = utils.LessThanFilter(logging.WARNING)
    record = logging.LogRecord("root", logging.DEBUG, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.INFO, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.WARNING, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)
    record = logging.LogRecord("root", logging.ERROR, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)
    record = logging.LogRecord("root", logging.CRITICAL, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)


def test_plot_dist():
    """
    Plot a given distribution, and check that output is as expected
    """
    values = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 10]
    limit = 3
    res_dir = os.path.join("test", "data")
    os.makedirs(res_dir, exist_ok=True)
    outfile = os.path.join(res_dir, "distrib.png")
    reffile = os.path.join("test", "data", "exp_files", "res_plot_distr.png")
    title = "Distribution test"
    text = "Max L90 ="
    utils.plot_distr(values, limit, outfile, title, text)
    assert util_tests.compare_files(outfile, reffile)
    os.remove(outfile)


def test_skipped_prokka(capsys):
    """
    Test that when the list of skipped genomes (because of prokka run) is not empty,
    it writes the right message.
    """
    skipped = ["toto", "genome", "genome2"]
    utils.write_warning_skipped(skipped)
    out, err = capsys.readouterr()
    assert ("Prokka had problems while annotating some genomes. Hence, they are not "
            "formatted, and absent from your output database. Please look at their "
            "Prokka logs (<output_directory>/tmp_files/<genome_name>-prokka.log) and "
            "to the current error log (<output_directory>/<input_filename>.log.err)"
            " to get more information, and run again to annotate and format them. "
            "Here are the genomes:") in err
    assert ("\\n\\t- toto\\n\\t- genome\\n\\t- genome2" in err or
            "\n\t- toto\n\t- genome\n\t- genome2" in err)


def test_skipped_format(capsys):
    """
    Test that when the list of skipped genomes (format step could not run) is not empty,
    it writes the right message.
    """
    skipped_format = ["toto", "genome", "genome2"]
    utils.write_warning_skipped(skipped_format, format=True)
    out, err = capsys.readouterr()
    assert ("Some genomes were annotated by prokka, but could not be formatted, "
            "and are hence absent from your output database. Please look at log "
            "filesto get more information about wh they could not be ") in err
    assert ("formatted.\n\t- toto\n\t- genome\n\t- genome2\n" in err or
            "formatted.\\n\\t- toto\\n\\t- genome\\n\\t- genome2" in err)


def test_write_discarded():
    """
    Test that the list of discarded genomes is written as expected.
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join("test", "data", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], 9876546, 6, 2]
              }
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3"]
    list_file = os.path.join("test", "data", "list_genomes.txt")
    res_path = os.path.join("test", "data")
    utils.write_discarded(genomes, kept_genomes, list_file, res_path)
    outfile = os.path.join("test", "data", "discarded-list_genomes.lst")
    exp_file = os.path.join("test", "data", "exp_files", "res_test_write_discard.lst")
    # There is no order in the discarded file. So, just check that the lines
    # written are as expected.
    with open (outfile, "r") as outf, open(exp_file, "r") as expf:
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
    list_file = os.path.join("test", "data", "list_genomes.txt")
    res_path = os.path.join("test", "data")
    utils.write_discarded(genomes, kept_genomes, list_file, res_path)
    outfile = os.path.join("test", "data", "discarded-list_genomes.lst")
    with open (outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "orig_name\tgsize\tnb_conts\tL90\n"
    os.remove(outfile)


def test_write_discarded_allKept():
    """
    Test that when all genomes are kept, the discarded lst file only contains the
    header line.
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join("test", "data", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], 9876546, 6, 2]
              }
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3", "genome2", "genome1"]
    list_file = os.path.join("test", "data", "list_genomes.txt")
    res_path = os.path.join("test", "data")
    utils.write_discarded(genomes, kept_genomes, list_file, res_path)
    outfile = os.path.join("test", "data", "discarded-list_genomes.lst")
    with open (outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "orig_name\tgsize\tnb_conts\tL90\n"
    os.remove(outfile)


def test_write_lstinfo():
    """
    Test that lstinfo file is written as expected.
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join("test", "data", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417.00010", gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417.00006", gpaths[1], 456464645, 5, 1],
               gnames[2]: ["genome.0417.00008", gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417.00008", gpaths[3], 6549, 16, 8],
               gnames[4]: ["genome.0417.00001", gpaths[4], 9876546, 6, 2]
              }
    list_file = os.path.join("test", "data", "list_genomes.txt")
    outdir = os.path.join("test", "data")
    utils.write_lstinfo(list_file, genomes, outdir)
    outfile = os.path.join(outdir, "LSTINFO-list_genomes.lst")
    exp_file = os.path.join("test", "data", "exp_files", "res_test_write_lstinfo.lst")
    with open (outfile, "r") as outf, open(exp_file, "r") as expf:
        for line_out, line_exp in zip(outf, expf):
            assert line_out == line_exp
    os.remove(outfile)


def test_write_lstinfo_nogenome():
    """
    Test that when there is no genome fully annotated, lstinfo contains
    only header.
    """
    genomes = {}
    list_file = os.path.join("test", "data", "list_genomes.txt")
    outdir = os.path.join("test", "data")
    utils.write_lstinfo(list_file, genomes, outdir)
    outfile = os.path.join(outdir, "LSTINFO-list_genomes.lst")
    exp_file = os.path.join("test", "data", "exp_files", "res_test_write_lstinfo.lst")
    with open (outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "gembase_name\torig_name\tgsize\tnb_conts\tL90\n"
    os.remove(outfile)


def test_sort_gene():
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
               "name8": ["toto.0417.00011", "path/to/genome", 123456, 50, 3],}
    sorted_genomes = sorted(genomes.items(), key=utils.sort_genomes)
    exp = [("name7", genomes["name7"]), ("name6", genomes["name6"]),
           ("name1", genomes["name1"]), ("name4", genomes["name4"]),
           ("name3", genomes["name3"]), ("name5", genomes["name5"]),
           ("name2", genomes["name2"]), ("name8", genomes["name8"])]
    assert sorted_genomes == exp


def test_read_genomes():
    """

    """
    name = "ESCO"
    date = "0417"
    list_file = os.path.join("test", "data", "test_files", "list_genomes.lst")
    genomes = utils.read_genomes(list_file, name, date)
    exp = {"mygenome_in-fasta.fasta": ["ESCO.0417"], "othergenome.fa": ["ESCO.0417"],
           "yet-another-genome": ["ABCD.0417"], "andAnother1.fna": ["TOTO.0417"],
           "species-not-given": ["ESCO.0417"], "another-genome-from-TOTO_species": ["TOTO.0417"]}
    assert exp == genomes



