#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of annote_pipeline.py
"""
import pytest
import argparse
import time

from genomeAPCAT.subcommands import annote as annot


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert "-d DB_PATH -r RES_PATH [-n NAME] [-Q] [--l90 L90]" in err
    assert "[--nbcont NBCONT] [--cutN CUTN] [--date DATE] [--tmp TMPDIR]" in err
    assert "[--prok PROKKADIR] [-F] [--threads THREADS] [-v] [-q] [-h]" in err
    assert "list_file" in err
    assert "the following arguments are required: list_file, -d, -r" in err


def test_parser_noname(capsys):
    """
    Test that when the script is called without any name for the genomes not -Q option,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "listfile -d dbpath -r respath".split())
    _, err = capsys.readouterr()
    assert ("You must specify your genomes dataset name in 4 characters with "
            "'-n name' option (type -h for more information). Or, if you do not want "
            "to annotate and format your genomes but just to run quality control, use "
            "option '-Q") in err


def test_parser_wrongname(capsys):
    """
    Test that when the script is called with a genome name with more than 4 characters,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n genome".split())
    _, err = capsys.readouterr()
    assert ("The genome name must contain 4 characters. For example, this name can "
            "correspond to the 2 first letters of genus, and 2 first letters of "
            "species, e.g. ESCO for Escherichia Coli.") in err


def test_parser_negative_cont(capsys):
    """
    Test that when the script is called with a limit of contig number higher than 9999,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --nbcont -5".split())
    _, err = capsys.readouterr()
    assert "The maximum number of contigs allowed must be a positive number." in err


def test_parser_high_cont(capsys):
    """
    Test that when the script is called with a negative limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --nbcont 10005".split())
    _, err = capsys.readouterr()
    assert "We do not support genomes with more than 9999 contigs." in err


def test_parser_wrong_cont(capsys):
    """
    Test that when the script is called with a non integer limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --nbcont 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --nbcont: invalid int value: 10.5" in err


def test_parser_wrongl90(capsys):
    """
    Test that when the user does not give an int for the l90 limit, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --l90 l90".split())
    _, err = capsys.readouterr()
    assert "argument --l90: invalid int value: 'l90'" in err


def test_parser_wrong_cut(capsys):
    """
    Test that when the user does not give an int for the cutN value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --cutN 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --cutN: invalid int value: '10.5'" in err


def test_parser_wrong_thread(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --threads 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --threads: invalid int value: '10.5'" in err


def test_parser_wrong_date(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 --date 417".split())
    _, err = capsys.readouterr()
    assert ("The date must contain 4 characters. Usually, it contains 4 digits, "
            "corresponding to the month (2 digits) and year (2 digits).") in err


def test_parser_q_and_v(capsys):
    """
    Test that when the user wants both quiet and verbose option, it gives an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 -q -v".split())
    _, err = capsys.readouterr()
    assert (("Choose between a verbose output (-v) or quiet output (-q)."
             " You cannot have both...")) in err


def test_parser_default():
    """
    Test that when run with the minimum required arguments, all default values are
    as expected.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, "list_file -d dbpath -r respath -n g123".split())
    assert options.list_file == "list_file"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "g123"
    assert options.l90 == 100
    assert options.nbcont == 999
    assert options.cutn == 5
    assert options.threads == 1
    assert options.date == time.strftime("%m%y")
    assert not options.force
    assert not options.qc_only


def test_parser_values():
    """
    Test that values for L90, nbcontig, cutn, threads, date are taken into account
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, ("list_file -d dbpath -r respath -n g123 --l90 2 "
                                   "--nbcont 10 --cutN 0 --threads 8 --date toto").split())
    assert options.list_file == "list_file"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "g123"
    assert options.l90 == 2
    assert options.nbcont == 10
    assert options.cutn == 0
    assert options.threads == 8
    assert options.date == "toto"
    assert not options.force
    assert not options.qc_only


def test_parser_force():
    """
    Test that when run with '-F' option, force is initialized to "--force".
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, "list_file -d dbpath -r respath -n g123 -F".split())
    assert options.list_file == "list_file"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "g123"
    assert options.l90 == 100
    assert options.nbcont == 999
    assert options.cutn == 5
    assert options.threads == 1
    assert options.date == time.strftime("%m%y")
    assert options.force
    assert not options.qc_only


def test_parser_wrongforce(capsys):
    """
    Test that when run with '-F' option + a value, it returns an error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "list_file -d dbpath -r respath -n g123 -F 10".split())
    _, err = capsys.readouterr()
    assert "unrecognized arguments: 10" in err


def test_parser_qc():
    """
    Test that when run with '-Q' option (for QC only) and no name given for the genome, it
    is set to "NONE"
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, "list_file -d dbpath -r respath -Q".split())
    assert options.list_file == "list_file"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "NONE"
    assert options.l90 == 100
    assert options.nbcont == 999
    assert options.cutn == 5
    assert options.threads == 1
    assert options.date == time.strftime("%m%y")
    assert not options.force
    assert options.qc_only
