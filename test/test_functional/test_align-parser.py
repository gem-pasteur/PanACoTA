#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of align subcommand
"""
import argparse
import pytest

from PanACoTA.subcommands import align


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Align families", add_help=False)
    align.build_parser(parser)
    with pytest.raises(SystemExit):
        align.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert "usage: " in err
    assert "-c COREPERS -l LIST_GENOMES -n DATASET_NAME -d DBPATH" in err
    assert "-o OUTDIR" in err
    assert "[--threads THREADS] [-F] [-v] [-q] [-h]" in err
    assert "[-h]" in err
    assert "the following arguments are required: -c, -l, -n, -d, -o" in err


def test_parser_thread_notint(capsys):
    """
    Test that when the number of threads given is not an int, it returns an error message
    """
    parser = argparse.ArgumentParser(description="Align families", add_help=False)
    align.build_parser(parser)
    with pytest.raises(SystemExit):
        align.parse(parser, "-l listgenome -n dname -d dbpath -o outdir --threads 1.5".split())
    _, err = capsys.readouterr()
    assert "argument --threads threads: invalid int value: 1.5" in err


def test_parser_thread_toomany(capsys):
    """
    Test that when the number of threads given is higher than the total number of threads,
    it returns the expected error message
    """
    import multiprocessing
    nb = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Align families", add_help=False)
    align.build_parser(parser)
    with pytest.raises(SystemExit):
        align.parse(parser,
                    "-l listgenome -n dname -d dbpath -o outdir "
                    "--threads {}".format(nb + 3).split())
    _, err = capsys.readouterr()
    assert ("You have {} threads on your computer, you cannot ask for more: "
            "invalid value: {}".format(nb, nb+3)) in err


def test_parser_thread_neg(capsys):
    """
    Test that when the number of threads given is a negative number, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Align families", add_help=False)
    align.build_parser(parser)
    with pytest.raises(SystemExit):
        align.parse(parser, "-l listgenome -n dname -d dbpath -o outdir --threads -5".split())
    _, err = capsys.readouterr()
    assert ("Please provide a positive number of threads (or 0 for all threads): "
            "Invalid value: -5") in err


def test_parser_default():
    """
    Test that when run with the minimum required arguments, all default values are
    as expected.
    """
    parser = argparse.ArgumentParser(description="Align families", add_help=False)
    align.build_parser(parser)
    options = align.parse(parser, "-c cp -l listgenome -n dname -d dbpath -o outdir".split())
    assert options.corepers == "cp"
    assert options.list_genomes == "listgenome"
    assert options.dataset_name == "dname"
    assert options.dbpath == "dbpath"
    assert options.outdir == "outdir"
    assert options.threads == 1
    assert options.force is False
    assert options.verbose == 0
    assert options.quiet is False


def test_parser_allthreads():
    """
    Test that when run with 0 for --threads option, it returns the total number of threads in
    computer
    """
    import multiprocessing
    nb = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Align families", add_help=False)
    align.build_parser(parser)
    options = align.parse(parser, "-c cp -l listgenome -n dname -d dbpath -o outdir "
                                  "--threads 0".split())
    assert options.corepers == "cp"
    assert options.list_genomes == "listgenome"
    assert options.dataset_name == "dname"
    assert options.dbpath == "dbpath"
    assert options.outdir == "outdir"
    assert options.threads == nb
    assert options.force is False
    assert options.verbose == 0
    assert options.quiet is False


def test_parser_all_ok():
    """
    Test that when all arguments given are ok, it builds expected parser
    """
    parser = argparse.ArgumentParser(description="Align families", add_help=False)
    align.build_parser(parser)
    options = align.parse(parser, "-c cp -l listgenome -n dname -d dbpath -o outdir "
                                  "--threads 1".split())
    assert options.corepers == "cp"
    assert options.list_genomes == "listgenome"
    assert options.dataset_name == "dname"
    assert options.dbpath == "dbpath"
    assert options.outdir == "outdir"
    assert options.threads == 1
    assert options.force is False
    assert options.verbose == 0
    assert options.quiet is False
