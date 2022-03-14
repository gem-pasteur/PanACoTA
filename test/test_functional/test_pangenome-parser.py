#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of 'pangenome' subcommand
"""

import pytest
import argparse

from PanACoTA.subcommands import pangenome


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert "usage: " in err
    assert "-l LSTINFO_FILE -n DATASET_NAME -d DBPATH" in err
    assert "[-i MIN_ID]" in err
    assert " -o OUTDIR" in err
    assert "[-f OUTFILE] [-c {0,1,2}]" in err
    assert "[-s SPEDIR] [--threads THREADS] [-v]" in err
    assert "[-q] [-h]" in err
    assert "the following arguments are required: -l, -n, -d, -o" in err


def test_perc_id_nonum(capsys):
    """
    Test that when the given percentage of identity is not a number, it returns the
    expected error.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -i toto -o outdir".split())
    _, err = capsys.readouterr()
    assert "argument -i percentage_id: invalid float value: toto" in err


def test_neg_perc_id(capsys):
    """
    Test that when the given percentage of identity is a negative number, it returns the
    expected error.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -i -0.2 -o outdir".split())
    _, err = capsys.readouterr()
    assert "The minimum %% of identity must be in [0, 1]. Invalid value: -0.2" in err


def test_big_perc_id(capsys):
    """
    Test that when the given percentage of identity is greater than 1, it returns the
    expected error.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -i 2 -o outdir".split())
    _, err = capsys.readouterr()
    assert "The minimum %% of identity must be in [0, 1]. Invalid value: 2" in err


def test_wrong_clust_mode(capsys):
    """
    Test that when the given cluster mode does not exist, it returns the expected error message
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -i 0.8 -o od -c 3".split())
    _, err = capsys.readouterr()
    assert "argument -c: invalid choice: 3 (choose from 0, 1, 2)" in err


def test_thread_no_int(capsys):
    """
    Test that when the given number of threads is not an int, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -o od --threads 1.5".split())
    _, err = capsys.readouterr()
    assert "argument --threads threads: invalid int value: 1.5" in err


def test_thread_too_many(capsys):
    """
    Test that when the given number of threads is greater than the number of available threads,
    it returns the expected error message.
    """
    import multiprocessing
    nb = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, f"-l lstinfo -n TEST4 -d dbpath -o od --threads {nb*10}".split())
    _, err = capsys.readouterr()
    assert ("You have {} threads on your computer, you cannot ask for more: "
            "invalid value: {}").format(nb, nb*10) in err


def test_thread_neg(capsys):
    """
    Test that when the given number of threads is negative,
    it returns the expected error message.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -o od --threads -1".split())
    _, err = capsys.readouterr()
    assert "Please provide a positive number of threads (or 0 for all threads)" in err


def test_parser_default():
    """
    Test that when run with the minimum required arguments, all default values are
    as expected.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    options = pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -o od".split())
    assert options.lstinfo_file == "lstinfo"
    assert options.dataset_name == "TEST4"
    assert options.dbpath == "dbpath"
    assert options.min_id == 0.8
    assert options.outdir == "od"
    assert options.clust_mode == 1
    assert not options.spedir
    assert options.threads == 1
    assert not options.outfile
    assert options.verbose == 0
    assert not options.quiet


def test_parser_all_threads():
    """
    Test that when run with the minimum required arguments, + '0' for the thread number,
    all default values are as expected, and the number of threads is equal to the total
    number of threads.
    """
    import multiprocessing
    nb = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    options = pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -o od --threads 0".split())
    assert options.lstinfo_file == "lstinfo"
    assert options.dataset_name == "TEST4"
    assert options.dbpath == "dbpath"
    assert options.min_id == 0.8
    assert options.outdir == "od"
    assert options.clust_mode == 1
    assert not options.spedir
    assert options.threads == nb
    assert not options.outfile
    assert options.verbose == 0
    assert not options.quiet


def test_parser_1thread():
    """
    Test that when run with the minimum required arguments, + '0' for the thread number,
    all default values are as expected, and the number of threads is equal to the total
    number of threads.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    options = pangenome.parse(parser, "-l lstinfo -n TEST4 -d dbpath -o od --threads 1".split())
    assert options.lstinfo_file == "lstinfo"
    assert options.dataset_name == "TEST4"
    assert options.dbpath == "dbpath"
    assert options.min_id == 0.8
    assert options.outdir == "od"
    assert options.clust_mode == 1
    assert not options.spedir
    assert options.threads == 1
    assert not options.outfile
    assert options.verbose == 0
    assert not options.quiet
