#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of 'all' subcommand
"""
import argparse
import pytest

from PanACoTA.subcommands import all_modules as allm


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    with pytest.raises(SystemExit):
        allm.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert "usage: " in err
    assert "-o OUTDIR [--threads THREADS]" in err
    assert "[-t NCBI_SPECIES_TAXID]" in err
    assert "-o OUTDIR" in err
    assert "[-s NCBI_SPECIES]" in err
    assert "[-l LEVELS]" in err
    assert "[--cutn CUTN] [--l90 L90]" in err
    assert "[--tol TOL]" in err
    assert "[-Mu] [-X]" in err
    assert "[--soft {fasttree,fastme,quicktree,iqtree,iqtree2}]" in err
    assert "the following arguments are required: -o, -n" in err


def test_parser_noconffile():
    """
    Test running with arguments ok for prepare, no configfile given.
    Only args given in line change, all others are by default
    """
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    options = allm.parse(parser, "-o out-all -n TEST -t 1234".split())
    assert options.outdir == "out-all"
    assert options.threads == 1
    assert options.ncbi_species_taxid == '1234'
    assert options.prodigal_only == False
    assert not options.norefseq
    assert options.l90 == 100
    assert options.clust_mode == 1
    assert options.min_id == 0.8
    assert options.tol == 1
    assert not options.boot


def test_parser_conffile():
    """
    Test running with arguments with config file.
    arguments given in cmd line are kept, others are those in configfile,
    and the ones not defined in configfile are by default
    """
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    options = allm.parse(parser,
                         "-c test/data/all/default-conffigfile.ini -o out-all -n TEST".split())
    assert options.outdir == "out-all"
    assert options.threads == 10
    assert not options.ncbi_species_taxid
    assert options.prodigal_only == False
    assert options.norefseq
    assert options.l90 == 99
    assert options.clust_mode == 1
    assert options.min_id == 0.85
    assert options.tol == 1
    assert not options.multi
    assert not options.boot


def test_parser_conffile_and_cmd():
    """
    Test that when some arguments given in config file are also given in cmd,
    the value kept is the one in cmd.
    """
    import multiprocessing
    nb = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    options = allm.parse(parser,
                         "-c test/data/all/default-conffigfile.ini -o out-all "
                         "-n TEST --threads 0 -i 0.99 -Mu".split())
    assert options.outdir == "out-all"
    assert options.threads == nb
    assert not options.ncbi_species_taxid
    assert options.prodigal_only == False
    assert options.norefseq
    assert options.l90 == 99
    assert options.clust_mode == 1
    assert options.min_id == 0.99
    assert options.tol == 1
    assert options.multi
    assert not options.boot


def test_parser_thread_quicktree(capsys):
    """
    Testing default that quicktree + threads != 1 gives an error
    """
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    with pytest.raises(SystemExit):
        allm.parse(parser,
                         "-c test/data/all/default-conffigfile.ini -o out-all "
                         "-n TEST --threads 2 --soft quicktree".split())
    _, err = capsys.readouterr()
    assert ("You cannot run quicktree with multiple threads. "
            "Choose another software, or remove the --threads option") in err


def test_duplicate_value(capsys):
    """
    test that when a parameter is given twice in the configfile, it raises an error
    """
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    with pytest.raises(SystemExit):
        allm.parse(parser,
                         "-c test/data/all/conffile-duplicate.ini -o out-all "
                         "-n TEST --threads 2 --soft quicktree".split())
    out, err = capsys.readouterr()
    assert ("option 'verbose' in section 'prepare' already exists") in out


def test_wrong_conffile_name(capsys):
    """
    test that when the given config file does not exist, it raises an error
    """
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    with pytest.raises(SystemExit):
        allm.parse(parser, "-c conf.ini -o out-all -t 1234 -n TEST".split())
    out, err = capsys.readouterr()
    assert ("Error: config file conf.ini not found.") in out


def test_annot_defvalue(capsys):
    """
    Test when a value is given in 'default' section of config file, but redefined in
    another section, it keeps the redefined value
    """
    parser = argparse.ArgumentParser(description="Run all modules", add_help=False)
    allm.build_parser(parser)
    options = allm.parse(parser,
                         "-c test/data/all/default-conffigfile.ini -o out-all -n TEST".split())
    assert options.outdir == "out-all"
    assert options.threads == 10
    assert not options.ncbi_species_taxid
    assert options.prodigal_only == False
    assert options.min_dist == '5'
    assert options.norefseq
    assert options.l90 == 99
    assert options.clust_mode == 1
    assert options.min_id == 0.85
    assert options.tol == 1
    assert not options.multi
    assert not options.boot


