#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of 'corepers' subcommand
"""

import pytest
import argparse

from PanACoTA.subcommands import corepers


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    with pytest.raises(SystemExit):
        corepers.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert "usage: " in err
    assert " -p PANGENOME -o OUTPUTDIR [-t TOL] [-M] [-X] [-F]" in err
    assert "[-v] [-q]" in err
    assert "[-h]" in err
    assert "the following arguments are required: -p, -o" in err


def test_parser_tol_nofloat(capsys):
    """
    Test that when the '-t tol' parameter given is not a number, it raises
    the expected error message
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    with pytest.raises(SystemExit):
        corepers.parse(parser, "-p pangenome -t mytol".split())
    _, err = capsys.readouterr()
    assert "argument -t tol: invalid float value: mytol" in err


def test_parser_tol_neg(capsys):
    """
    Test that when the '-t tol' parameter given is a negative number, it raises
    the expected error message
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    with pytest.raises(SystemExit):
        corepers.parse(parser, "-p pangenome -t -0.5".split())
    _, err = capsys.readouterr()
    assert ("The minimum %% of genomes required in a family to be persistent must be in [0, "
            "1]. Invalid value: -0.5") in err


def test_parser_multi_mixed(capsys):
    """
    Test that when the user asks for multi (-M) and mixed (-X) persistent
    genome, it returns an error message
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    with pytest.raises(SystemExit):
        corepers.parse(parser, "-p pangenome -o toto -M -X".split())
    _, err = capsys.readouterr()
    assert ("-M and -X options cannot be activated together. Choose if you want to:\n"
            "- allow several members in any number of genomes of a family (-M)\n"
            "- allow several members in only '1-tol'% of the genomes of a family "
            "(other 'tol'% genomes must have exactly 1 member) (-X)") in err


def test_parser_mixed_tol1(capsys):
    """
    Test that when the users asks for a mixed persistent genome with 100% of
    strains with exactly 1 member, it returns an error message
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    with pytest.raises(SystemExit):
        corepers.parse(parser, "-p pangenome -t 1 -X -o toto".split())
    _, err = capsys.readouterr()
    assert ("You are asking for mixed families, while asking for 100% of the genomes of "
            "a family to have exactly one member, which is not compatible. Do you want "
            "to \n- lower the percentage of genomes required to have exactly "
            "1 member (-t tol)\n- not allow mixed families (remove -X option)") in err


def test_parser_floor_tol1(capsys):
    """
    Test that when the user asks for 100% of strains having exactly 1 member,
    + the option to use 'floor', it returns an error message.
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    with pytest.raises(SystemExit):
        corepers.parse(parser, "-p pangenome -t 1 -F -o outdir".split())
    _, err = capsys.readouterr()
    assert ("You are asking to use floor('tol'*N) as a minimum number of genomes "
            "present in a family, but with 'tol'=1: the minimum number of genomes "
            "will always be equal to N, using floor or the default ceil! Either "
            "use a 'tol' lower than 1, or remove the '-F' option.") in err


def test_parser_default():
    """
    Test that when run with the minimum required arguments, all default values are
    as expected.
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    options = corepers.parse(parser, "-p pangenome -o outdir".split())
    assert options.pangenome == "pangenome"
    assert options.tol == 1
    assert options.multi is False
    assert options.mixed is False
    assert options.outputdir == "outdir"
    assert not options.floor
    assert options.verbose == 0
    assert not options.quiet


def test_parser_mixed_floor():
    """
    Test that when run with tol of 0.99, mixed, and use floor, it returns
    the expected values for all arguments
    """
    parser = argparse.ArgumentParser(description="Do corepers", add_help=False)
    corepers.build_parser(parser)
    options = corepers.parse(parser, "-p pangenome --tol 0.99 -X -F -o outdir".split())
    assert options.pangenome == "pangenome"
    assert options.tol == 0.99
    assert options.multi is False
    assert options.mixed is True
    assert options.outputdir == "outdir"
    assert options.floor is True
    assert options.verbose == 0
    assert not options.quiet
