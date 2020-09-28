#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the functions of utils-argparse.py dealing with checking arguments given to argparse
"""

import logging
import os
import pytest
import shutil
import argparse

import PanACoTA.utils_argparse as autils


def test_gen_name(capsys):
    """
    Test that, when giving a gene name, if it does not have 4 characters (letters or num),
    it returns an error message
    """
    assert autils.gen_name("TOTO") == "TOTO"
    assert autils.gen_name("1234") == "1234"
    assert autils.gen_name("T1O2") == "T1O2"
    with pytest.raises(argparse.ArgumentTypeError):
        autils.gen_name("a long name")
        out, err = capsys.readouterr()
        assert ("The genome name must contain 4 characters. For example, this name can correspond "
                "to the 2 first letters of genus, and 2 first letters of species, e.g. "
                "ESCO for Escherichia Coli") in err
    with pytest.raises(argparse.ArgumentTypeError):
        autils.gen_name("-gdd")
        out, err = capsys.readouterr()
        assert ("The genome name must contain 4 characters. For example, this name can correspond "
                "to the 2 first letters of genus, and 2 first letters of species, e.g. "
                "ESCO for Escherichia Coli") in err


def test_date_name(capsys):
    """
    Test that when a date is given, it returns expected message
    """
    assert autils.date_name("0920") == "0920"
    assert autils.date_name("se20") == "se20"
    with pytest.raises(argparse.ArgumentTypeError):
        autils.gen_name("september 2020")
        out, err = capsys.readouterr()
        assert ("The date must contain 4 characers. Uually, it contains 4 digits, "
                "corresponding to the month (2 digits) and year (2 digits)") in err

def test_get_date():
    """
    test that it returns current date
    """
    import time
    t = time.strftime("%m%y")
    assert autils.get_date() == t


def test_cont_num():
    """
    Test that given value of contig number is valid
    """
    assert autils.cont_num(10) == 10
    with pytest.raises(argparse.ArgumentTypeError) as raised_err:
        autils.cont_num("-2")
    assert ("The maximum number of contigs allowed must be a "
            "positive number") in str(raised_err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.cont_num("10000")
    assert ("We do not support genomes with more than 9999 contigs") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.cont_num("a")
    assert ("argument --nbcont: invalid int value: a") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.cont_num("1.1")
    assert ("argument --nbcont: invalid int value: 1.") in str(err.value)

