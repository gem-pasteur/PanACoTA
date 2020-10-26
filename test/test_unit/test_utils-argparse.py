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
import multiprocessing

import PanACoTA.utils_argparse as autils


def test_gen_name():
    """
    Test that, when giving a gene name, if it does not have 4 characters (letters or num),
    it returns an error message
    """
    assert autils.gen_name("TOTO") == "TOTO"
    assert autils.gen_name("1234") == "1234"
    assert autils.gen_name("T1O2") == "T1O2"
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.gen_name("a long name")
    assert ("The genome name must contain 4 characters. For example, this name "
            "can correspond to the 2 first letters of genus, and 2 first letters of "
            "species, e.g. ESCO for Escherichia Coli") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.gen_name("-gdd")
    assert ("The genome name must contain 4 characters. For example, this name can correspond "
            "to the 2 first letters of genus, and 2 first letters of species, e.g. "
            "ESCO for Escherichia Coli") in str(err.value)


def test_date_name(capsys):
    """
    Test that when a date is given, it returns expected message
    """
    assert autils.date_name("0920") == "0920"
    assert autils.date_name("se20") == "se20"
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.date_name("september 2020")
    assert ("The date must contain 4 characters. Usually, it contains 4 digits, "
            "corresponding to the month (2 digits) and year (2 digits)") in str(err.value)


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


def test_thread_num():
    """
    Test that given number of threads is as expected
    """
    assert autils.thread_num("1") == 1
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num("1.1")
    assert ("argument --threads threads: invalid int value: 1.1") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num("a")
    assert ("argument --threads threads: invalid int value: a") in str(err.value)
    nb_cpu = multiprocessing.cpu_count()
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num(str(nb_cpu*2))
    assert (f"You have {nb_cpu} threads on your computer, you cannot ask for more: invalid value: "
            f"{nb_cpu*2}") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num("-1")
    assert ("Please provide a positive number of threads (or 0 for all threads): "
            "Invalid value: -1") in str(err.value)
    assert autils.thread_num(0) == nb_cpu


def test_positive_int():
    """
    Test checking that given argument is a positive integer
    """
    assert autils.positive_int("1") == 1
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.positive_int("1.1")
    assert ("argument --cutn: invalid int value: '1.1'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.positive_int("-1")
    assert ("error: argument --cutn must be a positive integer: "
            "invalid int value: '-1'") in str(err.value)


def test_mash_dist():
    """
    Test checking that given value is ok for a mash distance
    """
    assert autils.mash_dist("0.05") == 0.05
    assert autils.mash_dist("1e-4") == 0.0001
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("1.1.1")
    assert ("error: mash distance: invalid float value: '1.1.1'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("one")
    assert ("error: mash distance: invalid float value: 'one'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("1.000001")
    assert ("error: mash distance must be between 0 and 1: "
            "invalid value: '1.000001'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("-1e-4")
    assert ("error: mash distance must be between 0 and 1: "
            "invalid value: '-0.0001'") in str(err.value)
