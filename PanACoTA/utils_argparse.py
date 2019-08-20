#!/usr/bin/env python3
# coding: utf-8

"""
Functions to check argparse aguments given by user


@author gem
April 2017
"""
from PanACoTA import utils
import argparse

def gen_name(param):
    if not utils.check_format(param):
        msg = ("The genome name must contain 4 characters. For example, this name can "
               "correspond to the 2 first letters of genus, and 2 first letters of "
               "species, e.g. ESCO for Escherichia Coli.")
        raise argparse.ArgumentTypeError(msg)
    return param


def date_name(param):
    if not utils.check_format(param):
        msg = ("The date must contain 4 characters. Usually, it contains 4 digits, "
               "corresponding to the month (2 digits) and year (2 digits).")
        raise argparse.ArgumentTypeError(msg)
    return param


def get_date():
    import time
    return time.strftime("%m%y")


def cont_num(param):
    try:
        param = int(param)
    except Exception:
        msg = "argument --nbcont: invalid int value: {}".format(param)
        raise argparse.ArgumentTypeError(msg)
    if param < 0:
        msg = "The maximum number of contigs allowed must be a positive number."
        raise argparse.ArgumentTypeError(msg)
    if param >= 10000:
        msg = "We do not support genomes with more than 9999 contigs."
        raise argparse.ArgumentTypeError(msg)
    return param


def thread_num(param):
    try:
        param = int(param)
    except Exception:
        msg = "argument --threads threads: invalid int value: {}".format(param)
        raise argparse.ArgumentTypeError(msg)
    nb_cpu = multiprocessing.cpu_count()
    if param > nb_cpu:
        msg = ("You have {} threads on your computer, you cannot ask for more: "
               "invalid value: {}").format(nb_cpu, param)
        raise argparse.ArgumentTypeError(msg)
    elif param < 0:
        msg = ("Please provide a positive number of threads (or 0 for all threads): "
               "Invalid value: {}").format(param)
        raise argparse.ArgumentTypeError(msg)
    elif param == 0:
        return nb_cpu
    return param