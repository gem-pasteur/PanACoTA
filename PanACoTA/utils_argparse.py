#!/usr/bin/env python3
# coding: utf-8

# ###############################################################################
# This file is part of PanACOTA.                                                #
#                                                                               #
# Authors: Amandine Perrin                                                      #
# Copyright © 2018-2020 Institut Pasteur (Paris).                               #
# See the COPYRIGHT file for details.                                           #
#                                                                               #
# PanACOTA is a software providing tools for large scale bacterial comparative  #
# genomics. From a set of complete and/or draft genomes, you can:               #
#    -  Do a quality control of your strains, to eliminate poor quality         #
# genomes, which would not give any information for the comparative study       #
#    -  Uniformly annotate all genomes                                          #
#    -  Do a Pan-genome                                                         #
#    -  Do a Core or Persistent genome                                          #
#    -  Align all Core/Persistent families                                      #
#    -  Infer a phylogenetic tree from the Core/Persistent families             #
#                                                                               #
# PanACOTA is free software: you can redistribute it and/or modify it under the #
# terms of the Affero GNU General Public License as published by the Free       #
# Software Foundation, either version 3 of the License, or (at your option)     #
# any later version.                                                            #
#                                                                               #
# PanACOTA is distributed in the hope that it will be useful, but WITHOUT ANY   #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS     #
# FOR A PARTICULAR PURPOSE. See the Affero GNU General Public License           #
# for more details.                                                             #
#                                                                               #
# You should have received a copy of the Affero GNU General Public License      #
# along with PanACOTA (COPYING file).                                           #
# If not, see <https://www.gnu.org/licenses/>.                                  #
# ###############################################################################

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
    import multiprocessing
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


def positive_int(param):
    try:
        param = int(param)
    except ValueError:
        msg = f"error: argument --cutn: invalid int value: '{param}'"
        raise argparse.ArgumentTypeError(msg)
    if param < 0:
        msg = f"error: argument --cutn must be a positive integer: invalid int value: '{param}'"
        raise argparse.ArgumentTypeError(msg)
    return param


def mash_dist(param):
    try:
        param = float(param)
    except ValueError:
        msg = f"error: mash distance: invalid float value: '{param}'"
        raise argparse.ArgumentTypeError(msg)
    if param < 0 or param > 1:
        msg = f"error: mash distance must be between 0 and 1: invalid value: '{param}'"
        raise argparse.ArgumentTypeError(msg)
    return param
