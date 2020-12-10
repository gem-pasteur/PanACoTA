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
import configparser
import sys


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


def percentage(param):
        try:
            param = float(param)
        except Exception:
            msg = "argument -t tol: invalid float value: {}".format(param)
            raise argparse.ArgumentTypeError(msg)
        if param < 0 or param > 1:
            msg = ("The minimum %% of genomes required in a family to be persistent must "
                   "be in [0, 1]. Invalid value: {}".format(param))
            raise argparse.ArgumentTypeError(msg)
        return param


def perc_id(param):
    try:
        param = float(param)
    except Exception:
        msg = "argument -i percentage_id: invalid float value: {}".format(param)
        raise argparse.ArgumentTypeError(msg)
    if param < 0 or param > 1:
        msg = ("The minimum %% of identity must be in [0, 1]. Invalid value: {}".format(param))
        raise argparse.ArgumentTypeError(msg)
    return param


class Conf_all_parser(configparser.ConfigParser):
    """
    Read configfile and return arguments found, according to required type
    """
    def __init__(self, conffile, sections):
        super().__init__()
        self.read(conffile)
        self.sec_dicts = {}
        for sec in sections:
            self.sec_dicts[sec] = dict(self[sec])

    def get_section_dict(self, section):
        """
        get dictionary of values for 'section' section
        """
        return self.sec_dicts[section]

    def add_default(self, defargs, section):
        """
        Add all default arguments (defargs) in section dict.
        """
        for key, val in defargs.items():
            if key not in self.sec_dicts[section]:
                self[section][key] = str(val)
                self.sec_dicts[section][key] = val

    def update(self, args, section):
        """
        Add all arguments from args. If key already exists in self, overwrite it.
        Otherwise, create it.
        """
        self.sec_dicts[section].update(args)
        for key, val in self.sec_dicts[section].items():
            self[section][key] = str(val)

    def set_boolean(self, section, param):
        """
        Change param of section to boolean
        """
        try:
            bool_param = self.getboolean(section, param)
            self.sec_dicts[section][param] = bool_param
        except ValueError as err:
            val = self[section][param]
            print(f"ERROR: {param} must be a boolean. Wrong value: {val}.")
            sys.exit(1)

    def set_int(self, section, param):
        """
        Change param of section to boolean
        """
        try:
            int_param = self.getint(section, param)
            self.sec_dicts[section][param] = int_param
        except ValueError as err:
            val = self[section][param]
            print(f"ERROR: {param} must be an int. Wrong value: {val}.")
            sys.exit(1)

    def set_float(self, section, param):
        """
        Change param of section to boolean
        """
        try:
            float_param = self.getfloat(section, param)
            self.sec_dicts[section][param] = float_param
        except ValueError as err:
            val = self[section][param]
            print(f"ERROR: {param} must be a float. Wrong value: {val}.")
            sys.exit(1)
