#!/usr/bin/env python3
# coding: utf-8

"""
Util functions and classes.


@author gem
April 2017
"""

import os
import sys
import logging
import subprocess

logger = logging.getLogger()


class LessThanFilter(logging.Filter):
    """
    When using log, when a level is set to a handler, it is a minimum level. All
    levels higher than it will be printed. If you want to print only until
    a given level (no levels higher than the specified one), use this class like this:
    handler.addFilter(LessThanFilter(level))
    """
    def __init__(self, level):
        self._level = level
        logging.Filter.__init__(self)

    def filter(self, rec):
        return rec.levelno < self._level

def check_installed(cmd):
    """
    Check that the given command exists
    """
    FNULL = open(os.devnull, 'w')
    try:
        returncode = subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
        return returncode
    except Exception as err:
        logger.error(("{0} failed: {1}").format(cmd[0], err))
        sys.exit(1)

def plot_distr(values, limit, outfile, title, text):
    """ Plot histogram of given values, and add a vertical line corresponding to the choosen
     'limit' and saves the image into the 'outfile'

    :param values: list of values
    :type values: list
    :param limit: limit for which a vertical line must be drawn
    :type limit: int
    :param outfile: file in which the output image must be saved
    :type outfile: str
    """
    import math
    import numpy as np
    import matplotlib
    matplotlib.use('AGG')
    from matplotlib import pyplot as plt
    plt.figure(figsize=(10,7))
    max_x = max(values)
    # if too many values, group them to have less bins in the histogram.
    # Put 'group_values' values in each bin ->
    # if less than 300 values, 1 value per bin, otherwise more values per bin
    group_values = int(max_x/300) + 1
    dec_ax = math.exp(0.001 * max_x) - 1
    dec_text = 3 * dec_ax
    bins = np.arange(0, max_x + 2*group_values, group_values) - 0.5
    axes = plt.hist(values, bins = bins, edgecolor="black", color="blue")
    plt.xlim(0.5, max_x + 0.5*group_values)
    plt.axvline(x=limit + 0.5*group_values + dec_ax, color="r")
    plt.text(x=limit + 0.5*group_values + dec_text, y=plt.ylim()[1]/2,
             s=text + " " + str(limit), color="r", rotation=90)
    plt.title(title)
    plt.savefig(outfile)