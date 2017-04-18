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
    except Exception as err:
        logger.error(("{0} failed: {1}").format(cmd[0], err))
        sys.exit(1)
