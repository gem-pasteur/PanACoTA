#!/usr/bin/env python3
# coding: utf-8

"""
Util functions and classes.


@author gem
April 2017
"""

import logging


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