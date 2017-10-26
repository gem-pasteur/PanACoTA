#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the utils_pangenome submodule of genomeAPCAT
"""

import pytest
import os

from genomeAPCAT import utils_pangenome as upan


def test_read_gene():
    """
    Check that when reading a given gene name, it extracts expected information
    """
