#!/usr/bin/env python3
# coding: utf-8

"""
This script contains functions used by several tests.
"""

import hashlib
import os


def get_file_hash(filepath):
    """
    Returns the MD5 hash of a file 'filepath'. This is the signature of a file, used to check data
    integrity of this file.
    """
    fileobj = open(filepath, 'rb')
    m = hashlib.md5()
    while True:
        d = fileobj.read(2**20)
        if not d:
            break
        m.update(d)
    fileobj.close()
    return m.hexdigest()


def compare_files(file1, file2):
    """
    Checks if the two files are the same or not, by comparing their signatures.
    """
    sig_file1 = get_file_hash(file1)
    sig_file2 = get_file_hash(file2)
    return sig_file1 == sig_file2
