#!/usr/bin/env python3
# coding: utf-8

"""
This script contains functions used by several tests.
"""

import hashlib


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


def compare_files_bin(file1, file2):
    """
    Checks if the two files are the same or not, by comparing their signatures.
    """
    sig_file1 = get_file_hash(file1)
    sig_file2 = get_file_hash(file2)
    return sig_file1 == sig_file2


def compare_file_content(file1, file2):
    """
    Check that the 2 files have the same content, but no necessarily in the same order
    """
    with open(file1, "r") as f1:
        lines_f1 = f1.readlines()

    with open(file2, "r") as f2:
        lines_f2 = f2.readlines()
        if len(lines_f1) != len(lines_f2):
            print(f"not same number of lines in {file1} and {file2}")
            return False
        for line2 in lines_f2:
            if line2 not in lines_f1:
                print(f"'{line2}' not found")
                return False
    return True


def compare_order_content(file1, file2):
    """
    Check that the 2 files have exactly the same content, in the same order
    """
    with open(file1, "r") as f1:
        lines_f1 = f1.readlines()

    with open(file2, "r") as f2:
        lines_f2 = f2.readlines()
        if len(lines_f1) != len(lines_f2):
            print(f"not same number of lines in {file1} and {file2}")
            return False
        for line1, line2 in zip(lines_f1, lines_f2):
            if line1 != line2:
                print(f"'{line1}' vs {line2} do not correspond.")
                return False
    return True