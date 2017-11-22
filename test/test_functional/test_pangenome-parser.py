#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of pangenome.py
"""

import pytest
import argparse

from genomeAPCAT.subcommands import pangenome


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Do pangenome", add_help=False)
    pangenome.build_parser(parser)
    with pytest.raises(SystemExit):
        pangenome.parse(parser, "".split())
    _, err = capsys.readouterr()
    print(err)
    assert "usage: " in err
    assert "-l LSTINFO_FILE -n DATASET_NAME -d DBPATH" in err
    assert "-i MIN_ID -o OUTDIR\n" in err
    assert " [-f OUTFILE] [-c {0,1,2}]" in err
    assert "[-s SPEDIR] [--threads THREADS] [-v]" in err
    assert "[-q] [-h]" in err
    assert "the following arguments are required: -l, -n, -d, -i, -o" in err