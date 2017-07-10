#!/usr/bin/env python3
# coding: utf-8

"""
align is a subcommand of genomeAPCAT


@author gem
June 2017
"""

import os
import sys


def main_from_parse(args):
    """
    Call main function from the arguments given by parser
    """
    main()


def main():
    """
    Align given core genome families
    """
    print("alignment module")


def build_parser(parser):
    """
    Method to create a parser for command-line options
    """
    import argparse

    # Create command-line parser for all options and arguments to give
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-p", dest="pangenome", required=True,
                        help=("PanGenome file (1 line per family, first column is fam number)"))

    optional = parser.add_argument_group('Optional arguments')
    # optional.add_argument("-t", dest="tol", default=1, type=percentage,
    #                       help=("min %% of genomes having at least 1 member in a family to "
    #                             "consider the family as persistent (between 0 and 1, "
    #                             "default is 1 = 100%% of genomes = Core genome)"))

    helper = parser.add_argument_group('Others')
    helper.add_argument("-v", "--verbose", dest="verbose", action="count", default=0,
                        help=("Increase verbosity in stdout/stderr."))
    helper.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                        help=("Do not display anything to stdout/stderr. log files will "
                              "still be created."))
    helper.add_argument("-h", "--help", dest="help", action="help",
                        help="show this help message and exit")


def check_args(parser, args):
    """
    Check that arguments given to parser are as expected.
    """
    # if args.multi and args.mixed:
    #     parser.error("-M and -X options cannot be activated together. Choose if you want to:\n"
    #                  "- allow several members in any number of genomes of a family (-M)\n"
    #                  "- allow several members in only '1-tol'% of the genomes of a family "
    #                  "(other 'tol'% genomes must have exactly 1 member) (-X)")
    # if args.mixed and args.tol==1:
    #     parser.error("You are asking for mixed families, while asking for 100% of the genomes of "
    #                  "a family to have exactly one member, which is not compatible. Do you want "
    #                  "to \n- lower the percentage of genomes required to have exactly "
    #                  "1 member (-t tol)\n- not allow mixed families (remove -X option)")
    return args


def parse(parser, argu):
    """
    Parse arguments given to parser
    """
    args = parser.parse_args(argu)
    return check_args(parser, args)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=("Align Core/Persistent families"),
                                     add_help=False)
    build_parser(parser)
    OPTIONS = parse(parser, sys.argv[1:])
    main_from_parse(OPTIONS)

