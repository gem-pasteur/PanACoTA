#!/usr/bin/env python3
# coding: utf-8

"""
pan-genome is a subcommand of genomeAPCAT

@author gem
May 2017
"""

def main_from_parse(arguments):
    """
    Call main function from the arguments given by parser
    """
    main()


def main():
    """
    Main method, doing all steps:
    -

    """
    # import needed packages


def build_parser(parser):
    """
    Method to create a parser for command-line options
    """
    import argparse
    # Create command-line parser for all options and arguments to give
    parser.add_argument(dest="list_file",
                        help=("File containing the list of genome filenames to annotate (1 genome"
                              " per line). Each genome is in multi-fasta format. You can "
                              "specify the species name (4 characters) you want to give to each "
                              "genome by adding it after the genome filename(s), separated "
                              "by '::'. If not given, the species name will be the one given in "
                              "'species' argument. You can also specify the date (4 digits) "
                              "by adding '.' + your date choice after the genome "
                              "filename(s), '::' and, if given, the species name."))


def parse(parser, argu):
    """
    Parse arguments given to parser
    """
    args = parser.parse_args(argu)
    return check_args(parser, args)


def check_args(parser, args):
  """
  Check that arguments given to parser are as expected.
  """
  return args


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=("Compute pan-genome"))
    build_parser(parser)
    OPTIONS = parse(parser, sys.argv[1:])
    main_from_parse(OPTIONS)
