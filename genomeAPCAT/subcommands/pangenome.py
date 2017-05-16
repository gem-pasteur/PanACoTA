#!/usr/bin/env python3
# coding: utf-8

"""
pan-genome is a subcommand of genomeAPCAT

@author gem
May 2017
"""

import sys


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
    def perc_id():
        try:
            param = float(param)
        except Exception:
            msg = "argument -i percentage_id: invalid float value: {}".format(param)
            raise argparse.ArgumentTypeError(msg)
        if param < 0 or param > 1:
            msg = ("The minimum %% of identity must be in [0, 1]. Invalid value: {}".format(param))
            raise argparse.ArgumentTypeError(msg)
        return param
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-l", dest="lstinfo_file", required=True,
                        help=("File containing the list of all genomes to include in "
                              "the pan-genome, 1 genome per line: it can be the "
                              "LSTINFO-<list_file>.lst file of 'genomeAPCAT annotate' module."
                              "Here, only the first column (genome name without extension) "
                              "will be used. All proteins of all these genomes will be "
                              "concatenated in a file called <dataset_name>.All.prt."))
    parser.add_argument("-n", dest="dataset_name", required=True,
                        help=("Name of the dataset which will be clustered (for exemple, "
                              "SAEN1234 for 1234 Salmonella enterica genomes). This name will "
                              "be used to name the protein databank, a well as the "
                              "pangenome files."))
    parser.add_argument("-d", dest="dbpath", required=True,
                        help=("Path to the folder containing all protein files corresponding "
                              "to the genomes of the dataset (output directory 'Proteins' "
                              "of 'genomeAPCAT annotate' module)."))
    parser.add_argument("-i", dest="min_id", required=True, type=perc_id,
                        help=("Minimum sequence identity to be considered in the same "
                              "cluster (number between 0 and 1)"))
    parser.add_argument("-o", dest="outdir", required=True,
                        help=("Output directory, where all results must be saved "
                              "(including tmp folder)"))
    parser.add_argument("-c", dest="clust_mode", choices=[0, 1, 2], default=1,
                        help=("Choose the clustering mode: 0 for 'set cover', 1 for "
                              "'single-linkage', 2 for 'CD-Hit'. Default "
                              "is 'single-linkage' (1)"))
    parser.add_argument("-s", dest="spedir",
                        help=("use this option if you want to save the concatenated protein "
                              "databank in another directory than the one containing all "
                              "individual protein files ('Proteins' folder)."))
    parser.add_argument("--threads", dest="threads",
                        help=("add this option if you want to parallelize on several threads. "
                              "Indicate on how many threads you want to parallelize. "
                              "By default, it uses all threads in the computer."))


def parse(parser, argu):
    """
    Parse arguments given to parser
    """
    args = parser.parse_args(argu)
    return args


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=("Compute pan-genome"))
    build_parser(parser)
    OPTIONS = parse(parser, sys.argv[1:])
    main_from_parse(OPTIONS)
