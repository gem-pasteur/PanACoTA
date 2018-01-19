#!/usr/bin/env python3
# coding: utf-8

"""
prepare is a subcommand of genomeAPCAT

It is a downloading a dataset if requested, and/or doing quality control on the dataset.
Steps are:

- if dataset to download:

    - download assembly_summary
    - download refseq NCBI genomes for the given species
    - add complete genomes (for GEM use)
- optional: find stretches of at least 'n' N (default 5), and cut into a new contig at this stretch
- for each genome, calc L90 and number of contigs (after cut at N stretches if used)
- order genomes by:

    - complete genomes first (if gembase use)
    - order by L90 and number of contigs
- iteratively run mash, to remove too close and too different genomes
- QC: keep only genomes with:

    - L90 <= x (default = 100)
    - #contig <= y (default = 999)
- rename those genomes and their contigs, with strain name increasing with quality (L90 and
  #contig)

Input:

- NCBI species taxid
- output directory where genomes must be downloaded
- number of threads to use
- the minimum mash distance required to keep a genome in the dataset (default 1e-4).

-- for GEM:
- gembase species (like ESCO001)
- gembase version (like 'Microbial_B_1116)

Output:

- In your given ``outdir``, you will find 5 folders:



@author gem
January 2018
"""
import sys


def main_from_parse(arguments):
    """
    Call main function from the arguments given by parser

    Parameters
    ----------
    arguments : argparse.Namespace
        result of argparse parsing of all arguments in command line

    """
    main(arguments.NCBI_species_taxid, arguments.outdir, arguments.min_dist, arguments.parallel)


def main(ncbi_taxid, outdir, min_dist, parallel):
    print(ncbi_taxid, outdir, min_dist, parallel)


def build_parser(parser):
    """
    Method to create a parser for command-line options

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser to configure

    """
    # TODO: add check that min distance is lower than 0.06 (and positive)
    import argparse
    import multiprocessing

    def thread_num(param):
        try:
            param = int(param)
        except Exception:
            msg = "argument --threads threads: invalid int value: {}".format(param)
            raise argparse.ArgumentTypeError(msg)
        nb_cpu = multiprocessing.cpu_count()
        if param > nb_cpu:
            msg = ("You have {} threads on your computer, you cannot ask for more: "
                   "invalid value: {}").format(nb_cpu, param)
            raise argparse.ArgumentTypeError(msg)
        elif param < 0:
            msg = ("Please provide a positive number of threads (or 0 for all threads): "
                   "Invalid value: {}").format(param)
            raise argparse.ArgumentTypeError(msg)
        elif param == 0:
            return nb_cpu
        return param

    # Create command-line parser for all options and arguments to give
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-t", dest="NCBI_species_taxid", required=True,
                          help=("Species taxid to download, corresponding to the "
                                "'species taxid' provided by the NCBI")
                          )

    # # gem options (to download complete genomes from gembase on abgfour)
    # gem = parser.add_argument_group('Specific arguments for GEM')
    # gem.add_argument("-n", dest="gembase_species",
    #                     help=("Species name in gembases, to add the complete genomes. "
    #                           "For example ESCO001 for Escherichia coli."))
    # gem.add_argument("-g", dest="gembase_version",
    #                     help=("Gembase version from which you want to get the complete "
    #                           "genomes. For example, 'Microbial_B_1116'."))

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-o", dest="outdir", default=".",
                          help=("Give the path to the directory where you want to save the "
                                "database. In the given diretory, it will create a folder with "
                                "the gembase species name. Inside this folder, you will find a "
                                "folder 'Database_init' containing all fasta files, as well as a "
                                "folder 'refseq' with files downloaded by ncbi_genome_download."))
    optional.add_argument("-p", dest="parallel", type=thread_num, default=1,
                          help=("Run 'N' downloads in parallel (default=1). Put 0 if "
                                "you want to use all cores of your computer."))
    # parser.add_argument("-R", dest="no_refseq", action="store_true",
    #                     help=("If you already downloaded refseq genomes and do not want to "
    #                           "check them, add this option to directly go to the next steps."))
    # parser.add_argument('--mash', dest="mash", action="store_true",
    #                     help=("Add this option if you already downloaded complete and refseq "
    #                           "genomes, and ran quality control (you have, in your result "
    #                           "folder, a file called "
    #                           "'info-genomes-list-<gembase_species>.lst'. "
    #                           "It will then get information on genomes quality from this "
    #                           "file, and run mash steps."))
    optional.add_argument("-m", dest="min_dist", default=1e-4, type=float,
                          help="By default, genomes whose distance to the reference is not "
                               "between 1e-4 and 0.06 are discarded. You can specify your own "
                               "lower limit (instead of 1e-4) with this option.")


def check_args(parser, args):
    """
    Check that arguments given to parser are as expected.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser used to parse command-line
    args : argparse.Namespace
        Parsed arguments

    Returns
    -------
    argparse.Namespace or None
        The arguments parsed, updated according to some rules. Exit program
        with error message if error occurs with arguments given.

    """
    # TODO: Check that if a gembase species is given, gembase version also, and vice versa
    return args


def parse(parser, argu):
    """
    arse arguments given to parser

    Parameters
    ----------
    parser : argparse.ArgumentParser
        the parser used
    argu : [str]
        command-line given by user, to parse using parser

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    args = parser.parse_args(argu)
    return check_args(parser, args)


if __name__ == '__main__':
    import argparse

    my_parser = argparse.ArgumentParser(description="Prepare the dataset", add_help=False)
    build_parser(my_parser)
    OPTIONS = parse(my_parser, sys.argv[1:])
    main_from_parse(OPTIONS)
