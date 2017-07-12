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
    main(args.corepers, args.list_genomes, args.dataset_name, args.dbpath, args.outdir,
         args.verbose, args.quiet)


def main(corepers, list_genomes, dname, dbpath, outdir, verbose=0, quiet=False):
    """
    Align given core genome families
    """
    # import needed packages
    import logging
    from genomeAPCAT import utils
    from genomeAPCAT.align_module import pan_to_pergenome as p2g
    from genomeAPCAT.align_module import get_seqs as gseqs
    from genomeAPCAT.align_module import alignment as ali

    os.makedirs(outdir, exist_ok=True)
    # name logfile, add timestamp if already existing
    logfile_base = os.path.join(outdir, "genomeAPCAT-align_" + dname)
    level = logging.DEBUG
    utils.init_logger(logfile_base, level, '', verbose=verbose, quiet=quiet)
    logger = logging.getLogger()

    all_genomes, aldir, listdir, fam_nums = p2g.get_per_genome(corepers, list_genomes,
                                                               dname, outdir)
    gseqs.get_all_seqs(all_genomes, dname, dbpath, listdir, quiet)
    prefix = os.path.join(aldir, dname)
    ali.align_all_families(prefix, fam_nums, len(all_genomes))


def build_parser(parser):
    """
    Method to create a parser for command-line options
    """
    import argparse

    # Create command-line parser for all options and arguments to give
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-c", dest="corepers", required=True,
                        help=("Core or persistent genome whose families must be aligned."))
    required.add_argument("-l", dest="list_genomes", required=True,
                        help=("File containing the list of all the genomes you want "
                              "to align from their core/persistent families. "
                              "1 genome per line: it can be the "
                              "LSTINFO-<list_file>.lst file of 'genomeAPCAT annotate' module. "
                              "Here, only the first column (genome name without extension) "
                              "will be used. The final alignment file will contain "
                              "1 alignment per genome in this file."))
    required.add_argument("-n", dest="dataset_name", required=True,
                        help=("Name of the dataset which will be aligned (for example, "
                              "SAEN1234 for 1234 Salmonella enterica genomes). This name will "
                              "be used to name the alignment file."))
    required.add_argument("-d", dest="dbpath", required=True,
                        help=("Path to the folder containing the directories 'Proteins' "
                              "and 'Genes', created by 'genomeAPCAT annotate'."))
    required.add_argument("-o", dest="outdir", required=True,
                          help=("Output directory, where all results must be saved "))

    # optional = parser.add_argument_group('Optional arguments')
    # optional.add_argument("-m", dest="mafft_options",
    #                       help=("If you want to add options to mafft, such as --quiet..."))

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

