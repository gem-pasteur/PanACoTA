#!/usr/bin/env python3
# coding: utf-8

"""
tree is a subcommand of genomeAPCAT

@author gem
June 2017
"""

import sys


def main_from_parse(args):
    """
    Call main function from the arguments given by parser
    """
    main(args.alignment, args.boot, args.outfile, args.threads, args.verbose, args.quiet)


def main(align, boot, outfile, threads, verbose, quiet):
    """
    Inferring a phylogenetic tree from an alignment file.
    """
    # import needed packages
    import logging
    import os
    from genomeAPCAT import utils
    from genomeAPCAT.tree_module import fasttree_func as ft

    # test if prokka is installed and in the path
    if not utils.check_installed("FastTreeMP"):  # pragma: no cover
        logger.error("FastTreeMP is not installed. 'genomeAPCAT tree' cannot run.")
        sys.exit(1)

    outdir = os.path.dirname(align)
    # name logfile, add timestamp if already existing
    logfile_base = os.path.join(outdir, "genomeAPCAT-tree")
    level = logging.DEBUG
    utils.init_logger(logfile_base, level, '', verbose=verbose, quiet=quiet)
    logger = logging.getLogger()

    ft.define_nb_threads(threads)
    ft.run_fasttree(alignfile, boot, outfile)


def build_parser(parser):
    """
    Method to create a parser for command-line options
    """
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
    required.add_argument("-a", dest="alignment", required=True,
                          help=("Alignment file in multi-fasta: each header will be a "
                                "leaf of the inferred tree."))

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-b", "--boot", dest="boot",
                          help=("Indicate how many bootstraps you want to compute. By "
                                "default, no bootstrap is calculated."))
    optional.add_argument("-o", dest="outfile",
                          help=("By default, the output tree file will be called "
                                "'<input_alignment_filename>.fasttree_tree.nwk'. You "
                                "can give a custom output name with this option."))
    optional.add_argument("--threads", dest="threads", default=1, type=thread_num,
                        help=("add this option if you want to parallelize on several threads. "
                              "Indicate on how many threads you want to parallelize. "
                              "By default, it uses 1 thread. Put 0 if you want to use "
                              "all threads of your computer."))

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
    parser = argparse.ArgumentParser(description=(("Infer phylogenetic tree based on "
                                                   "core/persistent genome")),
                                     add_help=False)
    build_parser(parser)
    OPTIONS = parse(parser, sys.argv[1:])
    main_from_parse(OPTIONS)
