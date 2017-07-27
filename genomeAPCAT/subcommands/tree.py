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
    main(args.alignment, args.boot, args.outfile, args.soft, args.model,
         args.write_boot, args.threads, args.verbose, args.quiet)


def main(align, boot, outfile, soft, model, write_boot, threads, verbose, quiet):
    """
    Inferring a phylogenetic tree from an alignment file.
    """
    # import needed packages
    import logging
    import os
    from genomeAPCAT import utils
    if soft == "fasttree":
        # test if fasttree is installed and in the path
        if not utils.check_installed("FastTreeMP"):  # pragma: no cover
            logger.error("FastTreeMP is not installed. 'genomeAPCAT tree' cannot run.")
            sys.exit(1)
        from genomeAPCAT.tree_module import fasttree_func as tree
    elif soft == "fastme":
        from genomeAPCAT.tree_module import fastme_func as tree
        # test if fastME is installed and in the path
        if not utils.check_installed("fastme"):  # pragma: no cover
            logger.error("fastme is not installed. 'genomeAPCAT tree' cannot run.")
            sys.exit(1)
    outdir = os.path.dirname(align)
    # name logfile, add timestamp if already existing
    logfile_base = os.path.join(outdir, "genomeAPCAT-tree-" + soft)
    level = logging.DEBUG
    utils.init_logger(logfile_base, level, '', verbose=verbose, quiet=quiet)
    logger = logging.getLogger()

    tree.run_tree(align, boot, threads, outfile, quiet, model, write_boot)


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

    # Choose with which soft inferring phylogenetic tree
    softparse = parser.add_argument_group('Choose soft to use (default is fasttree)')
    softs = ["fasttree", "fastme"]
    softparse.add_argument("-s", "--soft", dest="soft", choices=softs, default="fasttree",
                           help=("Choose with which software you want to infer the "
                                 "phylogenetic tree. Default is FastTree"))

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-b", "--boot", dest="boot",
                          help=("Indicate how many bootstraps you want to compute. By "
                                "default, no bootstrap is calculated."))
    optional.add_argument("-o", dest="outfile",
                          help=("By default, the output tree file will be called "
                                "'<input_alignment_filename>.<software_used>_tree.nwk'. You "
                                "can give a custom output name with this option."))
    optional.add_argument("--threads", dest="threads", default=1, type=thread_num,
                        help=("add this option if you want to parallelize on several threads. "
                              "Indicate on how many threads you want to parallelize. "
                              "By default, it uses 1 thread. Put 0 if you want to use "
                              "all threads of your computer."))
    optional.add_argument("-m", "--model", dest="model",
                          help=("Choose your DNA substitution model. "
                                "Default for FastTree: GTR. Default for FastME: TN93. "
                                "For FastTree, the choices are 'GTR' and 'JC'. "
                                "For FastME, choices are: 'p-distance' "
                                "(or 'p'), 'RY symmetric' (or 'Y'), 'RY' (or 'R'), "
                                "'JC69' (or 'J'), 'K2P' (or 'K'), 'F81' (or '1'), "
                                "'F84' (or '4'), 'TN93' (or 'T'), 'LogDet' (or 'L')."))
    optional.add_argument("-B", dest="write_boot", action="store_true",
                          help=("Add this option if you want to write all bootstrap "
                                "pseudo-trees. Only available with FastME."))

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
    models_fastme = {"p-distance": "p", "RY-symetric": "Y", "RY": "R",
                     "JC69": "J", "K2P": "K", "F81": "1", "F84": "4",
                     "TN93": "T", "LogDet": "L"}
    models_fasttree = {"GTR": "-gtr", "JC": ""}

    def check_model(models, choice):
        if choice in models.keys():
            return models[choice]
        elif choice in models.values():
            return choice
        msg = ("{} is not an available model for {}. Please choose an available DNA model "
               "(see -h for more details)").format(choice, args.soft)
        parser.error(msg)

    if args.soft == "fastme":
        if args.model:
            args.model = check_model(models_fastme, args.model)
        else:
            args.model = "T"
    elif args.soft == "fasttree":
        if args.write_boot:
            msg = ("'-B' option is only available with FastME, not with FastTree")
            parser.error(msg)
        if args.model:
            args.model = check_model(models_fasttree, args.model)
        else:
            args.model = "-gtr"

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
