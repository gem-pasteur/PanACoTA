#!/usr/bin/env python3
# coding: utf-8

# ###############################################################################
# This file is part of PanACOTA.                                                #
#                                                                               #
# Authors: Amandine Perrin                                                      #
# Copyright © 2018-2020 Institut Pasteur (Paris).                               #
# See the COPYRIGHT file for details.                                           #
#                                                                               #
# PanACOTA is a software providing tools for large scale bacterial comparative  #
# genomics. From a set of complete and/or draft genomes, you can:               #
#    -  Do a quality control of your strains, to eliminate poor quality         #
# genomes, which would not give any information for the comparative study       #
#    -  Uniformly annotate all genomes                                          #
#    -  Do a Pan-genome                                                         #
#    -  Do a Core or Persistent genome                                          #
#    -  Align all Core/Persistent families                                      #
#    -  Infer a phylogenetic tree from the Core/Persistent families             #
#                                                                               #
# PanACOTA is free software: you can redistribute it and/or modify it under the #
# terms of the Affero GNU General Public License as published by the Free       #
# Software Foundation, either version 3 of the License, or (at your option)     #
# any later version.                                                            #
#                                                                               #
# PanACOTA is distributed in the hope that it will be useful, but WITHOUT ANY   #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS     #
# FOR A PARTICULAR PURPOSE. See the Affero GNU General Public License           #
# for more details.                                                             #
#                                                                               #
# You should have received a copy of the Affero GNU General Public License      #
# along with PanACOTA (COPYING file).                                           #
# If not, see <https://www.gnu.org/licenses/>.                                  #
# ###############################################################################

"""
tree is a subcommand of PanACoTA

@author gem
June 2017
"""
import sys


def main_from_parse(args):
    """
    Call main function from the arguments given by parser

    Parameters
    ----------
    args : argparse.Namespace
        result of argparse parsing of all arguments in command line
    """
    cmd = "PanACoTA " + ' '.join(args.argv)
    main(cmd, args.alignment, args.outdir, args.soft, args.model, args.threads,
         args.boot, args.write_boot, args.memory, args.fast, args.verbose, args.quiet)


def main(cmd, align, outdir, soft, model, threads, boot=False, write_boot=False,
         memory=False, fast=False, verbose=0, quiet=False):
    """
    Inferring a phylogenetic tree from an alignment file, with the given software.

    Parameters
    ----------
    cmd: str
        command used to launch tree module
    align: str
        Path to file containing alignments of persistent families grouped by genome
    outdir: str or None
        Path to file which will contain the tree inferred
    soft: str
        Soft to use to infer the phylogenetic tree: 1 of quicktree, fasttree or fastme
    model: str or None
        DNA substitution model chosen by user, None if quicktree used
    threads: int
        Maximum number of threads to use
    boot: int or None
        Number of bootstraps to compute. None if no bootstrap asked
    write_boot: bool
        True if all bootstrap pseudo-trees must be saved into a file, False otherwise
    memory: str
        Maximal RAM usage in GB | MB | % - Only for iqtree
    fast: boolean
        use -fast option with IQtree
    verbose : int
        verbosity:
        - defaut 0 : stdout contains INFO, stderr contains ERROR.
        - 1: stdout contains INFO, stderr contains WARNING and ERROR
        - 2: stdout contains (DEBUG), DETAIL and INFO, stderr contains WARNING and ERROR
        - >=15: Add DEBUG in stdout
    quiet: bool
        True if nothing must be sent to stdout/stderr, False otherwise
    """
    # import needed packages
    import logging
    import os
    from PanACoTA import utils
    from PanACoTA import __version__ as version
    tree = None
    if soft == "fasttree":
        # test if fasttree is installed and in the path
        if not utils.check_installed("FastTreeMP"): # pragma: no cover
            print("FastTreeMP is not installed. 'PanACoTA tree' cannot run.")
            sys.exit(1)
        from PanACoTA.tree_module import fasttree_func as tree
    elif soft == "fastme":
        # test if fastME is installed and in the path
        if not utils.check_installed("fastme"): # pragma: no cover
            print("fastme is not installed. 'PanACoTA tree' cannot run.")
            sys.exit(1)
        from PanACoTA.tree_module import fastme_func as tree
    elif soft == "quicktree":
        # test if fastME is installed and in the path
        if not utils.check_installed("quicktree"):  # pragma: no cover
            print("quicktree is not installed. 'PanACoTA tree' cannot run.")
            sys.exit(1)
        from PanACoTA.tree_module import quicktree_func as tree
    elif soft == "iqtree2":
        # by default, iqtree2 (not iqtree).
        # So, if user did not specify, it means iqtree2. But if 'iqtree2' command
        # does not exist, use iqtree command instead.
        # test if iqtree2 is installed and in the path
        if not utils.check_installed("iqtree2"): # pragma: no cover
            if not utils.check_installed("iqtree"):
                print("IQtree2 is not installed. 'PanACoTA tree' cannot run.")
                sys.exit(1)
            else:
                soft == "iqtree"
        from PanACoTA.tree_module import iqtree_func as tree
    elif soft == "iqtree":
        # user specifically asked for iqtree (version 1)
        if not utils.check_installed("iqtree"): # pragma: no cover
            print("IQtree is not installed. 'PanACoTA tree' cannot run.")
            sys.exit(1)
        from PanACoTA.tree_module import iqtree_func as tree

    # If outdir does not already exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # name logfile, add timestamp if already existing
    logfile_base = os.path.join(outdir, "PanACoTA-tree-" + soft)
    # level is the minimum level that will be considered.
    # for verbose = 0 or 1, ignore details and debug, start from info
    if verbose <= 1:
        level = logging.INFO
    # for verbose = 2, ignore only debug
    if verbose >= 2 and verbose < 15:
        level = 15 # int corresponding to detail level
    # for verbose >= 15, write everything
    if verbose >= 15:
        level = logging.DEBUG

    utils.init_logger(logfile_base, level, 'tree', verbose=verbose, quiet=quiet, log_details=True)
    logger = logging.getLogger("tree")
    logger.info(f'PanACoTA version {version}')
    logger.info("Command used\n \t > " + cmd)
    tree.run_tree(align, boot, outdir, quiet, threads, model=model, wb=write_boot,
                  mem=memory, s=soft, f=fast)

    logger.info("END")


def build_parser(parser):
    """
    Method to create a parser for command-line options

    Parameters
    ----------
    parser : argparse.ArgumentParser
        parser to configure in order to extract command-line arguments
    """
    import argparse
    from PanACoTA import utils_argparse


    # Create command-line parser for all options and arguments to give
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-a", dest="alignment", required=True,
                          help=("Alignment file in multi-fasta: each header will be a "
                                "leaf of the inferred tree."))
    required.add_argument("-o", dest="outdir", required=True,
                          help=("Directory where tree results will be saved. "))

    # Choose with which soft inferring phylogenetic tree
    softparse = parser.add_argument_group('Choose soft to use (default is IQtree2)')
    softs = ["fasttree", "fastme", "quicktree", "iqtree", "iqtree2"]
    softparse.add_argument("-s", "--soft", dest="soft", choices=softs, default="iqtree2",
                           help=("Choose with which software you want to infer the "
                                 "phylogenetic tree. Default is IQtree2 "
                                 "(versions 2.x of IQtree). If you want version 1.x of "
                                 "IQtree, use '-s iqtree'"))

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-b", "--boot", dest="boot", type=int,
                          help=("Indicate how many bootstraps you want to compute. By "
                                "default, no bootstrap is calculated. For IQtree, it "
                                "will use ultrafast bootstrap (>=1000)."))

    optional.add_argument("--threads", dest="threads", default=1, type=utils_argparse.thread_num,
                          help=("add this option if you want to parallelize on several threads. "
                                "Indicate on how many threads you want to parallelize. "
                                "By default, it uses 1 thread. Put 0 if you want to use "
                                "all threads of your computer. Not available with quicktree."))
    optional.add_argument("-m", "--model", dest="model",
                          help=("Choose your DNA substitution model.\n"
                                "Default for FastTree and IQtree: GTR. Default for FastME: TN93.\n"
                                "For FastTree, the choices are 'GTR' and 'JC'.\n"
                                "For FastME, choices are: 'p-distance' "
                                "(or 'p'), 'RY symmetric' (or 'Y'), 'RY' (or 'R'), "
                                "'JC69' (or 'J'), 'K2P' (or 'K'), 'F81' (or '1'), "
                                "'F84' (or '4'), 'TN93' (or 'T'), 'LogDet' (or 'L').\n"
                                "For IQtree, choices are HKY, JC, F81, K2P, K3P, K81uf,"
                                " TNef, TIM, TIMef, TVM, TVMef, SYM, GTR."))
    optional.add_argument("-B", dest="write_boot", action="store_true",
                          help=("Add this option if you want to write all bootstrap "
                                "pseudo-trees. Only available with FastME and IQtree."))
    optional.add_argument("--mem", dest="memory",
                          help=("Maximal RAM usage in GB | MB. Only available with iqtree."))
    optional.add_argument("-fast", dest="fast", action="store_true",
                          help=("Use -fast option with iqtree."))

    helper = parser.add_argument_group('Others')
    helper.add_argument("-v", "--verbose", dest="verbose", action="count", default=0,
                        help="Increase verbosity in stdout/stderr.")
    helper.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                        help=("Do not display anything to stdout/stderr. log files will "
                              "still be created."))
    helper.add_argument("-h", "--help", dest="help", action="help",
                        help="show this help message and exit")


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
    models_fastme = {"p-distance": "p", "RY-symetric": "Y", "RY": "R",
                     "JC69": "J", "K2P": "K", "F81": "1", "F84": "4",
                     "TN93": "T", "LogDet": "L"}
    models_fasttree = {"GTR": "-gtr", "JC": ""}
    models_iqtree = set(["HKY", "JC", "F81", "K2P", "K3P", "K81uf",
                         "TNef", "TIM", "TIMef", "TVM", "TVMef", "SYM", "GTR"])
    models_iqtree = {mod: mod for mod in models_iqtree}

    def check_model(models, choice):
        if choice in models.keys():
            return models[choice]
        elif choice in models.values():
            return choice
        mmsg = ("{} is not an available model for {}. Please choose an available DNA model "
                "(see -h for more details)").format(choice, args.soft)
        parser.error(mmsg)

    if args.soft == "quicktree" and args.threads != 1:
        msg = ("You cannot run quicktree with multiple threads. Choose another software, "
               "or remove the --threads option.")
        parser.error(msg)

    if args.soft == "quicktree" and args.model:
        msg = "Quicktree only runs the NJ algorithm. You cannot choose a DNA substitution model."
        parser.error(msg)

    # Memory option only available with iqtree
    if args.soft != "iqtree" and args.soft != "iqtree2" and args.memory:
        msg = "'--mem' option is only available for IQtree."
        parser.error(msg)

    # If bootstraps are asked with iqtree, check the number is >= 1000
    if (args.soft == "iqtree" or args.soft == "iqtree2") and args.boot and int(args.boot) < 1000:
        msg = "With IQtree, number of replicates for bootstraps must be >= 1000."
        parser.error(msg)

    # Write bootstrap option only available for fastme and iqtree
    if (args.soft != "iqtree" and args.soft != "iqtree2" and args.soft != "fastme"
        and args.write_boot):
        msg = "'-B' option is only available with FastME and IQtree."
        parser.error(msg)

    # Fast option only available for iqtree
    if (args.fast and ((args.soft != "iqtree" and args.soft != "iqtree2")
        or (args.boot or args.write_boot))):
        msg = ("-fast option is available only for IQtree, and not compatible "
               "with '-B' and '-b' options (bootstraps).")
        parser.error(msg)

    # Check model name is valid for the chosen soft
    if args.soft == "fastme":
        if args.model:
            args.model = check_model(models_fastme, args.model)
        else:
            args.model = "T"
    elif args.soft == "fasttree":
        if args.model:
            args.model = check_model(models_fasttree, args.model)
        else:
            args.model = "-gtr"
    elif args.soft == "iqtree" or args.soft == "iqtree2":
        if args.model:
            args.model = check_model(models_iqtree, args.model)
        else:
            args.model = "GTR"
    return args


def parse(parser, argu):
    """
    Parse arguments given to parser

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

    myparser = argparse.ArgumentParser(description=(("Infer phylogenetic tree based on "
                                                     "core/persistent genome")),
                                       add_help=False)
    build_parser(myparser)
    OPTIONS = parse(myparser, sys.argv[1:])
    main_from_parse(OPTIONS)
