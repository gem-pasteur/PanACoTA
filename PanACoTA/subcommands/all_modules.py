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
'all' is a module of PanACoTA, allowing to run the whole pipeline at once.


@author gem
October 2020
"""

import os
import sys
from termcolor import colored
import sys

from PanACoTA import utils
from PanACoTA.subcommands import prepare
from PanACoTA.subcommands import annotate
from PanACoTA.subcommands import pangenome
from PanACoTA.subcommands import corepers
from PanACoTA.subcommands import align
from PanACoTA.subcommands import tree
from PanACoTA import __version__ as version


def main_from_parse(args):
    """
    Call main function from the arguments given by parser

    Parameters
    ----------
    args : argparse.Namespace
        result of argparse parsing of all arguments in command line
    """
    cmd = "PanACoTA " + ' '.join(args.argv)
    main(cmd, args.outdir, args.threads, args.NCBI_species_taxid, args.NCBI_species,
         args.levels, args.cutn, args.l90, args.nbcont, args.name, args.prodigal_only, args.min_id,
         args.tol, args.multi, args.mixed, args.soft, verbose=args.verbose, quiet=args.quiet)


def main(cmd, outdir, threads, NCBI_species_taxid, NCBI_species, levels, cutn, l90, nbcont,
         name, prodigal_only, min_id, tol, multi, mixed, soft, verbose=0, quiet=False):
    """
    Call all modules, one by one, using output of one as input for the next one
    """
    os.makedirs(outdir, exist_ok=True)
    # Initialize logger
    import logging
    # set level of logger: level is the minimum level that will be considered.
    if verbose <= 1:
        level = logging.INFO
    # for verbose = 2, ignore only debug
    if verbose >= 2 and verbose < 15:
        level = utils.detail_lvl() # int corresponding to detail level
    # for verbose >= 15, write everything
    if verbose >= 15:
        level = logging.DEBUG
    logfile_base = os.path.join(outdir, "PanACoTA-all_modules")
    logfile_base = utils.init_logger(logfile_base, level, name='all_modules',
                                     verbose=verbose, quiet=quiet)
    logger = logging.getLogger('all_modules')
    logger.info(f'PanACoTA version {version}')
    logger.info("Command used\n \t > " + cmd)

    # Run prepare module
    outdir_prepare = os.path.join(outdir, "1-prepare_module")
    tmp_dir = ""
    no_refseq = False
    db_dir = ""
    only_mash = False
    info_file = ""
    min_dist = 1e-4
    max_dist = 0.06

    info_file = prepare.main("PanACoTA prepare", NCBI_species, NCBI_species_taxid, levels,
                             outdir_prepare, tmp_dir, threads, no_refseq, db_dir, only_mash,
                             info_file, l90, nbcont, cutn, min_dist, max_dist, verbose, quiet)

    # Run annotate module
    list_file = ""
    db_path = ""
    outdir_annotate = os.path.join(outdir, "2-annotate_module")
    import time
    date = time.strftime("%m%y")
    force = False
    qc_only = False
    tmp_dir = ""
    res_annot_dir = None
    small = False

    lstinfo, nbgenomes = annotate.main("PanACoTA annotate", list_file, db_path, outdir_annotate,
                                       name, date, l90, nbcont, cutn, threads, force, qc_only,
                                       info_file, tmp_dir, res_annot_dir, verbose, quiet,
                                       prodigal_only, small)

    # Pangenome step
    name_pan = f"{name}_{nbgenomes}"
    outdir_pan = os.path.join(outdir, "3-pangenome_module")
    clust_mode = 1
    spe_dir = ""
    dbpath = os.path.join(outdir_annotate, "Proteins")
    panfile = pangenome.main("PanACoTA pangenome", lstinfo, name, dbpath, min_id, outdir_pan,
                             clust_mode, spe_dir, threads, outfile=None, verbose=verbose,
                             quiet=quiet)

    # Coregenome step
    outdir_corpers = os.path.join(outdir, "4-corepers_module")
    corepers_file = corepers.main("PanACoTA corepers", panfile, tol, multi, mixed, outdir_corpers,
                                  verbose=verbose, quiet=quiet)

    # Align step
    outdir_align = os.path.join(outdir, "5-align_module")
    force = False
    align_file = align.main("PanACoTA align", corepers_file, lstinfo, name_pan, outdir_annotate,
                            outdir_align, threads, force, verbose, quiet)
    print(align_file)

    # Tree step
    # models_fastme = {"p-distance": "p", "RY-symetric": "Y", "RY": "R",
    #                  "JC69": "J", "K2P": "K", "F81": "1", "F84": "4",
    #                  "TN93": "T", "LogDet": "L"}
    # models_fasttree = {"GTR": "-gtr", "JC": ""}
    # models_iqtree = set(["HKY", "JC", "F81", "K2P", "K3P", "K81uf",
    #                      "TNef", "TIM", "TIMef", "TVM", "TVMef", "SYM", "GTR"])
    # models_iqtree = {mod: mod for mod in models_iqtree}
    if soft == "fasttree":
        model = "-gtr"
    elif soft =="iqtree" or soft == "iqtree2":
        model = "GTR"
    elif soft == "quicktree":
        model = ""
    elif soft == "fastme":
        model = "T"
    else:
        logger.error(f"Soft {soft} is not possible.")
        sys.exit(1)
    outdir_tree = os.path.join(outdir, "6-tree_module")
    tree.main("PanACoTA tree", align_file, outdir_tree, soft, model, threads, boot=False,
              write_boot=False, memory=False, fast=False, verbose=verbose, quiet=quiet)

    logger.info("All modules of PanACOTA are finished.")


def build_parser(parser):
    """
    Method to create a parser for command-line options

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser to configure

    """
    import argparse
    from PanACoTA import utils_argparse

    # Create command-line parser for all options and arguments to give
    general = parser.add_argument_group("General arguments")
    general.add_argument("-o", dest="outdir", required=True,
                          help=("Path to your output folder, where all results "
                                "from all 6 modules will be saved.")
                          )
    general.add_argument("--threads", dest="threads", type=utils_argparse.thread_num, default=1,
                          help="Specify how many threads can be used (default=1)")

    prepare = parser.add_argument_group("'prepare' module arguments")
    prepare.add_argument("-t", dest="NCBI_species_taxid", default="",
                          help=("Species taxid to download, corresponding to the "
                                "'species taxid' provided by the NCBI. A comma-separated "
                                "list of taxid can also be provided.")
                         )
    prepare.add_argument("-s", dest="NCBI_species", default="",
                          help=("Species to download, corresponding to the "
                                "'organism name' provided by the NCBI. Give name between "
                                "quotes (for example \"escherichia coli\")")
                        )
    prepare.add_argument("-l", "--assembly_level", dest="levels", default="",
                          help=("Assembly levels of genomes to download (default: all). "
                                "Possible levels are: 'all', 'complete', 'chromosome', "
                                "'scaffold', 'contig'."
                                "You can also provide a comma-separated list of assembly "
                                "levels. For ex: 'complete,chromosome'")
                          )
    prepare_annote = parser.add_argument_group("Common arguments to 'prepare' "
                                               "and 'annotate' modules")
    prepare_annote.add_argument("--cutn", dest="cutn", type=utils_argparse.positive_int, default=5,
                          help=("By default, each genome will be cut into new contigs when "
                                "at least 5 'N' in a row are found in its sequence. "
                                "If you don't want to "
                                "cut genomes into new contigs when there are rows of 'N', "
                                "put 0 to this option. If you want to cut from a different number "
                                "of 'N' in a row, put this value to this option.")
                          )
    prepare_annote.add_argument("--l90", dest="l90", type=int, default=100,
                                help=("Maximum value of L90 allowed to keep a genome. "
                                      "Default is 100.")
                                )
    prepare_annote.add_argument("--nbcont", dest="nbcont", type=utils_argparse.cont_num,
                                default=999, help=("Maximum number of contigs allowed to "
                                                   "keep a genome. Default is 999."))

    # build_parser(prepare)
    # args_prepare = parser.parse_args(argu)
    # prepare.check_args(prepare_annote, args_prepare)
    # OPTIONS = parse(my_parser, sys.argv[1:])

    annote = parser.add_argument_group("'annotate' module arguments")
    annote.add_argument("--prodigal", dest="prodigal_only", action="store_true", default=False,
                        help="Add this option if you only want syntactical annotation, given "
                             "by prodigal, and not functional annotation requiring prokka and "
                             "is slower.")
    annote.add_argument("-n", dest="name", required=True, type=utils_argparse.gen_name,
                        help=("Choose a name for your annotated genomes. This name should "
                              "contain 4 alphanumeric characters. Generally, they correspond "
                              "to the 2 first letters of genus, and 2 first letters of "
                              "species, e.g. ESCO for Escherichia Coli."))

    pangenome = parser.add_argument_group("'pangenome' module arguments")
    pangenome.add_argument("-i", dest="min_id", type=utils_argparse.perc_id, default=0.8,
                           help=("Minimum sequence identity to be considered in the same "
                                 "cluster (float between 0 and 1). Default is 0.8."))

    corepers = parser.add_argument_group("'corepers' module arguments")
    corepers.add_argument("--tol", dest="tol", type=utils_argparse.percentage, default=1,
                          help=("min %% of genomes having at least 1 member in a family to "
                                "consider the family as persistent (between 0 and 1, "
                                "default is 1 = 100%% of genomes = Core genome)."
                                "By default, the minimum number of genomes will be "
                                "ceil('tol'*N) (N being the total number of genomes). If "
                                "you want to use floor('tol'*N) instead, add the '-F' option."))
    corepers.add_argument("-M", dest="multi", action='store_true',
                          help=("Add this option if you allow several members in any genome "
                                "of a family. By default, only 1 (or 0 if tol<1) member "
                                "per genome are allowed in all genomes. If you want to allow "
                                "exactly 1 member in 'tol'%% of the genomes, and 0, 1 "
                                "or several members in the '1-tol'%% left, use the option -X "
                                "instead of this one: -M and -X options are not compatible."))
    corepers.add_argument("-X", dest="mixed", action='store_true',
                          help="Add this option if you want to allow families having several "
                               "members only in '1-tol'%% of the genomes. In the other genomes, "
                               "only 1 member exactly is allowed. This option is not compatible "
                               "with -M (which is allowing multigenic families: having several "
                               "members in any number of genomes).")

    tree = parser.add_argument_group("'tree' module arguments")
    softs = ["fasttree", "fastme", "quicktree", "iqtree", "iqtree2"]
    tree.add_argument("--soft", dest="soft", choices=softs, default="iqtree",
                      help=("Choose with which software you want to infer the "
                            "phylogenetic tree. Default is IQtree."))

    helper = parser.add_argument_group('Others')
    helper.add_argument("-v", "--verbose", dest="verbose", action="count", default=0,
                        help="Increase verbosity in stdout/stderr.")
    helper.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                        help=("Do not display anything to stdout/stderr. log files will "
                              "still be created."))
    helper.add_argument("-h", "--help", dest="help", action="help",
                        help="show this help message and exit")


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
    import argparse

    args = parser.parse_args(argu)
    return args
    # return check_args(parser, args)


if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    header = '''
     ___                 _____  ___         _____  _____
    (  _`\              (  _  )(  _`\      (_   _)(  _  )
    | |_) )  _ _   ___  | (_) || ( (_)   _   | |  | (_) |
    | ,__/'/'_` )/' _ `\|  _  || |  _  /'_`\ | |  |  _  |
    | |   ( (_| || ( ) || | | || (_( )( (_) )| |  | | | |
    (_)   `\__,_)(_) (_)(_) (_)(____/'`\___/'(_)  (_) (_)


       Large scale comparative genomics tools

     -------------------------------------------
     '''
    my_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                        description=dedent(header), add_help=False)
    build_parser(my_parser)
    OPTIONS = parse(my_parser, sys.argv[1:])
    main_from_parse(OPTIONS)
