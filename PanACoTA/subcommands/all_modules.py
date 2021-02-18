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
from PanACoTA import utils_argparse
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

        verbosity:

    - defaut 0 : stdout contains INFO, stderr contains ERROR.
    - 1: stdout contains INFO, stderr contains WARNING and ERROR
    - 2: stdout contains (DEBUG), DETAIL and INFO, stderr contains WARNING and ERROR
    - >=15: Add DEBUG in stdout

    Parameters
    ----------
    args_all : tuple
        arguments common to all modules: output directory (str),
        threads (int), verbose (int), quiet (bool)
    args_prepare : tuple
        arguments for prepare module (see subcommands.prepare.py): NCBI_species_taxid (int),
        NCBI_species (str), levels (str), tmp_dir (str), norefseq (bool), db_dir (str),
        only_mash (bool), info_file (str), l90 (int), nbcont (int), cutn (int),
        min_dist (float), max_dist (float)
    args_annot : tuple
        arguments for annotate module (see subcommands/annotate.py): name (str), qc_only (bool),
        date (str), prodigal_only (bool), small (bool)
    args_pan : tuple
        arguments for pangenome module (see subcommands/pangenome.py): min_id (float),
        clust_mode (int), spe_dir (str), outfile (str)
    args_corepers : tuple
        arguments for corepers module (see subcommands.corepers.py): tol (float), mixed (bool),
        multi (bool), floor (bool)
    args_tree : tuple
        arguments for tree module (see subcommands.tree.py): soft (str), model (str), boot (bool),
        write_boot (bool), memory (str), fast (bool)
    """
    cmd = "PanACoTA " + ' '.join(args.argv)
    args_all = (args.outdir, args.threads, args.verbose, args.quiet)
    args_prepare = (args.ncbi_species_taxid, args.ncbi_species_name, args.ncbi_taxid, args.levels,
                    args.ncbi_section, args.tmp_dir, args.norefseq, args.db_dir, args.only_mash, 
                    args.info_file, args.l90, args.nbcont, args.cutn, args.min_dist, args.max_dist)
    args_annot = (args.name, args.qc_only, args.date, args.prodigal_only, args.small)
    args_pan = (args.min_id, args.clust_mode, args.spedir, args.outfile)
    args_cp = (args.tol, args.mixed, args.multi, args.floor)
    args_tree = (args.soft, args.model, args.boot, args.write_boot, args.memory, args.fast)
    main(cmd, args_all, args_prepare, args_annot, args_pan, args_cp, args_tree)


def main(cmd, args_all, args_prepare, args_annot, args_pan, args_corepers, args_tree):
    """
    Call all modules, one by one, using output of one as input for the next one

    Parameters
    ----------
    cmd : str
        command line used to launch the program
    args_all : tuple
        arguments common to all modules: output directory (str),
        threads (int), verbose (int), quiet (bool)
    args_prepare : tuple
        arguments for prepare module (see subcommands.prepare.py): NCBI_species_taxid (int),
        NCBI_species_name (str), NCBI_taxid (int), levels (str), NCBI_section (str),
        tmp_dir (str), norefseq (bool), db_dir (str),
        only_mash (bool), info_file (str), l90 (int), nbcont (int), cutn (int),
        min_dist (float), max_dist (float)
    args_annot : tuple
        arguments for annotate module (see subcommands/annotate.py): name (str), qc_only (bool),
        date (str), prodigal_only (bool), small (bool)
    args_pan : tuple
        arguments for pangenome module (see subcommands/pangenome.py): min_id (float),
        clust_mode (int), spe_dir (str), outfile (str)
    args_corepers : tuple
        arguments for corepers module (see subcommands.corepers.py): tol (float), mixed (bool),
        multi (bool), floor (bool)
    args_tree : tuple
        arguments for tree module (see subcommands.tree.py): soft (str), model (str), boot (bool),
        write_boot (bool), memory (str), fast (bool)
    """
    outdir, threads, verbose, quiet = args_all
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
    (NCBI_species_taxid, NCBI_species_name, NCBI_taxid, levels, NCBI_section,
     tmp_dir, norefseq, db_dir, only_mash, info_file,
     l90, nbcont, cutn, min_dist, max_dist) = args_prepare
    logger.info("prepare step")
    info_file = prepare.main("PanACoTA prepare", NCBI_species_name, NCBI_species_taxid,
                             NCBI_taxid, levels, NCBI_section,
                             outdir_prepare, tmp_dir, threads, norefseq, db_dir, only_mash,
                             info_file, l90, nbcont, cutn, min_dist, max_dist, verbose, quiet)

    # Run annotate module
    list_file = ""
    db_path = ""
    tmp_dir = ""
    force = False
    outdir_annotate = os.path.join(outdir, "2-annotate_module")
    (name, qc_only, date, prodigal_only, small) = args_annot
    res_annot_dir = None

    logger.info("annotate step")
    lstinfo, nbgenomes = annotate.main("PanACoTA annotate", list_file, db_path, outdir_annotate,
                                       name, date, l90, nbcont, cutn, threads, force, qc_only,
                                       info_file, tmp_dir, res_annot_dir, verbose, quiet,
                                       prodigal_only=prodigal_only, small=small)
    if qc_only:
        return "QC_only done"

    # Pangenome step
    name_pan = f"{name}_{nbgenomes}"
    outdir_pan = os.path.join(outdir, "3-pangenome_module")
    dbpath = os.path.join(outdir_annotate, "Proteins")
    (min_id, clust_mode, spe_dir, outfile) = args_pan
    logger.info("pangenome step")
    panfile = pangenome.main("PanACoTA pangenome", lstinfo, name_pan, dbpath, min_id, outdir_pan,
                             clust_mode, spe_dir, threads, outfile, verbose=verbose,
                             quiet=quiet)

    # Coregenome step
    outdir_corpers = os.path.join(outdir, "4-corepers_module")
    logger.info("corepers step")
    (tol, mixed, multi, floor) = args_corepers
    corepers_file = corepers.main("PanACoTA corepers", panfile, tol, multi, mixed, outdir_corpers,
                                  floor, verbose, quiet)
    # Align step
    outdir_align = os.path.join(outdir, "5-align_module")
    force = False
    logger.info("align step")
    align_file = align.main("PanACoTA align", corepers_file, lstinfo, name_pan, outdir_annotate,
                            outdir_align, threads, force, verbose=verbose, quiet=quiet)


    # Tree step
    (soft, model, boot, write_boot, memory, fast) = args_tree
    outdir_tree = os.path.join(outdir, "6-tree_module")
    logger.info("tree step")
    tree.main("PanACoTA tree", align_file, outdir_tree, soft, model, threads, boot,
              write_boot, memory, fast, verbose=verbose, quiet=quiet)
    logger.info("All modules of PanACOTA are finished.")
    return 0


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
    general.add_argument("-c", dest="configfile",
                          help=("Path to your configuration file, defining values of parameters.")
                          )
    general.add_argument("-o", dest="outdir", required=True,
                          help=("Path to your output folder, where all results "
                                "from all 6 modules will be saved.")
                          )
    general.add_argument("--threads", dest="threads", type=utils_argparse.thread_num,
                          help="Specify how many threads can be used (default=1)")
    # prepare arguments
    pprepare = parser.add_argument_group("'prepare' module arguments")
    pprepare.add_argument("-t", dest="ncbi_species_taxid",
                          help=("Species taxid to download, corresponding to the "
                                "'species taxid' provided by the NCBI. A comma-separated "
                                "list of taxid can also be provided.")
                         )
    pprepare.add_argument("-s", dest="ncbi_species",
                          help=("Species to download, corresponding to the "
                                "'organism name' provided by the NCBI. Give name between "
                                "quotes (for example \"escherichia coli\")")
                        )
    pprepare.add_argument("-l", "--assembly_level", dest="levels",
                          help=("Assembly levels of genomes to download (default: all). "
                                "Possible levels are: 'all', 'complete', 'chromosome', "
                                "'scaffold', 'contig'."
                                "You can also provide a comma-separated list of assembly "
                                "levels. For ex: 'complete,chromosome'")
                          )
    pprepare_annote = parser.add_argument_group("Common arguments to 'prepare' "
                                               "and 'annotate' modules")
    pprepare_annote.add_argument("--cutn", dest="cutn", type=utils_argparse.positive_int,

                          help=("By default, each genome will be cut into new contigs when "
                                "at least 5 'N' in a row are found in its sequence. "
                                "If you don't want to "
                                "cut genomes into new contigs when there are rows of 'N', "
                                "put 0 to this option. If you want to cut from a different number "
                                "of 'N' in a row, put this value to this option.")
                          )
    pprepare_annote.add_argument("--l90", dest="l90", type=int,
                                help=("Maximum value of L90 allowed to keep a genome. "
                                      "Default is 100.")
                                )
    pprepare_annote.add_argument("--nbcont", dest="nbcont", type=utils_argparse.cont_num,
                                help=("Maximum number of contigs allowed to "
                                                   "keep a genome. Default is 999."))

    pannote = parser.add_argument_group("'annotate' module arguments")
    pannote.add_argument("--prodigal", dest="prodigal_only", action="store_true",
                        help="Add this option if you only want syntactical annotation, given "
                             "by prodigal, and not functional annotation requiring prokka and "
                             "is slower.")
    pannote.add_argument("-n", dest="name", required=True, type=utils_argparse.gen_name,
                        help=("Choose a name for your annotated genomes. This name should "
                              "contain 4 alphanumeric characters. Generally, they correspond "
                              "to the 2 first letters of genus, and 2 first letters of "
                              "species, e.g. ESCO for Escherichia Coli."))

    ppangenome = parser.add_argument_group("'pangenome' module arguments")
    ppangenome.add_argument("-i", dest="min_id", type=utils_argparse.perc_id,
                           help=("Minimum sequence identity to be considered in the same "
                                 "cluster (float between 0 and 1). Default is 0.8."))

    pcorepers = parser.add_argument_group("'corepers' module arguments")
    pcorepers.add_argument("--tol", dest="tol", type=utils_argparse.percentage,
                          help=("min %% of genomes having at least 1 member in a family to "
                                "consider the family as persistent (between 0 and 1, "
                                "default is 1 = 100%% of genomes = Core genome)."
                                "By default, the minimum number of genomes will be "
                                "ceil('tol'*N) (N being the total number of genomes). If "
                                "you want to use floor('tol'*N) instead, add the '-F' option."))
    pcorepers.add_argument("-Mu", dest="multi", action='store_true',
                          help=("Add this option if you allow several members in any genome "
                                "of a family. By default, only 1 (or 0 if tol<1) member "
                                "per genome are allowed in all genomes. If you want to allow "
                                "exactly 1 member in 'tol'%% of the genomes, and 0, 1 "
                                "or several members in the '1-tol'%% left, use the option -X "
                                "instead of this one: -M and -X options are not compatible."))
    pcorepers.add_argument("-X", dest="mixed", action='store_true',
                          help="Add this option if you want to allow families having several "
                               "members only in '1-tol'%% of the genomes. In the other genomes, "
                               "only 1 member exactly is allowed. This option is not compatible "
                               "with -M (which is allowing multigenic families: having several "
                               "members in any number of genomes).")

    ptree = parser.add_argument_group("'tree' module arguments")
    softs = ["fasttree", "fastme", "quicktree", "iqtree", "iqtree2"]
    ptree.add_argument("--soft", dest="soft", choices=softs,
                      help=("Choose with which software you want to infer the "
                            "phylogenetic tree. Default is IQtree."))

    helper = parser.add_argument_group('Others')
    helper.add_argument("-v", "--verbose", dest="verbose", action="count",
                        help="Increase verbosity in stdout/stderr.")
    helper.add_argument("-q", "--quiet", dest="quiet", action="store_true",
                        help=("Do not display anything to stdout/stderr. log files will "
                              "still be created."))
    helper.add_argument("-h", "--help", dest="help", action="help",
                        help="show this help message and exit")


def parse(parser, argu):
    """
    parse arguments given to parser

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
    return check_args(parser, args)


def check_args(parser, argv):
    """
    Check arguments given by user

    Parameters
    ----------
    parser : argparse.ArgumentParser
        the parser used
    argv : argparse.Namespace
        parsed command-line given by user

    Returns
    -------
    argv : argparse.Namespace
        parsed command-line 
    """
    from argparse import Namespace

    # Get arguments given in command line
    # If an argument was not in the user command-line, it still exists in argv,
    # but is set to 'None' -> ignore it
    # If a bool argument (quiet, mixed, multi, prodigal_only) is not in the user command-line,
    # it still exists in argv but is set to False -> ignore it to replace by info from
    # config file if any, or default otherwise
    dict_argv = {key:val for key,val in vars(argv).items() if val is not None and val != False}
    final_dict = {}

    # PREPARE STEP
    prep_dict = get_prepare(dict_argv)
    a = Namespace(**prep_dict)
    prepare.check_args(parser, a)
    final_dict.update(prep_dict)

    # ANNOTATE STEP
    annot_dict = get_annotate(dict_argv)
    a = Namespace(**annot_dict)
    annotate.check_args(parser, a)
    final_dict.update(annot_dict)

    # PANGENOME STEP
    pan_dict = get_pangenome(dict_argv)
    final_dict.update(pan_dict)
    # Add default arguments if not found in commandline nor config file

    # COREPERS STEP
    cp_dict = get_corepers(dict_argv)
    a = Namespace(**cp_dict)
    corepers.check_args(parser, a)
    final_dict.update(cp_dict)

    # TREE STEP
    tree_dict = get_tree(dict_argv)
    # If we chose soft==quicktree, check_args will return error if threads > 1. So,
    # we put threads = 1 during check_args, and put back the value after
    th_save = tree_dict["threads"]
    a = Namespace(**tree_dict)
    checked_a = tree.check_args(parser, a)
    # Variable can be changed by check_args (like put default model if not given)
    tree_dict.update(vars(checked_a))
    tree_dict["threads"] = th_save
    final_dict.update(tree_dict)

    # Combine files
    # Put new args values in argv
    for key, val in final_dict.items():
        vars(argv)[key] = val
    return argv


def get_prepare(dict_argv):
    """
    Check that arguments given for prepare step are compatible
    """
    # Get arguments from config file and add them (overwrite if needed)
    if not "configfile" in dict_argv:
        conf_conffile = utils_argparse.Conf_all_parser("",["prepare"])
    else:
        conf_conffile = utils_argparse.Conf_all_parser(dict_argv['configfile'],
                                                       readsec=["prepare"])
    # Add arguments from commandline
    conf_conffile.update(dict_argv, "prepare")
    # Add default arguments if not found in comd line nor config file
    defaults = {"verbose": 0, "threads": 1, "cutn": 5, "l90": 100, "nbcont":999,
                "levels": "all", "quiet": False, "ncbi_species_name": "",
                "ncbi_species_taxid": "", "ncbi_taxid": "", "tmp_dir": "", "db_dir": "",
                "info_file": "", "min_dist": 1e-4, "max_dist": 0.06,
                "norefseq": False, "only_mash": False, "ncbi_section": "refseq"}
    conf_conffile.add_default(defaults, "prepare")
    # Change to expected types (boolean, int, float)
    conf_conffile.set_boolean("prepare", "quiet")
    conf_conffile.set_boolean("prepare", "only_mash")
    conf_conffile.set_boolean("prepare", "norefseq")
    conf_conffile.set_int("prepare", "threads")
    conf_conffile.set_int("prepare", "verbose")
    conf_conffile.set_int("prepare", "cutn")
    conf_conffile.set_int("prepare", "l90")
    conf_conffile.set_int("prepare", "nbcont")
    conf_conffile.set_float("prepare", "min_dist")
    conf_conffile.set_float("prepare", "max_dist")
    prep_dict = conf_conffile.get_section_dict("prepare")
    return prep_dict


def get_annotate(dict_argv):
    """
    Check that arguments given for annotate step are compatible
    """
    # Get arguments from config file and add them (overwrite if needed)
    if not "configfile" in dict_argv:
        conf_conffile = utils_argparse.Conf_all_parser("",["annotate"])
    else:
        conf_conffile = utils_argparse.Conf_all_parser(dict_argv['configfile'],
                                                       readsec=["annotate"])
    # Add arguments from commandline
    conf_conffile.update(dict_argv, "annotate")
    if "date" not in dict(conf_conffile["annotate"]):
        import time
        date = time.strftime("%m%y")
        conf_conffile.update({"date": date}, "annotate")
    # Add default arguments if not found in commandline nor config file
    defaults = {"verbose": 0, "threads": 1, "cutn": 5, "l90": 100, "nbcont":999,
                "quiet": False, "prodigal_only": False, "small": False, "qc_only": False,
                "list_file": "", "db_path": "", "from_info": True}
    conf_conffile.add_default(defaults, "annotate")
    conf_conffile.set_boolean("annotate", "quiet")
    conf_conffile.set_boolean("annotate", "prodigal_only")
    conf_conffile.set_boolean("annotate", "small")
    conf_conffile.set_boolean("annotate", "qc_only")
    conf_conffile.set_int("annotate", "verbose")
    conf_conffile.set_int("annotate", "threads")
    conf_conffile.set_int("annotate", "cutn")
    conf_conffile.set_int("annotate", "nbcont")
    conf_conffile.set_int("annotate", "l90")
    annot_dict = conf_conffile.get_section_dict("annotate")
    return annot_dict


def get_pangenome(dict_argv):
    """
    Check that arguments given for pangenome step are compatible
    """
    # Get arguments from config file and add them (overwrite if needed)
    if not "configfile" in dict_argv:
        conf_conffile = utils_argparse.Conf_all_parser("",["pangenome"])
    else:
        conf_conffile = utils_argparse.Conf_all_parser(dict_argv['configfile'],
                                                       readsec=["pangenome"])
    # Add arguments from commandline
    conf_conffile.update(dict_argv, "pangenome")
    # Add default arguments if not found in commandline nor config file
    defaults = {"verbose": 0, "threads": 1, "min_id": 0.8, "quiet": False, "clust_mode": 1,
                "outfile": "", "spedir": ""}
    conf_conffile.add_default(defaults, "pangenome")
    conf_conffile.set_boolean("pangenome", "quiet")
    conf_conffile.set_int("pangenome", "verbose")
    conf_conffile.set_int("pangenome", "threads")
    conf_conffile.set_float("pangenome", "min_id")
    pan_dict = conf_conffile.get_section_dict("pangenome")
    return pan_dict


def get_corepers(dict_argv):
    """
    Check that arguments given for corepers step are compatible
    """
    # Get arguments from config file and add them (overwrite if needed)
    if not "configfile" in dict_argv:
        conf_conffile = utils_argparse.Conf_all_parser("",["corepers"])
    else:
        conf_conffile = utils_argparse.Conf_all_parser(dict_argv['configfile'],
                                                       readsec=["corepers"])
    # Add arguments from commandline
    conf_conffile.update(dict_argv, "corepers")
    # Add default arguments if not found in commandline nor config file
    defaults = {"verbose": 0, "quiet": False, "tol": 1, "mixed": False, "multi": False,
                "floor": False, "threads": 1}
    conf_conffile.add_default(defaults, "corepers")
    conf_conffile.set_boolean("corepers", "quiet")
    conf_conffile.set_boolean("corepers", "floor")
    conf_conffile.set_boolean("corepers", "mixed")
    conf_conffile.set_boolean("corepers", "multi")
    conf_conffile.set_int("corepers", "verbose")
    conf_conffile.set_int("corepers", "tol")
    conf_conffile.set_int("corepers", "threads")
    corepers_dict = conf_conffile.get_section_dict("corepers")
    return corepers_dict


def get_tree(dict_argv):
    """
    Check that arguments given for tree step are compatible
    """
    # Get arguments from config file and add them (overwrite if needed)
    if not "configfile" in dict_argv:
        conf_conffile = utils_argparse.Conf_all_parser("",["tree"])
    else:
        conf_conffile = utils_argparse.Conf_all_parser(dict_argv['configfile'], readsec=["tree"])
    # Add arguments from commandline
    conf_conffile.update(dict_argv, "tree")
    # Add default arguments if not found in commandline nor config file
    defaults = {"verbose": 0, "quiet": False, "threads": 1,
                "soft": "iqtree", "model": None, "boot": False, "write_boot": False,
                "memory": None, "fast": False}
    conf_conffile.add_default(defaults, "tree")
    conf_conffile.set_boolean("tree", "quiet")
    conf_conffile.set_boolean("tree", "fast")
    conf_conffile.set_boolean("tree", "boot")
    conf_conffile.set_boolean("tree", "write_boot")
    conf_conffile.set_int("tree", "verbose")
    conf_conffile.set_int("tree", "threads")
    tree_dict = conf_conffile.get_section_dict("tree")
    return tree_dict


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
