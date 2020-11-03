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
pangenome is a subcommand of PanACoTA


@author gem
May 2017
"""

import sys
import os


def main_from_parse(args):
    """
    Call main function from the arguments given by parser

    Parameters
    ----------
    args : argparse.Namespace
        result of argparse parsing of all arguments in command line
    """
    cmd = "PanACoTA " + ' '.join(args.argv)
    main(cmd, args.lstinfo_file, args.dataset_name, args.dbpath, args.min_id, args.outdir,
         args.clust_mode, args.spedir, args.threads, args.outfile, args.verbose,
         args.quiet)


def main(cmd, lstinfo, name, dbpath, min_id, outdir, clust_mode, spe_dir, threads, outfile=None,
         verbose=0, quiet=False):
    """
    Main method, doing all steps:

    - concatenate all protein files
    - create database as ffindex
    - cluster all proteins
    - convert to pangenome file
    - creating summary and matrix of pangenome

    Parameters
    ----------
    lstinfo : str
        file with name of genomes to consider for pan in the first column, without extension.
        Other columns are ignored. The first column header must be 'gembase_name'
    name : str
        name given to the dataset. For example, ESCO44 for 44 *Escherichia coli* genomes.
    dbpath : str
        path to the folder containing all protein files (files called as the name of genome
        given in lstinfo + ".prt"
    min_id : float
        Minimum percentage of identity between 2 proteins to put them in the same family
    outdir : str
        path to folder which will contain pangenome results and tmp files
    clust_mode : [0, 1, 2]
        0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'
    spe_dir : str or None
        path to the folder where concatenated bank of proteins must be saved.
        None to use the same folder as protein files
    threads : int
        Max number of threads to use
    outfile : str or None
        Name of the pangenome. None to use the default name
    verbose : int
            verbosity:
        - defaut 0 : stdout contains INFO, stderr contains ERROR.
        - 1: stdout contains INFO, stderr contains WARNING and ERROR
        - 2: stdout contains (DEBUG), DETAIL and INFO, stderr contains WARNING and ERROR
        - >=15: Add DEBUG in stdout
    quiet : bool
        True if nothing must be sent to stdout/stderr, False otherwise
    """
    # import needed packages
    import logging
    from PanACoTA import utils
    from PanACoTA.pangenome_module import protein_seq_functions as protf
    from PanACoTA.pangenome_module import mmseqs_functions as mmf
    from PanACoTA.pangenome_module import post_treatment as pt
    from PanACoTA import __version__ as version

    # test if mmseqs is installed and in the path
    if not utils.check_installed("mmseqs"):  # pragma: no cover
        print("mmseqs is not installed. 'PanACoTA pangenome' cannot run.")
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    # level is the minimum level that will be considered.
    # for verbose = 0 or 1, ignore details and debug, start from info
    if verbose <= 1:
        level = logging.INFO
    # for verbose = 2, ignore only debug
    if verbose >= 2 and verbose < 15:
        level = utils.detail_lvl() # int corresponding to detail level
    # for verbose >= 15, write everything
    if verbose >= 15:
        level = logging.DEBUG
    # name logfile, add timestamp if already existing
    logfile_base = os.path.join(outdir, "PanACoTA-pangenome_" + name)
    utils.init_logger(logfile_base, level, '', verbose=verbose, quiet=quiet, log_details=True)
    logger = logging.getLogger("pangenome")
    logger.info(f'PanACoTA version {version}')
    logger.info("Command used\n \t > " + cmd)

    # Build bank with all proteins to include in the pangenome
    logger.info("build prt bank")
    prt_path = protf.build_prt_bank(lstinfo, dbpath, name, spe_dir, quiet)
    # Do pangenome
    families, panfile = mmf.run_all_pangenome(min_id, clust_mode, outdir,
                                              prt_path, threads, outfile, quiet)
    # Create matrix pan_quali, pan_quanti and summary file
    pt.post_treat(families, panfile)
    logger.info("DONE")
    return panfile


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
    required.add_argument("-l", dest="lstinfo_file", required=True,
                          help=("File containing the list of all genomes to include in "
                                "the pan-genome, 1 genome per line: it can be the "
                                "LSTINFO-<list_file>.lst file of 'PanACoTA annotate' module."
                                "Here, only the first column (genome name without extension) "
                                "will be used. All proteins of all these genomes will be "
                                "concatenated in a file called <dataset_name>.All.prt. The "
                                "column header must be 'gembase_name'."))
    required.add_argument("-n", dest="dataset_name", required=True,
                          help=("Name of the dataset which will be clustered (for example, "
                                "SAEN1234 for 1234 Salmonella enterica genomes). This name will "
                                "be used to name the protein databank, a well as the "
                                "pangenome files."))
    required.add_argument("-d", dest="dbpath", required=True,
                          help=("Path to the folder containing all protein files corresponding "
                                "to the genomes of the dataset (output directory 'Proteins' "
                                "of 'PanACoTA annotate' module)."))
    required.add_argument("-o", dest="outdir", required=True,
                          help=("Output directory, where all results must be saved "
                                "(including tmp folder)"))

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-i", dest="min_id", type=utils_argparse.perc_id, default=0.8,
                          help=("Minimum sequence identity to be considered in the same "
                                "cluster (float between 0 and 1). Default is 0.8."))
    optional.add_argument("-f", dest="outfile",
                          help=("Use this option if you want to give the name of the pangenome "
                                "output file (without path). Otherwise, by default, it is called "
                                "PanGenome-mmseq_<given_dataset_name>.All.prt_<"
                                "information_on_parameters>.lst"))
    optional.add_argument("-c", dest="clust_mode", type=int, choices=[0, 1, 2], default=1,
                          help=("Choose the clustering mode: 0 for 'set cover', 1 for "
                                "'single-linkage', 2 for 'CD-Hit'. Default "
                                "is 'single-linkage' (1)"))
    optional.add_argument("-s", dest="spedir",
                          help=("use this option if you want to save the concatenated protein "
                                "databank in another directory than the one containing all "
                                "individual protein files ('Proteins' folder)."))
    optional.add_argument("--threads", dest="threads", default=1, type=utils_argparse.thread_num,
                          help=("add this option if you want to parallelize on several threads. "
                                "Indicate on how many threads you want to parallelize. "
                                "By default, it uses 1 thread. Put 0 if you want to use "
                                "all threads of your computer."))

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
    Parse arguments given to parser

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to use to parse command-line arguments
    argu : [str]
        command-line given

    Returns
    -------
    argparse.Namespace or None
        The arguments parsed, updated according to some rules. Exit program
        with error message if error occurs with arguments given.
    """
    args = parser.parse_args(argu)
    return args


if __name__ == '__main__':
    import argparse

    my_parser = argparse.ArgumentParser(description="Compute pan-genome", add_help=False)
    build_parser(my_parser)
    OPTIONS = parse(my_parser, sys.argv[1:])
    main_from_parse(OPTIONS)
