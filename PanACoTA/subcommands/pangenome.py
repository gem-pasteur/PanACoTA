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
    main(cmd, args.lstinfo_file, args.dataset_name, args.dbpath, args.method, args.min_id, args.outdir,
         args.clust_mode, args.po_mode, args.evalue, args.conn, args.purity, args.minspec, args.spedir,
         args.threads, args.outfile, args.verbose, args.quiet)


def main(cmd, lstinfo, name, dbpath, method, min_id, outdir, clust_mode, po_mode, evalue, conn, purity, minspec,
         spe_dir, threads, outfile=None, verbose=0, quiet=False):
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
    method : str
        Method to construct the pangenome
    min_id : float
        Minimum percentage of identity between 2 proteins to put them in the same family
    outdir : str
        path to folder which will contain pangenome results and tmp files
    clust_mode : [0, 1, 2]
        0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'
    po_mode : {blastp+|blastn+|tblastx+|diamond|usearch|ublast|lastp|lastn|rapsearch|topaz|blatp|blatn|mmseqsp|mmseqsn}
        blast algorithm which is run on first steps of proteinortho
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
    from PanACoTA.pangenome_module import post_treatment as pt
    from PanACoTA import __version__ as version

    from PanACoTA.pangenome_module.MMseq import MMseq
    from PanACoTA.pangenome_module.ProteinOrtho import ProteinOrtho

    # check if method is valid
    if method not in ["mmseqs", "proteinortho"]:
        print("unknown method. 'PanACoTA pangenome' cannot run.")
        sys.exit(1)

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


    # Do pangenome
    if method == "mmseqs":
        # Build bank with all proteins to include in the pangenome
        prt_path = protf.build_prt_bank(lstinfo, dbpath, name, spe_dir, quiet)
        runner = MMseq(min_id, clust_mode, outdir, prt_path, threads, outfile, quiet)

    if method == "proteinortho":
        prt_path = protf.build_prt_bank(lstinfo, dbpath, name, spe_dir, quiet, dir=True)
        runner = ProteinOrtho(po_mode, evalue, conn, purity, minspec, name, outdir, prt_path, threads, outfile, quiet)

    # test if package required for pangenome building is installed and in the path
    if not runner.check_installed():  # pragma: no cover
        logger.error(f"{method} is not installed. 'PanACoTA pangenome' cannot run.")
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    families, panfile = runner.run()

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
    required.add_argument("-m", dest="method", choices=["mmseqs", "proteinortho"], required=True,
                          help=("Method to construct pangenome."))

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
    optional.add_argument("-p", dest="po_mode", choices="blastp+ blastn+ tblastx+ diamond usearch ublast lastp "
                                                        "lastn rapsearch topaz blatp blatn mmseqsp mmseqsn".split(),
                          default="diamond",
                          help=("For proteinortho algorithm. A blast algorithm which is used on initial steps."))
    optional.add_argument("--eval", dest="evalue", default=1e-5,
                          help=("Proteinortho search option. Evalue threshold for search algorithm."))
    optional.add_argument("--conn", dest="conn", default=0.1,
                          help=("Proteinortho clustering option. "
                                "Minimal required algebraic connectivity for each connected "
                                "component/group. This is the main parameter for the clustering algorithm. "
                                "The higher this value, the more splits are made"))
    optional.add_argument("--purity", dest="purity", default=1e-7,
                          help=("Proteinortho clustering option. Avoid spurious graph assignments, "
                                "the higher the more uncertain edges are cut"))
    optional.add_argument("--minspecies", dest="minspec", default=1,
                          help=("Proteinortho clustering option. "
                                "minimal number of genes per species for each group. If a group is found "
                                "with up to (minspecies) genes/species, it wont be split again "
                                "(regardless of the connectivity). A value of 0 is always satisfied"))
    optional.add_argument("-s", dest="spedir",
                          help=("use this option if you want to save the concatenated protein "
                                "databank in another directory than the one containing all "
                                "individual protein files ('Proteins' folder)."))
    optional.add_argument("--threads", dest="threads", default=1, type=utils_argparse.thread_num,
                          help=("add this option if you want to parallelize on several threads. "
                                "Indicate on how many threads you want to parallelize. "
                                "By default, it uses 1 thread. Put 0 if you want to use "
                                "all threads of your computer."))
    # TODO : добавить параметры поиска

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
