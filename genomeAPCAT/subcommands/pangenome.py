#!/usr/bin/env python3
# coding: utf-8

"""
pan-genome is a subcommand of genomeAPCAT


@author gem
May 2017
"""

import sys
import os


def main_from_parse(args):
    """
    Call main function from the arguments given by parser
    """
    main(args.lstinfo_file, args.dataset_name, args.dbpath, args.min_id, args.outdir,
         args.clust_mode, args.spedir, args.threads, args.outfile)


def main(lstinfo, name, dbpath, min_id, outdir, clust_mode, spe_dir, threads, outfile=None):
    """
    Main method, doing all steps:
    - concatenate all protein files
    - create database as ffindex
    - cluster all proteins
    - convert to pangenome file
    - summary/matrix?
    """
    # import needed packages
    import logging
    from genomeAPCAT import utils
    from genomeAPCAT.pangenome_module import protein_seq_functions as protf
    from genomeAPCAT.pangenome_module import mmseqs_functions as mmf

    os.makedirs(outdir, exist_ok=True)
    # name logfile, add timestamp if already existing
    logfile = os.path.join(outdir, "genomeAPCAT-pangenome_" + name + ".log")
    # if os.path.isfile(logfile):
    #     import time
    #     logfile = os.path.splitext(logfile)[0] + time.strftime("_%y-%m-%d_%H-%m-%S.log")
    # set level of logger (here debug to show everything during development)
    level = logging.DEBUG
    utils.init_logger(logfile, level)
    logger = logging.getLogger()

    # test if mmseqs is installed and in the path
    if not utils.check_installed("mmseqs"):  # pragma: no cover
        logger.error("mmseqs is not installed. 'genomeAPCAT pangenome' cannot run.")
        sys.exit(1)

    # Build bank with all proteins to include in the pangenome
    prt_path = protf.build_prt_bank(lstinfo, dbpath, name, spe_dir)
    # Do pangenome
    mmf.run_all_pangenome(min_id, clust_mode, outdir, prt_path, threads, outfile)


def build_parser(parser):
    """
    Method to create a parser for command-line options
    """
    import argparse
    import multiprocessing
    def perc_id(param):
        try:
            param = float(param)
        except Exception:
            msg = "argument -i percentage_id: invalid float value: {}".format(param)
            raise argparse.ArgumentTypeError(msg)
        if param < 0 or param > 1:
            msg = ("The minimum %% of identity must be in [0, 1]. Invalid value: {}".format(param))
            raise argparse.ArgumentTypeError(msg)
        return param
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
            return multiprocessing.cpu_count()
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
    parser.add_argument("-f", dest="outfile",
                        help=("Use this option if you want to give the name of the pangenome "
                              "output file (without path). Otherwise, by default, it is called "
                              "PanGenome-mmseq_<given_dataset_name>.All.prt_<"
                              "information_on_parameters>.lst"))
    parser.add_argument("-c", dest="clust_mode", type=int, choices=[0, 1, 2], default=1,
                        help=("Choose the clustering mode: 0 for 'set cover', 1 for "
                              "'single-linkage', 2 for 'CD-Hit'. Default "
                              "is 'single-linkage' (1)"))
    parser.add_argument("-s", dest="spedir",
                        help=("use this option if you want to save the concatenated protein "
                              "databank in another directory than the one containing all "
                              "individual protein files ('Proteins' folder)."))
    parser.add_argument("--threads", dest="threads", default=1, type=thread_num,
                        help=("add this option if you want to parallelize on several threads. "
                              "Indicate on how many threads you want to parallelize. "
                              "By default, it uses 1 thread. Put 0 if you want to use "
                              "all threads of your computer."))


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
