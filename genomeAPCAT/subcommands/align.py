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

    Parameters
    ----------
    args : argparse.Namespace
        result of argparse parsing of all arguments in command line
    """
    main(args.corepers, args.list_genomes, args.dataset_name, args.dbpath, args.outdir,
         args.threads, args.force, args.verbose, args.quiet)


def main(corepers, list_genomes, dname, dbpath, outdir, threads, force, verbose=0, quiet=False):
    """
    Align given core genome families

    Parameters
    ----------
    corepers : str
        File containing persistent genome families
    list_genomes : str
        File containing the list of all genomes in the dataset. Only first column is
        considered.
    dname : str
        Dataset name, used to name output files
    dbpath : str
        path to the directory containing 'Proteins' and 'Genes' folders
    outdir : str
        path to the directory where output files must be saved
    threads : int
        Max number of threads to use
    force : bool
        Remove existing output files and rerun everything if True.
    verbose : int
        verbosity:

        - defaut 0 : stdout contains DEBUG and INFO, stderr contains ERROR.
        - 1: stdout contains (DEBUG) and INFO, stderr contains WARNING and ERROR
        - 2: stdout contains (DEBUG), DETAIL and INFO, stderr contains WARNING and ERROR
    quiet : bool
        True if nothing must be sent to stdout/stderr, False otherwise
    """
    # import needed packages
    import logging
    import shutil
    from genomeAPCAT import utils
    from genomeAPCAT.align_module import pan_to_pergenome as p2g
    from genomeAPCAT.align_module import get_seqs as gseqs
    from genomeAPCAT.align_module import alignment as ali
    from genomeAPCAT.align_module import post_align as post

    # test if prokka is installed and in the path
    if not utils.check_installed("fftns"):  # pragma: no cover
        print("fftns (from mafft) is not installed. 'genomeAPCAT align' cannot run.")
        sys.exit(1)

    if force and os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # name logfile, add timestamp if already existing
    logfile_base = os.path.join(outdir, "genomeAPCAT-align_" + dname)
    level = logging.DEBUG
    utils.init_logger(logfile_base, level, '', verbose=verbose, quiet=quiet)
    logger = logging.getLogger()
    all_genomes, aldir, listdir, fam_nums = p2g.get_per_genome(corepers, list_genomes,
                                                               dname, outdir)
    gseqs.get_all_seqs(all_genomes, dname, dbpath, listdir, aldir, fam_nums, quiet)
    prefix = os.path.join(aldir, dname)
    status = ali.align_all_families(prefix, fam_nums, len(all_genomes), dname, quiet, threads)
    if not status:
        logger.error(("At least one alignment did not run well. See detailed log file for "
                      "more information. Program will stop here, alignments won't be "
                      "grouped by genome."))
        sys.exit(1)
    post.post_alignment(fam_nums, all_genomes, prefix, outdir, dname, quiet)
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
    required.add_argument("-c", dest="corepers", required=True,
                          help="Core or persistent genome whose families must be aligned.")
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
                          help="Output directory, where all results must be saved ")

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("--threads", dest="threads", default=1, type=thread_num,
                          help=("add this option if you want to parallelize on several threads. "
                                "Indicate on how many threads you want to parallelize. "
                                "By default, it uses 1 thread. Put 0 if you want to use "
                                "all threads of your computer."))
    optional.add_argument("-F", "--force", dest="force", action="store_true",
                          help=("Force run: Add this option if you want to redo all alignments "
                                "for all families, even if their result file already exists. "
                                "Without this option, if an alignment file already exists, "
                                "it will be used for the next step. If you want to redo only "
                                "a given alignment, just delete its file, without using "
                                "this option."))
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
        the parser used
    argu : [str]
        command-line given by user, to parse using parser

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    return parser.parse_args(argu)


if __name__ == '__main__':
    import argparse

    myparser = argparse.ArgumentParser(description="Align Core/Persistent families",
                                       add_help=False)
    build_parser(myparser)
    OPTIONS = parse(myparser, sys.argv[1:])
    main_from_parse(OPTIONS)
