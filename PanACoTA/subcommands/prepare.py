#!/usr/bin/env python3
# coding: utf-8

"""
Subcommand to prepare a dataset:

- Download all genomes of a given species from refseq
- Filter them with L90 and number of contigs thresholds
- Remove too close/far genomes using Mash


@author Amandine PERRIN
August 2019
"""
import os
import sys
import logging

from PanACoTA import utils
from PanACoTA.prepare_module import download_genomes_func as dgf

def main_from_parse(arguments):
    """
    Call main function from the arguments given by parser

    Parameters
    ----------
    arguments : argparse.Namespace
        result of argparse parsing of all arguments in command line

    """
    cmd = "PanACoTA " + ' '.join(arguments.argv)
    main(cmd, arguments.NCBI_species, arguments.NCBI_species_taxid, arguments.outdir,
         arguments.parallel, arguments.no_refseq, arguments.only_mash,
         arguments.min_dist, arguments.verbose, arguments.quiet)


def main(cmd, NCBI_species, NCBI_species_taxid, outdir, threads, no_refseq, only_mash, min_dist,
         verbose, quiet):
    """
    Main method, constructing the draft dataset for the given species

    verbosity:
    - defaut 0 : stdout contains INFO, stderr contains ERROR.
    - 1: stdout contains INFO, stderr contains WARNING and ERROR
    - 2: stdout contains (DEBUG), DETAIL and INFO, stderr contains WARNING and ERROR
    - >=15: Add DEBUG in stdout

    Parameters
    ----------
    cmd : str
        command line used to launch this program
    NCBI_species : str
        name of species to download, as given by NCBI
    NCBI_species_taxid : int
        species taxid given in NCBI
    outdir : str
        path to output directory (where created database will be saved).
    threads : int
        max number of threads to use
    no_refseq : bool
        True if user does not want to download again the database
    only_mash : bool
        True if user user already has the database and quality of each genome (L90, #contigs etc.)
    min_dist : int
        lower limit of distance between 2 genomes to keep them
    verbose : int
        verbosity:
        default (0): info in stdout, error and more in stderr
        1 = add warnings in stderr
        2 = like 1 + add DETAIL to stdout (by default only INFO)
        >15: add debug to stdout
    quiet : bool
        True if nothing must be sent to stdout/stderr, False otherwise

    """
    # Fixed limits. For now, we do not propose to user to give its own limits
    max_dist = 0.06
    max_l90 = 100
    max_cont = 999

    # get species name in NCBI format
    if NCBI_species:
        NCBI_species = NCBI_species.capitalize()
        species_linked = "_".join(NCBI_species.split())
        species_linked = "_".join(species_linked.split("/"))

    # if species name not given by user, use taxID instead to name files
    else:
        species_linked = str(NCBI_species_taxid)
    outdir = os.path.join(outdir, species_linked)
    db_dir = None
    os.makedirs(outdir, exist_ok=True)

    # Initialize logger
    # set level of logger: level is the minimum level that will be considered.
    if verbose <= 1:
        level = logging.INFO
    # for verbose = 2, ignore only debug
    if verbose >= 2 and verbose < 15:
        level = utils.detail_lvl() # int corresponding to detail level
    # for verbose >= 15, write everything
    if verbose >= 15:
        level = logging.DEBUG
    logfile_base = os.path.join(outdir, "PanACoTA_prepare_{}").format(species_linked)
    logfile_base = utils.init_logger(logfile_base, level, name='', details=True,
                                     verbose=verbose, quiet=quiet)

    # Start prepare step
    logger = logging.getLogger('')
    logger.info("Command used\n \t > " + cmd)
    logger.info("'PanACoTA prepare' Will run on {} cores".format(threads))
    sys.exit(1)

    # Run more than only mash filter
    if not only_mash:
        if norefseq:   # Do not download genomes, just do QC and mash filter on given genomes
            logger.info('You asked to skip refseq downloads.')

        else:  # Do all steps: download, QC, mash filter
            sum_file = dgf.download_summary(species_linked, outdir)
            db_dir = dgf.download_from_refseq(sum_file, NCBI_species, NCBI_taxid, outdir, threads)

        # If refseq genomes already downloaded but not in the database_init folder,
        # put them in it.
        if not db_dir:  # if norefseq and not already created before, db_dir can be missing
            db_dir = os.path.join(outdir, "Database_init")
            if norefseq and not os.path.exists(db_dir):
                # add genomes from refseq/bacteria folder to Database_init
                nb_gen, _ = dgf.to_database(outdir)
                logger.info("{} refseq genomes downloaded".format(nb_gen))
        genomes = fg.check_quality(outdir, species_linked, db_dir, max_l90, max_cont)
    # # Do only mash filter. Genomes must be alreday downloaded, and there must be a file with
    # # all information on these genomes (L90 etc.)
    else:
        info_file = os.path.join(outdir, "info-genomes-list-{}.lst".format(species_linked))
        if not os.path.exists(info_file):  # info-file missing -> error and exit
            logger.error(("You do not have the file called {} with all information about "
                          "genomes. Provide it with the right name, or remove the '--mash' "
                          "option to rerun quality control.".format(info_file)))
            sys.exit(1)
        logger.info(("You want to rerun only mash steps. Getting information "
                     "from {}").format(info_file))
        genomes = utils.get_info_genomes(info_file, species_linked)

    # # Run Mash
    sorted_genomes = fg.sort_genomes_minhash(genomes, max_l90, max_cont)
    removed = fg.iterative_mash(sorted_genomes, genomes, outdir, species_linked,
                                 min_dist, max_dist, threads)

    # Write list of genomes kept, and list of genomes removed
    fg.write_outputfiles(sorted_genomes, removed, outdir, species_linked, min_dist)
    logger.info("End")


def build_parser(parser):
    """
    Method to create a parser for command-line options

    Parameters
    ----------
    parser : argparse.ArgumentParser
        parser to configure in order to extract command-line arguments
    """
    import multiprocessing
    import argparse
    from PanACoTA import utils_argparse

    required = parser.add_argument_group('Required arguments')
    required.add_argument("-t", dest="NCBI_species_taxid", required=True,
                        help=("Species taxid to download, corresponding to the "
                               "'species taxid' provided by the NCBI")
                        )

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-s", dest="NCBI_species",
                        help=("Species to download, corresponding to the "
                               "'organism name' provided by the NCBI. Give name given between "
                               "quotes (for example \"escherichia coli\")")
                        )
    optional.add_argument("-o", dest="outdir", default=".",
                        help=("Give the path to the directory where you want to save the "
                              "database. In the given diretory, it will create a folder with "
                              "the gembase species name. Inside this folder, you will find a "
                              "folder 'Database_init' containing all fasta files, as well as a "
                              "folder 'refseq' with files downloaded by ncbi_genome_download."))
    optional.add_argument("-p", dest="parallel", type=utils_argparse.thread_num, default=1,
                        help=("Run 'N' downloads in parallel (default=1). Put 0 if "
                              "you want to use all cores of your computer."))
    optional.add_argument("-r", dest="no_refseq", action="store_true",
                        help=("If you already downloaded refseq genomes and do not want to "
                              "check them, add this option to directly go to the next steps:"
                              "quality control (L90, number of contigs...) and mash filter."))
    optional.add_argument('-M', '--only-mash', dest="only_mash", action="store_true",
                        help=("Add this option if you already downloaded complete and refseq "
                              "genomes, and ran quality control (you have, in your result "
                              "folder, a file called 'info-genomes-list-<gembase_species>.lst', "
                              "contaning all genome names, as well as their genome size, "
                              "number of contigs and L90 values). "
                              "It will then get information on genomes quality from this "
                              "file, and run mash steps."))
    optional.add_argument("-m", dest="min_dist", default=1e-4, type=float,
                        help="By default, genomes whose distance to the reference is not "
                             "between 1e-4 and 0.06 are discarded. You can specify your own "
                             "lower limit (instead of 1e-4) with this option.")
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
    return check_args(parser, args)


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
    from termcolor import colored
    if not args.NCBI_species:
        print(colored("WARNING: you did not provide a species name. All files will "
                      "be downloaded in a folder called with the NCBI species taxid, and "
                      "you won't be able to get the summary file (showing a summary of all "
                      "strains downloaded)", "yellow"))
        yes = ["y", "yes", "Y"]
        no = ["n", "no", "N"]
        while True:
            print("Do you still want to continue? [Y/n] ", end ="")
            choice = input().lower()
            if choice in no:
                sys.exit(0)
            elif choice in yes or choice == '':
                return args
            else:
                sys.stdout.write("Please answer with 'y' or 'n': ")
    return True


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
