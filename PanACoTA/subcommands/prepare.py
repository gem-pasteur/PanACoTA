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
from termcolor import colored

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
         arguments.parallel, arguments.no_refseq, arguments.only_mash, arguments.l90,
         arguments.nbcont, arguments.cutn, arguments.min_dist, arguments.verbose, arguments.quiet)



def main(cmd, NCBI_species, NCBI_taxid, outdir, threads, no_refseq, only_mash, l90,
         nbcont, cutn, min_dist, verbose, quiet):
    """
    Main method, constructing the draft dataset for the given species

    verbosity:
    - defaut 0 : stdout contains INFO, stderr contains ERROR, .log contains INFO and more, .log.err contains warning and more
    - 1: same as 0 + WARNING in stderr
    - 2: same as 1 + DETAILS in stdout + DETAILS in .log.details
    - >=15: same as 2 + Add DEBUG in stdout + create .log.debug with everything from info to debug


    Parameters
    ----------
    cmd : str
        command line used to launch this program
    NCBI_species : str
        name of species to download, as given by NCBI
    NCBI_taxid : int
        species taxid given in NCBI
    outdir : str
        path to output directory (where created database will be saved).
    threads : int
        max number of threads to use
    no_refseq : bool
        True if user does not want to download again the database
    only_mash : bool
        True if user user already has the database and quality of each genome (L90, #contigs etc.)
    l90 : int
        Max L90 allowed to keep a genome
    nbcont : int
        Max number of contigs allowed to keep a genome
    cutn : int
        cut at each stretch of this number of 'N'. Don't cut if equal to 0
    min_dist : int
        lower limit of distance between 2 genomes to keep them
    verbose : int
        verbosity:
        - defaut 0 : stdout contains INFO, stderr contains ERROR, .log contains INFO and more, .log.err contains warning and more
        - 1: same as 0 + WARNING in stderr
        - 2: same as 1 + DETAILS in stdout + DETAILS in .log.details
        - >=15: same as 2 + Add DEBUG in stdout + create .log.debug with everything from info to debug
    quiet : bool
        True if nothing must be sent to stdout/stderr, False otherwise

    """
    # Fixed limits. For now, we do not propose to user to give its own limits
    max_dist = 0.06


    # get species name in NCBI format
    # -> will be used to name output directory
    # -> will be used to download summary file if given species corresponds to NCBI name
    if NCBI_species:
        NCBI_species = NCBI_species.capitalize()
        species_linked = "_".join(NCBI_species.split())
        species_linked = "_".join(species_linked.split("/"))

    # if species name not given by user, use taxID instead to name output directory
    else:
        species_linked = str(NCBI_taxid)
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
    logfile_base, logger = utils.init_logger(logfile_base, level, 'prepare', details=True,
                                             verbose=verbose, quiet=quiet)

    # Message on what will be done (cmd, cores used)
    logger.info("Command used\n \t > " + cmd)
    message = f"'PanACoTA prepare' will run on {threads} "
    message += f"cores" if threads>1 else "core"
    logger.info(message)

    # Start prepare step
    # Run more than only mash filter (!only_mash):
    # - start from QC and mash (norefseq)
    # - start from genome download (!norefseq))
    if not only_mash:
        if no_refseq:   # Do not download genomes, just do QC and mash filter on given genomes
            logger.info('You asked to skip refseq downloads.')
        else:  # Do all steps: download, QC, mash filter
            # Download all genomes of the given taxID
            db_dir = dgf.download_from_refseq(species_linked, NCBI_species, NCBI_taxid,
                                              outdir, threads)
            sys.exit(1)

        # If refseq genomes already downloaded but not in the database_init folder,
        # put them in it.
        if not db_dir:  # if norefseq and not already created before, db_dir can be missing
            db_dir = os.path.join(outdir, "Database_init")
            if norefseq and not os.path.exists(db_dir):
                # add genomes from refseq/bacteria folder to Database_init
                nb_gen, _ = dgf.to_database(outdir)
                logger.info("{} refseq genomes downloaded".format(nb_gen))
        genomes = fg.check_quality(outdir, species_linked, db_dir, max_l90, max_cont)
    # # Do only mash filter. Genomes must be already downloaded, and there must be a file with
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
    import argparse
    from PanACoTA import utils_argparse


    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-t", dest="NCBI_species_taxid",
                          help=("Species taxid to download, corresponding to the "
                                "'species taxid' provided by the NCBI")
                         )
    optional.add_argument("-s", dest="NCBI_species",
                          help=("Species to download, corresponding to the "
                                "'organism name' provided by the NCBI. Give name between "
                                "quotes (for example \"escherichia coli\")")
                        )
    optional.add_argument("-o", dest="outdir", default=".",
                          help=("Give the path to the directory where you want to save the "
                               "database. In the given directory, it will create a folder with "
                               "the gembase species name. Inside this folder, you will find a "
                               "folder 'Database_init' containing all fasta files, as well as a "
                               "folder 'refseq' with files downloaded from refseq."))
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
    optional.add_argument("--l90", dest="l90", type=int, default=100,
                          help="Maximum value of L90 allowed to keep a genome. Default is 100.")
    optional.add_argument("--nbcont", dest="nbcont", type=utils_argparse.cont_num, default=999,
                          help=("Maximum number of contigs allowed to keep a genome. "
                                "Default is 999."))
    optional.add_argument("--cutn", dest="cutn", type=int, default=5,
                          help=("By default, each genome will be cut into new contigs when "
                                "at least 5 'N' at a stretch are found in its sequence. "
                                "If you don't want to "
                                "cut genomes into new contigs when there are stretches of 'N', "
                                "put 0 to this option. If you want to cut from a different number "
                                "of 'N' stretches, put this value to this option."))
    helper = parser.add_argument_group('Others')
    helper.add_argument("-v", "--verbose", dest="verbose", action="count", default=0,
                        help="Increase verbosity in stdout/stderr and log files.")
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
    # print(f"before check: {args}")
    # toto = check_args(parser, args)
    # print(f"res: {toto}")
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

    # Message if user kept default thresholds for L90 and nbcont. Just to warn him, to be sure
    # it was on purpose
    def thresholds_message(l90, nbcont):
        return ("  !! Your genomes will be filtered, and only the ones with 'L90' <= {} "
                "and 'number of contigs' < {} will be kept. If you want to change those "
                "thresholds, use '--l90' and '--nbcont' "
                "options.".format(args.l90, args.nbcont))

    #  ERRORS
    # Check that at least taxid or species name was given
    if not args.NCBI_species_taxid and not args.NCBI_species:
        parser.error("Give at least an NCBI species name or taxID.")

    # WARNINGS
    # User did not specify a species name
    if not args.NCBI_species:
        print(colored("WARNING: you did not provide a species name ('-s species' option'). "
                      "All files will be downloaded in a folder called with the NCBI species "
                      f"taxid {args.NCBI_species_taxid} instead of the species name.", "yellow"))
    # If user wants to cut genomes, warn him to check that it is on purpose (because default is cut at each 5'N')
    if args.cutn == 0 or args.cutn == 5:
        message = ("  !! Your genomes will be split when sequence contains at "
                   "least {}'N' at a stretch. If you want to change this threshold, use "
                   "'--cutn n' option (n=0 if you do not want to cut)").format(args.cutn)
        print(colored(message, "yellow"))

    # Warn user about selection of genomes thresholds
    if args.l90 == 100 or args.nbcont == 999:
        print(colored(thresholds_message(args.l90, args.nbcont), "yellow"))


    return args


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
