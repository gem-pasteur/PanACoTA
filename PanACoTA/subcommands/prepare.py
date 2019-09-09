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
from PanACoTA.prepare_module import filter_genomes as fg


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
         arguments.tmp_dir, arguments.parallel, arguments.no_refseq, arguments.only_mash,
         arguments.from_info, arguments.l90, arguments.nbcont, arguments.cutn, arguments.min_dist,
         arguments.verbose, arguments.quiet)


def main(cmd, NCBI_species, NCBI_taxid, outdir, tmp_dir, threads, no_refseq, only_mash, info_file,
         l90, nbcont, cutn, min_dist, verbose, quiet):
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
    tmp_dir : str
        Path to directory where tmp files are saved (sequences split at each stretch of 'N')
    threads : int
        max number of threads to use
    no_refseq : bool
        True if user does not want to download again the database
    only_mash : bool
        True if user user already has the database and quality of each genome (L90, #contigs etc.)
    info_file : str
        File containing information on QC if it was already ran before (columns to_annotate,
        gsize, nb_conts and L90).
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
        - defaut 0 : stdout contains INFO, stderr contains ERROR, .log contains INFO and more,
          .log.err contains warning and more
        - 1: same as 0 + WARNING in stderr
        - 2: same as 1 + DETAILS in stdout + DETAILS in .log.details
        - >=15: same as 2 + Add DEBUG in stdout + create .log.debug with everything
          from info to debug
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
    # Default outdir is species name if given, or species taxID
    if not outdir:
        outdir = species_linked
    # Default tmp_dir is outdir/tmp_files
    if not tmp_dir:
        tmp_dir = os.path.join(outdir, "tmp_files")
    db_dir = None
    # directory that will be created by ncbi_genome_download
    refseqdir = os.path.join(outdir, "refseq", "bacteria")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

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
        # Not only mash, so a new info file will be created. If the user still gave an info
        # file (he will be warned that it will be ignored), rename it with '.bak'
        # to avoid erasing it
        if info_file:
            os.rename(info_file, info_file + ".back")
        # 'no_refseq = True" : Do not download genomes, just do QC and mash filter on given genomes
        # (sequences must, at least, be in outdir/refeq/bacteria/<genome_name>.fna.gz)
        # (they can also be in Database_init/<genome_name>.fna)
        if no_refseq:
            logger.warning('You asked to skip refseq downloads.')
            # Check that db_dir exists (folder with uncompressed fna files)
            db_dir = os.path.join(outdir, "Database_init")
            # If it does not exist, check that original dir exists
            if not os.path.exists(db_dir):
                logger.warning(f"Database folder {db_dir} containing fna sequences does not "
                               "exist. We will check that refseq donwload directory exists, "
                               "to be able to do next steps (genomes filter)")
                db_dir = None

                if not os.path.exists(refseqdir):
                    logger.error(f"{refseqdir} does not exist. You do not have any "
                                 "genome to analyse")
                    sys.exit(1)
        # No sequence: Do all steps -> download, QC, mash filter
        else:
            # Download all genomes of the given taxID
            db_dir = dgf.download_from_refseq(species_linked, NCBI_species, NCBI_taxid,
                                              outdir, threads)
        # if norefseq: refseq/bacteria must exist, but Database_init can be absent (if
        # genomes were downloaded but not uncompressed)
        # -> create and fill it
        if not db_dir:
            db_dir = os.path.join(outdir, "Database_init")
            # If db_dir does not exist: create and fill it
            if no_refseq and not os.path.exists(db_dir):
                # add genomes from refseq/bacteria folder to Database_init
                nb_gen, _ = dgf.to_database(outdir)
                # If no genome found, error -> nothing to analyse
                if nb_gen == 0:
                    logger.error(f"There is no genome in {refseqdir}.")
                    sys.exit(1)
                logger.info("{} refseq genome(s) downloaded".format(nb_gen))
        # Now that genomes are downloaded and uncompressed, check their quality to remove bad ones
        genomes = fg.check_quality(outdir, species_linked, db_dir, tmp_dir, l90, nbcont, cutn)
    # Do only mash filter. Genomes must be already downloaded, and there must be a file with
    # all information on these genomes (L90 etc.)
    else:
        logger.warning('You asked to run only mash steps.')
        if not os.path.exists(info_file):  # info-file missing -> error and exit
            logger.error(f"Your info file {info_file} does not exist. Please Provide the  "
                          "right name/path, or remove the '--mash-only option to rerun "
                          "quality control.")
            sys.exit(1)
        logger.info(("You want to run only mash steps. Getting information "
                     "from {}").format(info_file))
        genomes = utils.read_genomes_info(info_file, species_linked, )
    # Run Mash
    # genomes : {genome_file: [genome_name, orig_name, path_to_seq_to_annotate, size,
                             # nbcont, l90]}
    # sorted_genome : [genome_file] ordered by L90/nbcont (keys of genomes)
    sorted_genomes = fg.sort_genomes_minhash(genomes, l90, nbcont)
    removed = fg.iterative_mash(sorted_genomes, genomes, outdir, species_linked,
                                min_dist, max_dist, threads)
    # Write list of genomes kept, and list of genomes removed
    fg.write_outputfiles(genomes, sorted_genomes, removed, outdir, species_linked, min_dist)
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
    optional.add_argument("-o", dest="outdir",
                          help=("Give the path to the directory where you want to save the "
                               "database. In the given directory, it will create a folder with "
                               "the gembase species name. Inside this folder, you will find a "
                               "folder 'Database_init' containing all fasta files, as well as a "
                               "folder 'refseq' with files downloaded from refseq."))
    optional.add_argument("--tmp", dest="tmp_dir",
                          help=("Specify where the temporary files (sequence split by stretches "
                                "of 'N', sequence with new contig names etc.) must be saved. "
                                "By default, it will be saved in your "
                                "result_directory/tmp_files."))
    optional.add_argument("-p", dest="parallel", type=utils_argparse.thread_num, default=1,
                          help=("Run 'N' downloads in parallel (default=1). Put 0 if "
                                "you want to use all cores of your computer."))
    optional.add_argument("--norefseq", dest="no_refseq", action="store_true",
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
    optional.add_argument("--info", dest="from_info",
                          help="If you already ran the 'prepare' data module, or already "
                               "calculated yourself the size, L90 and number of contigs for each "
                               "genome, you can give this information, to go directly to "
                               "Mash filtering step. This file contains at "
                               "least 4 columns, tab separated, with the following headers: "
                               "'to_annotate', 'gsize', 'nb_conts', 'L90'. Any other column "
                               "will be ignored.")
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

    # If user wants only mash steps, check that he gave info file
    if args.only_mash and not args.from_info:
        parser.error("If you want to run only Mash filtering steps, please give the "
                     "info file with the required information (see '--info' option)")

    # WARNINGS
    # User did not specify a species name
    if not args.NCBI_species:
        print(colored("WARNING: you did not provide a species name ('-s species' option'). "
                      "All files will be downloaded in a folder called with the NCBI species "
                      f"taxid {args.NCBI_species_taxid} instead of the species name.", "yellow"))
    # If user wants to cut genomes, warn him to check that it is on purpose (because default is cut at each 5'N')
    if args.cutn == 5:
        message = ("  !! Your genomes will be split when sequence contains at "
                   "least {}'N' at a stretch. If you want to change this threshold, use "
                   "'--cutn n' option (n=0 if you do not want to cut)").format(args.cutn)
        print(colored(message, "yellow"))

    # Warn user about selection of genomes thresholds
    if args.l90 == 100 or args.nbcont == 999:
        print(colored(thresholds_message(args.l90, args.nbcont), "yellow"))

    # Warn if user gave info file, but does not ask to run only Mash -> info file will be ignored
    if args.from_info and not args.only_mash:
        message = ("  !! You gave an info file (--info option), but did not ask to run only Mash "
                   "step (-M option). Your info file will be ignored (and renamed with '.back' "
                   "at the end), and another one will "
                   "be created with the new calculated values.")
        print(colored(message, "yellow"))

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
