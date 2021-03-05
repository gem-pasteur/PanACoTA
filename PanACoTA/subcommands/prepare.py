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

from PanACoTA import __version__ as version
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
    main(cmd, arguments.ncbi_species_name, arguments.ncbi_species_taxid, arguments.ncbi_taxid, arguments.levels,
         arguments.ncbi_section, arguments.outdir, arguments.tmp_dir, arguments.parallel, arguments.norefseq,
         arguments.db_dir, arguments.only_mash,
         arguments.info_file, arguments.l90, arguments.nbcont, arguments.cutn, arguments.min_dist,
         arguments.max_dist, arguments.verbose, arguments.quiet)


def main(cmd, ncbi_species_name, ncbi_species_taxid, ncbi_taxid, levels, ncbi_section,
         outdir, tmp_dir, threads, norefseq, db_dir,
         only_mash, info_file, l90, nbcont, cutn, min_dist, max_dist, verbose, quiet):
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
    ncbi_species_name : str
        name of species to download, as given by NCBI
    ncbi_species_taxid : int
        species taxid given in NCBI
    ncbi_taxid : int
        NCBI taxid of strain
    levels: str
        Level of assembly to download. Choice between 'all', 'complete', 'chromosome',
        'scaffold', 'contig'. Default is 'all'
    outdir : str
        path to output directory (where created database will be saved).
    tmp_dir : str
        Path to directory where tmp files are saved (sequences split at each row of 5 'N')
    threads : int
        max number of threads to use
    norefseq : bool
        True if user does not want to download again the database
    db_dir : str
        Name of the folder where already downloaded fasta files are saved.
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
        cut at each when there are 'cutn' N in a row. Don't cut if equal to 0
    min_dist : int
        lower limit of distance between 2 genomes to keep them
    max_dist : int
        upper limit of distance between 2 genomes to keep them (default is 0.06)
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

    # get species name in NCBI format
    # -> will be used to name output directory
    # -> will be used to download summary file if given species corresponds to NCBI name
    if ncbi_species_name:
        species_linked = "_".join(ncbi_species_name.split())
        species_linked = "_".join(species_linked.split("/"))

    # if species name not given by user, use species taxID (if given) to name output directory
    elif ncbi_species_taxid:
        species_linked = str(ncbi_species_taxid)
    # if species name not species taxid by user, use taxID (if given) to name output directory
    elif ncbi_taxid:
        species_linked = str(ncbi_taxid)
    # if neither speName, speID nor taxID given (--norefseq, mashonly), name is NA
    else:
        species_linked = "NA"
    # Default outdir is species name if given, or species taxID
    if not outdir:
        outdir = species_linked
    # Default tmp_dir is outdir/tmp_files
    if not tmp_dir:
        tmp_dir = os.path.join(outdir, "tmp_files")
    # directory that will be created by ncbi_genome_download
    ncbidir = os.path.join(outdir, ncbi_section, "bacteria")
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
    logfile_base, logger = utils.init_logger(logfile_base, level, 'prepare', log_details=True,
                                             verbose=verbose, quiet=quiet)

    # Message on what will be done (cmd, cores used)
    logger.info(f'PanACoTA version {version}')
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
        if info_file and os.path.isfile(info_file):
            os.rename(info_file, info_file + ".back")

        # 'norefseq = True" : Do not download genomes, just do QC and mash filter on given genomes
        # -> if not, error and exit
        if norefseq:
            logger.warning(f'You asked to skip {ncbi_section} downloads.')

            # -> if db_dir given, watch for sequences there. If does not exist, error and exit
            # (user gave a directory (even if it does not exist), so we won't look for
            # the sequences in other folders)
            if db_dir:
                if not os.path.exists(db_dir):
                    logger.error(f"Database folder {db_dir} supposed to contain fasta "
                                 "sequences does not "
                                 "exist. Please give a valid folder, or leave the default "
                                 "directory (no '-d' option).")
                    sys.exit(1)
            # -> If user did not give db_dir, genomes could be in
            # outdir/Database_init/<genome_name>.fna
            else:
                db_dir = os.path.join(outdir, "Database_init")
                # If it does not exist, check if default compressed files folder exists.
                if not os.path.exists(db_dir):
                    logger.warning(f"Database folder {db_dir} supposed to contain fasta "
                                   "sequences does not "
                                   "exist. We will check if the download folder (with compressed "
                                   "sequences) exists.")
                    # -> if not in database_init, genomes must be in
                    # outdir/refeq/bacteria/<genome_name>.fna.gz. In that case,
                    # uncompress and add them to Database_init
                    if not os.path.exists(ncbidir):
                        logger.error(f"Folder {ncbidir} does not exist. You do not have any "
                                     "genome to analyse. Possible reasons:\n"
                                     "- if you want to rerun analysis in the same folder as "
                                     "sequences were downloaded (my_outdir/Database_init or "
                                     f"my_outdir/{ncbi_section}), make sure you have '-o my_outdir' "
                                     "option\n"
                                     "- if you want to rerun analysis and save them in a new "
                                     "output folder called 'new_outdir', make sure you have "
                                     "'-o new_outdir' option, "
                                     "and you specified where the uncompressed sequences to "
                                     "use are ('-d sequence_database_path'). ")
                        sys.exit(1)
                    # add genomes from refseq/bacteria folder to Database_init
                    nb_gen, _ = dgf.to_database(outdir, ncbi_section)
        # No sequence: Do all steps -> download, QC, mash filter
        else:
            # Download all genomes of the given taxID
            db_dir, nb_gen = dgf.download_from_ncbi(species_linked, ncbi_section, ncbi_species_name, ncbi_species_taxid,
                                                      ncbi_taxid, levels, outdir, threads)
            logger.info(f"{nb_gen} {ncbi_section} genome(s) downloaded")

        # Now that genomes are downloaded and uncompressed, check their quality to remove bad ones
        genomes = fg.check_quality(species_linked, db_dir, tmp_dir, l90, nbcont, cutn)

    # Do only mash filter. Genomes must be already downloaded, and there must be a file with
    # all information on these genomes (L90 etc.)
    else:
        logger.warning('You asked to run only mash steps.')
        if not os.path.exists(info_file):  # info-file missing -> error and exit
            logger.error(f"Your info file {info_file} does not exist. Please provide the  "
                          "right name/path, or remove the '--mash-only option to rerun "
                          "quality control.")
            sys.exit(1)
        logger.info(("You want to run only mash steps. Getting information "
                     "from {}").format(info_file))
        genomes = utils.read_genomes_info(info_file, species_linked, )

    # Run Mash
    # genomes : {genome_file: [genome_name, orig_name, path_to_seq_to_annotate, size, nbcont, l90]}
    # sorted_genome : [genome_file] ordered by L90/nbcont (keys of genomes)
    sorted_genomes = fg.sort_genomes_minhash(genomes, l90, nbcont)

    # Write discarded genomes to a file -> orig_name, to_annotate, gsize, nb_conts, L90
    discQC = f"by-L90_nbcont-{species_linked}.txt"
    utils.write_genomes_info(genomes, sorted_genomes, discQC, outdir)

    # Remove genomes not corresponding to mash filters
    removed = fg.iterative_mash(sorted_genomes, genomes, outdir, species_linked,
                                min_dist, max_dist, threads, quiet)
    # Write list of genomes kept, and list of genomes discarded by mash step
    info_file = fg.write_outputfiles(genomes, sorted_genomes, removed, outdir, species_linked,
                                     min_dist, max_dist)
    logger.info("End")
    return info_file


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

    general = parser.add_argument_group('General arguments')
    general.add_argument("-T", dest="ncbi_species_taxid", default="",
                          help=("Species taxid to download, corresponding to the "
                                "'species taxid' provided by the NCBI. "
                                "This will download all sequences of this species and all its sub-species."
                                "A comma-separated list of species taxids can also be provided. "
                                "(Ex: -T 573 for Klebsiella pneumoniae)")
                         )
    general.add_argument("-t", dest="ncbi_taxid", default="",
                          help=("Taxid to download. "
                                "This can be the taxid of a sub-species, or of a specific strain. "
                                "A taxid of a subspecies will download all strains in this subspecies "
                                "EXCEPT the ones which have a specific taxid."
                                "A comma-separated list of taxids can also be provided."
                                "Ex: '-t 72407' will download all 'general' K. pneumoniae subsp. pneumoniae strains, "
                                "and '-t 1123862' will download the strain K. pneumoniae subsp. pneumoniae Kp13 "
                                "(not included in -t 72407, as it is a strain of the subspecies with a specific taxid).")
                         )
    general.add_argument("-g", dest="ncbi_species_name", default="",
                          help=("Species to download, corresponding to the "
                                "'organism name' provided by the NCBI. Give name between "
                                "quotes (for example \"escherichia coli\")")
                        )
    general.add_argument("-s", dest="ncbi_section", default="refseq", choices = ["refseq", "genbank"],
                          help=("NCBI section to download: all genbank, or only refseq (default)")
                        )
    general.add_argument("-l", "--assembly_level", dest="levels", default="",
                          help=("Assembly levels of genomes to download (default: all). "
                                "Possible levels are: 'all', 'complete', 'chromosome', "
                                "'scaffold', 'contig'."
                                "You can also provide a comma-separated list of assembly "
                                "levels. For ex: 'complete,chromosome'")
                          )
    general.add_argument("-o", dest="outdir",
                          help=("Give the path to the directory where you want to save the "
                               "downloaded database. In the given directory, it will create "
                               "a folder 'Database_init' containing all fasta "
                               "files that will be sent to the control procedure, as well as "
                               "a folder 'refseq' with all original compressed files "
                               "downloaded from refseq. By default, this output dir name is the "
                               "ncbi_species name if given, or ncbi_species_taxid or ncbi_taxid otherwise.")
                          )
    general.add_argument("--tmp", dest="tmp_dir",
                          help=("Specify where the temporary files (sequence split by rows "
                                "of 'N', sequence with new contig names etc.) must be saved. "
                                "By default, it will be saved in your "
                                "out_dir/tmp_files.")
                          )
    general.add_argument("--cutn", dest="cutn", type=utils_argparse.positive_int, default=5,
                          help=("By default, each genome will be cut into new contigs when "
                                "at least 5 'N' in a row are found in its sequence. "
                                "If you don't want to "
                                "cut genomes into new contigs when there are rows of 'N', "
                                "put 0 to this option. If you want to cut from a different number "
                                "of 'N' in a row, put this value to this option.")
                          )
    general.add_argument("--l90", dest="l90", type=int, default=100,
                          help="Maximum value of L90 allowed to keep a genome. Default is 100.")
    general.add_argument("--nbcont", dest="nbcont", type=utils_argparse.cont_num, default=999,
                          help=("Maximum number of contigs allowed to keep a genome. "
                                "Default is 999."))
    general.add_argument("--min_dist", dest="min_dist", default=1e-4,
                        type=utils_argparse.mash_dist,
                        help="By default, genomes whose distance to the reference is not "
                             "between 1e-4 and 0.06 are discarded. You can specify your own "
                             "lower limit (instead of 1e-4) with this option.")
    general.add_argument("--max_dist", dest="max_dist", default=0.06,
                         type=utils_argparse.mash_dist,
                         help="By default, genomes whose distance to the reference is not "
                              "between 1e-4 and 0.06 are discarded. You can specify your own "
                              "lower limit (instead of 0.06) with this option.")
    general.add_argument("-p", "--threads", dest="parallel", type=utils_argparse.thread_num,
                         default=1, help=("Run 'N' downloads in parallel (default=1). Put 0 if "
                                "you want to use all cores of your computer."))

    optional = parser.add_argument_group('Alternatives')
    optional.add_argument("--norefseq", dest="norefseq", action="store_true",
                          help=("If you already downloaded refseq genomes and do not want to "
                                "check them, add this option to directly go to the next steps:"
                                "quality control (L90, number of contigs...) and mash filter. "
                                "Don't forget to specify the db_dir (-d option) where you "
                                "already have your genome sequences."))
    optional.add_argument("-d", dest="db_dir",
                          help=("If your already downloaded sequences are not in the default "
                                "directory (outdir/Database_init), you can specify here the "
                                "path to those fasta files."))
    optional.add_argument('-M', '--only-mash', dest="only_mash", action="store_true",
                          help=("Add this option if you already downloaded complete and refseq "
                                "genomes, and ran quality control (you have, a file "
                                "containing all genome names, as well as their genome size, "
                                "number of contigs and L90 values). "
                                "It will then get information on genomes quality from this "
                                "file, and run mash steps."))
    optional.add_argument("--info", dest="info_file",
                          help=("If you already ran the quality control, specify from which "
                                "file PanACoTA can read this information, in order to proceed "
                                "to the mash step. This file must contain at "
                                "least 4 columns, tab separated, with the following headers: "
                                "'to_annotate', 'gsize', 'nb_conts', 'L90'. Any other column "
                                "will be ignored."))

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

    # We don't want to run only mash, nor only quality control, but don't give a NCBI taxID.
    # -> Give at least 1!
    if (not args.only_mash and not args.norefseq and
        not args.ncbi_species_taxid and not args.ncbi_species_name and not args.ncbi_taxid):
        parser.error("As you did not put the '--norefseq' nor the '-M' option, it means that "
                     "you want to download refseq (or genbank) genomes. But you did not provide any "
                     "information, so PanACoTA cannot guess which species you want to download. "
                     "Specify NCBI_taxid (-t), and/or NCBI species taxid (-T) and/or NCBI_species (-g) to download, or add one of "
                     "the 2 options (--norefseq or -M) if you want to skip the 'download step'.")

    # If norefseq, give output directory
    #  - folder containing Database_init, with all sequences
    #  - or new folder where you want to put the new results
    if args.norefseq and not args.outdir:
        parser.error("You must provide an output directory, where your results will be saved.")

    # If user wants only mash steps, check that he gave info file, and outdir
    if args.only_mash and not args.info_file:
        parser.error("If you want to run only Mash filtering steps, please give the "
                     "info file with the required information (see '--info' option)")
    if args.only_mash and not args.outdir:
        parser.error("If you want to run only Mash filtering steps, please give the "
                     "output directory where you want to save your results (see '-o' option)")

    # Cannot be verbose and quiet at the same time
    if int(args.verbose) > 0 and args.quiet:
        parser.error("Choose between a verbose output (-v) or a quiet output (-q)."
                     " You cannot have both.")

    # min_dist must be higher than max_dist
    if float(args.min_dist) >= float(args.max_dist):
        parser.error(f"min_dist ({args.min_dist}) cannot be higher "
                     f"than max_dist ({args.max_dist})")

    # Check that levels, if given, are among possible ones
    possible = ["all", "complete", "chromosome", "scaffold", "contig"]
    if args.levels:
        for level in args.levels.split(","):
            if level not in possible:
                parser.error("Please choose between available assembly levels: 'all', 'complete', "
                             "'chromosome', 'scaffold', 'contig'. If several levels, provide a "
                             f"comma-separated list. Invalid value: '{args.levels}'")

    # WARNINGS
    # User did not specify a species name
    if not args.ncbi_species_name and not args.outdir:
        if args.ncbi_species_taxid:
            print(colored("WARNING: you did not provide a species name ('-g species' option) "
                          "nor an output directory ('-o outdir'). "
                          "All files will be downloaded in a folder called with the NCBI species "
                          f"taxid {args.ncbi_species_taxid} instead of the species name.", "yellow"))
        else:
            print(colored("WARNING: you did not provide a species name ('-g species' option) "
                          "nor a species taxid ('-T spetaxid') nor an output directory ('-o outdir'). "
                          "All files will be downloaded in a folder called with the NCBI "
                          f"taxid {args.ncbi_taxid}.", "yellow"))
    # If user wants to cut genomes, warn him to check that it is on purpose (because default is cut at each 5'N')
    if args.cutn == 5:
        message = ("  !! Your genomes will be split when sequence contains at "
                   "least {}'N' in a row. If you want to change this threshold, use "
                   "'--cutn n' option (n=0 if you do not want to cut)").format(args.cutn)
        print(colored(message, "yellow"))

    # Warn user about selection of genomes thresholds
    if args.l90 == 100 or args.nbcont == 999:
        print(colored(thresholds_message(args.l90, args.nbcont), "yellow"))

    # Warn if user gave info file, but does not ask to run only Mash -> info file will be ignored
    if (args.info_file and not args.only_mash) or (args.info_file and not args.norefseq):
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
