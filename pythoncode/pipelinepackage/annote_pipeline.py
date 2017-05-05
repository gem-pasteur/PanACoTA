#!/usr/bin/env python3
# coding: utf-8

"""
Pipeline to annotate genomes. Steps are:
- optional: find stretches of at least 'n' N (default 5), and cut into a new contig at this stretch
- for each genome, calc L90 and number of contigs (after cut at N stretches if used)
- keep only genomes with:
    - L90 <= x (default = 100)
    - #contig <= y (default = 999)
- rename those genomes and their contigs, with strain name increasing with quality (L90 and
 #contig)
- annotate kept genomes with prokka
- gembase format
- find essential genes and give a distribution graph/table so that the user can choose which
genomes he wants to remove on the 'essential genes' criteria

Input:
- list_file: list of genome filenames to annotate. Each genome is in a multi-fasta file. This file
contains 1 genome filename per line. It can contain a second column, with the species in 4 letters.
Genomes without this 2nd column will be renamed with the name given in species
- species: with 4 letters to rename genomes (except those whose species name is precised
in 2nd column of list file)
- dbpath: path to folder containing all multi-fasta sequences of genomes
- respath: path to folder where outputs must be saved (folders Genes, Replicons, Proteins,
LSTINFO and LSTINFO_dataset.lst file)
- threads: number of threads that can be used (default 1)

Output:
- In your given respath, you will find 4 folders: LSTINFO (information on each genome, with gene annotations), Genes (nuc. gene sequences), Proteins (aa proteins sequences), Replicons (input sequences but with formatted headers).
- In your given respath, you will find a "tmp_files" folder, where folders with prokka results will be created for each input genome (<genome_name>-prokkaRes). If errors are generated during prokka step, you can look at the log file to see what was wrong (<genome_name>-prokka.log).
- In your given respath, a file called `annote-genomes-<list_file>.log` will be generated. You can find there all logs: problems during annotation (hence no formatting step ran), and problems during formatting step.
- In your given respath, a file called `annote-genomes-<list_file>.log.err` will be generated, containing information on errors and warnings that occured. If this file is empty, then annotation and formatting steps finished without any problem for all genomes.
- In your given respath, you will find a file called `LSTINFO-<list_file>.lst` with information on all genomes: gembase_name, original_name, genome_size, L90, nb_contigs
- In your given respath, you will find a file called `discarded-<list_file>.lst` with information on genomes that were discarded (and hence not annotated) because of the L90 and/or nb_contig threshold: original_name, genome_size, L90, nb_contigs
- In your given respath, you will find 2 png files: `QC_L90-<list_file>.png` and `QC_nb-contigs-<list_file>.png`, containing the histograms of L90 and nb_contigs values for all genomes, with a vertical red line representing the limit applied here.

Requested:
- in prokka results, all genes are called <whatever>_<number> -> the number will be kept.
- The number of the genes annotated by prokka are in increasing order in tbl, faa and ffn files
- genome names given to prokka should not end with '_<number>'. Ideally, they should always have
the same format: <spegenus>.<date>.<strain_number> but they can have another format, as long as
they don't end by '_<number>', which is the format of a gene name.

@author gem
April 2017
"""


import os
import sys
import subprocess
import shutil
import logging
from logging.handlers import RotatingFileHandler

from pipelinepackage import genome_seq_functions as gfunc
from pipelinepackage import prokka_functions as pfunc
from pipelinepackage import format_functions as ffunc
from pipelinepackage import utils


def main(list_file, db_path, res_dir, name, date, l90=100, nbcont=999, cutn=5,
         threads=1, force=False, qc_only=False, tmp_dir=None, prok_dir=None):
    """
    Main method, doing all steps:
    - analyse genomes (nb contigs, L90, stretches of N...)
    - keep only genomes with 'good' (according to user thresholds) L90 and nb_contigs
    - rename genomes with strain number in decreasing quality
    - annotate genome with prokka
    - format annotated genomes

    """
    # By default, all tmp files (split sequences, renamed sequences, prokka results) will
    # be saved in the given <res_dir>/tmp_files.
    # Create output (results, tmp...) directories if not already existing
    if not tmp_dir:
        tmp_dir = os.path.join(res_dir, "tmp_files")
    if not prok_dir:
        prok_dir = tmp_dir
    os.makedirs(res_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(prok_dir, exist_ok=True)
    # If force was set, remove result folders (Proteins, Replicons, Genes, LSTINFO)
    if force:
        shutil.rmtree(os.path.join(res_dir, "LSTINFO"), ignore_errors=True)
        shutil.rmtree(os.path.join(res_dir, "Proteins"), ignore_errors=True)
        shutil.rmtree(os.path.join(res_dir, "Genes"), ignore_errors=True)
        shutil.rmtree(os.path.join(res_dir, "Replicons"), ignore_errors=True)
    else:
        # Check that resultdir is not already used
        utils.check_out_dirs(res_dir)

    # get only filename of list_file, without extension
    listfile_base = os.path.basename(os.path.splitext(list_file)[0])
    # name logfile, add timestamp if already existing
    logfile = os.path.join(res_dir, "annote-genomes-" + listfile_base + ".log")
    # if os.path.isfile(logfile):
    #     import time
    #     logfile = os.path.splitext(logfile)[0] + time.strftime("_%y-%m-%d_%H-%m-%S.log")
    # set level of logger (here debug to show everything during development)
    level = logging.DEBUG
    init_logger(logfile, level)
    logger = logging.getLogger()
    if not qc_only:
        # test if prokka is installed and in the path
        prokka_cmd = ["prokka", "-h"]
        utils.check_installed(prokka_cmd)

    # Read genome names.
    # genomes = {genome: [spegenus.date]}
    genomes = utils.read_genomes(list_file, name, date, db_path)
    # Get L90, nbcontig, size for all genomes, and cut at stretches of 'N' if asked
    # genomes = {genome: [spegenus.date, path_to_splitSequence, size, nbcont, l90]}
    gfunc.analyse_all_genomes(genomes, db_path, tmp_dir, cutn)
    # Plot L90 and nb_contigs distributions
    gfunc.plot_distributions(genomes, res_dir, listfile_base, l90, nbcont)
    # Get list of genomes kept (according to L90 and nbcont thresholds)
    kept_genomes = {genome: info for genome, info in genomes.items()
                    if info[-2] <= nbcont and info[-1] <= l90}
    # Write discarded genomes to a file
    utils.write_discarded(genomes, list(kept_genomes.keys()), list_file, res_dir)
    # If only QC, stop here.
    if qc_only:
        return genomes, kept_genomes
    # Rename genomes kept, ordered by quality
    # kept_genomes = {genome: [gembase_name, path_split_gembase, gsize, nbcont, L90]}
    gfunc.rename_all_genomes(kept_genomes, tmp_dir)
    # Write lstinfo file (list of genomes kept with info on L90 etc.)
    utils.write_lstinfo(list_file, kept_genomes, res_dir)
    # Annotate all kept genomes
    results = pfunc.run_prokka_all(kept_genomes, threads, force, prok_dir)
    # Generate database (folders Proteins, Genes, Replicons, LSTINFO)
    skipped, skipped_format = ffunc.format_genomes(genomes, results, res_dir, prok_dir)
    if skipped:
        utils.write_warning_skipped(skipped)
    if skipped_format:
        utils.write_warning_skipped(skipped_format, format=True)
    return genomes, kept_genomes, skipped, skipped_format


def init_logger(logfile, level, name= None):
    """
    Create logger and its handlers, and set them to the given level

    level hierarchy:
    CRITICAL > ERROR > WARNING > INFO > DEBUG

    Messages from all levels are written in 'logfile'
    Messages for levels less than WARNING (only INFO and DEBUG) written to stdout
    Messages for levels equal or higher than WARNING written to stderr

    level: minimum level that must be considered.
    """
    # create logger
    if name:
        logger = logging.getLogger(name)
    else:
        logger = logging.getLogger()
    # set level of logger (here debug to show everything during development)
    logger.setLevel(level)
    # create formatter for log messages: "timestamp :: level :: message"
    formatterFile = logging.Formatter('[%(asctime)s] :: %(levelname)s :: %(message)s')
    formatterStream = logging.Formatter('  * %(message)s')

    # Create handler 1: writing to 'logfile'. mode 'write', max size = 1Mo. If logfile is 1Mo, it is renamed to logfile.1, and next logs are still written to logfile. Then, logfile.1 is renamed to logfile.2, logfile to logfile.1 etc. We allow maximum 5 log files.
    open(logfile, "w").close()  # empty logfile if already existing
    errfile = logfile + ".err"
    open(errfile, "w").close()
    logfile_handler = RotatingFileHandler(logfile, 'w', 1000000, 5)
    # set level to the same as the logger level
    logfile_handler.setLevel(level)
    logfile_handler.setFormatter(formatterFile)  # add formatter
    logger.addHandler(logfile_handler)  # add handler to logger

    # Create handler 2: errfile
    errfile_handler = RotatingFileHandler(errfile, 'w', 1000000, 5)
    # set level to the same as the logger level
    errfile_handler.setLevel(logging.WARNING)
    errfile_handler.setFormatter(formatterFile)  # add formatter
    logger.addHandler(errfile_handler)  # add handler to logger

    # Create handler 3: write to stdout
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)  # write any message
    stream_handler.addFilter(utils.LessThanFilter(logging.WARNING)) # don't write messages >= WARNING
    stream_handler.setFormatter(formatterStream)
    logger.addHandler(stream_handler)  # add handler to logger

    # Create handler 4: write to stderr
    err_handler = logging.StreamHandler(sys.stderr)
    err_handler.setLevel(logging.WARNING)  # write any message >= WARNING
    err_handler.setFormatter(formatterStream)
    logger.addHandler(err_handler)  # add handler to logger


def parse(argu=sys.argv[1:]):
    """
    Method to create a parser for command-line options
    """
    import argparse
    def gen_name(param):
        if not utils.check_format(param):
            msg = ("The genome name must contain 4 characters. For example, this name can "
                   " correspond to the 2 first letters of genus, and 2 first letters of "
                   "species, e.g. ESCO for Escherichia Coli.")
            raise argparse.ArgumentTypeError(msg)
        return param

    def date_name(param):
        if not utils.check_format(param):
            msg = ("The date must contain 4 characters. Usually, it contains 4 digits, "
                   "corresponding to the month (2 digits) and year (2 digits).")
            raise argparse.ArgumentTypeError(msg)
        return param

    def get_date():
        import time
        return time.strftime("%m%y")

    def cont_num(param):
        try:
            param = int(param)
        except Exception:
            msg = "argument --nbcont: invalid int value: {}".format(param)
            raise argparse.ArgumentTypeError(msg)
        if param < 0 :
            msg = ("The maximum number of contigs allowed must be a positive number.")
            raise argparse.ArgumentTypeError(msg)
        if param >= 10000:
            msg = ("We do not support genomes with more than 9999 contigs.")
            raise argparse.ArgumentTypeError(msg)
        return param

    parser = argparse.ArgumentParser(description=("Annotate all genomes"))
    # Create command-line parser for all options and arguments to give
    parser.add_argument(dest="list_file",
                        help=("File containing the list of genome filenames to annotate (1 genome"
                              " per line). Each genome is in multi-fasta format. You can "
                              "specify the species name (4 letters) you want to give to each "
                              "genome by adding it after the genome filename, separated by a "
                              "space. If not given, the species name will be the one given in "
                              "'species' argument. "))
    parser.add_argument("-d", dest="db_path", required=True,
                        help=("Path to folder containing all multifasta genome files"))
    parser.add_argument("-r", dest="res_path", required=True,
                        help=("Path to folder where output annotated genomes must be saved"))
    parser.add_argument("-n", dest="name", type=gen_name,
                        help=("Choose a name for your annotated genomes. This name should contain 4 letters. Generally, they correspond to the 2 first letters of genus, and 2 first letters of species, e.g. ESCO for Escherichia Coli."))
    parser.add_argument("--l90", dest="l90", type=int, default=100,
                        help=("Maximum value of L90 allowed to keep a genome. Default is 100."))
    parser.add_argument("--nbcont", dest="nbcont", type=cont_num, default=999,
                        help=("Maximum number of contigs allowed to keep a genome. "
                              "Default is 999."))
    parser.add_argument("--cutN", dest="cutn", type=int, default=5,
                        help=("By default, each genome will be cut into new contigs at each "
                              "stretch of at least 5 'N' in its sequence. If you don't want to "
                              "cut genomes into new contigs when there are stretches of 'N', "
                              "put 0 to this option. If you want to cut from a different number "
                              "of 'N' stretches, put this value to this option."))
    parser.add_argument("--threads", dest="threads", type=int, default=1,
                        help=("Specify how many threads can be used (default=1)"))
    parser.add_argument("--date", dest="date", default=get_date(), type=date_name,
                        help=("Specify the date (MMYY) to give to your annotated genomes. "
                              "By default, will give today's date. The only requirement on the"
                              " given date is that it is 4 characters long. You can use letters"
                              " if you want. But the common way is to use 4 digits, "
                              "corresponding to MMYY."))
    parser.add_argument("--tmp", dest="tmpdir",
                        help=("Specify where the temporary files (sequence split by stretches "
                              "of 'N', sequence with new contig names etc.) must be saved. "
                              "By default, it will be saved in your result_directory/tmp_files."))
    parser.add_argument("--prok", dest="prokkadir",
                        help=("Specify in which directory the prokka output files "
                              "(1 folder per genome, called <genome_name>-prokkaRes) must be "
                              "saved. By default, they are saved in the same directory as "
                              "your temporary files (see --tmp option to change it)."))
    parser.add_argument("-F", "--force", dest="force", action="store_true",
                        help=("Force run: Add this option if you want to run prokka and "
                              "formatting steps for all genomes "
                              "even if their result folder (for prokka step) or files (for "
                              "format step) already exist: override "
                              "existing results.\n"
                              "Without this option, if there already are results in the given "
                              "result folder, the program stops. If there are no results, but "
                              "prokka folder already exists, prokka won't run again, and the "
                              "formating step will use the already existing folder if correct, "
                              "or skip the genome if there are problems in prokka folder."))
    parser.add_argument("-Q", dest="qc_only", action="store_true", default=False,
                        help=("Add this option if you want only to do quality control on your "
                              "genomes (cut at 5N if asked, calculate L90 and number of contigs "
                              "and plot their distributions). This allows you to check which "
                              "genomes would be annotated with the given parameters, and to "
                              "modify those parameters if you want, before you launch the "
                              "annotation and formatting steps."))
    args = parser.parse_args(argu)
    if not args.qc_only and not args.name:
        parser.error("You must specify your genomes dataset name in 4 characters with "
                     "'-n name' option (type -h for more information). Or, if you do not want "
                     "to annotate and format your genomes but just to run quality control, use "
                     "option '-Q")
    if args.qc_only and not args.name:
        args.name = "NONE"
    return args


if __name__ == '__main__':
    OPTIONS = parse(sys.argv[1:])
    main(OPTIONS.list_file, OPTIONS.db_path, OPTIONS.res_path, OPTIONS.name, OPTIONS.date,
         OPTIONS.l90, OPTIONS.nbcont, OPTIONS.cutn, OPTIONS.threads,
         OPTIONS.force, OPTIONS.qc_only, OPTIONS.tmpdir, OPTIONS.prokkadir)
