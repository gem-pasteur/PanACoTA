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


@author gem
April 2017
"""


import os
import sys
import subprocess
import logging
from logging.handlers import RotatingFileHandler

from pipelinepackage import genome_seq_functions as gfunc
from pipelinepackage import prokka_functions as pfunc
from pipelinepackage import format_functions as ffunc
from pipelinepackage import utils


def main(list_file, db_path, res_path, name, l90, nbcont, cutn, threads, date, force):
    """
    Main method, doing all steps:
    - analyse genomes (nb contigs, L90, stretches of N...)
    - keep only genomes with 'good' (according to user thresholds) L90 and nb_contigs
    - rename genomes with strain number in decreasing quality
    - annotate genome with prokka
    - format annotated genomes

    """
    # get only filename of list_file, without extension
    listfile_base = os.path.basename(os.path.splitext(list_file)[0])
    # name logfile, add timestamp if already existing
    logfile = os.path.join(res_path, "annote-genomes-" + listfile_base + ".log")
    # if os.path.isfile(logfile):
    #     import time
    #     logfile = os.path.splitext(logfile)[0] + time.strftime("_%y-%m-%d_%H-%m-%S.log")
    # set level of logger (here debug to show everything during development)
    level = logging.DEBUG
    init_logger(logfile, level)
    logger = logging.getLogger()

    # test if prokka is installed and in the path
    prokka_cmd = ["prokka", "-h"]
    utils.check_installed(prokka_cmd)

    # Read genome names
    genomes = read_genomes(list_file, name, date)
    # Get L90, nbcontig, size for all genomes, and cut at stretches of 'N' if asked
    gfunc.analyse_all_genomes(genomes, db_path, res_path, cutn)
    # Get list of genomes kept (according to L90 and nbcont thresholds)
    kept_genomes = {genome: info for genome, info in genomes.items()
                    if info[-2] <= nbcont and info[-1] <= l90}
    # Write discarded genomes to a file
    write_discarded(genomes, list(kept_genomes.keys()), list_file, res_path)
    # Rename genomes kept, ordered by quality
    gfunc.rename_all_genomes(kept_genomes, res_path)
    # Write lstinfo file (list of genomes kept with info on L90 etc.)
    write_lstinfo(list_file, kept_genomes, res_path)
    # Annotate all kept genomes
    results = pfunc.run_prokka_all(kept_genomes, threads, force)
    # Generate database (folders Proteins, Genes, Replicons, LSTINFO)
    skipped = ffunc.format_genomes(genomes, results, res_path)
    if skipped:
        write_warning_skipped(skipped)


def init_logger(logfile, level):
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


def write_warning_skipped(skipped):
    """
    At the end of the script, write a warning to the user with the names of the genomes
    which had problems with prokka.

    skipped: list of genomes with problems
    """
    logger = logging.getLogger()
    list_to_write = "\n".join(["\t- " + genome for genome in skipped])
    logger.warning(("Prokka had problems while annotating some genomes. Hence, they are not "
                    "formatted, and absent from your output database. Please look at their "
                    "Prokka logs (<output_directory>/tmp_files/<genome_name>-prokka.log) and "
                    "to the current error log (<output_directory>/<input_filename>.log.err) to "
                    "get more information, and run again to annotate and format them. "
                    "Here are the genomes: \n{}").format(list_to_write))


def write_discarded(genomes, kept_genomes, list_file, res_path):
    """
    Write the list of genomes discarded to a file, so that users can keep a trace of them,
    with their information (nb contigs, L90 etc.)

    genomes: {genome: [gembase_start_name, seq_file, genome_size, nb_contigs, L90]}
    kept_genomes: list of genomes kept
    list_file: input file containing the list of genomes
    res_path: folder where results must be saved
    """
    _, name_lst = os.path.split(list_file)
    outdisc = os.path.join(res_path, "discarded-" + ".".join(name_lst.split(".")[:-1]) +".lst")
    with open(outdisc, "w") as outdf:
        outdf.write("\t".join(["orig_name", "gsize", "nb_conts", "L90"]) + "\n")
        for genome, values in genomes.items():
            if genome in kept_genomes:
                continue
            _, _, gsize, nbcont, l90 = [str(x) for x in values]
            outdf.write("\t".join([genome, gsize, nbcont, l90]) + "\n")


def write_lstinfo(list_file, genomes, outdir):
    """
    Write lstinfo file, with following columns:
    gembase_name, orig_name, size, nbcontigs, l90

    list_file: input file containing the list of genomes
    genomes: {genome: [gembase_start_name, seq_file, genome_size, nb_contigs, L90]}
    outdir: folder where results must be saved
    """
    _, name_lst = os.path.split(list_file)

    outlst = os.path.join(outdir, "LSTINFO-" + ".".join(name_lst.split(".")[:-1]) +".lst")
    with open(outlst, "w") as outf:
        outf.write("\t".join(["gembase_name", "orig_name", "gsize", "nb_conts", "L90"]) + "\n")
        for genome, values in sorted(genomes.items(), key=sort_genomes):
            gembase, _, gsize, nbcont, l90 = [str(x) for x in values]
            outf.write("\t".join([gembase, genome, gsize, nbcont, l90]) + "\n")


def sort_genomes(x):
    """
    x = (genome_orig, [gembase, path, gsize, nbcont, L90]}
    with gembase = species.date.strain

    order by:
    - species
    - in each species, by strain number
    """
    return (x[1][0], int(x[1][0].split(".")[-1]))


def read_genomes(list_file, name, date):
    """
    Read list of genomes, and return them.
    If a genome has a name, also return it. Otherwise, return the name given by user.

    genomes = {genome: name or None}
    """
    logger = logging.getLogger()
    logger.info("Reading genomes")
    genomes = {}
    with open(list_file, "r") as lff:
        for line in lff:
            # empty line: go to the next one
            if line.strip() == "":
                continue
            elems = line.strip().split()
            if len(elems) == 1:
                genomes[elems[0]] = [name + "." + date]
            else:
                genomes[elems[0]] = [elems[1] + "." + date]
    return genomes


def parse():
    """
    Method to create a parser for command-line options
    """
    import argparse
    def gen_name(param):
        if len(param) != 4:
            msg = ("The genome name must contain 4 characters. For example, this name can "
                   " correspond to the 2 first letters of genus, and 2 first letters of "
                   "species, e.g. ESCO for Escherichia Coli.")
            raise argparse.ArgumentTypeError(msg)
        return param

    def get_date():
        import time
        return time.strftime("%m%y")

    def cont_num(param):
        param = int(param)
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
    parser.add_argument("-s", dest="name", required=True, type=gen_name,
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
    parser.add_argument("--date", dest="date", default=get_date(),
                        help=("Specify the date (MMYY) to give to your annotated genomes. "
                              "By default, will give today's date."))
    parser.add_argument("-F", "--force", dest="force", const="--force", action="store_const",
                        help=("Add this option if you want to run prokka and formatting steps "
                              "even if their result folder (for prokka step) or files (for "
                              "format step) already exist for the given genome: override "
                              "existing results.\n"
                              "Otherwise, without this option, if the prokka folder exists, the "
                              "pipeline will run the formatting step with the already existing "
                              "results in prokka folder. If the Genes, Proteins, Replicons and "
                              "LSTINFO files already exist, the pipeline will skip the "
                              "formatting step too.\n"
                              "Note that this will be applied to all genomes: all prokka result "
                              "folders and format result files will be overridden. If you want "
                              "to rerun prokka (resp. rerun formatting step) only on a specific "
                              "genome, remove its prokka result folder (resp. Genes, Proteins, "
                              "LSTINFO, Replicons result files) before running this script "
                              "without this '-F' option."))
    args = parser.parse_args()
    # if args.multi and args.mixed:
    #     parser.error("-M and -X options cannot be activated together. Choose if you want to:\n"
    #                  "- allow several members in any number of genomes of a family (-M)\n"
    #                  "- allow several members in only '1-tol'% of the genomes of a family "
    #                  "(other 'tol'% genomes must have exactely 1 member)")
    return args

if __name__ == '__main__':
    OPTIONS = parse()
    main(OPTIONS.list_file, OPTIONS.db_path, OPTIONS.res_path, OPTIONS.name, OPTIONS.l90,
         OPTIONS.nbcont, OPTIONS.cutn, OPTIONS.threads, OPTIONS.date, OPTIONS.force)
