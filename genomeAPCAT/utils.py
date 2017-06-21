#!/usr/bin/env python3
# coding: utf-8

"""
Util functions and classes.


@author gem
April 2017
"""

import os
import sys
import glob
import logging
from logging.handlers import RotatingFileHandler
import subprocess
import shutil
import shlex
import progressbar

logger = logging.getLogger()


def init_logger(logfile_base, level, name= None, verbose=0, quiet=False):
    """
    Create logger and its handlers, and set them to the given level

    level hierarchy:
    CRITICAL > ERROR > WARNING > INFO > DEBUG

    Messages from all levels are written in 'logfile'.log
    Messages for levels less than WARNING (only INFO and DEBUG) written to stdout
    Messages for levels equal or higher than WARNING written to stderr
    Messages for levels equal or higher than WARNING written in `logfile`.log.err

    level: minimum level that must be considered.
    name: if we need to name the logger (used for tests)
    verbose: be more verbose: 1 = add warnings in stderr (by default, only error),
    2 = like 1 + add DETAIL to stdout (by default only INFO)
    quiet: do not print anything to stderr/stdout
    """
    import time
    time_start = time.strftime("_%y-%m-%d_%H-%m-%S.log")
    # create logger
    if name:
        logger = logging.getLogger(name)
    else:  # pragma: no cover
        logger = logging.getLogger()

    # Determine logfile names
    logfile = logfile_base + ".log"
    if os.path.isfile(logfile):
        logfile = logfile_base + "-" + time_start + ".log"
    errfile = logfile_base + ".log.err"
    if os.path.isfile(errfile):
        errfile = logfile_base + "-" + time_start + ".log.err"
    detailfile = logfile_base + ".log.details"
    if os.path.isfile(detailfile):
        detailfile = logfile_base + "-" + time_start + ".log.details"

    # Create a new logging level: details (between info and debug)
    # Used to add details to the log file, but not to stdout, while still having
    # the possibility to put debug messages, used only for development.
    logging.addLevelName(15, "DETAIL")
    logging.DETAIL = 15
    def details(self, message, *args, **kws):
        if self.isEnabledFor(logging.DETAIL):
            self._log(logging.DETAIL, message, args, **kws)
    logging.Logger.details = details

    # set level of logger
    logger.setLevel(level)

    # create formatter for log messages: "timestamp :: level :: message"
    formatterFile = logging.Formatter('[%(asctime)s] :: %(levelname)s :: %(message)s',
                                      '%Y-%m-%d %H:%M:%S')
    formatterStream = logging.Formatter('  * [%(asctime)s] %(message)s', '%Y-%m-%d %H:%M:%S')

    # Create handler 1: writing to 'logfile'. mode 'write', max size = 1Mo.
    # If logfile is 1Mo, it is renamed to logfile.1, and next logs are still
    # written to logfile. Then, logfile.1 is renamed to logfile.2, logfile to
    # logfile.1 etc. We allow maximum 5 log files.
    # logfile contains everything from INFO level (INFO, WARNING, ERROR)
    logfile_handler = RotatingFileHandler(logfile, 'w', 1000000, 5)
    # set level to the same as the logger level
    logfile_handler.setLevel(logging.INFO)
    logfile_handler.setFormatter(formatterFile)  # add formatter
    logger.addHandler(logfile_handler)  # add handler to logger

    # Create handler 2: errfile. Write only warnings and errors
    errfile_handler = RotatingFileHandler(errfile, 'w', 1000000, 5)
    errfile_handler.setLevel(logging.WARNING)
    errfile_handler.setFormatter(formatterFile)  # add formatter
    logger.addHandler(errfile_handler)  # add handler to logger

    # Create handler 3: detailsfile. Write everything to this file.
    # Create it only if level is less than INFO. otherwise, it is the
    # same file as .log
    if level < logging.INFO:
        detfile_handler = RotatingFileHandler(detailfile, 'w', 1000000, 5)
        detfile_handler.setLevel(level)
        detfile_handler.setFormatter(formatterFile)  # add formatter
        logger.addHandler(detfile_handler)  # add handler to logger

    # If not quiet, add handlers for stdout and stderr
    if not quiet:
        # Create handler 4: write to stdout
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(logging.DEBUG)  # write any message
        # don't write messages >= WARNING
        stream_handler.addFilter(LessThanFilter(logging.WARNING))
        if verbose < 2:
            stream_handler.addFilter(NoLevelFilter(logging.DETAIL))
        stream_handler.setFormatter(formatterStream)
        logger.addHandler(stream_handler)  # add handler to logger

        # Create handler 5: write to stderr
        err_handler = logging.StreamHandler(sys.stderr)
        if verbose > 0:
            err_handler.setLevel(logging.WARNING)  # write all messages >= WARNING
        else:
            err_handler.setLevel(logging.ERROR)  # write all messages >= ERROR
        err_handler.setFormatter(formatterStream)
        logger.addHandler(err_handler)  # add handler to logger


class LessThanFilter(logging.Filter):
    """
    When using log, when a level is set to a handler, it is a minimum level. All
    levels higher than it will be printed. If you want to print only until
    a given level (no levels higher than the specified one), use this class like this:
    handler.addFilter(LessThanFilter(level))
    """
    def __init__(self, level):
        self._level = level
        logging.Filter.__init__(self)

    def filter(self, rec):
        return rec.levelno < self._level


class NoLevelFilter(logging.Filter):
    """
    When using log, specify a given level that must not be taken into account by the handler.
    This is used for the stdout handler. We want to print, by default,
    DEBUG (for development use) and INFO levels, but not DETAILS level (which is between
    DEBUG and INFO). We want to print DETAIL only if verbose option was set
    """
    def __init__(self, level):
        self._level = level
        logging.Filter.__init__(self)

    def filter(self, rec):
        return rec.levelno != self._level


def check_installed(cmd):
    """
    Check if the command 'cmd' is in $PATH and can then be executed

    return True or False
    """
    torun = "which " + cmd
    trying = subprocess.Popen(shlex.split(torun), stdout=subprocess.PIPE)
    out, _ = trying.communicate()
    if trying.returncode == 0:
        if os.path.isfile(out.strip()):
            return True
    return False


def run_cmd(cmd, error, eof=False, **kwargs):
    """
    Run the given command line. If the return code is not 0, print error message.
    if eof (exit on fail) is True, exit program if error code is not 0.
    """
    if not "stdout" in kwargs:
        kwargs["stdout"] = None
    if not "stderr" in kwargs:
        kwargs["stderr"] = None
    try:
        retcode = subprocess.call(shlex.split(cmd), stdout=kwargs["stdout"],
                                  stderr=kwargs["stderr"])
    except OSError:
        logger.error(error + ": " + "{} does not exist".format(cmd))
        if eof:
            sys.exit(-1)
        else:
            return -1
    if retcode != 0:
        logger.error(error)
        if eof:
            sys.exit(retcode)
    return retcode


def plot_distr(values, limit, outfile, title, text):
    """ Plot histogram of given 'values', and add a vertical line corresponding to the chosen
     'limit' and saves the image into the 'outfile'

    :param values: list of values
    :type values: list
    :param limit: limit for which a vertical line must be drawn
    :type limit: int
    :param outfile: file in which the output image must be saved
    :type outfile: str
    """
    import math
    import numpy as np
    import matplotlib
    matplotlib.use('AGG')
    from matplotlib import pyplot as plt
    plt.figure(figsize=(10,7))
    max_x = max(values)
    # if too many values, group them to have less bins in the histogram.
    # Put 'group_values' values in each bin ->
    # if less than 300 values, 1 value per bin, otherwise more values per bin
    group_values = int(max_x/300) + 1
    dec_ax = math.exp(0.001 * max_x) - 1
    dec_text = 3 * dec_ax
    bins = np.arange(0, max_x + 2*group_values, group_values) - 0.5
    axes = plt.hist(values, bins = bins, edgecolor="black", color="blue")
    plt.xlim(0.5, max_x + 0.5*group_values)
    plt.axvline(x=limit + 0.5*group_values + dec_ax, color="r")
    plt.text(x=limit + 0.5*group_values + dec_text, y=plt.ylim()[1]/2,
             s=text + " " + str(limit), color="r", rotation=90)
    plt.title(title)
    plt.savefig(outfile)
    plt.close("all")


def write_warning_skipped(skipped, format=False):
    """
    At the end of the script, write a warning to the user with the names of the genomes
    which had problems with prokka.

    skipped: list of genomes with problems

    format: if False, genomes were not skipped because of format step, but before that.
    if True, they were skipped because of format
    """
    logger = logging.getLogger()
    list_to_write = "\n".join(["\t- " + genome for genome in skipped])
    if not format:
        logger.warning(("Prokka had problems while annotating some genomes, or did not "
                        "find any gene. Hence, they are not "
                        "formatted, and absent from your output database. Please look at their "
                        "Prokka logs (<output_directory>/tmp_files/<genome_name>-prokka.log) and "
                        "to the current error log (<output_directory>/<input_filename>.log.err)"
                        " to get more information, and run again to annotate and format them. "
                        "Here are the genomes (problem with prokka or no "
                        "gene found): \n{}").format(list_to_write))
    else:
        logger.warning(("Some genomes were annotated by prokka, but could not be formatted, "
                        "and are hence absent from your output database. Please look at log "
                        "files to get more information about why they could not be "
                        "formatted.\n{}").format(list_to_write))


def write_discarded(genomes, kept_genomes, list_file, res_path, qc=False):
    """
    Write the list of genomes discarded to a file, so that users can keep a trace of them,
    with their information (nb contigs, L90 etc.)

    genomes: {genome: [gembase_start_name, seq_file, genome_size, nb_contigs, L90]}
    kept_genomes: list of genomes kept
    list_file: input file containing the list of genomes
    res_path: folder where results must be saved
    qc: if it is the file written after qc only, call it info-genomes-<list_file>.txt
    otherwise, call it discarded-<list_file>.txt
    """
    _, name_lst = os.path.split(list_file)
    if not qc:
        outdisc = os.path.join(res_path,
                               "discarded-" + ".".join(name_lst.split(".")[:-1]) + ".lst")
        logger.info("Writing discarded genomes to {}".format(outdisc))
    else:
        outdisc = os.path.join(res_path,
                               "info-genomes-" + ".".join(name_lst.split(".")[:-1]) + ".lst")
        logger.info("Writting information on genomes in {}".format(outdisc))
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
    return (x[1][0].split(".")[0], int(x[1][0].split(".")[-1]))


def read_genomes(list_file, name, date, dbpath):
    """
    Read list of genomes, and return them.
    If a genome has a name, also return it. Otherwise, return the name given by user.

    Check that the given genome file exists in dbpath. Otherwise, put an error message,
    and ignore this file.

    genomes = {genome: spegenus.date} spegenus.date = name.date
    """
    logger = logging.getLogger()
    logger.info("Reading genomes")
    genomes = {}
    if not os.path.isfile(list_file):
        logger.error(("ERROR: Your list file '{}' does not exist. "
                      "Please provide a list file.\n Ending program.").format(list_file))
        sys.exit(1)
    with open(list_file, "r") as lff:
        for line in lff:
            line = line.strip()
            # empty line: go to the next one
            if line == "":
                continue
            # If separator ::, look for species and/or date
            if "::" in line:
                genomes_inf, name_inf = line.split("::")
                genomes_inf = genomes_inf.strip()
                cur_name, cur_date = read_info(name_inf, name, date, genomes_inf)
            else:
                genomes_inf = line.strip()
                cur_name = name
                cur_date = date
            # If several file names, check that each one exists, and concatenate the existing files
            genomes_inf = genomes_inf.split()
            if len(genomes_inf) > 1:
                to_concat = []
                for file in genomes_inf:
                    if os.path.isfile(os.path.join(dbpath, file)):
                        to_concat.append(file)
                    else:
                        logger.warning(("{} genome file does not exist. Its file will be "
                                        "ignored when concatenating {}").format(file, genomes_inf))
                # If there are files to concatenate, concatenate them
                if to_concat != []:
                    genome_name = to_concat[0] + "-all.fna"
                    concat_file = os.path.join(dbpath, genome_name)
                    to_concat = [os.path.join(dbpath, gname) for gname in to_concat]
                    cat(to_concat, concat_file)
                else:
                    logger.warning(("None of the genome files in {} exist. "
                                    "This genome will be ignored.").format(genomes_inf))
                    genome_name = ""
            # If only 1 sequence file, check that it exists, and take its name
            else:
                if not os.path.isfile(os.path.join(dbpath, genomes_inf[0])):
                    logger.warning(("{} genome file does not exist. "
                                    "It will be ignored.").format(genomes_inf[0]))
                    genome_name = ""
                else:
                    genome_name = genomes_inf[0]
            if genome_name != "":
                genomes[genome_name] = [cur_name + "." + cur_date]
    return genomes


def read_info(name_inf, name, date, genomes_inf):
    """
    From the given information in 'name_inf', check if there is a name (and if its
    format is ok) and if there is a date (and if its format is ok).
    If no name (resp. no date), return default name (resp. default date).

    Return name and date.
    """
    name_inf = name_inf.strip().split(".")
    # if only species provided
    if len(name_inf) == 1:
        if name_inf[0] == "":
            cur_name = name
        elif check_format(name_inf[0]):
            cur_name = name_inf[0]
        else:
            logger.warning(("Invalid name {} given for genome {}. Only put "
                            "4 alphanumeric characters in your date and name. "
                            "For this genome, the default name ({}) will be "
                            "used.").format(name_inf[0], genomes_inf, name))
            cur_name = name
        cur_date = date
    elif len(name_inf) > 2:
        logger.warning(("Invalid name/date given for genome {}. Only put "
                        "4 alphanumeric characters in your date and name. For "
                        "this genome, the default name ({}) and date ({}) will "
                        "be used.").format(genomes_inf, name, date))
        cur_name = name
        cur_date = date
    else:
        cur_name, cur_date = name_inf
        if cur_name == "":
            cur_name = name
        if cur_date == "":
            cur_date = date
        if not check_format(cur_name):
            logger.warning(("Invalid name {} given for genome {}. Only put "
                            "4 alphanumeric characters in your date and name. "
                            "For this genome, the default name ({}) "
                            "will be used.").format(cur_name, genomes_inf, name))
            cur_name = name
        if not check_format(cur_date):
            logger.warning(("Invalid date {} given for genome {}. Only put "
                            "4 alphanumeric characters in your date and name. "
                            "For this genome, the default date ({}) "
                            "will be used.").format(cur_date, genomes_inf, date))
            cur_date = date
    return cur_name, cur_date


def cat(list_files, output, title = None):
    """
    Concatenate all files in 'list_files' and save result in 'output'
    Concat using shutil.copyfileobj, in order to copy by chunks, to
    avoid memory problems if files are big.
    """
    if title:
        nbfiles = len(list_files)
        widgets = [title + ': ', progressbar.Bar(marker='█', left='', right='', fill=' '),
               ' ', progressbar.Counter(), "/{}".format(nbfiles), ' (',
               progressbar.Percentage(), ") - ", progressbar.Timer()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbfiles, term_width=100).start()
        curnum = 1
    with open(output, "w") as outf:
        for file in list_files:
            if title:
                bar.update(curnum)
                curnum += 1
            with open(file, "r") as inf:
                shutil.copyfileobj(inf, outf)
    if title:
        bar.finish()


def check_format(info):
    """
    Check that the given information (can be the genomes name or the date) is in the right
    format: it should have 4 characters, all alphanumeric.
    """
    if len(info) != 4:
        return False
    return info.isalnum()


def check_out_dirs(resdir):
    """
    Check that there is no file in:
    - resdir/LSTINFO
    - resdir/Genes
    - resdir/Proteins
    - resdir/Replicons
    """
    if glob.glob(os.path.join(resdir, "LSTINFO", "*.lst")) != []:
        logger.error("ERROR: Your output directory already has .lst files in the "
                     "LSTINFO folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)
    if glob.glob(os.path.join(resdir, "Proteins", "*.prt")) != []:
        logger.error("ERROR: Your output directory already has .prt files in the "
                     "Proteins folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)
    if glob.glob(os.path.join(resdir, "Genes", "*.gen")) != []:
        logger.error("ERROR: Your output directory already has .gen files in the "
                     "Genes folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)
    if glob.glob(os.path.join(resdir, "Replicons", "*.fna")) != []:
        logger.error("ERROR: Your output directory already has .fna files in the "
                     "Replicons folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)


def rename_genome_contigs(gembase_name, gpath, outfile):
    """
    For the given genome (sequence in gpath), rename all its contigs
    with the new name: 'gembase_name', and save the output sequence in outfile
    """
    contig_num = 1
    with open(gpath, "r") as gpf, open(outfile, "w") as grf:
        for line in gpf:
            if line.startswith(">"):
                new_cont = ">" + gembase_name + "." + str(contig_num).zfill(4)
                contig_num += 1
                grf.write(new_cont + "\n")
            else:
                grf.write(line)

