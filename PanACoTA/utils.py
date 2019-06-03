#!/usr/bin/env python3
# coding: utf-8

"""
Util functions and classes.


@author gem
April 2017
"""

import os
import sys
import re
import glob
import logging
from logging.handlers import RotatingFileHandler
import subprocess
import shutil
import shlex
import progressbar

try:
    import cPickle as pickle
except:
    try:
        import _pickle as pickle
    except:  # pragma: no cover
        import pickle


def init_logger(logfile_base, level, name, details=False, verbose=0, quiet=False):
    """
    Create logger and its handlers, and set them to the given level

    level hierarchy: ``CRITICAL > ERROR > WARNING > INFO > DETAILS > DEBUG``

    Messages from all levels are written in 'logfile'.log

    Messages for levels less than WARNING (only INFO and DEBUG) written to stdout

    Messages for levels equal or higher than WARNING written to stderr

    Messages for levels equal or higher than WARNING written in `logfile`.log.err


    Parameters
    ----------
    logfile_base : str
        base of filename to use for logs. Will add '.log', '.log.details' and '.log.err' for\
        the 3 log files created
    level : int
        minimum level that must be considered.
    name : str or None
        if we need to name the logger (used for tests)
    verbose : int
        be more verbose:
        default (0): info in stdout, error and more in stderr
        1 = add warnings in stderr
        2 = like 1 + add DETAIL to stdout (by default only INFO)
        >15: add debug to stdout
    quiet : bool
        True if nothing must be sent to stdout/stderr, False otherwise
    """
    import time
    time_start = time.strftime("_%y-%m-%d_%H-%m-%S")
    # create logger
    logger = logging.getLogger(name)

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
    debugfile = logfile_base + ".log.debug"
    if os.path.isfile(debugfile):
        debugfile = logfile_base + "-" + time_start + ".log.debug"

    # Create a new logging level: details (between info and debug)
    # Used to add details to the log file, but not to stdout, while still having
    # the possibility to put debug messages, used only for development.
    logging.addLevelName(detail_lvl(), "DETAIL")
    logging.DETAIL = detail_lvl()

    def details(self, message, *args, **kws):
        """
        Define a new log level: details
        """
        if self.isEnabledFor(logging.DETAIL):
            self._log(logging.DETAIL, message, args, **kws)

    logging.Logger.details = details

    # set level of logger
    logger.setLevel(level)

    # create formatter for log messages: "timestamp :: level :: message"
    # :: %(name)s  to add the logger name
    # my_format = '[%(asctime)s] :: from %(name)s %(levelname)s :: %(message)s'
    my_format = '[%(asctime)s] :: %(levelname)s :: %(message)s'
    formatter_file = logging.Formatter(my_format,
                                       '%Y-%m-%d %H:%M:%S')
    formatter_stream = logging.Formatter('  * [%(asctime)s] %(message)s', '%Y-%m-%d %H:%M:%S')

    # Create handler 1: writing to 'logfile'. mode 'write', max size = 1Mo.
    # If logfile is 1Mo, it is renamed to logfile.1, and next logs are still
    # written to logfile. Then, logfile.1 is renamed to logfile.2, logfile to
    # logfile.1 etc. We allow maximum 5 log files.
    # logfile contains everything from INFO level (INFO, WARNING, ERROR)
    logfile_handler = RotatingFileHandler(logfile, 'w', 10000000, 5)
    # set level to the same as the logger level
    logfile_handler.setLevel(logging.INFO)
    logfile_handler.setFormatter(formatter_file)  # add formatter
    logger.addHandler(logfile_handler)  # add handler to logger

    # Create handler 2: errfile. Write only warnings and errors
    errfile_handler = RotatingFileHandler(errfile, 'w', 10000000, 5)
    errfile_handler.setLevel(logging.WARNING)
    errfile_handler.setFormatter(formatter_file)  # add formatter
    logger.addHandler(errfile_handler)  # add handler to logger

    # Create handler 3: detailsfile. Write everything to this file, except debug
    # Create it only if:
    # - level is <= info (for modules which have no details, so detailsfile is the same as
    # logfile)
    # - details==True force creation of detailsfile
    # - quiet==True nothing in stdout, put all log files so that user can check
    if level < logging.INFO or quiet or details:
        detfile_handler = RotatingFileHandler(detailfile, 'w', 10000000, 5)
        detfile_handler.setLevel(logging.DETAIL)
        detfile_handler.setFormatter(formatter_file)  # add formatter
        logger.addHandler(detfile_handler)  # add handler to logger

    # Create handler 3: debug file. Write everything
    if level < logging.DETAIL:
        debugfile_handler = RotatingFileHandler(debugfile, 'w', 10000000, 5)
        debugfile_handler.setLevel(logging.DEBUG)
        debugfile_handler.setFormatter(formatter_file)  # add formatter
        logger.addHandler(debugfile_handler)  # add handler to logger

    # If not quiet, add handlers for stdout and stderr
    if not quiet:
        # Create handler 4: write to stdout
        stream_handler = logging.StreamHandler(sys.stdout)
        # By default, write everything
        stream_handler.setLevel(logging.DEBUG)
        # BUT: don't write messages >= WARNING (warning, error, critical)
        stream_handler.addFilter(LessThanFilter(logging.WARNING))
        # if not verbose (level 0 or 1): only put info in stdout (remove details and debug)
        if verbose < 2:
            stream_handler.addFilter(NoLevelFilter(logging.DETAIL))
            stream_handler.addFilter(NoLevelFilter(logging.DEBUG))
        # if verbose (level 2): put info and details in stdout: only remove debug
        if verbose < 2:
            stream_handler.addFilter(NoLevelFilter(logging.DEBUG))
        stream_handler.setFormatter(formatter_stream)
        logger.addHandler(stream_handler)  # add handler to logger

        # Create handler 5: write to stderr
        err_handler = logging.StreamHandler(sys.stderr)

        if verbose > 0:
            err_handler.setLevel(logging.WARNING)  # write all messages >= WARNING
        else:
            err_handler.setLevel(logging.ERROR)  # write all messages >= ERROR
        err_handler.setFormatter(formatter_stream)
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
        """
        Function to decide if given log has to be logged or not, according to its level

        Parameters
        ----------
        rec : current record handled by logger

        Returns
        -------
        bool
            True if level of current log is less than the defined limit, False otherwise
        """
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
        """
        Function to decide if given log has to be logged or not, according to its level

        Parameters
        ----------
        rec : current record handled by logger

        Returns
        -------
        bool
            True if level of current log is different from forbidden level, False if it is the same
        """
        return rec.levelno != self._level


def check_installed(cmd):
    """
    Check if the command 'cmd' is in $PATH and can then be executed

    Parameters
    ----------
    cmd : str
        command to run

    Returns
    -------
    bool
        True if installed, False otherwise
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

    Parameters
    ----------
    cmd : str
        command to run
    error : str
        error message to print if error while running command
    eof : bool
        True: exit program if command failed, False: do not exit even if command fails
    kwargs : Object
        Can provide a logger, stdout and/or stderr streams

    Returns
    -------
    subprocess.Popen
        returns object of subprocess call (has attributes returncode, pid, communicate etc.)

    """
    if "logger" not in kwargs:
        logger = logging.getLogger("utils.run_cmd")
    else:
        logger = kwargs["logger"]
    if "stdout" not in kwargs:
        kwargs["stdout"] = None
    if "stderr" not in kwargs:
        kwargs["stderr"] = None
    try:
        call = subprocess.Popen(shlex.split(cmd), stdout=kwargs["stdout"],
                                stderr=kwargs["stderr"])
        call.wait()
        retcode = call.returncode
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
    return call


def plot_distr(values, limit, title, text):
    """
    Plot histogram of given 'values', and add a vertical line corresponding to the chosen
    'limit' and return the mpl figure

    Parameters
    ----------
    values : list
        list of values
    limit : int
        limit for which a vertical line must be drawn
    title : str
        Title to give to plot
    text : str
        text to write near the vertical line representing the limit

    Returns
    -------
    matplotlib.figure.Figure
        figure generated
    """
    import math
    import numpy as np
    import matplotlib
    matplotlib.use('AGG')
    from matplotlib import pyplot as plt
    plt.close("all")
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(1, 1, 1)
    max_x = max(values)
    # if too many values, group them to have less bins in the histogram.
    # Put 'group_values' values in each bin ->
    # if less than 300 values, 1 value per bin, otherwise more values per bin
    group_values = int(max_x / 300) + 1
    dec_ax = math.exp(0.001 * max_x) - 1
    dec_text = 3 * dec_ax
    bins = np.arange(0, max_x + 2 * group_values, group_values) - 0.5
    ax.hist(values, bins=bins, edgecolor="black", color="blue")
    ax.set_xlim(0.5, max_x + 0.5 * group_values)
    ax.axvline(x=limit + 0.5 * group_values + dec_ax, color="r")
    ax.text(x=limit + 0.5 * group_values + dec_text, y=plt.ylim()[1] / 2,
            s=text + " " + str(limit), color="r", rotation=90)
    ax.set_title(title)
    return fig


def write_warning_skipped(skipped, do_format=False, prodigal_only=False, logfile=""):
    """
    At the end of the script, write a warning to the user with the names of the genomes
    which had problems with prokka.

    Parameters
    ----------
    skipped : list
        list of genomes with problems
    do_format : bool
        if False, genomes were not skipped because of format step, but before that.\
        if True, they were skipped because of format
    prodigal_only : bool
        if False: used prokka to annotate
        if True: used prodigal to annotate
    """
    if not prodigal_only:
        soft = "Prokka"
    else:
        soft = "Prodigal"
    logger = logging.getLogger("utils")
    list_to_write = "\n".join(["\t- " + genome for genome in skipped])
    if not do_format:
        logger.warning(("{0} had problems while annotating some genomes, or did not "
                        "find any gene. Hence, they are not "
                        "formatted, and absent from your output database. Please look at their "
                        "{0} logs (<output_directory>/tmp_files/<genome_name>-{0}.log and "
                        ".log.err) "
                        " and to the current error log "
                        "<output_directory>/<input_filename>.log.err)"
                        " to get more information, and run again to annotate and format them. "
                        "Here are the genomes (problem with {1} or no "
                        "gene found): \n{0}").format(soft, list_to_write))
    else:
        logger.info("WARNING: Some genomes could not be formatted. See {}".format(logfile))
        logger.warning(("Some genomes were annotated by {}, but could not be formatted, "
                        "and are hence absent from your output database. Please look at "
                        "log.details files to get more information about why they could not be "
                        "formatted.\n{}").format(soft, list_to_write))


def write_discarded(genomes, kept_genomes, list_file, res_path, qc=False):
    """
    Write the list of genomes discarded to a file, so that users can keep a trace of them,
    with their information (nb contigs, L90 etc.)

    Parameters
    ----------
    genomes : dict
        {genome: [gembase_start_name, seq_file, genome_size, nb_contigs, L90]}
    kept_genomes : list
        list of genomes kept
    list_file : str
        input file containing the list of genomes
    res_path : str
        folder where results must be saved
    qc : bool
        if it is the file written after qc only (True), call it info-genomes-<list_file>.txt\
        otherwise (False), call it discarded-<list_file>.txt

    """
    logger = logging.getLogger("utils")
    nb_disc = len(genomes) - len(kept_genomes)
    if not qc and nb_disc < 2:
        logger.info("{} genome was discarded.".format(nb_disc))
    elif not qc:
        logger.info("{} genomes were discarded.".format(nb_disc))
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

    Parameters
    ----------
    list_file : str
        input file containing the list of genomes
    genomes : dict
        {genome: [gembase_start_name, seq_file, genome_size, nb_contigs, L90]}
    outdir : str
        folder where results must be saved

    """
    _, name_lst = os.path.split(list_file)

    outlst = os.path.join(outdir, "LSTINFO-" + ".".join(name_lst.split(".")[:-1]) + ".lst")
    with open(outlst, "w") as outf:
        outf.write("\t".join(["gembase_name", "orig_name", "gsize", "nb_conts", "L90"]) + "\n")
        for genome, values in sorted(genomes.items(), key=sort_genomes):
            gembase, _, gsize, nbcont, l90 = [str(x) for x in values]
            outf.write("\t".join([gembase, genome, gsize, nbcont, l90]) + "\n")


def sort_genomes(x):
    """
    order by:

        - species
        - in each species, by strain number

    Parameters
    ----------
    x : tuple or str
        [genome_orig, [gembase, path, gsize, nbcont, L90]] with gembase = species.date.strain

    Returns
    -------
    str
        variable to take into account for sorting. If format is ESCO.1512.00001 return\
        ESCO and 00001. Otherwise, just return x itself (sort by alphabetical order)
    """
    if isinstance(x, tuple):
        x = x[1][0]

    # if format is ESCO.1512.00001 sort by ESCO, then 00001
    if "." in x and len(x.split(".")) >= 3:
        return x.split(".")[0], int(x.split(".")[-1])
    # if format is not like this, just return alphabetical order
    return x,


def sort_proteins(x):
    """
    order by:

    - species
    - in each species, strain number
    - in each species and strain number, by protein number

    Parameters
    ----------
    x : str
        species.date.strain.contig_protnum

    Returns
    -------
    str
        variable to take into account for sorting. If format is ESCO.1512.00001.i0002_12124,\
        return ESCO, 00001 and 12124. If not, it must be something_00001:\
        return something and 00001.
    """
    try:
        # if format is ESCO.1512.00001.i0002_12124, sort by ESCO, then 00001, then 12124
        if "." in x and len(x.split(".")) >= 3:
            return x.split(".")[0], int(x.split(".")[2].split("_")[0]), int(x.split("_")[-1])
        # if format is not like this, it must be something_00001:
        # sort by 'something' and then 00001
        return "_".join(x.split("_")[:-1]), int(x.split("_")[-1])
    except (IndexError, ValueError):
        logger = logging.getLogger("utils")
        logger.error(("ERROR: Protein {} does not have the required format. "
                      "It must contain, at least <alpha-num>_<num_only>, and at best "
                      "<name>.<date>.<strain_num>.<contig_info>_<prot_num>. "
                      "Please change its name.").format(x))
        sys.exit(1)


def read_genomes(list_file, name, date, dbpath, tmp_path):
    """
    Read list of genomes, and return them.
    If a genome has a name, also return it. Otherwise, return the name given by user.

    Check that the given genome file exists in dbpath. Otherwise, put an error message,
    and ignore this file.

    Parameters
    ----------
    list_file : str
        input file containing the list of genomes
    name : str
        Default species name
    date : str
        Default date
    dbpath : str
        path to folder containing original genome files
    tmp_path : str
        path to folder which will contain the genome files to use before annotation, if\
        needed to change them from original file (for example, merging several contig files\
        in one file, split at each stretch of 5 'N', etc.).

    Returns
    -------
    dict
        {genome: spegenus.date} spegenus.date = name.date
    """
    logger = logging.getLogger("utils")
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
                if to_concat:
                    genome_name = to_concat[0] + "-all.fna"
                    concat_file = os.path.join(tmp_path, genome_name)
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

    Parameters
    ----------
    name_inf : str
        information on current genome, which could contain a species name and a date
    name : str
        default species name
    date : str
        default date
    genomes_inf : str
        current genome filename. Used to complete information when there is a warning (species\
        name or date given not in the right format...)

    Returns
    -------
    (cur_name, cur_date) : tuple
        with:

        - curname: name to use for this genome (can be the default one, or the one read from\
        'name_inf'
        - curdate: date to use for this genome (default or read from 'name_inf')
    """
    logger = logging.getLogger("utils")
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


def cat(list_files, output, title=None):
    """
    Equivalent of 'cat' unix command.

    Concatenate all files in 'list_files' and save result in 'output'
    Concat using shutil.copyfileobj, in order to copy by chunks, to
    avoid memory problems if files are big.

    Parameters
    ----------
    list_files : list
        list of filenames to concatenate
    output : str
        output filename, where all concatenated files will be written
    title : str or None
        if you want to show a progressbar while concatenating files, add a title for this\
        progressbar here. If no title, nothing will be shown during concatenation.

    """
    bar = None
    curnum = None
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


def grep(filein, pattern, counts=False):
    """
    Equivalent of 'grep' unix command

    By default, returns all the lines containing the given pattern.
    If counts = True, returns the number of lines containing the pattern.

    Parameters
    ----------
    filein : str
        path to the file in which pattern must be searched
    pattern : str
        pattern to search
    counts : bool
        True if you want to count how many lines have the pattern and return this number,\
        False if you want to return all lines containing the pattern.

    Returns
    -------
    list or int
        list of lines if counts=False; number of lines if counts=True
    """
    num = 0
    lines = []
    with open(filein, "r") as inf:
        for line in inf:
            if re.search(pattern, line):
                lines.append(line.strip())
                num += 1
    if counts:
        return num
    else:
        return lines


def count(filein, get="lines"):
    """
    Similar to 'wc' unix command.

    Count the number of what is given in 'get'. It can be:

    - lines (default)
    - words

    Parameters
    ----------
    filein : str
        path to the file for which we want to count lines or words
    get : ["lines", "words"]
        either lines to count the number of lines in the file, or words to count the number\
        of words.

    Returns
    -------
    int
        Number of lines or words according to value of 'get' parameter.
    """
    gets = ["lines", "words"]
    if get not in gets:
        logger = logging.getLogger("utils")
        logger.error("Choose what you want to count among {}.".format(gets))
        sys.exit(1)
    num = 0
    with open(filein, "r") as inf:
        for line in inf:
            if get == "lines":
                num += 1
            elif get == "words":
                num += len(line.split())
    return num


def check_format(info):
    """
    Check that the given information (can be the genomes name or the date) is in the right
    format: it should have 4 characters, all alphanumeric.

    Parameters
    ----------
    info : str
        information to check

    Returns
    -------
    bool
        True if right format, False otherwise
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
    - resdir/gff3

    Parameters
    ----------
    resdir : str
        path to result directory

    """
    logger = logging.getLogger("utils")
    if glob.glob(os.path.join(resdir, "LSTINFO", "*.lst")):
        logger.error("ERROR: Your output directory already has .lst files in the "
                     "LSTINFO folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)
    if glob.glob(os.path.join(resdir, "Proteins", "*.prt")):
        logger.error("ERROR: Your output directory already has .prt files in the "
                     "Proteins folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)
    if glob.glob(os.path.join(resdir, "Genes", "*.gen")):
        logger.error("ERROR: Your output directory already has .gen files in the "
                     "Genes folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)
    if glob.glob(os.path.join(resdir, "Replicons", "*.fna")):
        logger.error("ERROR: Your output directory already has .fna files in the "
                     "Replicons folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)
    if glob.glob(os.path.join(resdir, "gff3", "*.gff")):
        logger.error("ERROR: Your output directory already has .gff files in the "
                     "gff3 folder. Provide another result directory, or remove the "
                     "files in this one.\nEnding program.")
        sys.exit(1)


def rename_genome_contigs(gembase_name, gpath, outfile):
    """
    For the given genome (sequence in gpath), rename all its contigs
    with the new name: 'gembase_name', and save the output sequence in outfile

    Parameters
    ----------
    gembase_name : str
        genome name to use (species.date.strain)
    gpath : str
        path to the genome sequence
    outfile : str
        path to the new file, containing 'gpath' sequence, but with 'gembase_name' in headers

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


def logger_thread(q):
    """
    Queue listener used in a thread to handle the logs put to a QueueHandler
    by several processes (multiprocessing.pool.map_async for example)

    Parameters
    ----------
    q : multiprocessing.managers.AutoProxy[Queue]
        queue to listen

    """
    while True:
        record = q.get()
        if record is None:
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)


def detail_lvl():
    """
    Get the int level corresponding to "DETAIL"

    Returns
    -------
    int
        int corresponding to the level "DETAIL"
    """
    return 15


def save_bin(objects, fileout):
    """
    Save python 'objects' in a binary file called 'fileout'

    Parameters
    ----------
    objects : Object
        python object to save
    fileout : str
        path to binary file where objects must be saved

    """
    with open(fileout, "wb") as binf:
        pickle.dump(objects, binf)


def load_bin(binfile):
    """
    Unpickle python objects from the binary file 'binfile'

    Parameters
    ----------
    binfile : str
        path to binary file containing python object

    Returns
    -------
    Object
        The python objects unpickled

    """
    with open(binfile, "rb") as binf:
        objects = pickle.load(binf)
    return objects


def write_list(liste):
    """
    Return a string corresponding to the given liste, with all elements separated
    by a space. Used to write a list into a file. Ex::

        [1, 2, "toto"] -> "1 2 toto"

    Parameters
    ----------
    liste : list
        list of elements that we would like to write

    Returns
    -------
    str
        the string to write
    """
    list_write = [str(l) for l in liste]
    return " ".join(list_write) + "\n"


def remove(infile):
    """
    Remove the given file if it exists

    Parameters
    ----------
    infile : str
        path to file to remove

    """
    if os.path.isfile(infile):
        os.remove(infile)
