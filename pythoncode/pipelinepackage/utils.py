#!/usr/bin/env python3
# coding: utf-8

"""
Util functions and classes.


@author gem
April 2017
"""

import os
import sys
import logging
import subprocess

logger = logging.getLogger()


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


def check_installed(cmd):
    """
    Check that the given command exists
    """
    FNULL = open(os.devnull, 'w')
    try:
        returncode = subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
        return returncode
    except Exception as err:
        logger.error(("{0} failed: {1}").format(cmd[0], err))
        sys.exit(1)


def plot_distr(values, limit, outfile, title, text):
    """ Plot histogram of given values, and add a vertical line corresponding to the choosen
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
        logger.warning(("Prokka had problems while annotating some genomes. Hence, they are not "
                        "formatted, and absent from your output database. Please look at their "
                        "Prokka logs (<output_directory>/tmp_files/<genome_name>-prokka.log) and "
                        "to the current error log (<output_directory>/<input_filename>.log.err)"
                        " to get more information, and run again to annotate and format them. "
                        "Here are the genomes: \n{}").format(list_to_write))
    else:
        logger.warning(("Some genomes were annotated by prokka, but could not be formatted, "
                        "and are hence absent from your output database. Please look at log "
                        "files to get more information about wh they could not be "
                        "formatted.\n{}").format(list_to_write))


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
    return (x[1][0].split(".")[0], int(x[1][0].split(".")[-1]))


def read_genomes(list_file, name, date):
    """
    Read list of genomes, and return them.
    If a genome has a name, also return it. Otherwise, return the name given by user.

    genomes = {genome: spegenus.date} spegenus.date = name.date
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
