#!/usr/bin/env python3
# coding: utf-8

"""
Functions to:
- analyse a genome (cut at stretch of N if asked, calc L90, nb cont, size...)
- if genome cut by stretch of N, write the new sequence to new file
- rename genome contigs and write new sequence to new file


@author gem
April 2017
"""
import os
import re
import numpy as np
import logging
import progressbar

from genomeAPCAT import utils


logger = logging.getLogger("qc_annote.gseq")


def analyse_all_genomes(genomes, dbpath, tmp_path, nbn, quiet=False):
    """
    genomes: {genome: spegenus.date}
    dbpath: path to folder containing genomes
    tmp_path: path to put out files
    nbn: minimum number of 'N' required to cut into a new contig

    returns: genomes {genome: [spegenus.date, path, size, nbcont, l90]}
    """
    cut = nbn > 0
    pat = 'N' * nbn + "+"
    nbgen = len(genomes)
    if cut:
        logger.info(("Cutting genomes at each stretch of at least {} 'N', "
                     "and then, calculating genome size, number of contigs and L90.").format(nbn))
    else:
        logger.info("Calculating genome size, number of contigs, L90")
    if not quiet:
        # Create progressbar
        widgets = ['Analysis: ', progressbar.Bar(marker='█', left='', right='', fill=' '),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ') - ', progressbar.Timer(), ' - ',
                   progressbar.ETA()
                  ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=100).start()
        curnum = 1
    toremove = []
    for genome, name in genomes.items():
        if not quiet:
            bar.update(curnum)
            curnum += 1
        res = analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes)
        if not res:
            toremove.append(genome)
    if toremove:
        for gen in toremove:
            del genomes[gen]
    if not quiet:
        bar.finish()


def analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes):
    """
    Analyse given genome:
    - if cut is asked:
        - cut its contigs at each time that 'pat' is seen
        - save cut genome in new file
    - calculate genome size, L90, nb contigs and save it into genomes

    * genome: given genome to analyse
    * dbpath: path to the folder containing the given genome sequence
    * tmp_path: path to folder where output files must be saved.
    * cut: True if contigs must be cut, False otherwise
    * pat: pattern on which contigs must be cut
    * genomes: {genome: [spegenus.date]} -> {genome: [spegenus.date, path, gsize, nbcont, L90]}
    """
    gpath = os.path.join(dbpath, genome)  # genome file is in dbpath
    # except if it was in several files in dbpath, in which case it has been concatenated to
    # a file in tmp_path
    if not os.path.isfile(gpath):
        gpath = os.path.join(tmp_path, genome)
    if cut:
        grespath = os.path.join(tmp_path, genome + "-split{}N.fna".format(len(pat) - 1))
    else:
        grespath = os.path.join(tmp_path, genome + "-short-contig.fna")
    with open(gpath, 'r') as genf, open(grespath, "w") as gresf:
        contig_sizes = {}  # {contig_name: size}
        cur_contig_name = ""
        cur_contig = ""
        num = 0
        for line in genf:
            if line.startswith(">"):
                # If it is not the first contig, find stretch(es) of N in previous contig
                if cur_contig != "":
                    if cut:
                        num = save_contig(pat, cur_contig, cur_contig_name,
                                          contig_sizes, gresf, num)
                    else:
                        gresf.write(">" + cur_contig_name[:15] + "_" + str(num) + "\n")
                        gresf.write(cur_contig + "\n")
                        num += 1
                        contig_sizes[cur_contig_name] = len(cur_contig)
                cur_contig_name = "_".join(line.strip().split())
                cur_contig = ""
            else:
                cur_contig += line.strip().upper()
        if cur_contig != "":
            if cut:
                num = save_contig(pat, cur_contig, cur_contig_name, contig_sizes, gresf, num)
            else:
                gresf.write(">" + cur_contig_name[:15] + "_" + str(num) + "\n")
                gresf.write(cur_contig + "\n")
                contig_sizes[cur_contig_name] = len(cur_contig)
    nbcont = len(contig_sizes)
    gsize = sum(contig_sizes.values())
    l90 = calc_l90(contig_sizes)
    if not l90:
        logger.error("Problem with genome {}. Not possible to get L90".format(genome))
        return False
    else:
        genomes[genome] += [grespath, gsize, nbcont, l90]
    return True


def save_contig(pat, cur_contig, cur_contig_name, contig_sizes, gresf, num):
    """
    Save the contig read just before into dicts and write it to sequence file.
    Contig name must be at most 20 characters (required by prokka)

    pat: pattern to split a contig
    cur_contig: sequence of current contig, to save once split according to pat
    cur_contig_name: name of current contig to save once split according to pat
    contig_sizes: {name: size} save cur_contig once split according to pat
    gresf: file open in w mode to save the split sequence
    """
    # split contig each time a stretch of at least nbn 'N' is found
    cont_parts = re.split(pat, cur_contig)
    # save contig parts
    for seq in cont_parts:
        # Only save non empty contigs (with some patterns, it could arrive that
        # we get empty contigs, if 2 occurrences of the pattern are side by side).
        if len(seq) == 0:
            continue
        num += 1
        cur_name = cur_contig_name[:15] + "_" + str(num)
        contig_sizes[cur_name] = len(seq)
        gresf.write(cur_name + "\n")
        gresf.write(seq + "\n")
    return num


def calc_l90(contig_sizes):
    """
    Calc L90 of a given genome
    """
    gsize = sum(contig_sizes.values())
    sizes = [contig_sizes[cont] for cont in contig_sizes]
    cum_sizes = np.cumsum(sorted(sizes, reverse=True))
    lim = 0.9 * gsize
    for num, val in enumerate(cum_sizes):
        if val >= lim:
            return num + 1


def rename_all_genomes(genomes, tmp_path):
    """
    Sort kept genomes by L90 and then nb contigs.
    For each genome, assign a strain number, and rename all its contigs.

    genomes: {genome: [name, path, gsize, nbcont, L90]} ->
    genomes: {genome: [gembase_name, path, gsize, nbcont, L90]}
    """
    logger.info("Renaming kept genomes according to their quality.")
    last_name = ""
    last_strain = 0
    #"SAEN.1015.{}".format(str(last_strain).zfill(5))
    for genome, [name, gpath, size, nbcont, l90] in sorted(genomes.items(), key=sort_genomes):
        if last_name != name.split(".")[0]:
            last_strain = 1
            last_name = name.split(".")[0]
        else:
            last_strain += 1
        gembase_name = ".".join([name, str(last_strain).zfill(5)])
        genomes[genome][0] = gembase_name
        outfile = os.path.join(tmp_path, os.path.basename(gpath) + "-gembase.fna")


def sort_genomes(x):
    """
    Sort all genomes with the following criteria:
    - sort by species (x[1][0] is species.date)
    - for each species, sort by l90
    - for same l90, sort by nb contigs
    """
    return (x[1][0].split(".")[0], x[1][-1], x[1][-2])


def plot_distributions(genomes, res_path, listfile_base, l90, nbconts):
    """
    genomes: {genome: [name, path, size, nbcont, l90]}
    res_path: path to put all output files
    listfile_base: name of list file
    l90: max value of l90 allowed
    nbconts: max value of nb contigs allowed
    """
    logger.info("Generating distribution of L90 and #contigs graphs.")
    L90_vals = [val for _, (_, _, _, _, val) in genomes.items()]
    outl90 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
    nbcont_vals = [val for _, (_, _, _, val, _) in genomes.items()]
    outnbcont = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
    dist1 = utils.plot_distr(L90_vals, l90, "L90 distribution for all genomes",
                             "max L90 =")
    dist2 = utils.plot_distr(nbcont_vals, nbconts,
                             "Distribution of number of contigs among all genomes",
                             "max #contigs =")
    dist1.savefig(outl90)
    dist2.savefig(outnbcont)
    return L90_vals, nbcont_vals, dist1, dist2

