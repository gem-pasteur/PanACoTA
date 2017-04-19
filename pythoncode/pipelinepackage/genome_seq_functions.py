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

from pipelinepackage import utils


logger = logging.getLogger()


def analyse_all_genomes(genomes, dbpath, res_path, nbn):
    """
    genomes: {genome: name}
    dbpath: path to folder containing genomes
    res_path: path to put out files
    nbn: minimum number of 'N' required to cut into a new contig

    returns: genomes {genome: [name, path, size, nbcont, l90]}
    """
    cut = nbn > 0
    pat = 'N' * nbn + "+"
    if cut:
        logger.info(("Cutting genomes at each stretch of at least {} 'N', "
                     "and then, calculating genome size, number of contigs and L90.").format(nbn))
    else:
        logger.info("Calculating genome size, number of contigs, L90")
    for genome, name in genomes.items():
        analyse_genome(genome, dbpath, res_path, cut, pat, genomes)


def analyse_genome(genome, dbpath, res_path, cut, pat, genomes):
    """
    Anayse given genome:
    - if cut is asked:
        - cut its contigs at each time that 'pat' is seen
        - save cut genome in new file
    - calculate genome size, L90, nb contigs and save it into genomes

    * genome: given genome to analyse
    * dbpath: path to the folder containing the given genome sequence
    * res_path: path to folder where output files must be saved. cut genome will be saved
    in res_path/tmp_files
    * cut: True if contigs must be cut, False otherwise
    * pat: pattern on which contigs must be cut
    * genomes: {genome: [name]} -> {genome: [name, path, gsize, nbcont, L90]}
    """
    gpath = os.path.join(dbpath, genome)
    if cut:
        tmp_res = os.path.join(res_path, "tmp_files")
        if not os.path.isdir(tmp_res):
            os.makedirs(tmp_res)
        grespath = os.path.join(tmp_res, genome + "-split{}N.fna".format(len(pat) - 1))
        gresf = open(grespath, "w")
    else:
        grespath = gpath
    with open(gpath, 'r') as genf:
        contig_sizes = {}  # {name: size}
        cur_contig_name = ""
        cur_contig = ""
        for line in genf:
            if line.startswith(">"):
                # If it is not the first contig, find stretch(es) of N in previous contig
                if cur_contig != "":
                    if cut:
                        save_contig(pat, cur_contig, cur_contig_name, contig_sizes, gresf)
                    else:
                        contig_sizes[cur_contig_name] = len(cur_contig)
                cur_contig_name = line.strip().split()[0]
                cur_contig = ""
            else:
                cur_contig += line.strip()
        if cut:
            save_contig(pat, cur_contig, cur_contig_name, contig_sizes, gresf)
        else:
            contig_sizes[cur_contig_name] = len(cur_contig)
        nbcont = len(contig_sizes)
        gsize = sum(contig_sizes.values())
        l90 = calc_l90(gsize, contig_sizes)
        genomes[genome] += [grespath, gsize, nbcont, l90]


def save_contig(pat, cur_contig, cur_contig_name, contig_sizes, gresf):
    """
    Save the contig read just before into dicts and write it to sequence file.

    pat: pattern to split a contig
    cur_contig: sequence of current contig, to save once split according to pat
    cur_contig_name: name of current contig to save once split according to pat
    contig_sizes: {name: size} save cur_contig once split according to pat
    gresf: file open in w mode to save the split sequence
    """
    # split contig each time a stretch of at least nbn 'N' is found
    cont_parts = re.split(pat, cur_contig)
    # save contig parts
    num = 0  # contig number
    for seq in cont_parts:
        # Only save non empty contigs (with some patterns, it could arrive that
        # we get empty contigs, if 2 occurrences of the pattern are side by side).
        if len(seq) == 0:
            continue
        cur_name = cur_contig_name + "_" + str(num)
        contig_sizes[cur_name] = len(seq)
        gresf.write(cur_name + "\n")
        gresf.write(seq + "\n")
        num += 1


def calc_l90(gsize, contig_sizes):
    """
    Calc L90 of a given genome
    """
    sizes = [contig_sizes[cont] for cont in contig_sizes]
    cum_sizes = np.cumsum(sorted(sizes, reverse=True))
    lim = 0.9 * gsize
    for num, val in enumerate(cum_sizes):
        if val >= lim:
            return num + 1

def rename_all_genomes(genomes, res_path):
    """
    Sort kept genomes by L90 and then nb contigs.
    For each genome, assign a strain number, and rename all its contigs.

    genomes: {genome: [name, path, gsize, nbcont, L90]} ->
    {genome: [gembase_name, path, gsize, nbcont, L90]}
    """
    last_name = ""
    last_strain = 0
    #"SAEN.1015.{}".format(str(last_strain).zfill(5))
    for genome, [name, gpath, size, nbcont, l90] in sorted(genomes.items(), key=sort_genomes):
        if last_name != name:
            last_strain = 1
            last_name = name
        else:
            last_strain += 1
        gembase_name = ".".join([name, str(last_strain).zfill(5)])
        genomes[genome][0] = gembase_name
        genomes[genome][1] = rename_genome_contigs(gembase_name, gpath, res_path)


def rename_genome_contigs(gembase_name, gpath, res_path):
    """
    For the given genome (sequence in gpath), rename all its contigs
    with the new name: 'gembase_name'.
    """
    contig_num = 1
    outfile = os.path.join(res_path, "tmp_files", os.path.basename(gpath) + "-gembase.fna")
    with open(gpath, "r") as gpf, open(outfile, "w") as grf:
        for line in gpf:
            if line.startswith(">"):
                new_cont = ">" + gembase_name + "." + str(contig_num).zfill(4)
                contig_num += 1
                grf.write(new_cont + "\n")
            else:
                grf.write(line)
    return outfile


def sort_genomes(x):
    """
    Sort all genomes with the following criteria:
    - sort by species
    - for each species, sort by l90
    - for same l90, sort by nb contigs
    """
    return (x[1][0], x[1][-1], x[1][-2])


def plot_distributions(genomes, res_path, listfile_base, l90, nbconts):
    """
    genomes: {genome: [name, path, size, nbcont, l90]}
    res_patj: path to put all output files
    listfile_base: name of list file
    l90: max value of l90
    nbconts: max value of nb contigs
    """
    L90_vals = [val for _, (_, _, _, _, val) in genomes.items()]
    outl90 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
    nbcont_vals = [val for _, (_, _, _, val, _) in genomes.items()]
    outnbcont = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
    utils.plot_distr(L90_vals, l90, outl90, "L90 distribution for all genomes",
                     "max L90 =")
    utils.plot_distr(nbcont_vals, nbconts, outnbcont,
                     "Distribution of number of contigs among all genomes", "max #contigs =")
