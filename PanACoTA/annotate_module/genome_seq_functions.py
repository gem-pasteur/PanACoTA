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
Functions to:

- analyse a genome (cut at stretch of N if asked, calc L90, nb cont, size...)
- if genome cut by stretch of N, write the new sequence to new file
- rename genome contigs and write new sequence to new file


@author gem
April 2017
"""
import os
import re
import sys
import numpy as np
import logging
import progressbar

from PanACoTA import utils

logger = logging.getLogger("annotate.gseq_functions")

def analyse_all_genomes(genomes, dbpath, tmp_path, nbn, soft, logger, quiet=False):
    """

    Parameters
    ----------
    genomes : dict
        {genome: spegenus.date}
    dbpath : str
        path to folder containing genomes
    tmp_path : str
        path to put out files
    nbn : int
        minimum number of 'N' required to cut into a new contig
    soft : str
        soft used (prokka, prodigal, or None if called by prepare module)
    logger : logging.Logger
        logger object to write log information. Because this function can be called from
        prepare module, where sub logger name is different
    quiet : bool
        True if nothing must be written to stdout/stderr, False otherwise

    Returns
    -------
    dict
        {genome: [spegenus.date, orig_name, path_to_seq_to_annotate, size, nbcont, l90]}

    """
    cut = nbn > 0
    pat = None  ## To put pattern with which sequence must be cut
    if cut:
        pat = 'N' * nbn + "+"
    nbgen = len(genomes)
    bar = None
    curnum = None
    if cut:
        logger.info(("Cutting genomes at each time there are at least {} 'N' in a row, "
                     "and then, calculating genome size, number of contigs and L90.").format(nbn))
    else:
        logger.info("Calculating genome size, number of contigs, L90")
    if not quiet:
        # Create progressbar
        widgets = ['Analysis: ', progressbar.Bar(marker='█', left='', right=''),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ') - ', progressbar.Timer(), ' - ',
                   progressbar.ETA()
                   ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=79).start()
        curnum = 1
    toremove = []
    # Analyse genomes 1 by 1
    for genome, name in genomes.items():
        # If not quiet option, show progress bar
        if not quiet:
            bar.update(curnum)
            curnum += 1
        # analyse genome, and check everything went well.
        # exception if binary file
        try:
            res = analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes, soft, logger=logger)
        except UnicodeDecodeError:
            logger.warning(f"'{genome}' does not seem to be a fasta file. It will be ignored.")
            res = False
        # Problem while analysing genome -> genome ignored
        if not res:
            toremove.append(genome)
    # If there are some genomes to remove (analysis failed), remove them from genomes dict.
    if toremove:
        for gen in toremove:
            del genomes[gen]
    if not genomes:
        logger.error(f"No genome was found in the database folder {dbpath}. See logfile "
                     "for more information.")
        sys.exit(1)
    if not quiet:
        bar.finish()
    return 0


def analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes, soft, logger):
    """
    Analyse given genome:

    - if cut is asked:

        - cut its contigs at each time that 'pat' is seen
        - save cut genome in new file

    - calculate genome size, L90, nb contigs and save it into genomes

    Parameters
    ----------
    genome : str
        given genome to analyse
    dbpath : str
        path to the folder containing the given genome sequence
    tmp_path : str
        path to folder where output files must be saved.
    cut : bool
        True if contigs must be cut, False otherwise
    pat : str
        pattern on which contigs must be cut. ex: "NNNNN"
    genomes : dict
        {genome_file: [genome_name]} as input, and will be changed to\
        {genome_file: [genome_name, path, path_annotate, gsize, nbcont, L90]}
    soft : str
        soft used (prokka, prodigal, or None if called by prepare module)

    Returns
    -------
    bool
        True if genome analysis went well, False otherwise
        Modifies 'genomes' for the analysed genome: -> {genome_file: [genome_name, path,
        path_annotate, gsize, nbcont, L90]}
    """
    gpath, grespath = get_output_dir(soft, dbpath, tmp_path, genome, cut, pat)
    if not os.path.exists(gpath):
        logger.error(f"The file {gpath} does not exist")
        return False
    # Open original sequence file
    with open(gpath, "r") as genf:
        # If a new file must be created (sequences cut), open it
        gresf = None
        if grespath:
            gresf = open(grespath, "w")

        # Initialize variables
        cur_contig_name = "" # header text
        contig_sizes = {}  # {header text: size}
        cur_seq = "" # sequence
        num = 1 # Used to get unique contig names

        # Read each line of original sequence
        for line in genf:
            #### NEW CONTIG
            # Line corresponding to a new contig
            if line.startswith(">"):
                # If not first contig, write info to output file (if needed)
                if cur_seq != "":
                    num = format_contig(cut, pat, cur_seq, cur_contig_name, genome, contig_sizes,
                                        gresf, num, logger)
                    # If problem while formatting contig, return False -> genome ignored
                    if num == -1:
                        return False
                # Get contig name for next turn, and reset sequence
                cur_contig_name = line.strip()
                # Initialize for next contig
                cur_seq = ""
            # #### SEQUENCE LINE
            # If sequence line, keep it, all in upper case
            else:
                # Add this line without \n to sequence (1 sequence per line)
                cur_seq += line.strip().upper()

        # LAST CONTIG
        if cur_contig_name != "":
            num = format_contig(cut, pat, cur_seq, cur_contig_name, genome, contig_sizes, gresf,
                                num, logger)
            if num == -1:
                return False
    # GLOBAL INFORMATION
    nbcont = len(contig_sizes)
    gsize = sum(contig_sizes.values())
    if nbcont == 0 or gsize == 0:
        logger.warning(f"Your file {gpath} does not contain any gene. Please check that you "
                       "really gave a fasta sequence file")
        if grespath:
            gresf.close()
            os.remove(grespath)
        return False
    l90 = calc_l90(contig_sizes)
    # Everything ok for this genome -> complete its list of information in genomes dict
    # [spegenus.date] -> [spegenus.date, gpath, gpath_to_annotate, gsize, nbcont, L90]}
    if grespath:
        genomes[genome] += [gpath, grespath, gsize, nbcont, l90]
    else:
        genomes[genome] += [gpath, gpath, gsize, nbcont, l90]
    # If we wrote a new sequence file, close it
    if grespath:
        gresf.close()
    return True


def get_output_dir(soft, dbpath, tmp_path, genome, cut, pat):
    """
    Get output file to put sequence cut and/or sequence with shorter contigs (prokka)

    Parameters
    ----------
    soft : str
        soft used (prokka, prodigal, or None if called by prepare module)
    dbpath : str
        path to the folder containing the given genome sequence
    tmp_path : str
        path to folder where output files must be saved.
    genome : str
        genome name
    cut : bool
        True if contigs must be cut, False otherwise
    pat : str
        pattern on which contigs must be cut. ex: "NNNNN"

    Return
    ------
    grespath : str
        path to ouput file. None if no need to create new sequence file
    """
    # Path to sequence to analyze
    gpath = os.path.join(dbpath, genome)
    # genome file is in dbpath except if it was in several files in dbpath,
    # in which case it has been concatenated to a file in tmp_path
    if not os.path.isfile(gpath):
        gpath = os.path.join(tmp_path, genome)
    # New file create if needed. If not (prodigal and not cut), empty filename
    grespath = None
    # If user asks to cut at each 'pat', need to create a new sequence file,
    # whatever the annotation soft used
    if cut:
        new_file = genome + "_{}-split{}N.fna".format(soft, len(pat) - 1)
        grespath = os.path.join(tmp_path, new_file)
    # If no cutl, just keep original sequence, no need to create new file.
    # Just check that contigs have different names
    return gpath, grespath


def format_contig(cut, pat, cur_seq, cur_contig_name, genome, contig_sizes, gresf, num, logger):
    """
    Format given contig, and save to output file if needed

    - if cut: cut it and write each subsequence
    - write new contig just check that contig names are different

    Parameters
    ----------
    cut : bool
        True if contigs must be cut, False otherwise
    pat : str
        pattern on which contigs must be cut. ex: "NNNNN"
    cur_seq : str
        current sequence (aa)
    cur_contig_name : str
        name of current contig
    genome : str
        name of current genome
    cont_sizes : dict
        {contig_name : sequence length}
    gresf : io.TextIOWrappe
        open file to write new sequence. If we are annotating with prodigal and not cutting,
        there is no new sequence -> gref is None
    num : int
        current contig number
    logger : logging.Logger
        logger object to write log information

    Returns
    -------
    bool
        True if contig has been written without problem, False if any problem
    """
    # "CUT" if cut: need to cut at each 'pat' -> write new header + seq in new file
    if cut:
        # Cut sequence and write header + sequence to res file
        num = split_contig(pat, cur_seq, cur_contig_name, contig_sizes, gresf, num)
    # No cut -> no new file created, but check contig unique names
    else:
        if cur_contig_name in contig_sizes.keys():
            logger.error(f"In genome {genome}, '{cur_contig_name}' contig name is used for "
                         "several contigs. Please put different names for each contig. "
                         "This genome will be ignored.")
            return -1
        else:
            contig_sizes[cur_contig_name] = len(cur_seq)
    return num


def split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num):
    """
    Save the contig read just before into dicts and write it to sequence file.
    Unique ID of contig must be in the first field of header, before the first space
    (required by prokka)

    Parameters
    ----------
    pat : str
        pattern to split a contig. None if we do not want to split
    whole_seq : str
        sequence of current contig, to save once split according to pat
    cur_contig_name : str
        name of current contig to save once split according to pat
    contig_sizes : dict
        {name: size} save cur_contig once split according to pat
    gresf : _io.TextIOWrapper
        file open in w mode to save the split sequence
    num : int
        current contig number.

    Returns
    -------
    int
        new contig number, after giving number(s) to the current contig
    """
    # split contig each time a stretch of at least nbn 'N' is found (pattern pat)
    if not pat:
        cont_parts = [whole_seq]
    else:
        cont_parts = re.split(pat, whole_seq)

    # save contig parts
    for seq in cont_parts:
        # Only save non empty contigs (with some patterns, it could arrive that
        # we get empty contigs, if 2 occurrences of the pattern are side by side).
        if len(seq) == 0:
            continue
        new_contig_name = ">{}_{}\n".format(num, cur_contig_name.split(">")[1])
        contig_sizes[new_contig_name] = len(seq)
        gresf.write(new_contig_name)
        gresf.write(seq + "\n")
        num += 1
    return num


def calc_l90(contig_sizes):
    """
    Calc L90 of a given genome

    Parameters
    ----------
    contig_sizes : dict
        {name: size}

    Returns
    -------
    None or int
        if L90 found, returns L90. Otherwise, returns nothing
    """
    gsize = sum(contig_sizes.values())
    sizes = [contig_sizes[cont] for cont in contig_sizes]
    cum_sizes = np.cumsum(sorted(sizes, reverse=True))
    lim = 0.9 * gsize
    for num, val in enumerate(cum_sizes):
        if val >= lim:
            return num + 1


def rename_all_genomes(genomes):
    """
    FUNCTION DIRECTLY CALLED FROM MAIN ANNOTATE MODULE (step 3)
    Sort kept genomes by L90 and then nb contigs.
    For each genome, assign a strain number, and rename all its contigs.

    Parameters
    ----------
    genomes : dict
        {genome: [name, path, path_to_seq, gsize, nbcont, L90]} as input, and will become\
        {genome: [gembase_name, path, path_to_seq, gsize, nbcont, L90]} at the end

    Return
    ------
        change 1st field of genomes dict. name -> gembase_name (with strain number)

    """
    logger.info(f"Renaming kept genomes according to their quality ({len(genomes)} genomes)")
    # Keep first genome name to give to prodigal for training
    first_gname = None
    # Keep previous genome name (ESCO.0109 -> ESCO)
    last_name = ""
    # Keep last strain number
    last_strain = 0
    # "SAEN.1015.{}".format(str(last_strain).zfill(5))
    # Sort genomes by species, L90 and nb_contigs
    for genome, [name, _, _, _, _, _] in sorted(genomes.items(),
                                                key=utils.sort_genomes_byname_l90_nbcont):
        if not first_gname:
            first_gname = genome
        # first genome, or new strain name (ex: ESCO vs EXPL)
        # -> keep this new name, and add 1 to next strain number
        if last_name != name.split(".")[0]:
            last_strain = 1
            last_name = name.split(".")[0]
        # same strain name
        # -> write this new sequence, and go to next one (strain += 1)
        else:
            last_strain += 1
        # Write information to "genomes" dict.
        gembase_name = ".".join([name, str(last_strain).zfill(5)])
        genomes[genome][0] = gembase_name
    return first_gname


def plot_distributions(genomes, res_path, listfile_base, l90, nbconts):
    """
    FUNCTION DIRECTLY CALLED FROM MAIN ANNOTATE MODULE (step2)
    Plot distributions of L90 and nbcontig values.

    Parameters
    ----------
    genomes : dict
        {genome: [name, orig_path, to_annotate_path, size, nbcont, l90]}
    res_path : str
        path to put all output files
    listfile_base : str
        name of list file
    l90 : int
        L90 threshold
    nbconts : int
        nb contigs threshold

    Returns
    -------
    (l90_vals, nbcont_vals, dist1, dist2) :

    - l90_vals : list of l90 values for all genomes
    - nbcont_vals : list of nbcontigs for all genomes
    - dist1 : matplotlib figure of distribution of L90 values
    - dist2 : matplotlib figure of distribution of nb contigs values

    """
    logger.info("Generating distribution of L90 and #contigs graphs.")
    l90_vals = [val for _, (_, _, _, _, _, val) in genomes.items()]
    outl90 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
    nbcont_vals = [val for _, (_, _, _, _, val, _) in genomes.items()]
    outnbcont = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
    dist1 = utils.plot_distr(l90_vals, l90, "L90 distribution for all genomes",
                             "max L90 =", logger)
    dist2 = utils.plot_distr(nbcont_vals, nbconts,
                             "Distribution of number of contigs among all genomes",
                             "max #contigs =", logger)
    dist1.savefig(outl90)
    dist2.savefig(outnbcont)
    return l90_vals, nbcont_vals, dist1, dist2
