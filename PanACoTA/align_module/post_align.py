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
Concatenate all alignment files of all families. Then,
group alignments by genome.

@author: GEM, Institut Pasteur
March 2017
"""

import os
import sys
import logging
import progressbar
import multiprocessing
from PanACoTA import utils

logger = logging.getLogger("align.post")


def post_alignment(fam_nums, all_genomes, prefix, outdir, dname, quiet):
    """
    After the alignment of all proteins by family:

    - concatenate all alignment files
    - group the alignment by genome

    Parameters
    ----------
    fam_nums : []
        list of family numbers
    all_genomes : []
        list of all genomes in dataset
    prefix : str
        path to ``aldir/<name of dataset>`` (used to get extraction, alignment and btr files easily)
    outdir : str
        path to output directory, containing Aldir and Listdir, and that will also contain Treedir
    dname : str
        name of dataset (used to name concat and grouped files, as well as tree folder)
    quiet : bool
        True if nothing must be sent to sdtout/stderr, False otherwise
    """
    all_alns, status = concat_alignments(fam_nums, prefix, quiet)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(treedir, exist_ok=True)
    outfile = os.path.join(treedir, dname + ".grp.aln")
    res = launch_group_by_genome(all_genomes, all_alns, status, outfile, dname, quiet)
    if not res:
        logger.error("An error occurred. We could not group sequences by genome.")
    return outfile


def concat_alignments(fam_nums, prefix, quiet):
    """
    Concatenate all family alignment files to a unique file

    Parameters
    ----------
    fam_nums : []
        list of family numbers
    prefix : str
        path to ``aldir/<name of dataset>`` (used to get extraction, alignment and btr files easily)
    quiet : bool
        True if nothing must be sent to sdtout/stderr, False otherwise

    Returns
    -------
    tuple
        (output, str) with:

        - output: path to file containing concatenation of all alignments
        - str: "OK" if concatenation file already exists, "Done" if just did concatenation
    """
    output = "{}-complete.cat.aln".format(prefix)
    if os.path.isfile(output):
        logger.info("Alignments already concatenated")
        logger.warning(("Alignments already concatenated in {}. Program will use "
                        "it for next steps. If you want to redo it, remove it before "
                        "running.").format(output))
        return output, "OK"
    logger.info("Concatenating all alignment files")
    list_files = ["{}-mafft-prt2nuc.{}.aln".format(prefix, num_fam) for num_fam in fam_nums]
    # Check that all files exist
    for f in list_files:
        if not os.path.isfile(f):
            logger.error("The alignment file {} does not exist. Please check the families you "
                         "want, and their corresponding alignment files".format(f))
            sys.exit(1)
    if quiet:
        utils.cat(list_files, output)
    else:
        utils.cat(list_files, output, title="Concatenation")
    return output, "Done"


def launch_group_by_genome(all_genomes, all_alns, status, outfile, dname, quiet):
    """
    Function calling group_by_genome in a pool, while giving information to user
    (time elapsed)

    Parameters
    ----------
    all_genomes : []
        list of all genomes in the dataset
    all_alns : str
        path to file containing all alignments concatenated
    status : str
        "OK" if concatenation file already existed before running, "Done" if just did concatenation
    outfile : str
        file containing all families align by genome
    dname : str
        name of dataset
    quiet : bool
        True if nothing must be sent to sdtout/stderr, False otherwise

    Returns
    -------
    bool
        - True if everything went well or was already done
        - False if error occurred in at least one step
    """
    # Status = Done means that we just did the concatenation. So, if grouped by genome
    # file already exists, remove it.
    if status == "Done":
        if os.path.isfile(outfile):
            utils.remove(outfile)
    # Status was not 'Done' (it was 'ok', concat file already existed). And by_genome file
    # also already exists. Warn user
    if os.path.isfile(outfile):
        logger.info("Alignments already grouped by genome")
        logger.warning((f"Alignments already grouped by genome in {outfile}. "
                        "Program will end. "))
        return True
    logger.info("Grouping alignments per genome")
    bar = None
    if not quiet:
        widgets = [progressbar.BouncingBar(marker=progressbar.RotatingMarker(markers="◐◓◑◒")),
                   "  -  ", progressbar.Timer()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=20, term_width=50)
    pool = multiprocessing.Pool(1)
    args = [all_genomes, all_alns, outfile]
    final = pool.map_async(group_by_genome, [args], chunksize=1)
    pool.close()
    if not quiet:
        while True:
            if final.ready():
                break
            bar.update()
        bar.finish()
    pool.join()
    return False not in final.get()


def group_by_genome(args):
    """
    From the alignment file 'all_alns' containing all proteins, group the alignments of
    proteins by their genome (listed in 'all_genomes'), and save the result
    in 'treedir'

    Parameters
    ----------
    args : tuple
        - all_genomes: list of all genomes
        - all_alns: path to file containing all alignments concatenated
        - outfile: path to file which will contain alignments grouped by genome

    Returns
    -------
    bool
        - True if everything went well
        - False if problem when trying to group by genomes
    """
    all_genomes, all_alns, outfile = args
    sequences = read_alignments(all_alns, all_genomes)
    if not sequences:
        return False
    write_groups(outfile, sequences)
    return True


def read_alignments(all_alns, all_genomes):
    """
    Read alignment file, and assign each sequence to a genome

    Parameters
    ----------
    all_alns : str
        path to file containing all alignments concatenated
    all_genomes : []
        list of all genomes

    Returns
    -------
    dict or None
        - {genome_name: [list of sequences for this genome]}
        - None if problem with a protein for which we don't find the genome
    """
    sequences = {}  # name: [ordered list of sequences]
    genome = None
    seq = ""
    with open(all_alns, 'r') as alnf:
        for line in alnf:
            if line.startswith(">"):
                # If new header, write previous protein name/sequence to 'sequences'
                if genome and seq:
                    sequences[genome].append(seq)
                    seq = ""
                # Get new genome header
                genome = get_genome(line, all_genomes)
                if not genome:
                    return None
                if genome not in sequences:
                    sequences[genome] = []
            else:
                seq += line.strip()
    if genome and seq:
        sequences[genome].append(seq)
    per_genome = [len(seq) for seq in sequences.values()]
    if len(set(per_genome)) != 1:
        logger.error("Problems occurred while grouping alignments by genome: all genomes "
                     "do not have the same number of sequences. Check that each protein "
                     "name contains the name of the genome from which it comes.")
        return None
    logger.details("{} sequences found per genome".format(per_genome[0]))
    return sequences


def write_groups(outfile, sequences):
    """
    Writing alignments per genome to output file.

    Parameters
    ----------
    outfile : str
        path to file that will contain alignments grouped by genome
    sequences : dict
        {genome_name: [list of sequences (DNA, prot...) for this genome]}
    """
    logger.details("Writing alignments per genome")
    with open(outfile, "w") as outf:
        for genome in sorted(sequences, key=utils.sort_genomes_by_name):
            # write header for genome
            outf.write(">" + genome + "\n")
            # Write all sequences
            outf.write("".join(sequences[genome]) + "\n")


def get_genome(header, all_genomes):
    """
    Find to which genome belongs 'header'

    Parameters
    ----------
    header : str
        header read in alignment file
    all_genomes : []
        list of all genomes

    Returns
    -------
    str or None
        name of genome from which the header is
        None if no genome found
    """
    header = header.split(">")[1].split()[0]
    for genome in all_genomes:
        if genome in header:
            return genome
    logger.error((f"Protein {header} does not correspond to any genome name "
                  f"given... {all_genomes}"))
    return None
