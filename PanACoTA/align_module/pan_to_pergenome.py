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
From the Persistent Genome file, group all persistent proteins per genome, in order to be
able to extract them faster after.

Input:

- Persistent genome (1 line per family, 1st column is family number, others are all members)
- list of all genomes (1 genome per line, only first column is considered)
- output directory
- dname: the name of the dataset, used to name output files

Output:

- for each genome: ``<outdir>/List-<dname>/<dname>-getEntry_gen_<genome>.txt``: list of all
  genes from the 'genome' which correspond to persistent proteins. 2 columns:
  the first one is the protein name, the second is the filename to which it must be extracted
  (corresponding to the family in the persistent genome).
- for each genome: ``<outdir>/List-<dname>/<dname>-getEntry_prt_<genome>.txt``: same as
  previous file but with the list of proteins to extract instead of genes.
- for each family: ``<outdir>/Align-<dname>/<dname>-current.<fam_num>.miss.lst``: list of
  genomes which do not have members in family 'fam_num'.

@author GEM
November 2016
"""

import os
import sys
import logging

logger = logging.getLogger("align.pan_to_pergenome")


def get_per_genome(persgen, list_gen, dname, outdir):
    """
    From persistent genome and list of all genomes, sort persistent proteins by genome

    For each genome, write all persistent proteins to a file, with the family from which they
    are, in order to extract those proteins after.
    For each family, also save the names of genomes which do not have any member. This will
    be used to complete the alignments by stretches of '-'.

    Parameters
    ----------
    persgen : str
        file containing persistent genome
    list_gen : str
        file containing the list of all genomes
    dname : str
        name of the dataset
    outdir : str
        Directory where files must be saved. Will create 2 subfolders: ``Align-<dname>``
        and ``List-<dname>``

    Returns
    -------
    (all_genomes, aldir, listdir, families) : tuple

        * all_genomes : [] list of all genome names
        * aldir : str, path to align directory
        * listdir : str, path to List directory
        * families : str, list of family numbers
    """
    # Define output directories
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    os.makedirs(aldir, exist_ok=True)
    os.makedirs(listdir, exist_ok=True)

    # Get list of all genomes
    all_genomes = get_all_genomes(list_gen)
    logger.info(f"Found {len(all_genomes)} genomes.")

    # Sort proteins by strain
    logger.info("Reading PersGenome and constructing lists of missing genomes in each family.")
    all_prots, fam_genomes, several = proteins_per_strain(persgen)
    # Write output files
    write_getentry_files(all_prots, several, listdir, aldir, dname, all_genomes)
    write_missing_genomes(fam_genomes, several, all_genomes, aldir, dname)
    return all_genomes, aldir, listdir, fam_genomes.keys()


def get_all_genomes(list_gen):
    """
    Read the file containing the genomes list and get all genome names

    Parameters
    ----------
    list_gen : str
        File containing the list of all genomes

    Returns
    -------
    list
        list of all genome names

    """
    all_genomes = []
    with open(list_gen, "r") as lgf:
        for line in lgf:
            if "gembase" not in line:
                genome = line.strip().split()[0]
                all_genomes.append(genome)
    return all_genomes


def proteins_per_strain(persgen):
    """
    From the persistentGenome file, get all persistent proteins, and classify them
    according to the strain from which they are.

    Parameters
    ----------
    persgen : str
        File containing persistent genome

    Returns
    -------
    (all_prots, fam_genomes, several) : tuple

        * all_prots: dict, {strain: {member: fam_num}}
        * fam_genomes: dict, {fam_num: [genomes having a member in fam]}
        * several: dict, {fam_num: [genomes having several members in fam]}
    """
    logger.info("Getting all persistent proteins and classify by strain.")
    all_prots = {}  # {strain: {member: fam_num}}
    fam_genomes = {}  # {fam_num: [genomes having a member in fam]}
    several = {}  # {fam_num: [genomes having several members in fam]}
    with open(persgen, "r") as pgf:
        for line in pgf:
            fam_num = line.split()[0]
            members = line.strip().split()[1:]
            fam_genomes[fam_num] = []
            several[fam_num] = []
            for mem in members:
                # if format is ESCO.1512.00001.i0002_12124, strain is 3 first fields
                # separated by '.'
                if "." in mem and len(mem.split(".")) >= 3:
                    strain = ".".join(mem.split(".")[:3])
                # if format is not like this, it must be something_00001:
                else:
                    strain = "_".join(mem.split("_")[:-1])
                # if strain not already in fam_genomes, add it
                if strain not in fam_genomes[fam_num]:
                    fam_genomes[fam_num].append(strain)
                # If strain already in fam_genomes, it has several members: add it to several
                elif strain not in several[fam_num]:
                    several[fam_num].append(strain)
                if strain not in all_prots:
                    all_prots[strain] = {}
                if mem in all_prots[strain]:
                    logger.warning((" problem: {} already exists, in family {}. Conflict with "
                                    "family {}.").format(mem, all_prots[strain][mem], fam_num))
                all_prots[strain][mem] = fam_num
    return all_prots, fam_genomes, several


def write_getentry_files(all_prots, several, listdir, aldir, dname, all_genomes):
    """
    For each species, write all its persistent proteins into a file, with the
    persistent family in which they are. Those files will be used to extract all
    proteins.

    Parameters
    ----------
    all_prots : dict
        {strain: {member: fam_num}}
    several : dict
        {fam_num: [genomes having several members in family]}
    listdir : str
        directory where lists of proteins per genome must be saved
    aldir : str
         directory where extracted proteins per family must be saved
    dname : str
        name of dataset
    all_genomes : list
        list of all genomes
    """
    for strain, member in all_prots.items():
        write_genome_file(listdir, aldir, dname, strain, member, several)
    error = []
    for strain in all_genomes:
        gegenfile = os.path.join(listdir, dname + "-getEntry_gen_" + strain + ".txt")
        geprtfile = os.path.join(listdir, dname + "-getEntry_prt_" + strain + ".txt")
        if not os.path.isfile(gegenfile) or not os.path.isfile(geprtfile):
            error.append(strain)
    if error:
        for gen in error:
            logger.error(("There is not any protein for genome {} in any family! "
                          "The program will close, please fix this problem to be able to "
                          "run the alignments").format(gen))
        sys.exit(1)


def write_genome_file(listdir, aldir, dname, strain, member, several):
    """
    For a given genome, write all the proteins and genes to extract to its file.
    If one of the 2 files (proteins and genes) already exists, overwrite it.
    If no file exists -> write them.
    If the 2 files exist -> warning saying that we use already existing files.

    Parameters
    ----------
    listdir : str
        path to List directory
    aldir : str
        path to Align directory
    dname : str
        name of dataset
    strain : str
        current genome name
    member : dict
        {member: fam_num}
    several : dict
        {fam_num: [genomes having several members in family]}
    """
    # If the 2 files exist, use them as they are
    gegenfile = os.path.join(listdir, dname + "-getEntry_gen_" + strain + ".txt")
    geprtfile = os.path.join(listdir, dname + "-getEntry_prt_" + strain + ".txt")
    if os.path.isfile(gegenfile) and os.path.isfile(geprtfile):
        logger.warning(f"For genome {strain}, {geprtfile} and {gegenfile} already exist. "
                       "The program will use them "
                       "to extract proteins and genes. If you prefer to rewrite them, use "
                       "option -F (or --force).")
        return

    # If at least one of the 2 files already exists, overwrite both files
    with open(gegenfile, "w") as gegf, open(geprtfile, "w") as gepf:
        for mem, fam in member.items():
            if strain not in several[fam]:
                genfile = os.path.join(aldir, dname + "-current." + fam + ".gen")
                gegf.write(mem + " " + genfile + "\n")
                prtfile = os.path.join(aldir, dname + "-current." + fam + ".prt")
                gepf.write(mem + " " + prtfile + "\n")


def write_missing_genomes(fam_genomes, several, all_genomes, aldir, dname):
    """
    For each family, write the names of all genomes which do not have any member in
    the family.

    Parameters
    ----------
    fam_genomes : dict
        {fam_num: [genomes having a member in fam]}
    several : dict
        {fam_num: [genomes having several members in family]}
    all_genomes : list
        list of all genomes
    aldir : str
        directory where the lists of missing genomes per family must be saved
    dname : str
        name of dataset
    """
    for fam, genomes in fam_genomes.items():
        # File where missing genomes will be written
        missfile = os.path.join(aldir, f"{dname}-current.{fam}.miss.lst")
        with open(missfile, "w") as mff:
            # missing = missing or several members:
            # miss: all genomes - genomes in the family
            # several: add to 'miss' genomes with several members in the family
            missing = (set(all_genomes) - set(genomes)).union(set(several[fam]))
            if missing:
                mff.write("\n".join(missing) + "\n")
