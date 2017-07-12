#!/usr/bin/env python3
# coding: utf-8

"""
From the Persistent Genome file, group all persistent proteins per genome, in order to be
able to extract them faster after.

Input:
- Persistent genome (1 line per family, 1st column is family number, others are all members)
- list of all genomes (1 genome per line, only first column is considered)
- output directory
- dname: the name of the dataset, used to name output files

Output:
- for each genome: <outdir>/List-<dname>/<dname>-getEntry_gen_<genome>.txt: list of all
genes from the 'genome' which correspond to persistent proteins. 2 columns:
the first one is the protein name, the second is the filename to which it must be extracted
(corresponding to the family in the persistent genome).
- for each genome: <outdir>/List-<dname>/<dname>-getEntry_prt_<genome>.txt: same as
previous file but with the list of proteins to extract instead of genes.
- for each family: <outdir>/Align-<dname>/<dname>-current.<fam_num>.miss.lst: list of
 genomes which do not have members in family 'fam_num'.

@author GEM
November 2016
"""

import os
import sys
import datetime
import logging

logger = logging.getLogger("align.pan_to_pergenome")


def get_per_genome(persgen, list_gen, dname, outdir):
    """ From persistent genome and list of all genomes, sort persistent proteins by genome

    For each genome, write all persistent proteins to a file, with the family from which they
    are, in order to extract those proteins after.
    For each family, also save the names of genomes which do not have any member. This will
    be used to complete the alignments by stretches of '-'.
    persgen: file containing persistent genome
    list_gen: file containing the list of all genomes
    dname: name of the dataset
    outdir: Directory where files must be saved. Will create 2 subfolders: Align-<dname>
    and List-<dname>
    """
    logger.info("Reading PersGenome and constructing lists of missing genomes in each family.")
    # Define output directories
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    os.makedirs(aldir, exist_ok=True)
    os.makedirs(listdir, exist_ok=True)

    # Get list of all genomes
    all_genomes = get_all_genomes(list_gen)
    # Sort proteins by strain
    all_prots, fam_genomes, several = proteins_per_strain(persgen)
    # Write output files
    write_getentry_files(all_prots, several, listdir, aldir, dname, all_genomes)
    write_missing_genomes(fam_genomes, several, all_genomes, aldir, dname)
    return all_genomes, aldir, listdir, fam_genomes.keys()


def get_all_genomes(list_gen):
    """ Read the file containing the genomes list and get all genome names

    :param list_gen: file containing the list of genome names in 1st column
    :type list_gen: str
    """
    all_genomes = []
    with open(list_gen, "r") as lgf:
        for line in lgf:
            if "gembase" not in line:
                genome = line.strip().split()[0]
                all_genomes.append(genome)
    return all_genomes


def proteins_per_strain(persgen):
    """ From the persistentGenome file, get all persistent proteins, and classify them
    according to the strain from which they are.

    :param persgen: file containing persistent genome.
    :type persgen: str
    """
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
                strain = ".".join(mem.split(".")[:3])
                if strain not in fam_genomes[fam_num]:
                    fam_genomes[fam_num].append(strain)
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
    """ For each species, write all its persistent proteins into a file, with the
    persistent family in which they are. Those files will be used to extract all
    proteins.

    all_prots: {strain: {member: fam_num}}
    several: {fam_num: [genomes having several members in family]}
    listdir: directory where lists of proteins per genome must be saved
    aldir: directory where extracted proteins per family must be saved
    dname: name of dataset.
    all_genomes: list of all genomes
    """
    for strain, member in all_prots.items():
        gegenfile = os.path.join(listdir, dname + "-getEntry_gen_" + strain + ".txt")
        geprtfile = os.path.join(listdir, dname + "-getEntry_prt_" + strain + ".txt")
        with open(gegenfile, "w") as gegf, open(geprtfile, "w") as gepf:
            for mem, fam in member.items():
                if strain not in several[fam]:
                    genfile = os.path.join(aldir, dname + "-current." + fam + ".gen")
                    gegf.write(mem + " " + genfile + "\n")
                    prtfile = os.path.join(aldir, dname + "-current." + fam + ".prt")
                    gepf.write(mem + " " + prtfile + "\n")
    error = []
    for strain in all_genomes:
        gegenfile = os.path.join(listdir, dname + "-getEntry_gen_" + strain + ".txt")
        geprtfile = os.path.join(listdir, dname + "-getEntry_prt_" + strain + ".txt")
        if not os.path.isfile(gegenfile) or not os.path.isfile(geprtfile):
            error.append(strain)
    if error != []:
        for gen in error:
            logger.error(("There is not any protein for genome {} in any family! "
                          "The program will close, please fix this problem to be able to "
                          "run the alignments").format(gen))
        sys.exit(1)


def write_missing_genomes(fam_genomes, several, all_genomes, aldir, dname):
    """ For each family, write the names of all genomes which do not have any member in
    the family.

    fam_genomes: {fam_num: [genomes having a member in fam]}
    several: {fam_num: [genomes having several members in family]}
    all_genomes: list of all genomes
    aldir: directory where the lists of missing genomes per family must be saved
    dname: name of dataset
    """
    for fam, genomes in fam_genomes.items():
        missfile = os.path.join(aldir, dname + "-current." + fam + ".miss.lst")
        with open(missfile, "w") as mff:
            missing = (set(all_genomes) - set(genomes)).union(set(several[fam]))
            if missing:
                mff.write("\n".join(missing) + "\n")
