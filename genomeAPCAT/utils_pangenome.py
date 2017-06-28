#!/usr/bin/env python3
# coding: utf-8

"""
Functions used to deal with pangenome file

@author gem
April 2017
"""
import logging
import os
from genomeAPCAT import utils

logger = logging.getLogger("utils.pan")


def read_pangenome(pangenome, families=None):
    """ Read pangenome information

    Read pangenome according to what is available. First, check if python objects are available,
    then if not, search for the binary file, and if not, read the text file.
    pangenome: pangenome file
    families: {num: [members]}
    """
    if families:
        fams_by_strain, all_strains = get_fams_info(families)
        if not os.path.isfile(pangenome + ".bin"):
            logger.info("Saving all information to a binary file for later use")
            utils.save_bin([fams_by_strain, families, all_strains], pangenome + ".bin")
    elif os.path.isfile(pangenome + ".bin"):
        logger.info("Retrieving info from binary file")
        fams_by_strain, families, all_strains = utils.load_bin(pangenome + ".bin")
    else:
        fams_by_strain, families, all_strains = read_pan_file(pangenome)
        logger.info("Saving all information to a binary file for later use")
        utils.save_bin([fams_by_strain, families, all_strains], pangenome + ".bin")
    return fams_by_strain, families, all_strains


def get_fams_info(families):
    """
    input: families = {num: [members]}
    output:
    - fams_by_strain = {fam_num: {strain: [members], strain: [members]}}
    - all_strains = list of all strains found in the pangenome
    """
    logger.info("Calculating pan_summary, pan_quali and pan_quanti")
    fams_by_strain = {}
    all_strains = []
    for num, fam in families.items():
        fams_by_strain[num] = {}
        for gene in fam:
            read_gene(gene, num, fams_by_strain, all_strains)
    sort_all_strains = sorted(all_strains, key=utils.sort_genomes)
    return fams_by_strain, sort_all_strains


def read_pan_file(filein):
    """ Read PanGenome file in 'filein', and put it into Python objects

    Save all the python objects to a binary file, so that, if the script is ran several
    times, there is no need to parse again all the file.

    :param filein: name of the PanGenome file
    :type filein: str
    :returns: fams_by_strain = {fam_num: {strain: [members], strain: [members]}} and
    families = {fam_num: [all members]}
    all_strains = list of all strains found in the pangenome
    :rtype: tuple(dict, dict, list)
    """
    logger.info("Reading and getting information from pangenome file")
    fams_by_strain = {}
    families = {}
    all_strains = []
    nfam = 0
    with open(filein, 'r') as coref:
        for line in coref:
            genes = line.strip().split()
            fam_num = genes[0]
            fams_by_strain[fam_num] = {}
            genes_ok = genes[1:]
            for gene in genes_ok:
                read_gene(gene, fam_num, fams_by_strain, all_strains)
            families[fam_num] = genes_ok
            nfam += 1
    sort_all_strains = sorted(all_strains, key=utils.sort_genomes)
    return fams_by_strain, families, sort_all_strains


def read_gene(gene, num, fams_by_strain, all_strains):
    """
    Read information from a given gene name, and save it to appropriate dicts
    """
    strain = ".".join(gene.split(".")[:3])
    if strain in fams_by_strain[num]:
        fams_by_strain[num][strain].append(gene)
    else:
        fams_by_strain[num][strain] = [gene]
    if strain not in all_strains:
        all_strains.append(strain)

