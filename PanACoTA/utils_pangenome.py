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
Functions used to deal with pangenome file

@author gem
April 2017
"""
import logging
import os
from PanACoTA import utils

logger = logging.getLogger("utils.pan")


def read_pangenome(pangenome, logger, families=None):
    """
    Read pangenome information

    Read pangenome according to what is available. First, check if python objects are available,
    then if not, search for the binary file, and if not, read the text file.

    Parameters
    ----------
    pangenome : str
        path to pangenome file
    families : dict or None
        {num: [members]} if families are given. If not (must read them from binary file if exists\
        or pangenome file otherwise), None.

    Returns
    -------
    (fams_by_strain, families, all_strains) : tuple
        with:

        - fams_by_strain: {fam_num: {strain: [members]}}
        - families: {fam_num: [all members]}
        - all_strains: list of all genome names
    """
    if families:
        fams_by_strain, all_strains = get_fams_info(families, logger)
        if not os.path.isfile(pangenome + ".bin"):
            logger.details("Saving all information to a binary file for later use")
            utils.save_bin([fams_by_strain, families, all_strains], pangenome + ".bin")
    elif os.path.isfile(pangenome + ".bin"):
        logger.info("Retrieving info from binary file")
        fams_by_strain, families, all_strains = utils.load_bin(pangenome + ".bin")
    else:
        fams_by_strain, families, all_strains = read_pan_file(pangenome, logger)
        logger.info("Saving all information to a binary file for later use")
        utils.save_bin([fams_by_strain, families, all_strains], pangenome + ".bin")
    return fams_by_strain, families, all_strains


def get_fams_info(families, logger):
    """
    From all families as list of members, get more information:

    - all strains found, sorted by species name
    - for each family, sort members by strain

    Parameters
    ----------
    families : dict
        {num: [members]}
    logger : logging.Logger
        logger object to write log information

    Returns
    -------
    (fams_by_strain, sorted_all_strains) : tuple
        with:

        - fams_by_strain: {fam_num: {strain: [members], strain: [members]}}
        - sorted_all_strains: list of all strains found, sorted by species
    """
    logger.info("Retrieving information from pan families")
    fams_by_strain = {}
    all_strains = set()
    for num, fam in families.items():
        fams_by_strain[num] = {}
        for gene in fam:
            read_gene(gene, num, fams_by_strain, all_strains)
    sort_all_strains = sorted(list(all_strains), key=utils.sort_genomes_by_name)
    return fams_by_strain, sort_all_strains


def read_pan_file(filein, logger):
    """
    Read PanGenome file in 'filein', and put it into Python objects

    Parameters
    ----------
    filein : str
        path to pangenome file
    logger :

    Returns
    -------
    (fams_by_strain, families, sort_all_strains) : tuple
        with:

        - fams_by_strain: {fam_num: {strain: [members]}}
        - families: {fam_num: [all members]}
        - sort_all_strains: list of all genome names, sorted by species name
    """
    logger.info("Reading and getting information from pangenome file")
    fams_by_strain = {}
    families = {}
    all_strains = set()
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
    sort_all_strains = sorted(list(all_strains), key=utils.sort_genomes_by_name)
    return fams_by_strain, families, sort_all_strains


def read_gene(gene, num, fams_by_strain, all_strains):
    """
    Read information from a given gene name, and save it to appropriate dicts

    Parameters
    ----------
    gene : str
        gene name (species.date.strain.contig_number
    num : str
        num of family from which the given gene is
    fams_by_strain : dict
        {fam_num: {strain: [members]}}
    all_strains : set
        set of all strains

    """
    # if format is ESCO.1512.00001.i001_12313 genome name is ESCO.1512.00001
    if "." in gene and len(gene.split(".")) >= 3:
        strain = ".".join(gene.split("_")[0].split(".")[:3])
    # otherwise, genename is everything before the last "_"
    else:
        strain = "_".join(gene.split("_")[:-1])
    if strain in fams_by_strain[num]:
        fams_by_strain[num][strain].append(gene)
    else:
        fams_by_strain[num][strain] = [gene]
    if strain not in all_strains:
        all_strains.add(strain)
