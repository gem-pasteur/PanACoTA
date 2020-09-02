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
Functions to generate the matrix pan_quali, pan_quanti, as well
as a summary file for the pangenome.

@author gem
April 2017
"""
import logging
import numpy as np

from PanACoTA import utils
from PanACoTA import utils_pangenome as utilsp

logger = logging.getLogger("pangenome.post-treat")


def post_treat(families, pangenome):
    """
    From clusters = {num: [members]}, create:

    - a pan_quali matrix (lines = families, columns = genomes, 1 if genome present in\
     family, 0 otherwise)
    - a pan_quanti matrix (lines = families, columns = genomes, number of members from given\
     genome in the given family)
    - a summary file: lines = families. For each family:

        - nb_members: total number of members
        - sum_quanti: should be the same as nb_members!
        - sum_quali: number of different genomes in family
        - nb_0: number of missing genomes in family
        - nb_mono: number of genomes with exactly 1 member
        - nb_multi: number of genomes with more than 1 member
        - sum_0-mono-multi: should be equal to the total number of genomes in dataset
        - max_multi: maximum number of members from 1 genome

    Parameters
    ----------
    families : dict
        {num_fam: [list of members]}. Can be None, and then they will be retrieved from the\
        pangenome file
    pangenome : str
        file containing pangenome
    """
    fams_by_strain, families, all_strains = utilsp.read_pangenome(pangenome, logger, families)
    open_outputs_to_write(fams_by_strain, families, all_strains, pangenome)
    # result of open_outputs_to_write = (qualis, quantis, summaries)


def open_outputs_to_write(fams_by_strain, families, all_strains, pangenome):
    """
    Open output files, and call function to generate the matrix and summary file,
    and write it in those output files

    Parameters
    ----------
    fams_by_strain : dict
        {fam_num: {strain: [members]}}
    families : dict
        {fam_num: [all members]}
    all_strains : list
        list of all genome names
    pangenome : str
        filename containing pangenome. Will be extended for the 3 output files

    Returns
    -------
    (qualis, quantis, summaries) : tuple

        with:

        - qualis = {fam_num: [0 if no gene for species, 1 if at least 1 gene, for each\
          species in all_strains]}
        - quantis = {fam_num: [number of genes for each strain in all_strains]}
        - summaries = {fam_num: [nb_members, sum_quanti, sum_quali,\
          nb_0, nb_mono, nb_multi, sum_0-mono-multi, max_multi]}

    """
    panquali = pangenome + ".quali.txt"
    panquanti = pangenome + ".quanti.txt"
    pansum = pangenome + ".summary.txt"
    with open(pansum, "w") as psf:
        psf.write("num_fam,nb_members,sum_quanti,sum_quali,"
                  "nb_0,nb_mono,nb_multi,sum_0_mono_multi,max_multi\n")
        res = generate_and_write_outputs(fams_by_strain, families,
                                         all_strains, panquali, panquanti, psf)
    return res


def generate_and_write_outputs(fams_by_strain, families, all_strains, panquali, panquanti, psf):
    """
    From the python objects of pangenome, generate qualitative and quantitative matrix,
    as well as summary file.

    Parameters
    ----------
    fams_by_strain : dict
        {fam_num: {strain: [members]}}
    families : dict
        {fam_num: [all members]}
    all_strains : list
        list of all strains
    pqlf : _io.TextIOWrapper
        open file where qualitative matrix will be written
    pqtf : _io.TextIOWrapper
        open file where quantitative matrix will be written
    psf : _io.TextIOWrapper
        open file where summary will be written

    Returns
    -------
    (qualis, quantis, summaries) : tuple

        with:

        - qualis = {fam_num: [0 if no gene for species, 1 if at least 1 gene, for each\
        species in all_strains]}
        - quantis = {fam_num: [number of genes for each strain in all_strains]}
        - summaries = {fam_num: [nb_members, sum_quanti, sum_quali,\
         nb_0, nb_mono, nb_multi, sum_0-mono-multi, max_multi]}

    """
    logger.info("Generating qualitative and quantitative matrix, and summary file")

    # Matrix has:
    # - 1 row per family (header will be added after)
    # - 1 column for fam nums + 1 column per strain
    matrix_quali = np.empty((len(families), len(all_strains) + 1), dtype=int)
    matrix_quanti = np.empty((len(families), len(all_strains) + 1), dtype=int)

    # also save matrix as python objects
    qualis = {}
    quantis = {}
    summaries = {}
    row = 0
    for fam_num in sorted(fams_by_strain, key=lambda x: int(x)):
        strains = fams_by_strain[fam_num]
        quali = [1 if strain in strains else 0 for strain in all_strains]
        quanti = [len(strains[strain]) if strain in strains else 0 for strain in all_strains]
        nb_0 = quanti.count(0)
        nb_mono = quanti.count(1)
        nb_multi = len(quanti) - nb_0 - nb_mono
        max_multi = max(quanti)
        # Add line to quali and quanti matrices
        matrix_quali[row,:] = [fam_num] + quali
        matrix_quanti[row,:] = [fam_num] + quanti
        # Write summary line
        summ = [len(families[fam_num]), sum(quanti), sum(quali),
                nb_0, nb_mono, nb_multi, nb_0 + nb_mono + nb_multi, max_multi]
        psf.write(f"{fam_num},{utils.list_to_str(summ, sep=',')}")
        # Complete python objects with quali, quanti, sumary
        qualis[fam_num] = quali
        quantis[fam_num] = quanti
        summaries[fam_num] = summ
        row += 1
    # Add headers to quali and quanti matrix
    header = np.array(["fam_num"] + all_strains)
    matrix_quali = np.vstack((header, matrix_quali))
    matrix_quanti = np.vstack((header, matrix_quanti))
    # Transpose matrix: lines = genomes, columns = families
    tmatrix_quali = matrix_quali.transpose()
    np.savetxt(panquali, tmatrix_quali, delimiter=",", fmt="%s")
    tmatrix_quanti = matrix_quanti.transpose()
    np.savetxt(panquanti, tmatrix_quanti, delimiter=",", fmt="%s")
    return qualis, quantis, summaries
