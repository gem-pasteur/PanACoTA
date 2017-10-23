#!/usr/bin/env python3
# coding: utf-8

"""
Functions to generate the matrix pan_quali, pan_quanti, as well
as a summary file for the pangenome.

@author gem
April 2017
"""
import logging
import os
from genomeAPCAT import utils
from genomeAPCAT import utils_pangenome as utilsp

logger = logging.getLogger("pangenome.post-treat")


def post_treat(families, pangenome):
    """
    families: {num: [members]}
    pangenome: pangenome filename
    From clusters = {num: [members]}, create:
    - a pan_quali matrix (lines = families, columns = genomes, 1 if genome present in
    family, 0 otherwise)
    - a pan_quanti matrix (lines = families, columns = genomes, number of members from given
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
    """
    fams_by_strain, families, all_strains = utilsp.read_pangenome(pangenome, families)
    res = open_outputs_to_write(fams_by_strain, families, all_strains, pangenome)
    # res = (qualis, quantis, summaries)


def open_outputs_to_write(fams_by_strain, families, all_strains, pangenome):
    """ Open output files, and call function to generate the matrix and summary file,
    and write it in those output files

    :param fams_by_strain: {fam_num: {strain: [members]}}
    :type fams_by_strain: dict
    :param families: {fam_num: [all members]}
    :type families: dict
    :param all_strains: list of all genome names
    :type all_strains: list
    :param pangenome: filename containing pangenome. Will be extended for the 3 output files
    :type pangenome: str
    :returns: qualis = {fam_num: [0 if no gene for species, 1 if at least 1 gene, for each
    species in all_strains]}
    quantis = {fam_num: [number of genes for each strain in all_strains]}
    summaries = {fam_num: [nb_members, sum_quanti, sum_quali,
                           nb_0, nb_mono, nb_multi, sum_0-mono-multi, max_multi]}
    :rtype: tuple(dict, dict, dict)
    """
    panquali = pangenome + ".quali.txt"
    panquanti = pangenome + ".quanti.txt"
    pansum = pangenome + ".summary.txt"
    with open(panquali, "w") as pqlf, open(panquanti, "w") as pqtf, open(pansum, "w") as psf:
        pqlf.write("fam_num " + utils.write_list(all_strains))
        pqtf.write("fam_num " + utils.write_list(all_strains))
        psf.write("num_fam nb_members sum_quanti sum_quali "
                  "nb_0 nb_mono nb_multi sum_0_mono_multi max_multi\n")
        res = generate_and_write_outputs(fams_by_strain, families,
                                         all_strains, pqlf, pqtf, psf)
    return res


def generate_and_write_outputs(fams_by_strain, families, all_strains, pqlf, pqtf, psf):
    """ From the python objects of pangenome, generate qualitative and quantitative matrix,
    as well as summary file.

    :param fams_by_strain: {fam_num: {strain: [members]}}
    :type fams_by_strain: dict
    :param families: {fam_num: [all members]}
    :type families: dict
    :param all_strains: list of all strains
    :type all_strains: list
    :returns: qualis = {fam_num: [0 if no gene for species, 1 if at least 1 gene, for each
    species in all_strains]}
    quantis = {fam_num: [number of genes for each strain in all_strains]}
    summaries = {fam_num: [nb_members, sum_quanti, sum_quali,
                           nb_0, nb_mono, nb_multi, sum_0-mono-multi, max_multi]}
    :rtype: tuple(dict, dict, dict)
    """
    logger.info("Generating qualitative and quantitative matrix, and summary file")
    qualis = {}
    quantis = {}
    summaries = {}
    for fam_num in sorted(fams_by_strain, key=lambda x: int(x)):
        strains = fams_by_strain[fam_num]
        quali = [1 if strain in strains else 0 for strain in all_strains]
        quanti = [len(strains[strain]) if strain in strains else 0 for strain in all_strains]
        nb_0 = quanti.count(0)
        nb_mono = quanti.count(1)
        nb_multi = len(quanti) - nb_0 -nb_mono
        max_multi = max(quanti)
        summ = [len(families[fam_num]), sum(quanti), sum(quali),
                nb_0, nb_mono, nb_multi, nb_0 + nb_mono + nb_multi, max_multi]
        pqlf.write(str(fam_num) + " " + utils.write_list(quali))
        pqtf.write(str(fam_num) + " " + utils.write_list(quanti))
        psf.write(str(fam_num) + " " + utils.write_list(summ))
        qualis[fam_num] = quali
        quantis[fam_num] = quanti
        summaries[fam_num] = summ
    return qualis, quantis, summaries

