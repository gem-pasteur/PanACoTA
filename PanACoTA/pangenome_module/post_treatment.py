#!/usr/bin/env python3
# coding: utf-8

"""
Functions to generate the matrix pan_quali, pan_quanti, as well
as a summary file for the pangenome.

@author gem
April 2017
"""
import logging
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
    fams_by_strain, families, all_strains = utilsp.read_pangenome(pangenome, families)
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
    with open(panquali, "w") as pqlf, open(panquanti, "w") as pqtf, open(pansum, "w") as psf:
        pqlf.write("fam_num " + utils.list_to_str(all_strains))
        pqtf.write("fam_num " + utils.list_to_str(all_strains))
        psf.write("num_fam nb_members sum_quanti sum_quali "
                  "nb_0 nb_mono nb_multi sum_0_mono_multi max_multi\n")
        res = generate_and_write_outputs(fams_by_strain, families,
                                         all_strains, pqlf, pqtf, psf)
    return res


def generate_and_write_outputs(fams_by_strain, families, all_strains, pqlf, pqtf, psf):
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
    qualis = {}
    quantis = {}
    summaries = {}
    for fam_num in sorted(fams_by_strain, key=lambda x: int(x)):
        strains = fams_by_strain[fam_num]
        quali = [1 if strain in strains else 0 for strain in all_strains]
        quanti = [len(strains[strain]) if strain in strains else 0 for strain in all_strains]
        nb_0 = quanti.count(0)
        nb_mono = quanti.count(1)
        nb_multi = len(quanti) - nb_0 - nb_mono
        max_multi = max(quanti)
        summ = [len(families[fam_num]), sum(quanti), sum(quali),
                nb_0, nb_mono, nb_multi, nb_0 + nb_mono + nb_multi, max_multi]
        pqlf.write(str(fam_num) + " " + utils.list_to_str(quali))
        pqtf.write(str(fam_num) + " " + utils.list_to_str(quanti))
        psf.write(str(fam_num) + " " + utils.list_to_str(summ))
        qualis[fam_num] = quali
        quantis[fam_num] = quanti
        summaries[fam_num] = summ
    return qualis, quantis, summaries
