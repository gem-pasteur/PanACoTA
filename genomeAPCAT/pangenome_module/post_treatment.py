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

logger = logging.getLogger("pangenome.post-treat")

def post_treat(families, pangenome):
    """
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
    fams_by_strain, families, all_strains = read_pangenome(pangenome, families)
    res = open_outputs_to_write(fams_by_strain, families, all_strains, pangenome)
    # res = (qualis, quantis, summaries)


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
    sort_all_strains = sorted(all_strains, key=lambda x: (x.split(".")[0], int(x.split(".")[-1])))
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
    logger.info("Calculating pan_summary, pan_quali and pan_quanti from pangenome file")
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
    sort_all_strains = sorted(all_strains, key=lambda x: int(x.split(".")[-1]))
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

