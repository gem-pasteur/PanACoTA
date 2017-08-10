#!/usr/bin/env python3
# coding: utf-8

"""
Functions to generate a persistent genome from a pangenome.

@author gem
April 2017
"""
import logging
import math

from genomeAPCAT import utils

logger = logging.getLogger("corepers.pers")


def get_pers(fam_by_strain, fam_all_members, nb_strains, tol=1, multi=False, mixed=False):
    """
    From the list of families, get the Pers Genome families, that are families having at least
    tol% of 'nb_strain' members.

    Parameters:
    fam_by_strain: dict, {fam_num: {strain: [members], strain: [members]}}
    fam_all_members: dict, {fam_num: [all members]}
    nb_strains: int, number of strains in dataset
    tol: float, percentage of isolates a gene must be in to be persistent
    multi: bool, True if multiple genes from the same isolate are tolerated, False otherwise
    mixed: bool, True if mixed families are allowed (exactly 1 member per family
    for at least tol% of the genomes, 0 or several members allowed for other (1-tol)%)
    """
    logger.info(("Generating Persistent genome of a dataset "
                 "containing {} genomes").format(nb_strains))
    pers = {} # {fam_num: {strain1: [genes from strain1], strain2: [genes from strain2]}}
    fams = {} # {fam_num: [list of members]}
    min_members = math.ceil(tol * nb_strains)
    for fam_num, family in fam_by_strain.items():
        # If enough strains and multi accepted, or multi not accepted but 1 member per strain, add family to core
        if mixed:
            if mixed_family(family, min_members):
                pers[fam_num] = family
                fams[fam_num] = fam_all_members[fam_num]
        else:
            if len(family) >= min_members and (multi or uniq_members(family)):
                pers[fam_num] = family
                fams[fam_num] = fam_all_members[fam_num]
    logger.info("The persistent genome contains {} families.".format(len(pers)))
    return fams


def mixed_family(family, thres):
    """Returns True if at least 'thres' genomes of the family have exactly 1 member.

    :param family: dict, {strain1: [members in strain1]}
    :type family: dict
    :param thres: minimum number of genomes which must have 1 member
    :type thres: number
    """
    nb_uniq = 0
    for ref, members in family.items():
        if nb_uniq < thres:
            if len(members) == 0:
                logger.warning("problem, no members for {}!".format(ref))
            if len(members) == 1:
                nb_uniq += 1
        else:
            return True
    return nb_uniq >= thres


def uniq_members(family, num=1):
    """
    Returns True if, in the family, each genome has no more than 'num' member(s),
    False otherwise (multigenic family)

    Parameters:
    family: dict, {strain1: [members in strain1], strain2: [members in strain2]}
    num: int, max number of members allowed in each genome to return True
    """
    for ref, members in family.items():
        if len(members) == 0:
            logger.warning("problem, no members for {}!".format(ref))
        if len(members) > num:
            return False
    return True


def write_persistent(fams, outfile):
    """ Write persistent families into output file

    :param fams: {num_fam: [members]}
    :type fams: dict
    :param outfile: output file to write all families
    :type outfile: str
    """
    with open(outfile, "w") as outf:
        for num_fam in sorted(fams, key=lambda x: int(x)):
            outf.write(str(num_fam))
            fam = fams[num_fam]
            for mem in sorted(fam, key=utils.sort_proteins):
                outf.write(" " + mem)
            outf.write("\n")

