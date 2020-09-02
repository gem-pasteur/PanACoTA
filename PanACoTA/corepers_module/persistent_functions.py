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
Functions to generate a persistent genome from a pangenome.

@author gem
April 2017
"""
import logging
import math

from PanACoTA import utils

logger = logging.getLogger("corepers.pers")


def get_pers(fam_by_strain, fam_all_members, nb_strains, tol=1, multi=False, mixed=False,
             floor=False):
    """
    From the list of families, get the Pers Genome families, that are families having at least
    tol% of 'nb_strain' members.

    Parameters
    ----------
    fam_by_strain : dict
        {fam_num: {genome1: [members], genome2: [members]}, fam_num2: {genome1: [members]}}
    fam_all_members : dict
        {fam_num: [all members]}
    nb_strains : int
        total number of strains/genomes in dataset
    tol : float
        min percentage of different genomes present in a family
        ex: if tol=50%, and there are 8 genomes. If a family contains 3 genomes, it is not
        persistent. If it contains 7 genomes, it can be persistent (depends on multi and
        mixed parameters)
    multi : bool
        True if multiple genes from the same genome/strain in a family are tolerated. -> a family is considered as multi-persistent if it has members from at least 'tol%' genomes
        False otherwise
    mixed : bool
        True if mixed families are allowed (mixed family = exactly 1 member per genome
        for at least tol% of the genomes, 0 or several members allowed for other (1-tol)% genomes)
    floor : bool
        Use a minimum number of genomes containing a gene to consider the family
        persistent equal to: floor(nb_strains*tol) genomes if True, ceil(nb_strains*tol)
        if False.

    Returns
    -------
    dict
        {fam_num: [list of members]} for persistent families
    """
    logger.info("Generating Persistent genome of a dataset "
                f"containing {nb_strains} genomes")
    pers = {}  # {fam_num: {strain1: [genes from strain1], strain2: [genes from strain2]}}
    fams = {}  # {fam_num: [list of members]}
    if floor:
        min_members = math.floor(tol * nb_strains)
    else:
        min_members = math.ceil(tol * nb_strains)
    for fam_num, family in fam_by_strain.items():
        # If enough strains and multi accepted, or multi not accepted but 1 member per
        # strain, add family to core
        if mixed:
            if mixed_family(family, min_members):
                pers[fam_num] = family
                fams[fam_num] = fam_all_members[fam_num]
        else:
            if len(family) >= min_members and (multi or uniq_members(family)):
                pers[fam_num] = family
                fams[fam_num] = fam_all_members[fam_num]
    # coregenome computed
    if tol == 1 and not multi and not mixed:
        logger.info(f"The core genome contains {len(pers)} families, each one having "
                    f"exactly {int(min_members)} members, from the {nb_strains} different genomes.")
    # multi persistent genome with multigenic families allowed
    elif multi:
        logger.info(f"The persistent genome contains {len(pers)} families with members present "
                    f"in at least {min_members} different genomes ({tol*100}% of the total number of "
                    "genomes).")
    # mixed persistent genome, tol% families with exactly 1 member from each genome,
    # multigenic families allowed for the '1-tol'% remaining families
    elif mixed:
        logger.info(f"The persistent genome contains {len(pers)} families, "
                    f"each one having exactly 1 member from at least {tol*100}% of the genomes ({min_members} "
                    f"genomes). In the remaining {round((1-tol)*100,3)}% genomes, there can be 0, 1 or "
                    "several members.")
    # Strict persistent genome. tol% families with exactly one member in each genome
    else:
        logger.info(f"The persistent genome contains {len(pers)} families, each one having "
                    f"exactly 1 member from at least {tol*100}% of the {nb_strains} "
                    f"different genomes (that is {min_members} genomes). The other genomes are absent from "
                    "the family.")
    return fams


def mixed_family(family, thres):
    """
    1 family = several genomes (genome=strain), each containing x members
    Returns True if at least 'thres' genomes of the family have exactly 1 member.

    Parameters
    ----------
    family : dict
        {strain1: [members in strain1]}
    thres : float
         minimum number of genomes which must have exactly 1 member

    Returns
    -------
    bool
    """
    nb_uniq = 0
    for genome, members in family.items():
        if nb_uniq < thres:
            if len(members) == 0:
                logger.warning(f"Problem, no members for {genome}!")
            if len(members) == 1:
                nb_uniq += 1
        else:
            return True
    return nb_uniq >= thres


def uniq_members(family, num=1):
    """
    Returns True if, in the family, each genome has no more than 'num' member(s),
    False otherwise (multigenic family)

    Parameters
    ----------
    family : dict
        {strain1: [members in strain1], strain2: [members in strain2]}
    num : int
        max number of members allowed in each genome to return True

    Returns
    -------
    bool
    """
    for genome, members in family.items():
        if len(members) == 0:
            logger.warning(f"Problem, no members for {genome}!")
        if len(members) > num:
            return False
    return True


def write_persistent(fams, outfile):
    """
    Write persistent families into output file

    Parameters
    ----------
    fams : dict
        {num_fam: [members]}
    outfile : str
        output file to write all families
    """
    with open(outfile, "w") as outf:
        # Order by family number
        for num_fam in sorted(fams, key=lambda x: int(x)):
            outf.write(str(num_fam))  # Write family num
            fam = fams[num_fam]  # Get list of members of the family
            for mem in sorted(fam, key=utils.sort_proteins):
                outf.write(" " + mem)
            outf.write("\n")
