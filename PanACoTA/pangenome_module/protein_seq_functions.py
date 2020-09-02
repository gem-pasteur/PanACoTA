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
Functions to build a bank of all proteins to include in the pangenome

@author gem
April 2017
"""
from PanACoTA import utils
import logging
import os

logger = logging.getLogger('pangenome.bank')


def build_prt_bank(lstinfo, dbpath, name, spedir, quiet):
    """
    Build a file containing all proteins of all genomes contained in lstinfo.

    Parameters
    ----------
    lstinfo : str
        1 line per genome, only 1st column considered here, as the genome name
        without extension
    dbpath : str
        Proteins folder, containing all proteins for each genome. Each genome has
        its own protein file, called `<genome_name>.prt`.
    name : str
        dataset name, used to name the output databank: <outdir>/<name>.All.prt
    spedir : str or None
        By default, output file is saved in dbpath directory. If it must be saved somewhere
        else, it is specified here.
    quiet : bool
        True if nothing must be written in stdout/stderr, False otherwise

    Returns
    -------
    str
        name (with path) of the protein databank generated
    """
    if not spedir:
        outdir = dbpath
    else:
        os.makedirs(spedir, exist_ok=True)
        outdir = spedir
    outfile = os.path.join(outdir, name + ".All.prt")
    if os.path.isfile(outfile):
        logger.warning((f"Protein bank {outfile} already exists. "
                        "It will be used by mmseqs."))
        return outfile
    logger.info(f"Building bank with all proteins to {outfile}")
    genomes = []
    with open(lstinfo) as lstf:
        for line in lstf:
            # skip header
            if "_name" in line:
                continue
            genome = line.strip().split()[0].strip()
            genomes.append(genome)
    all_names = [os.path.join(dbpath, gen + ".prt") for gen in genomes]
    if quiet:
        utils.cat(all_names, outfile)
    else:
        utils.cat(all_names, outfile, title="Building bank")
    return outfile
