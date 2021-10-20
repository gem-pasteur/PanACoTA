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
from PanACoTA import utils_pangenome as utilsp
import logging
import os
import functools as fnc
import shutil

logger = logging.getLogger('pangenome.bank')


def build_prt_bank(lstinfo, dbpath, name, spedir, quiet, dir=False):
    """
    Build a file or directory containing all proteins of all genomes contained in lstinfo.

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
    dir : bool
        If True, create a folder with required genomes instead of pulling them to one file. Some pangenome building
        methods require such format.

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

    genomes = utilsp.read_lstinfo(lstinfo, logger)
    all_names = [os.path.join(dbpath, gen + ".prt") for gen in genomes]

    if not dir:
        outfile = os.path.join(outdir, name + ".All.prt")
        if os.path.isfile(outfile):
                logger.warning((f"Protein bank {outfile} already exists. "
                            "It will be used by pangenome builder."))
                return outfile

        logger.info(f"Building bank with all proteins to {outfile}")

        if quiet:
            utils.cat(all_names, outfile)
        else:
            utils.cat(all_names, outfile, title="Building bank")
        return outfile
    else:
        outbase = os.path.join(outdir, name + "-All")
        destinations = [os.path.join(outbase, gen + ".prt") for gen in genomes]

        all_exist = fnc.reduce(lambda a, b: a and b, map(os.path.isfile, destinations))

        if all_exist:
            logger.warning((f"Protein bank {outbase} already exists. "
                            "It will be used by pangenome builder."))
            return outbase

        if os.path.isdir(outbase):
            shutil.rmtree(outbase)
            logger.warning((f"Protein bank {outbase} exists. "
                            "But it is invalid, so it will be removed and re-created."))

        os.mkdir(outbase)
        logger.info(f"Building bank with all proteins to {outbase} directory")
        for ini, dest in zip(all_names, destinations):
            shutil.copyfile(ini, dest)

        return outbase
