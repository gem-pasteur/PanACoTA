#!/usr/bin/env python3
# coding: utf-8

"""
Functions to build a bank of all proteins to include in the pangenome

@author gem
April 2017
"""
from genomeAPCAT import utils
import logging
import os

logger = logging.getLogger()


def build_prt_bank(lstinfo, dbpath, name, spedir):
    """
    Build a file containing all proteins of all genomes contained in lstinfo.

    lstinfo: 1 line per genome, only 1st column considered here, as the genome name
    without extension
    dbpath: Proteins folder, containing all proteins for each genome. Each genome has
    its own protein file, called `<genome_name>.prt`.
    name: dataset name, used to name the output databank: <outdir>/<name>.All.prt
    spedir: By default, output file is saved in dbpath directory. If it must be saved somewhere
    else, it is specified here.
    """
    if not spedir:
        outdir = dbpath
    else:
        os.makedirs(spedir, exist_ok=True)
        outdir = spedir
    outfile = os.path.join(outdir, name + ".All.prt")
    if os.path.isfile(outfile):
        logger.warning(("Protein bank {} already exists. "
                        "It will be used by mmseqs.").format(outfile))
        return outfile
    logger.info("Building bank with all proteins to {}".format(name + ".All.prt"))
    genomes = []
    with open(lstinfo, 'r') as lstf:
        for line in lstf:
            # skip header
            if "_name" in line:
                continue
            genome = line.strip().split()[0]
            genomes.append(genome)
    all_names = [os.path.join(dbpath, gen + ".prt") for gen in genomes]
    utils.cat(all_names, outfile, title="Building bank")
    return outfile
