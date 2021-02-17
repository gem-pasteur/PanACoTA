#!/usr/bin/env python3

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
Functions helping for downloading refseq genomes of a species,
gunzip them, adding complete genomes...

@author gem
August 2017
"""
import os
import logging
import shutil
import sys
import glob
import urllib.request
import ncbi_genome_download as ngd

from PanACoTA import utils


logger = logging.getLogger("prepare.dds")


def download_from_ncbi(species_linked, section, ncbi_species_name, 
    ncbi_species_taxid, ncbi_taxid, levels, outdir, threads):
    """
    Download ncbi genomes of given species

    Parameters
    ----------
    species_linked : str
        given NCBI species with '_' instead of spaces, or NCBI taxID if species
        name not given
    section : str
        genbank or only refseq (default = refseq)
    ncbi_species_name : str or None
        name of species to download: user given NCBI species. None if
        no species name given
    ncbi_species_taxid : int
        species taxid given in NCBI (-T option)
    ncbi_taxid : int
        taxid given in NCBI (-t option)
    outdir : str
        Directory where downloaded sequences must be saved
    threads : int
        Number f threads to use to download genome sequences

    Returns
    -------
    str :
        Output filename of downloaded summary

    """
    # Name of summary file, with metadata for each strain:
    sumfile = os.path.join(outdir, f"assembly_summary-{species_linked}.txt")
    abs_sumfile = os.path.abspath(sumfile)

    # arguments needed to download all genomes of the given species
    abs_outdir = os.path.abspath(outdir)
    keyargs = {"section": section, "file_formats": "fasta", 
               "output": abs_outdir,
               "parallel": threads, "groups": "bacteria",
               "metadata_table":abs_sumfile}
    message = "Downloading all genomes for "
    # If NCBI species given, add it to arguments to download genomes, and write it to info message
    if ncbi_species_name:
        keyargs["genera"] = ncbi_species_name
        message += f"NCBI species = {ncbi_species_name}"
    # If NCBI species given, add it to arguments to download genomes, and write it to info message
    if ncbi_species_taxid:
        keyargs["species_taxids"] = ncbi_species_taxid
        if ncbi_species_name:
            message += f" (NCBI_species_taxid = {ncbi_species_taxid})."
        else:
            message += f" NCBI_species_taxid = {ncbi_species_taxid}"
    if ncbi_taxid:
        keyargs["taxids"] = ncbi_taxid
        if ncbi_species_name or ncbi_species_taxid:
            message += f" (and NCBI_taxid = {ncbi_taxid})."
        else:
            message += f" NCBI_taxid = {ncbi_taxid}"

    # If assembly level(s) given, add it to arguments, and write to info message
    if levels:
        keyargs["assembly_levels"] = levels
        message += f" (Only those assembly levels: {levels}). "
    logger.info(f"Metadata for all genomes will be saved in {sumfile}")
    logger.info(message)

    # Download genomes
    max_retries = 15 # If connection to NCBI fails, how many retry downloads must be done
    error_message = ("Could not download genomes. Check that you gave valid NCBI taxid and/or "
                     "NCBI species name. If you gave both, check that given taxID and name really "
                     "correspond to the same species.")
    # widgets = [progressbar.BouncingBar(marker=progressbar.RotatingMarker(markers="◐◓◑◒")),
    #            "  -  ", progressbar.Timer()]
    # bar = progressbar.ProgressBar(widgets=widgets, max_value=20, term_width=50)
    try:
        # Download genomes
        # ret = None
        # while True:
        #     if ret:
        #         break
        #     bar.update()
        ret = ngd.download(**keyargs)

    except: # pragma: no cover
        # Error message if crash during execution of ncbi_genome_download
        logger.error(error_message)
        # bar.finish()
        sys.exit(1)
    attempts = 0
    while ret == 75 and attempts < max_retries: # pragma: no cover
        # bar.update()
        attempts += 1
        logging.error(('Downloading from NCBI failed due to a connection error, '
                       'retrying. Already retried so far: %s'), attempts)
        ret = ngd.download(**keyargs)
    # bar.finish()
    # Message if NGD did not manage to download the genomes (wrong species name/taxid)
    if ret != 0:
        # Error message
        logger.error(error_message)
        sys.exit(1)
    nb_gen, db_dir = to_database(outdir, section)
    return db_dir, nb_gen


def to_database(outdir, section):
    """
    Move .fna.gz files to 'database_init' folder, and uncompress them.

    Parameters
    ----------
    outdir : str
        directory where all results are (for now, refseq/genbank folders, assembly summary and log
    section : str
        refseq (default) or genbank

    Returns
    -------
        nb_gen : number of genomes downloaded
        db_dir : directory where are all fna files downloaded from refseq/genbank
    """
    # Copy .gz files in a new folder, and Unzip them in this new folder
    logger.info("Uncompressing genome files.")
    # Folder where are .gz files
    download_dir = os.path.join(outdir, section, "bacteria")
    # If no folder output/refseq/bacteria: error, no genome found
    # (or output/genbank/bacteria)
    if not os.path.exists(download_dir):
        logger.error(f"The folder containing genomes downloaded from NCBI {section} "
                     f"({download_dir}) does not exist. Check that you really downloaded "
                     "sequences (fna.gz) and that they are in this folder.")
        sys.exit(1)
    # If folder output/<refseq or genbank>/bacteria empty: error, no genome found
    list_downloads = os.listdir(download_dir)
    if list_downloads == []:
        logger.error(f"The folder supposed to contain genomes downloaded from NCBI {section} "
                     f"({download_dir}) exists but is empty. Check that you really downloaded "
                     "sequences (fna.gz).")
        sys.exit(1)
    # Create directory to put uncompressed genomes
    db_dir = os.path.join(outdir, "Database_init")
    os.makedirs(db_dir, exist_ok=True)
    nb_gen = 0
    # For each subfolder of download dir, move the .gz file it contains (if possible)
    # to the new database folder
    for g_folder in os.listdir(download_dir):
        fasta = glob.glob(os.path.join(download_dir, g_folder, "*.fna.gz"))
        # No .gz file in folder
        if len(fasta) == 0:
            logger.warning("Problem with genome in {}: no compressed fasta file downloaded. "
                           "This genome will be ignored.".format(g_folder))
            continue
        # Several gz files in folder
        elif len(fasta) > 1:
            logger.warning("Problem with genome in {}: several compressed fasta files found. "
                           "This genome will be ignored.".format(g_folder))
            continue
        # Copy gz file to new folder
        fasta_file = os.path.basename(fasta[0])
        fasta_out = os.path.join(db_dir, fasta_file)
        shutil.copy(fasta[0], fasta_out)
        # Uncompress file copied
        cmd = f"gunzip {fasta_out} -f"
        error = f"Error while trying to uncompress {fasta_out}. This genome will be ignored."
        call = utils.run_cmd(cmd, error)
        # Problem with uncompressing: genome ignored (remove gz file from new folder)
        if call.returncode != 0:
            os.remove(fasta_out)
            continue
        nb_gen += 1
    return nb_gen, db_dir
