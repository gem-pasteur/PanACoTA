#!/usr/bin/env python3

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


def download_from_refseq(species_linked, NCBI_species, NCBI_taxid, outdir, threads):
    """
    Download refseq genomes of given species

    Parameters
    ----------
    species_linked : str
        given NCBI species with '_' instead of spaces, or NCBI taxID if species
        name not given
    NCBI_species : str
        name of species to download: user given NCBI species with '_' instead of spaces. None if
        no species name given
    NCBI_taxid : int
        species taxid given in NCBI
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
    sumfile = os.path.join(outdir, "assembly_summary-{}.txt".format(species_linked))
    abs_sumfile = os.path.abspath(sumfile)

    # arguments needed to download all genomes of the given species
    abs_outdir = os.path.abspath(outdir)
    keyargs = {"section": "refseq", "file_format": "fasta", "output": abs_outdir,
               "parallel": threads, "group": "bacteria",
               "species_taxid": NCBI_taxid, "metadata_table":abs_sumfile}
    message = "Downloading all genomes for "
    # If NCBI species given, add it to arguments to download genomes, and write it to info message
    if NCBI_species:
        keyargs["genus"] = NCBI_species
        message += f"NCBI species = {NCBI_species}"
    # If NCBI species given, add it to arguments to download genomes, and write it to info message
    if NCBI_taxid:
        keyargs["species_taxid"] = NCBI_taxid
        if NCBI_species:
            message += f" (NCBI_taxid = {NCBI_taxid})."
        else:
            message += f" NCBI_taxid = {NCBI_taxid}"
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

    except:
        # Error message if crash during execution of ncbi_genome_download
        logger.error(error_message)
        # bar.finish()
        sys.exit(1)
    attempts = 0
    while ret == 75 and attempts < max_retries:
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
    nb_gen, db_dir = to_database(outdir)
    logger.info("Downloaded {} genomes.".format(nb_gen))
    return db_dir


def to_database(outdir):
    """
    Move .fna.gz files to 'database_init' folder, and uncompress them.

    Parameters
    ----------
    outdir : str
        directory where all results are (for now, refseq folders, assembly summary and log

    Returns
    -------
        nb_gen : number of genomes downloaded
        db_dir : directory where are all fna files downloaded from refseq
    """
    # Unzip fasta files and put to a same folder
    logger.info("Uncompressing genome files.")
    download_dir = os.path.join(outdir, "refseq", "bacteria")
    # If no folder output/refseq/bacteria: error, no genome found
    if not os.path.exists(download_dir):
        logger.error(f"The folder containing genomes downloaded from NCBI refseq "
                     f"({download_dir}) does not exist. Check that you really downloaded "
                     "sequences (fna.gz) and that they are in this folder.")
        sys.exit(1)
    # If folder output/refseq/bacteria empty: error, no genome found
    list_downloads = os.listdir(download_dir)
    if list_downloads == []:
        logger.error(f"The folder supposed to contain genomes downloaded from NCBI refseq "
                     f"({download_dir}) exists but is empty. Check that you really downloaded "
                     "sequences (fna.gz).")
        sys.exit(1)
    # Create directory to put uncompressed genomes
    db_dir = os.path.join(outdir, "Database_init")
    os.makedirs(db_dir, exist_ok=True)
    nb_gen = 0
    for g_folder in os.listdir(download_dir):
        fasta = glob.glob(os.path.join(download_dir, g_folder, "*.fna.gz"))
        if len(fasta) == 0:
            logger.warning("Problem with genome {}: no fasta file downloaded.".format(g_folder))
            continue
        elif len(fasta) > 1:
            logger.warning("Problem with genome {}: several fasta files found.".format(g_folder))
            continue
        nb_gen += 1
        fasta_file = os.path.basename(fasta[0])
        fasta_out = os.path.join(db_dir, fasta_file)
        shutil.copy(fasta[0], fasta_out)
        cmd = "gunzip {} -f".format(fasta_out)
        error = "Error while trying to uncompress {}".format(fasta_out)
        utils.run_cmd(cmd, error)
    return nb_gen, db_dir
