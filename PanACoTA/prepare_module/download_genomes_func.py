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


logger = logging.getLogger("ddg.log_dds")


def download_from_refseq(sum_file, NCBI_species, NCBI_taxid, outdir, threads):
    """
    Download refseq genomes of given species
    """
    # arguments needed to download a species genomes
    keyargs = {"section": "refseq", "file_format": "fasta", "output": outdir,
               "parallel": threads, "group": "bacteria",
               "species_taxid": NCBI_taxid}
    # summary file could not be downloaded because given species does not match
    #  any NCBI species. Just download genomes with the given taxID
    if not sum_file:
        logger.info("Downloading refseq genomes for taxid={}".format(NCBI_taxid))
    else:
        with open(sum_file, "r") as sum_lines:
            for line in sum_lines:
                infos = line.split()
                if len(infos)>=6:
                    try:
                        number = int(infos[6])
                    except ValueError:
                        continue
                if number != int(NCBI_taxid):
                    logger.error("Your NCBI_taxid ({}) does not match with your provided NCBI "
                                 "species ({}). The NCBI_taxid for this species is "
                                 "{}".format(NCBI_taxid, NCBI_species, infos[6]))
                    sys.exit(1)
        keyargs["genus"] = NCBI_species
        logger.info("Downloading refseq genomes for {} (taxid={})".format(NCBI_species,
                                                                          NCBI_taxid))
    max_retries = 15

    try:
        ret = ngd.download(**keyargs)
    except:
        logger.error("Could not download species taxID {}. Check that you gave the good "
                     "one.".format(NCBI_taxid))
        sys.exit(1)
    attempts = 0
    while ret == 75 and attempts < max_retries:
        attempts += 1
        logging.error(('Downloading from NCBI failed due to a connection error, '
                       'retrying. Already retried so far: %s'), attempts)
        ret = ngd.download(**keyargs)
    nb_gen, db_dir = to_database(outdir)
    logger.info("Downloaded {} genomes.".format(nb_gen))
    return db_dir


def download_summary(species_linked, outdir):
    """
    Get assembly_summary file for the given species if it exists. To be able to download it,
    the given NCBI species name must be exalctly as the name given on NCBI website.

    species_linked : given NCBI species with '_' instead of spaces, or NCBI taxID if species
    name not given (then, assembly file won't be found

    outdir: directory where downloaded assembly file must be saved
    """
    logger.info("Retrieving assembly_summary file for {}".format(species_linked))
    url = ("ftp://ftp.ncbi.nih.gov/genomes/refseq/"
           "bacteria/{}/assembly_summary.txt").format(species_linked)
    outfile = os.path.join(outdir, "assembly_summary-{}.txt".format(species_linked))
    try:
        urllib.request.urlretrieve(url, outfile)
    except:
        logger.warning("assembly_summary file cannot be downloaded. Please check that you "
                       "provided the exact species name, as given in NCBI")
        return
    return outfile


def to_database(outdir):
    """
    Move .fna.gz files to 'database_init' folder, and uncompress them.

    outdir : directory where all results are (for now, refseq folders, assembly summary and log

    Return:
        nb_gen : number of genomes downloaded
        db_dir : directory where are all fna files downloaded from refseq
    """
    # Unzip fasta files and put to a same folder
    logger.info("Uncompressing genome files.")
    db_dir = os.path.join(outdir, "Database_init")
    download_dir = os.path.join(outdir, "refseq", "bacteria")
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


def get_by_genome(all_list):
    """
    Given the list of replicons to download, sort them by genome, in order to directly
    concatenate all replicons from a same genome in the same file.
    """
    to_download = {}
    for file in all_list:
        filename = os.path.basename(file)
        strain = ".".join(filename.split(".")[:3])
        if strain in to_download:
            to_download[strain].append(file)
        else:
            to_download[strain] = [file]
    return to_download



