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
Functions to use mmseqs to create a pangenome

@author gem
April 2017
"""
import logging
import os
import time
import multiprocessing
import progressbar
import copy

from PanACoTA import utils

logger = logging.getLogger("pangenome.mmseqs")


def run_all_pangenome(min_id, clust_mode, outdir, prt_path, threads, panfile=None, quiet=False):
    """
    Run all steps to build a pangenome:

    - create mmseqs database from protein bank
    - cluster proteins
    - convert to pangenome


    Parameters
    ----------
    min_id : float
        minimum percentage of identity to be in the same family
    clust_mode : [0, 1, 2]
        0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'
    outdir : str
        directory where output cluster file must be saved
    prt_path : str
        path to file containing all proteins to cluster.
    threads : int
        number of threads which can be used
    panfile : str or None
        name for output pangenome file. Otherwise, will use default name
    quiet : bool
        True if nothing must be written on stdout, False otherwise.

    Returns
    -------
    (families, outfile) : tuple

        - families : {fam_num: [all members]}
        - outfile : pangenome filename
    """
    # Get general information and file/directory names
    prt_bank = os.path.basename(prt_path)
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    mmseqdb = os.path.join(outdir, prt_bank + "-msDB")
    information = ("Will run MMseqs2 with:\n"
                   f"\t- minimum sequence identity = {min_id*100}%\n"
                   f"\t- cluster mode {clust_mode}")
    if threads > 1:
        information += "\n\t- {} threads".format(threads)
    logger.info(information)

    infoname = get_info(threads, min_id, clust_mode, start)
    logmmseq = get_logmmseq(outdir, prt_bank, infoname)
    # Create ffindex of DB if not already done
    create_mmseqs_db(mmseqdb, prt_path, logmmseq)

    # Cluster with mmseqs
    if panfile:
        panfile = os.path.join(outdir, panfile)
    families, outfile = do_pangenome(outdir, prt_bank, mmseqdb, min_id,
                                     clust_mode, threads, start, panfile, quiet)
    return families, outfile


def get_info(threads, min_id, clust_mode, start):
    """
    Get string containing all information on future run

    Parameters
    ----------
    threads : int
        max number of threads to use
    min_id : float
        min percentage of identity to consider 2 proteins in hte same family
    clust_mode : [0, 1, 2]
        0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'
    start : str
        string containing starting date and time

    Returns
    -------
    str
        string containing info on current run, to put in filenames
    """
    if threads != 1:
        threadinfo = "-th" + str(threads)
    else:
        threadinfo = ""
    infoname = str(min_id) + "-mode" + str(clust_mode) + threadinfo + "_" + start
    return infoname


def get_logmmseq(outdir, prt_bank, infoname):
    """
    Get filename of logfile, given information

    Parameters
    ----------
    outdir : str
        output directory
    prt_bank : str
        name of file (without path) containing all proteins to cluster
    infoname : str
        string containing information on the current run

    Returns
    -------
    str
        path to mmseq logfile
    """
    return os.path.join(outdir, "mmseq_" + prt_bank + "_" + infoname + ".log")


def do_pangenome(outdir, prt_bank, mmseqdb, min_id, clust_mode, threads, start, panfile=None,
                 quiet=False):
    """
    Use mmseqs to cluster proteins

    Parameters
    ----------
    outdir : str
        directory where output files are saved
    prt_bank : str
        name of the file containing all proteins to cluster, without path
    mmseqdb : str
        path to base filename for output mmseq db
    min_id : float
        min percentage of identity to be considered in the same family (between 0 and 1)
    clust_mode : [0, 1, 2]
        0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'
    threads : int
        max number of threads to use
    start : str
        string containing trat date and time
    panfile : str
        if a pangenome file is specified. Otherwise, default pangenome name will be used
    quiet : bool
        true if nothing must be print on stdout/stderr, false otherwise (show progress bar)

    Returns
    -------
    (families, outfile) : tuple

        - families : {fam_num: [all members]}
        - outfile : pangenome filename
    """
    infoname = get_info(threads, min_id, clust_mode, start)
    logmmseq = get_logmmseq(outdir, prt_bank, infoname)
    mmseqclust = os.path.join(outdir, prt_bank + "-clust-" + infoname)
    tmpdir = os.path.join(outdir, "tmp_" + prt_bank + "_" + infoname)
    os.makedirs(tmpdir, exist_ok=True)
    bar = None
    logger.debug(mmseqclust)
    if os.path.isfile(mmseqclust):
        logger.warning(("mmseqs clustering {} already exists. The program will now convert "
                        "it to a pangenome file.").format(mmseqclust))
    else:
        logger.info("Clustering proteins...")
        if not quiet:
            widgets = [progressbar.BouncingBar(marker=progressbar.RotatingMarker(markers="◐◓◑◒")),
                       "  -  ", progressbar.Timer()]
            bar = progressbar.ProgressBar(widgets=widgets, max_value=20, term_width=50)
        pool = multiprocessing.Pool(1)
        args = [mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode]
        final = pool.map_async(run_mmseqs_clust, [args], chunksize=1)
        pool.close()
        if not quiet:
            while True:
                if final.ready():
                    break
                bar.update()
            bar.finish()
        pool.join()
    # Convert output to tsv file (one line per comparison done)
    #  # Convert output to tsv file (one line per comparison done)
    # -> returns (families, outfile)
    return mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, start, panfile)


def run_mmseqs_clust(args):
    """
    Run mmseqs clustering

    Parameters
    ----------
    args : tuple
         (mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode), with:

            * mmseqdb: path to base filename for output mmseq db
            * mmseqclust: path to base filename for output of mmseq clustering
            * tmpdir : path to folder which will contain mmseq temporary files
            * logmmseq : path to file where logs must be written
            * min_id : min percentage of identity to be considered in the same family
            *         (between 0 and 1)
            * threads : max number of threads to use
            * clust_mode : [0, 1, 2], 0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'

    """
    mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode = args
    cmd = ("mmseqs cluster {} {} {} --min-seq-id {} --threads {} --cluster-mode "
           "{}").format(mmseqdb, mmseqclust, tmpdir, min_id, threads, clust_mode)
    logger.details(f"MMseqs command: {cmd}")
    msg = f"Problem while clustering proteins with mmseqs. See log in {logmmseq}"
    with open(logmmseq, "a") as logm:
        utils.run_cmd(cmd, msg, eof=False, stdout=logm, stderr=logm)


def mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, start, outfile=None):
    """
    Convert mmseqs clustering to a pangenome file:

    - convert mmseqs results to tsv file
    - convert tsv file to pangenome

    Parameters
    ----------
    mmseqdb : str
         path to base filename for output mmseq db
    mmseqclust : str
        path to base filename for output of mmseq clustering
    logmmseq : str
         path to file where logs must be written
    start : str
        string containing start date and time
    outfile : str
        pangenome filename, or None if default one must be used

    Returns
    -------
    (families, outfile) : tuple

        - families : {fam_num: [all members]}
        - outfile : pangenome filename
    """
    cmd = f"mmseqs createtsv {mmseqdb} {mmseqdb} {mmseqclust} {mmseqclust}.tsv"
    msg = "Problem while trying to convert mmseq result file to tsv file"
    logger.details(f"MMseqs command: {cmd}")
    with open(logmmseq, "a") as logf:
        utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
    # Convert the tsv file to a 'pangenome' file: one line per family
    families, outfile = mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start, outfile)
    return families, outfile


def mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start, outfile=None):
    """
    Convert the tsv output file of mmseqs to the pangenome file

    Parameters
    ----------
    mmseqclust : str
        path to base filename for output of mmseq clustering
    logmmseq : str
        path to file where logs must be written
    start : str
        string containing start date and time
    outfile : str
        pangenome filename, or None if default one must be used

    Returns
    -------
    (families, outfile) : tuple

        - families : {fam_num: [all members]}
        - outfile : pangenome filename
    """
    logger.info("Converting mmseqs results to pangenome file")
    tsvfile = mmseqclust + ".tsv"
    if not outfile:
        outpath = os.path.dirname(tsvfile)
        base = os.path.basename(tsvfile)
        outfile = os.path.join(outpath, "PanGenome-" + base + ".lst")
    clusters = mmseq_tsv_to_clusters(tsvfile)
    families = clusters_to_file(clusters, outfile)
    end = time.strftime('%Y-%m-%d_%H-%M-%S')
    with open(logmmseq, "a") as logm:
        logm.write(f"\n------------\n\nStart: {start} \n")
        logm.write(f"End: {end}")
    return families, outfile


def mmseq_tsv_to_clusters(mmseq):
    """
    Reads the output of mmseq as a tsv file, and converts it to a python dict

    Parameters
    ----------
    mmseq : str
        filename of mmseq clustering output in tsv format

    Returns
    -------
    dict
        {representative_of_cluster: [list of members]}

    """
    clusters = {}  # {representative: [all members]}
    with open(mmseq) as mmsf:
        for line in mmsf:
            repres, other = line.strip().split()
            if repres in clusters:
                clusters[repres].append(other)
            else:
                clusters[repres] = [repres]
    logger.info("Pangenome has {} families.".format(len(clusters)))
    return clusters


def clusters_to_file(clust, fileout):
    """
    Write all clusters to a file

    Parameters
    ----------
    clust : {first_member: [all members = protein names]}
    fileout : filename of pangenome where families must be written

    Returns
    -------
    dict
        families : {famnum: [members]}
    """
    families = {}  # {famnum: [members]}
    with open(fileout, "w") as fout:
        num = 1
        for _, fam in clust.items():
            families[num] = []
            fout.write(str(num))
            for mem in sorted(fam, key=utils.sort_proteins):
                families[num].append(mem)
                fout.write(" " + mem)
            fout.write("\n")
            num += 1
    return families


def create_mmseqs_db(mmseqdb, prt_path, logmmseq):
    """
    Create ffindex of protein bank (prt_path) if not already done. If done, just write a message
    to tell the user that the current existing file will be used.

    Parameters
    ----------
    mmseqdb : str
         path to base filename for output mmseq db
    prt_path : str
        path to the file containing all proteins to cluster
    logmmseq : str
         path to file where logs must be written
    """
    outext = ["", ".index", ".dbtype", ".lookup", "_h", "_h.index", "_h.dbtype"]
    files_existing = []
    if os.path.isfile(mmseqdb):
        for file in [mmseqdb + ext for ext in outext]:
            if not os.path.isfile(file):
                continue
            files_existing.append(file)
        if len(files_existing) != len(outext):
            logger.warning(f"mmseqs database {mmseqdb} already exists, but at least 1 associated "
                            "file (.dbtype, .index etc). is missing. The program will "
                            "remove existing files and recreate the database.")
            files_remaining = copy.deepcopy(files_existing)
            for file in files_existing:
                os.remove(file)  # Delete file
                files_remaining.remove(file)  # Remove file from list of existing files
                logger.details(f"Removing '{file}'.")
            files_existing = copy.deepcopy(files_remaining)
        else:
            logger.warning(f"mmseqs database {mmseqdb} already exists. The program will "
                           "use it.")
            return 0
    if len(files_existing) != len(outext):
        logger.info("Creating database")
        logger.debug("Existing files: {}".format(len(files_existing)))
        logger.debug("Expected extensions: {}".format(len(outext)))
        cmd = f"mmseqs createdb {prt_path} {mmseqdb}"
        msg = (f"Problem while trying to convert database {prt_path} to mmseqs "
               "database format.")
        logger.details(f"MMseqs command: {cmd}")
        with open(logmmseq, "w") as logf:
            utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
