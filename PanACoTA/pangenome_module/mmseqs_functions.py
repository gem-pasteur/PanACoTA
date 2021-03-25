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
import sys
import time
import threading
import progressbar
import copy

from PanACoTA import utils
from PanACoTA import utils_pangenome as utils_pan

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
    # start = time.strftime('%Y-%m-%d_%H-%M-%S')
    information = ("Will run MMseqs2 with:\n"
                   f"\t- minimum sequence identity = {min_id*100}%\n"
                   f"\t- cluster mode {clust_mode}")
    if threads > 1:
        information += f"\n\t- {threads} threads"
    logger.info(information)
    infoname = get_info(threads, min_id, clust_mode)
    logmmseq = get_logmmseq(outdir, prt_bank, infoname)

    tmpdir = os.path.join(outdir, "tmp_" + prt_bank + "_" + infoname)
    mmseqdb = os.path.join(tmpdir, prt_bank + "-msDB")
    mmseqclust = os.path.join(tmpdir, prt_bank + "-clust-" + infoname)
    mmseqstsv = mmseqclust + ".tsv"
    # Get pangenome filename
    if not panfile:
        base = os.path.basename(mmseqstsv)
        panfile = os.path.join(outdir, f"PanGenome-{prt_bank}-clust-{infoname}.lst")
    else:
        panfile = os.path.join(outdir, panfile)
    # If pangenome file already exists, read it to get families
    if os.path.isfile(panfile):
        logger.warning(f"Pangenome file {panfile} already exists. PanACoTA will read it to get families.")
        _, families, _ = utils_pan.read_pan_file(panfile, logger)
    else:
        os.makedirs(tmpdir, exist_ok=True)
        # Create ffindex of DB if not already done
        status = do_mmseqs_db(mmseqdb, prt_path, logmmseq, quiet)
        # status = create_mmseqs_db(mmseqdb, prt_path, logmmseq)
        # Status = ok means that mmseqs_db files already existed and were not re-done
        # If they were redone (or just done), remove any existing following file (mmseqs clust, tsv, csv)
        # Cluster with mmseqs
        families, panfile = do_pangenome(outdir, prt_bank, mmseqdb, mmseqclust, tmpdir, logmmseq, min_id,
                                         clust_mode, status, threads, panfile, quiet)
    return families, panfile


def get_info(threads, min_id, clust_mode):
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

    Returns
    -------
    str
        string containing info on current run, to put in filenames
    """
    if threads != 1:
        threadinfo = "-th" + str(threads)
    else:
        threadinfo = ""
    infoname = str(min_id) + "-mode" + str(clust_mode) + threadinfo
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


def do_mmseqs_db(mmseqdb, prt_path, logmmseq, quiet):
    """
    Runs create_mmseqs_db with an "infinite progress bar" in the background.
    
    create_mmseqs_db does :
    Create ffindex of protein bank (prt_path) if not already done. If done, just write a message
    to tell the user that the current existing file will be used.

    Parameters
    ----------
    mmseqdb : str
         path to base filename for output of mmseqs createdb
    prt_path : str
        path to the file containing all proteins to cluster
    logmmseq : str
         path to file where logs must be written
    quiet : bool
        True if no output in stderr/stdout, False otherwise

    Returns
    -------
    bool
        True if mmseqs db just created, False if already existed
    """
    logger.info("Creating database")
    try:
        stop_bar = False
        if quiet:
            widgets = []
        # If not quiet, start a progress bar while clustering proteins. We cannot guess
        # how many time it will take, so we start an "infinite" bar, and send it a signal
        # when it has to stop. If quiet, we start a thread that will immediatly stop
        else:
            widgets = [progressbar.BouncingBar(marker=progressbar.RotatingMarker(markers="◐◓◑◒")),
                       "  -  ", progressbar.Timer()]
        x = threading.Thread(target=utils.thread_progressbar, args=(widgets, lambda : stop_bar,))
        x.start()
        res = create_mmseqs_db(mmseqdb, prt_path, logmmseq)
    # except KeyboardInterrupt: # pragma: no cover
    except: # pragma: no cover
        stop_bar = True
        x.join()
        sys.exit(1)
    # Clustering done, stop bar and join (if quiet, it was already finished, so we just join it)
    stop_bar = True
    x.join()
    return res


def do_pangenome(outdir, prt_bank, mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, clust_mode, 
                just_done, threads, panfile, quiet=False):
    """
    Use mmseqs to cluster proteins

    Parameters
    ----------
    outdir : str
        directory where output files are saved
    prt_bank : str
        name of the file containing all proteins to cluster, without path
    mmseqdb : str
        path to base filename of output of mmseqs createdb
    mmseqclust : str
        mmseqs clust
    tmp_dir : str
        path to tmp directory
    logmmseq : str
        path to file for mmseqs logs
    min_id : float
        min percentage of identity to be considered in the same family (between 0 and 1)
    clust_mode : [0, 1, 2]
        0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'
    just_done : str
        True if mmseqs db was just (re)created -> remove mmseqs clust. 
        False if mmseqs db was kept from previous run -> no need to rerun mmseqs clust if already exists
    threads : int
        max number of threads to use
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
    mmseqstsv = mmseqclust + ".tsv"
    # If we just made the database, we must redo all next steps
    # -> if existing, remove
    # mmseqsclust (created by run_mmseqs_clust)
    # mmseqstsv (created by mmseqs_to_pangenome)
    # pangenome file
    if just_done and os.path.isfile(mmseqclust) or os.path.isfile(mmseqstsv) or os.path.isfile(panfile):
        logger.details("Removing existing clustering and/or pangenome files.")
        utils.remove(mmseqclust)
        utils.remove(mmseqstsv)
        utils.remove(panfile)
    bar = None
    logger.debug(mmseqclust)
    if os.path.isfile(mmseqclust):
        logger.warning((f"mmseqs clustering {mmseqclust} already exists. The program will now convert "
                        "it to a pangenome file."))
    else:
        logger.info("Clustering proteins...")
        try:
            stop_bar = False
            if quiet:
                widgets = []
            # If not quiet, start a progress bar while clustering proteins. We cannot guess
            # how many time it will take, so we start an "infinite" bar, and send it a signal
            # when it has to stop. If quiet, we start a thread that will immediatly stop
            else:
                widgets = [progressbar.BouncingBar(marker=progressbar.RotatingMarker(markers="◐◓◑◒")),
                           "  -  ", progressbar.Timer()]
            x = threading.Thread(target=utils.thread_progressbar, args=(widgets, lambda : stop_bar,))
            x.start()
            args = (mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode)
            run_mmseqs_clust(args)
        # except KeyboardInterrupt: # pragma: no cover
        except: # pragma: no cover
            stop_bar = True
            x.join()
            sys.exit(1)
        # Clustering done, stop bar and join (if quiet, it was already finished, so we just join it)
        stop_bar = True
        x.join()
    # Convert output to tsv file (one line per comparison done)
    #  # Convert output to tsv file (one line per comparison done)
    # -> returns (families, outfile)
    families = mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, panfile)
    return families, panfile


def run_mmseqs_clust(args):
    """
    Run mmseqs clustering

    Parameters
    ----------
    args : tuple
         (mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode), with:

            * mmseqdb: path to base filename (output created by mmseq db)
            * mmseqclust: path to base filename for output of mmseq clustering
            * tmpdir : path to folder which will contain mmseq temporary files
            * logmmseq : path to file where logs must be written
            * min_id : min percentage of identity to be considered in the same family
            *         (between 0 and 1)
            * threads : max number of threads to use
            * clust_mode : [0, 1, 2], 0 for 'set cover', 1 for 'single-linkage', 2 for 'CD-Hit'

    """
    mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode = args
    cmd = (f"mmseqs cluster {mmseqdb} {mmseqclust} {tmpdir} --min-seq-id {min_id} --threads {threads} --cluster-mode "
           f"{clust_mode}")
    logger.details(f"MMseqs command: {cmd}")
    msg = f"Problem while clustering proteins with mmseqs. See log in {logmmseq}"
    with open(logmmseq, "a") as logm:
        utils.run_cmd(cmd, msg, eof=False, stdout=logm, stderr=logm)


def mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, outfile):
    """
    Convert mmseqs clustering to a pangenome file:

    - convert mmseqs results to tsv file
    - convert tsv file to pangenome

    Parameters
    ----------
    mmseqdb : str
         path to base filename of output of mmseqs createdb
    mmseqclust : str
        path to base filename of output of mmseqs cluster
    logmmseq : str
         path to file where logs must be written
    outfile : str
        pangenome filename

    Returns
    -------
    dict
        - families : {fam_num: [all members]}
    """
    cmd = f"mmseqs createtsv {mmseqdb} {mmseqdb} {mmseqclust} {mmseqclust}.tsv"
    msg = "Problem while trying to convert mmseq result file to tsv file"
    logger.details(f"MMseqs command: {cmd}")
    with open(logmmseq, "a") as logf:
        utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
    # Convert the tsv file to a 'pangenome' file: one line per family
    families = mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, outfile)
    return families


def mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, outfile):
    """
    Convert the tsv output file of mmseqs to the pangenome file

    Parameters
    ----------
    mmseqclust : str
        path to base filename for output of mmseq clustering
    logmmseq : str
        path to file where logs must be written
    outfile : str
        pangenome filename, or None if default one must be used

    Returns
    -------
    dict

        - families : {fam_num: [all members]}
    """
    logger.info("Converting mmseqs results to pangenome file")
    tsvfile = mmseqclust + ".tsv"
    clusters = mmseq_tsv_to_clusters(tsvfile)
    families = clusters_to_file(clusters, outfile)
    end = time.strftime('%Y-%m-%d_%H-%M-%S')
    with open(logmmseq, "a") as logm:
        logm.write(f"End: {end}")
    return families


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
         path to base filename for output of mmseqs createdb
    prt_path : str
        path to the file containing all proteins to cluster
    logmmseq : str
         path to file where logs must be written


    Returns
    -------
    bool
        True if mmseqs db just created, False if already existed
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
            return False
    logger.debug("Existing files: {}".format(len(files_existing)))
    logger.debug("Expected extensions: {}".format(len(outext)))
    cmd = f"mmseqs createdb {prt_path} {mmseqdb}"
    msg = (f"Problem while trying to convert database {prt_path} to mmseqs "
           "database format.")
    logger.details(f"MMseqs command: {cmd}")
    with open(logmmseq, "w") as logf:
        utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
    return True
