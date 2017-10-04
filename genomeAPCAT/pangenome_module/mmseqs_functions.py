#!/usr/bin/env python3
# coding: utf-8

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
from genomeAPCAT import utils

logger = logging.getLogger("pangnome.mmseqs")


def run_all_pangenome(min_id, clust_mode, outdir, prt_path, threads, panfile=None, quiet=False):
    """
    Run all steps to build a pangenome:
    - create mmseqs database from protein bank
    - cluster proteins
    - convert to pangenome
    """
    # Get general information and file/directory names
    prt_bank = os.path.basename(prt_path)
    start = time.strftime('%Y-%m-%d_%H-%M-%S')
    mmseqdb = os.path.join(outdir, prt_bank + "-msDB")
    information = ("Will run MMseqs2 with:\n"
                   "\t- minimum sequence identity = {}\n"
                   "\t- cluster mode {}").format(min_id, clust_mode)
    if threads > 1:
        information += "\n\t- {} threads".format(threads)
    logger.info(information)

    infoname = get_info(prt_bank, threads, min_id, clust_mode, start)
    logmmseq = get_logmmseq(outdir, prt_bank, infoname)
    # Create ffindex of DB if not already done
    create_mmseqs_db(mmseqdb, prt_path, logmmseq)

    # Cluster with mmseqs
    if panfile:
        panfile = os.path.join(outdir, panfile)
    families, outfile = do_pangenome(outdir, prt_bank, mmseqdb, min_id,
                            clust_mode, threads, start, panfile, quiet)
    return families, outfile


def get_info(prt_bank, threads, min_id, clust_mode, start):
    """
    Get string containing all information on future run
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
    """
    return os.path.join(outdir, "mmseq_" + prt_bank + "_" + infoname + ".log")


def do_pangenome(outdir, prt_bank, mmseqdb, min_id, clust_mode, threads, start, panfile=None,
                 quiet=False):
    """
    Use mmseqs to cluster proteins
    """
    infoname = get_info(prt_bank, threads, min_id, clust_mode, start)
    logmmseq = get_logmmseq(outdir, prt_bank, infoname)
    mmseqclust = os.path.join(outdir, prt_bank + "-clust-" + infoname)
    tmpdir = os.path.join(outdir, "tmp_" + prt_bank + "_" + infoname)
    os.makedirs(tmpdir, exist_ok = True)
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
            while(True):
                if final.ready():
                    break
                remaining = final._number_left
                bar.update()
            bar.finish()
        pool.join()
    # Convert output to tsv file (one line per comparison done)
    families, outfile = mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, start, panfile)
    return families, outfile


def run_mmseqs_clust(args):
    """
    Run mmseqs clustering
    """
    mmseqdb, mmseqclust, tmpdir, logmmseq, min_id, threads, clust_mode = args
    cmd = ("mmseqs cluster {} {} {} --min-seq-id {} --threads {} --cluster-mode "
           "{}").format(mmseqdb, mmseqclust, tmpdir, min_id, threads, clust_mode)
    msg = "Problem while clustering proteins with mmseqs. See log in {}".format(logmmseq)
    with open(logmmseq, "a") as logm:
        utils.run_cmd(cmd, msg, eof=True, stdout=logm, stderr=logm)


def mmseqs_to_pangenome(mmseqdb, mmseqclust, logmmseq, start, outfile=None):
    """
    Convert mmseqs clustering to a pangenome file:
    - convert mmseqs results to tsv file
    - convert tsv file to pangenome
    """
    cmd = "mmseqs createtsv {0} {0} {1} {1}.tsv".format(mmseqdb, mmseqclust)
    msg = "Problem while trying to convert mmseq result file to tsv file"
    with open(logmmseq, "a") as logf:
        utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
    # Convert the tsv file to a 'pangenome' file: one line per family
    families, outfile = mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start, outfile)
    return families, outfile


def mmseqs_tsv_to_pangenome(mmseqclust, logmmseq, start, outfile=None):
    """
    Convert the tsv output file of mmseqs to the pangenome file
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
        logm.write("\n------------\n\nStart: {}".format(start) + "\n")
        logm.write("End: {}".format(end) + "\n")
    return families, outfile


def mmseq_tsv_to_clusters(mmseq):
    """
    Reads the output of mmseq as a tsv file, and converts it to a python dist
    """
    clusters = {}  # {representative: [all members]}
    with open(mmseq, "r") as mmsf:
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
    clust: {first_member: [all members = protein names]}
    """
    families = {}  # {famnum: [members]}
    with open(fileout, "w") as fout:
        num = 1
        for _, fam in clust.items():
            families[num] = fam
            fout.write(str(num))
            for mem in sorted(fam, key=utils.sort_proteins):
                fout.write(" " + mem)
            fout.write("\n")
            num += 1
    return families


def create_mmseqs_db(mmseqdb, prt_path, logmmseq):
    """
    Create ffindex of protein bank if not already done. If done, just write a message
    to tell the user that the current existing file will be used.
    """
    if not os.path.isfile(mmseqdb):
        logger.info("Creating database")
        cmd = "mmseqs createdb {} {}".format(prt_path, mmseqdb)
        msg = ("Problem while trying to convert database {} to mmseqs "
               "database format.").format(prt_path)
        with open(logmmseq, "w") as logf:
            utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
    else:
        logger.warning(("mmseq database {} already exists. The program will "
                        "use it.").format(mmseqdb))

