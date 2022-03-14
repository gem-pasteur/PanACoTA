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
Functions to deal with prokka or prodigal only, according to user request


@author gem
April 2017
"""

import os
import sys
import shutil
import glob
import logging
import subprocess
import shlex
import multiprocessing
import progressbar
import threading

import PanACoTA.utils as utils

logger = logging.getLogger('annotate.run_annotation_all')


def run_annotation_all(genomes, threads, force, annot_folder, fgn, prodigal_only=False,
                       small=False, quiet=False):
    """
    For each genome in genomes, run prokka (or only prodigal) to annotate the genome.

    Parameters
    ----------
    genomes : dict
        {genome: [gembase_name, path_to_origfile, path_split_gembase, gsize, nbcont, L90]}
    threads : int
        max number of threads that can be used
    force : bool
        if False, do not override prokka/prodigal outdir and result dir if they exist. If\
        True, rerun prokka and override existing results, for all genomes.
    annot_folder : str
        folder where prokka/prodigal results must be written: for each genome,
        a directory <genome_name>-prokkaRes or <genome_name>-prodigalRes> will be created
        in this folder, and all the results
        of prokka/prodigal for the genome will be written inside
    fgn : str
        name (key in genomes dict) of the fist genome, which will be used for prodigal training
    prodigal_only : bool
        True if only prodigal must run, False if prokka must run
    small : bool
        True -> use -p meta option with prodigal. Do not use training
    quiet : bool
        True if nothing must be written to stderr/stdout, False otherwise

    Returns
    -------
    dict
        {genome: boolean} -> with True if prokka/prodigal ran well, False otherwise.
    """

    # Update information according to annotation soft used and write message
    if prodigal_only:
        message = "Annotating all genomes with prodigal"
        run_annot = run_prodigal
        main_logger = logging.getLogger("annotate.prodigal")
    else:
        message = "Annotating all genomes with prokka"
        run_annot = run_prokka
        main_logger = logging.getLogger("annotate.prokka")
    main_logger.info(message)
    # Get total number of genomes to annotate, used to show annotation progress
    nbgen = len(genomes)
    bar = None
    # If user did not ask for quiet, create progressbar
    if not quiet:
        # Create progress bar
        widgets = ['Annotation: ', progressbar.Bar(marker='█', left='', right=''),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ') - ', progressbar.Timer(), ' - '
                  ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen,
                                      term_width=79).start()
    # Get resource availability:
    # - number of threads used by prokka/prodigal (cores_annot)
    # - how many genomes can be annotated at the same time (pool_size)
    # - prodigal does not run with several threads: with prodigal, always cores_annot == 1
    # and pool_size == threads
    gpath_train = ""  # by default, no training genome
    if prodigal_only:
        cores_annot = 1
        pool_size = threads
        # If prodigal, train on the first genome
        # fgn is key of genomes, genomes[fgn] = [_,_,annote_file,_,_,_]
        gtrain = genomes[fgn][2]
        # If problem, gpath_train will be empty, but this will be checked while
        # trying to run prodigal, because we also need to check that genomes are not simply
        # already annotated
        if not small:
            gpath_train = prodigal_train(gtrain, annot_folder)
        else:
            gpath_train = "small option"
    elif threads <= 3:
        # less than 3 threads: run prokka 1 by 1 with all threads
        cores_annot = threads
        pool_size = 1
    else:
        # use multiprocessing
        # if there are more threads than genomes, use as many threads as possible per genome
        if len(genomes) <= threads:
            cores_annot = int(threads / len(genomes))
        # otherwise, use 2 threads per genome (and nb_genome/2 genomes at the same time)
        else:
            cores_annot = 2
        pool_size = int(threads / cores_annot)
    #  Create pool with a given size (=number of tasks to be launched in parallel)
    pool = multiprocessing.Pool(pool_size)
    # Create a Queue to put logs from processes, and handle them after from a single thread
    m = multiprocessing.Manager()
    q = m.Queue()
    # {genome: [gembase_name, path_to_origfile, path_toannotate_file, gsize, nbcont, L90]}
    # arguments: gpath, prok_folder, threads, name, force, nbcont, small(for prodigal), q
    arguments = [(genomes[g][2], annot_folder, cores_annot, genomes[g][0],
                  force, genomes[g][4], gpath_train, q)
                 for g in sorted(genomes)]
    try:
        # Start pool (run 'run_annot' n each set of arguments)
        final = pool.map_async(run_annot, arguments, chunksize=1)
        # Close pool: no more data will be put on this pool
        pool.close()
        # Listen for logs in processes
        lp = threading.Thread(target=utils.logger_thread, args=(q,))
        lp.start()
        if not quiet:
            while True:
                # Final is ready when all pool elements are done
                if final.ready():
                    break
                # If not done, get number of genomes left
                remaining = final._number_left
                # Add this to start progressbar with 0% instead of N/A%
                if remaining == nbgen:
                    bar.update(0.0000001)
                else:
                    # Update progress bar
                    bar.update(nbgen - remaining)
            # End progress bar
            bar.finish()
        pool.join()
        # Put None to tell 'q' that everything is finished. It can stopped and be joined.
        q.put(None)
        # join lp (tell to stop once all log processes are done, which is the case here)
        lp.join()
        final = final.get()
    # # If user stops programm (ctrl+C), end it
    # except KeyboardInterrupt as ki:
    #     print("error")
    #     for worker in pool._pool:
    #         print(worker.is_alive())
    #     pool.join()
    #     print("closed")
    #     pool.terminate()
    #     print("--------------terminate ok----------------")
    #     lp.join()
    #     print("thread stopped")
    #     # run_event.clear()
    #     # lp.terminate()
    #     # print("--------------JOIN--------------")
    #     # pool.terminate()
    #     main_logger.error("Process killed by CTRL+C")
    #     return "coucou"
    # If an error occurs, terminate pool, write error and exit
    except Exception as excp:  # pragma: no cover
        pool.terminate()
        main_logger.error(excp)
        sys.exit(1)
    final = {genome: res for genome, res in zip(sorted(genomes), final)}
    return final


def prodigal_train(gpath, annot_folder):
    """
    Use prodigal training mode.
    First, train prodigal on the first genome ('gpath'), and write it to 'genome'.trn,
    file which will be used for the annotation of all next sequence
    Parameters
    ----------
    gpath : str
        path to genome to train on
    annot_folder : str
        path to folder where the log files and train file will be saved

    Returns
    -------
    str
        path and name of train file (will be used to annotate all next genomes)
        If problem, returns empty string
    """
    logger.info(f"Prodigal will train using {gpath}")
    gname = os.path.basename(gpath)             # path/to/original/genome.fasta -> genome.fasta
    gpath_train = os.path.join(annot_folder, gname + ".trn") # path/to/prodiRes/genome.fasta.trn
    if os.path.isfile(gpath_train):
        logger.info(f"A training file already exists ({gpath_train}). "
                     "It will be used to annotate all genomes.")
        return gpath_train
    prodigal_logfile = gpath_train + "-prodigal-train.log"  # path/to/genome-prodigal-train.log
    prodigal_logfile_err = gpath_train + "-prodigal-train.log.err"
    cmd = (f"prodigal -i {gpath} -t {gpath_train}")
    error = (f"Error while trying to train prodigal on {gname}. See {prodigal_logfile_err}.")
    logger.log(utils.detail_lvl(), "prodigal command: " + cmd)
    prodigalf = open(prodigal_logfile, "w")
    prodigalferr = open(prodigal_logfile_err, "w")
    ret = utils.run_cmd(cmd, error, eof=False, stderr=prodigalferr, stdout=prodigalf,
                        logger=logger)
    prodigalf.close()
    prodigalferr.close()
    if ret.returncode == 0:
        logger.log(utils.detail_lvl(), f"End training on {gpath}")
        return gpath_train
    else:
        return ""


def run_prokka(arguments):
    """
    Run prokka for the given genome.

    Parameters
    ----------
    arguments : tuple
        (gpath, prok_folder, cores_annot, name, force, nbcont, small, q) with:

        * gpath: path and filename of genome to annotate
        * prok_folder: path to folder where all prokka folders for all genomes are saved
        * cores_annot: how many cores can use prokka
        * name: output name of annotated genome
        * force: True if force run (override existing files), False otherwise
        * nbcont: number of contigs in the input genome, to check prokka results
        * small: used for prodigal, if sequences to annotate are small. Not used here
        * q : queue where logs are put

    Returns
    -------
    boolean
        True if eveything went well (all needed output files present,
        corresponding numbers of proteins, genes etc.). False otherwise.
    """
    gpath, prok_folder, threads, name, force, nbcont, _, q = arguments
    # Set logger for this process
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('annotate.run_prokka')
    logger.log(utils.detail_lvl(), f"Start annotating {name} from {gpath} with Prokka")

    # Define prokka directory and logfile, and check their existence
    prok_dir = os.path.join(prok_folder, os.path.basename(gpath) + "-prokkaRes")
    fnull = open(os.devnull, 'w')
    prok_logfile = os.path.join(prok_folder, os.path.basename(gpath) + "-prokka.log")
    # import sys
    # sys.exit(1)
    # If result dir already exists, check if we can use it or next step or not
    if os.path.isdir(prok_dir) and not force:
        logger.warning(f"Prokka results folder {prok_dir} already exists.")
        ok = check_prokka(prok_dir, prok_logfile, name, gpath, nbcont, logger)
        # If everything ok in the result dir, do not rerun prokka,
        # use those results for next step (formatting)
        if ok:
            logger.log(utils.detail_lvl(), "Prokka did not run again, "
                       "formatting step used already generated results of "
                       f"Prokka in {prok_dir}. If you want to re-run prokka, first "
                       "remove this result folder, or use '-F' or '--force' "
                       "option if you want to rerun prokka for all genomes.")
            logger.log(utils.detail_lvl(), f"End annotating {name} {gpath}")
        # If missing files, or other problems in result dir, error message,
        # ask user to force or remove this folder.
        else:
            logger.warning("Problems in the files contained in your already existing output dir "
                           "({}). Please check it, or remove it to "
                           "re-annotate.".format(prok_dir))
        # If everything was ok -> everything is ready for next step -> return True
        # If something is wrong -> cannot use those results, genome won't be annotated
        # -> return False
        return ok
    # If result dir exists but user wants to force, remove this result dir
    elif os.path.isdir(prok_dir) and force:
        shutil.rmtree(prok_dir)
        logger.warning("Prokka results folder already exists, but removed because --force option "
                       "used")
    # Now that we checked and solved those cases:
    #     - outdir exists (problems or not, we returned appropriate boolean)
    #     - if outdir exists exists but force, remove this outdir.
    # So, outdir does not exist -> run prokka
    cmd = (f"prokka --outdir {prok_dir} --cpus {threads} "
           f"--prefix {name} --centre prokka {gpath}")
    error = (f"Error while trying to run prokka on {name} from {gpath}")
    logger.log(utils.detail_lvl(), "Prokka command: " + cmd)
    prokf = open(prok_logfile, "w")
    ret = utils.run_cmd(cmd, error, eof=False, stderr=prokf, logger=logger)
    prokf.close()
    if ret.returncode != 0:
        return False
    ok = check_prokka(prok_dir, prok_logfile, name, gpath, nbcont, logger)
    logger.log(utils.detail_lvl(), f"End annotating {name} from {gpath}.")
    return ok


def run_prodigal(arguments):
    """
    Run prodigal for the given genome.

    Parameters
    ----------
    arguments : tuple
        (gpath, prodigal_folder, cores_annot, name, force, nbcont, q) with:

        * gpath: path and filename of genome to annotate
        * prodigal_folder: path to folder where all prodigal folders for all genomes are saved
        * cores_annot: how many cores can use prodigal
        * name: output name of annotated genome
        * force: True if force run (override existing files), False otherwise
        * nbcont: number of contigs in the input genome, to check prodigal results
        * small: ifcontigs are too small (<20000bp), use -p meta option
        * q : queue where logs are put

    Returns
    -------
    boolean
        True if eveything went well (all needed output files present,
        corresponding numbers of proteins, genes etc.). False otherwise.
    """
    gpath, prodigal_folder, threads, name, force, nbcont, gpath_train, q = arguments
    # Set logger for this process, which will be given to all subprocess
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('annotate.run_prodigal')
    # Define prodigal directory and logfile, and check their existence
    # By default, prodigal is in tmp_folder -> resdir/tmp_files/genome-prodigalRes
    g_ori_name = os.path.basename(gpath)
    prodigal_dir = os.path.join(prodigal_folder, g_ori_name + "-prodigalRes")
    prodigal_logfile = os.path.join(prodigal_folder, g_ori_name + "-prodigal.log")
    prodigal_logfile_err = os.path.join(prodigal_folder, g_ori_name + "-prodigal.log.err")

    # If result dir exists but user wants to force, remove this result dir
    if os.path.isdir(prodigal_dir) and force:
        shutil.rmtree(prodigal_dir)
        logger.warning("Prodigal results folder already exists, but is removed because "
                       "--force option was used.")

    # Training file can be "small option", meaning that we did not use the training mode.
    # If not "small option", we used the training mode. If training file does not exist 
    # and prodigal result directory neither, return False
    # We cannot annotate using nothing.
    # Happens if there was a problem while training
    if (gpath_train != "small option" and not os.path.isfile(gpath_train) 
        and not os.path.isdir(prodigal_dir)):
        return False

    logger.log(utils.detail_lvl(), f"Start annotating {name} (from {gpath} sequence) "
                                     "with Prodigal")
    # If prodigal results dir already exists (meaning user did not want to force,
    # otherwise it would have been deleted just before),
    # can we use it for next step ? -> check content.
    if os.path.isdir(prodigal_dir):
        logger.warning(f"Prodigal results folder {prodigal_dir} already exists.")
        ok = check_prodigal(gpath, name, prodigal_dir, logger)
        # If everything ok in the result dir, do not rerun prodigal,
        # use those results for next step (formatting)
        if ok:
            logger.log(utils.detail_lvl(), "Prodigal did not run again. "
                                           "Formatting step will use already generated results of "
                                           "Prodigal in {}. If you want to re-run Prodigal, first "
                                           "remove this result folder, or use '-F' or '--force' "
                                           "option.".format(prodigal_dir))

            logger.log(utils.detail_lvl(), f"End annotating {name} (from {gpath})")
        # If missing files, or other problems in result dir, error message,
        # ask user to force or remove this folder.
        else:
            logger.warning("Problems in the files contained in your already existing output dir "
                           f"({prodigal_dir}). Please check it, or remove it to "
                           "re-annotate.")
        # If everything was ok -> everything is ready for next step -> return True
        # If something is wrong -> cannot use those results, genome won't be annotated
        # -> return False
        return ok
    else:
        # We are sure prodigal result dir does not exist yet, because either:
        #     - never existed
        #     - removed because user asked to force
        #     - exists but left function, so does not go until this line
        #        -> either if files inside are ok or not
        # So make prodigal_dir (not automatically created by prodigal)
        os.makedirs(prodigal_dir)

    # Prodigal_directory is empty and ready to get prodigal results
    basic_outname = os.path.join(prodigal_dir, name)
    # Define cmd, stderr and stdout files, and error to write if problem.
    error = (f"Error while trying to run prodigal. See {prodigal_logfile_err}.")
    prodigalf = open(prodigal_logfile, "w")
    prodigalferr = open(prodigal_logfile_err, "w")
    if gpath_train == "small option":
        training = "-p meta"
    else:
        training = f"-t {gpath_train}"
    cmd = (f"prodigal -i {gpath} -d {basic_outname + '.ffn'} -a {basic_outname + '.faa'} "
           f"-f gff -o {basic_outname + '.gff'} {training} -q")
    logger.log(utils.detail_lvl(), "Prodigal command: " + cmd)

    ret = utils.run_cmd(cmd, error, eof=False, stderr=prodigalferr, stdout=prodigalf,
                        logger=logger)
    prodigalf.close()
    prodigalferr.close()
    if ret.returncode == 0:
        logger.log(utils.detail_lvl(), f"End annotating {name} (from {gpath})")
        return True
    else:
        return False


def check_prokka(outdir, logf, name, gpath, nbcont, logger):
    """
    Prokka writes everything to stderr, and always returns a non-zero return code. So, we
    check if it ran well by checking the content of output directory.
    This function is also used when prokka files already exist (prokka was run previously), to
    check if everything is ok before going to next step.

    Parameters
    ----------
    outdir : str
        output directory, where all files are written by prokka
    logf : str
        prokka/prodigal logfile, containing stderr of prokka
    name : str
        genome name in gembase format
    gpath : str
        path to fasta file given as input for prokka
    nbcont : int
        number of contigs in fasta file given to prokka
    logger : logging.Logger
        logger object to get logs

    Returns
    -------
    bool
        True if everything went well, False otherwise
    """
    missing_file = False
    problem = False
    if not os.path.isdir(outdir):
        logger.error(("Previous annotation could not run properly. "
                      "Look at {} for more information.").format(logf))
        missing_file = True
    else:
        oriname = os.path.basename(gpath)
        fnafile = glob.glob(os.path.join(outdir, "*.fna"))
        tblfile = glob.glob(os.path.join(outdir, "*.tbl"))
        faafile = glob.glob(os.path.join(outdir, "*.faa"))
        ffnfile = glob.glob(os.path.join(outdir, "*.ffn"))
        gfffile = glob.glob(os.path.join(outdir, "*.gff"))
        if len(fnafile) == 0:
            logger.error("{} {}: no .fna file".format(name, oriname))
            missing_file = True
        elif len(fnafile) > 1:
            logger.error("{} {}: several .fna files".format(name, oriname))
            missing_file = True
        if len(tblfile) == 0:
            logger.error("{} {}: no .tbl file".format(name, oriname))
            missing_file = True
        elif len(tblfile) > 1:
            logger.error("{} {}: several .tbl files".format(name, oriname))
            missing_file = True
        if len(faafile) == 0:
            logger.error("{} {}: no .faa file".format(name, oriname))
            missing_file = True
        elif len(faafile) > 1:
            logger.error("{} {}: several .faa files".format(name, oriname))
            missing_file = True
        if len(ffnfile) == 0:
            logger.error("{} {}: no .ffn file".format(name, oriname))
            missing_file = True
        elif len(ffnfile) > 1:
            logger.error("{} {}: several .ffn files".format(name, oriname))
            missing_file = True
        if len(gfffile) == 0:
            logger.error("{} {}: no .gff file".format(name, oriname))
            missing_file = True
        elif len(gfffile) > 1:
            logger.error("{} {}: several .gff files".format(name, oriname))
            missing_file = True
        if not missing_file:
            tblfile = tblfile[0]
            faafile = faafile[0]
            ffnfile = ffnfile[0]
            fnbcont, tnb_cds, nb_gene = count_tbl(tblfile)
            faaprot = count_headers(faafile)
            ffngene = count_headers(ffnfile)
            if nbcont != fnbcont:
                logger.error(("{} {}: no matching number of contigs; "
                              "nbcontig={}; in tbl ={}").format(name, oriname, nbcont, fnbcont))
                problem = True
            if tnb_cds != faaprot:
                logger.error(("{} {}: no matching number of proteins between tbl and faa; "
                              "faa={}; in tbl ={}").format(name, oriname, faaprot, tnb_cds))
                problem = True
    return not problem and not missing_file


def check_prodigal(gpath, name, prodigal_dir, logger):
    """
    When prodigal result folder already exists, check that the ouput files exist.
    We cannot check all content, but check that they are present.

    Parameters
    ----------
    gpath : str
        path to fasta file given as input for prokka
    name : str
        genome name in gembase format
    prodigal_dir : str
        output directory, where all files are written by prodigal
    logger : logging.Logger
        logger object to get logs

    Returns
    -------
    bool
        True if everything went well, False otherwise
    """
    oriname = os.path.basename(gpath)
    faafile = glob.glob(os.path.join(prodigal_dir, "*.faa"))
    ffnfile = glob.glob(os.path.join(prodigal_dir, "*.ffn"))
    gfffile = glob.glob(os.path.join(prodigal_dir, "*.gff"))
    missing_file = False

    if len(faafile) != 1:
        logger.error("{} {}: no or several .faa file(s)".format(name, oriname))
        missing_file = True
    if len(ffnfile) !=  1:
        logger.error("{} {}: no or several .ffn file(s)".format(name, oriname))
        missing_file = True
    if len(gfffile) != 1:
        logger.error("{} {}: no or several .gff file(s)".format(name, oriname))
        missing_file = True

    # If we have all result files, check they are not empty
    if not missing_file:
        if (os.path.getsize(faafile[0]) == 0 or os.path.getsize(ffnfile[0]) == 0
            or os.path.getsize(gfffile[0]) == 0):
            logger.error("Genome {} (from {}): At least one of your Prodigal result file "
                         "is empty.".format(name, oriname))
            return False
    return not missing_file


def count_tbl(tblfile):
    """
    Count the different features found in the tbl file:

    - number of contigs
    - number of proteins (CDS)
    - number of genes (locus_tag)
    - number of CRISPR arrays (repeat_region) -> ignore crisprs

    Parameters
    ----------
    tblfile : str
        tbl file generated by prokka

    Returns
    -------
    (nbcont, nb_cds, nb_gene, nb_crispr)
        information on features found in the tbl file.
    """
    nbcont = 0
    nb_cds = 0
    nb_gene = 0
    # nb_crispr = 0
    with open(tblfile) as tblf:
        for line in tblf:
            if line.startswith(">"):
                nbcont += 1
            if "CDS" in line:
                nb_cds += 1
            if "locus_tag" in line:
                nb_gene += 1
            # if "repeat_region" in line or (len(line.split()) == 3 and "CRISPR" in line):
            #     nb_crispr += 1
    return nbcont, nb_cds, nb_gene  #, nb_crispr


def count_headers(seqfile):
    """
    Count how many sequences there are in the given file

    Parameters
    ----------
    seqfile : str
        file containing a sequence in multi-fasta format

    Returns
    -------
    int
        number of contigs in the given multi-fasta file
    """
    nbseq = 0
    with open(seqfile) as faaf:
        for line in faaf:
            if line.startswith(">"):
                nbseq += 1
    return nbseq
