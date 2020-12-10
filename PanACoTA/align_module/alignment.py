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
For a given family:

- align its proteins using mafft
- back translate to nucleotides
- add "-" of same size as alignment for genomes not having members in the family

@author: GEM, Institut Pasteur
March 2017
"""

import os
import sys
import logging
import multiprocessing
import progressbar
import threading

from PanACoTA import utils

main_logger = logging.getLogger("align.alignment")


def align_all_families(prefix, all_fams, ngenomes, dname, quiet, threads):
    """
    For each family:

    - align all its proteins with mafft
    - back-translate to nucleotides
    - add missing genomes

    Parameters
    ----------
    prefix :  str
        path to ``aldir/<name of dataset>`` (used to get extraction, alignment and btr files
        easily)
    all_fams : []
        list of all family numbers
    ngenomes : int
        total number of genomes in dataset
    dname : str
        name of dataset (used to name concat and grouped files, as well as tree folder)
    quiet : bool
        True if nothing must be written in stdout/stderr, False otherwise
    threads : int
        max number of threads that can be used by mafft

    Returns
    -------
    bool
        True if everything went well, False if there was a problem in at least 1 family.
    """
    main_logger.info(("Starting alignment of all families: protein alignment, "
                      "back-translation to nucleotides, and add missing genomes in the family"))
    nbfam = len(all_fams)
    bar = None
    if not quiet:
        # Create progressbar
        widgets = ['Alignment: ', progressbar.Bar(marker='█', left='', right='', fill=' '),
                   ' ', progressbar.Counter(), "/{}".format(nbfam), ' (',
                   progressbar.Percentage(), ') - ', progressbar.Timer(), ' - '
                   ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbfam,
                                      term_width=79).start()
    final = []
    if threads == 1:
        update_bar = 1
        for num_fam in all_fams:
            f = handle_family_1thread((prefix, num_fam, ngenomes))
            final.append(f)
            bar.update(update_bar)
            update_bar+=1

    else:
        pool = multiprocessing.Pool(threads)

        # Create a Queue to put logs from processes, and handle them after from a single thread
        m = multiprocessing.Manager()
        q = m.Queue()
        # arguments : (gpath, cores_prokka, name, force, nbcont, q) for each genome
        arguments = [(prefix, num_fam, ngenomes, q) for num_fam in all_fams]
        try:
            final = pool.map_async(handle_family, arguments, chunksize=1)
            pool.close()
            # Listen for logs in processes
            lp = threading.Thread(target=utils.logger_thread, args=(q,))
            lp.start()
            if not quiet:
                while True:
                    if final.ready():
                        break
                    remaining = final._number_left
                    bar.update(nbfam - remaining)
                bar.finish()
            pool.join()
            q.put(None)
            lp.join()
            final = final.get()
        # If an error occurs (or user kills with keybord), terminate pool and exit
        except Exception as excp:  # pragma: no cover
            pool.terminate()
            main_logger.error(excp)
            sys.exit(1)
        # We re-aligned at least one family -> remove concatenated files and groupby
    if set(final) != {"OK"}:
        aldir = os.path.split(prefix)[0]
        concat = os.path.join(aldir, f"{dname}-complete.cat.aln")
        outdir = os.path.split(aldir)[0]
        treedir = os.path.join(outdir, "Phylo-" + dname)
        grpfile = os.path.join(treedir, dname + ".grp.aln")
        utils.remove(concat)
        utils.remove(grpfile)
    return False not in final


def handle_family_1thread(args):
    """
    For the given family:

    - align its proteins with mafft
    - back-translate to nucleotides
    - add missing genomes

    Parameters
    ----------
    args : ()
         (prefix, num_fam, ngenomes, q) with:

         - prefix: path to ``aldir/<name of dataset>``
         - num_fam: the current family number
         - ngenomes: the total number of genomes in dataset

    Returns
    -------
    bool
        - "OK" if the files were not re-created, and have the expected format. This is used by
          ``align_all_families`` function, to know if something was regenerated, or if everything
          already existed with the expected format. If something was redone and concat/group files
          exist, it removes them.
        - False if any problem (extractions, alignment, btr, add missing genomes...)
        - True if just generated all files, and everything is ok
    """
    prefix, num_fam, ngenomes = args
    logger = logging.getLogger('align.align_family')
    # Get file names
    prt_file = "{}-current.{}.prt".format(prefix, num_fam)
    gen_file = "{}-current.{}.gen".format(prefix, num_fam)
    miss_file = "{}-current.{}.miss.lst".format(prefix, num_fam)
    mafft_file = "{}-mafft-align.{}.aln".format(prefix, num_fam)
    btr_file = "{}-mafft-prt2nuc.{}.aln".format(prefix, num_fam)
    status1 = family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file,
                               num_fam, ngenomes, logger)
    # If it returned true, Add missing genomes
    if status1:
        return add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, status1, logger)
    return False


def handle_family(args):
    """
    For the given family:

    - align its proteins with mafft
    - back-translate to nucleotides
    - add missing genomes

    Parameters
    ----------
    args : ()
         (prefix, num_fam, ngenomes, q) with:

         - prefix: path to ``aldir/<name of dataset>``
         - num_fam: the current family number
         - ngenomes: the total number of genomes in dataset
         - q: a queue, which will be used by logger to put logs while in other process

    Returns
    -------
    bool
        - "OK" if the files were not re-created, and have the expected format. This is used by
          ``align_all_families`` function, to know if something was regenerated, or if everything
          already existed with the expected format. If something was redone and concat/group files
          exist, it removes them.
        - False if any problem (extractions, alignment, btr, add missing genomes...)
        - True if just generated all files, and everything is ok
    """
    prefix, num_fam, ngenomes, q = args
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('align.align_family')
    return handle_family_1thread((prefix, num_fam, ngenomes))


def add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, status1, logger):
    """
    Once all family proteins are aligned, and back-translated to nucleotides,
    add missing genomes for the family to the alignment with '-'.

    Parameters
    ----------
    btr_file : str
        path to file containing protein alignments back-translated to nucleotides
    miss_file : str
        path to file containing the list of missing genomes in this family
    num_fam : int
        family number
    ngenomes : int
        total number of genomes in dataset
    status1 : bool or str
        - "OK" if we did not redo the alignments as they already were as expected. In that case,
          if missing genomes are already present, just add a warning message saying that we
          used the already existing btr file.
        - True if we just did the alignments and backtranslate. So no warning message needed.
        - False if problem with extraction, alignment or backtranslation (will never happen as
          this function is not called if status1 == False)
    logger : logging.Logger
        the logger, having a queue Handler, to give logs to the main logger

    Returns
    -------
    bool or str
        - "OK" if btr file was not recreated, and already has the right number of sequences,
          and all with the same length.
        - False if problem in btr file alignment, so missing genomes not added
        - True if alignment + adding missing genomes is ok. Can happen if there is no missing
          genome for this family (in that case, btr generated already has the right number of
          sequences), or if we just added the missing genomes.

    """
    # btr_file should always exist.
    # Sometimes it comes from previous step ('missing genomes' are missing)
    # Sometimes it comes from a previous run (all genomes should be here)
    status = check_add_missing(btr_file, num_fam, ngenomes, logger, prev=True)
    # If btr_file has the correct number of sequences, all the same length, return True
    if status is True:
        if status1 == "OK":
            logger.warning(("Alignment already done for family {}. The program will use "
                            "it for next steps").format(num_fam))
            return "OK"
        else:
            return True
    # If btr_files has problem in alignment (not all sequences with same size)
    elif status is False:
        return False
    # All sequences have same length but some genomes are missing -> Add missing genomes
    logger.log(utils.detail_lvl(), "Adding missing genomes for family {}".format(num_fam))
    len_aln = status
    with open(miss_file, "r") as missf, open(btr_file, "a") as btrf:
        for genome in missf:
            genome = genome.strip()
            toadd = ">" + genome + "\n" + "-" * len_aln + "\n"
            btrf.write(toadd)
    ret = check_add_missing(btr_file, num_fam, ngenomes, logger)
    # Return True only if all genomes found with all same length (check_add_missing = True,
    # and not check_add_missing = len_alignment)
    if ret is True:
        return True
    else:
        return False


def check_add_missing(btr_file, num_fam, ngenomes, logger, prev=False):
    """
    Check back-translated alignment while missing genomes have been added

    Parameters
    ----------
    btr_file : str
        path to file containing back-translated alignment
    num_fam : int
        current family number
    ngenomes : int
        total number of genomes in dataset
    logger : logging.Logger
        logger with queueHandler to give logs to main logger
    prev : bool
        True if we are checking alignments before adding missing genomes (in case it comes
        from a previous run and is already done, or in case there are no missing genomes for
        this family, so nothing to do to btr file). In that case, do not write error message
        if the number of sequences does not correspond to the total number of genomes.
        False if we just added missing genomes. In that case, nb sequences should be equal to
        the total number of genomes. If not, write error message.

    Returns
    -------
    bool or int
        - True if right number of sequences in btr file and all the same length (everything is ok)
        - False if problem in alignment (sequences do not all have the same size)
        - alignment length if sequences aligned all have the same length, but missing
          genomes are not added yet (so, will have to add lines with this number of '-')
    """
    res = check_lens(btr_file, num_fam, logger)
    # If all sequences have the same length, get number of sequences
    if res:
        len_aln, nb = res
    # otherwise, return False: problem while aligning or back-translating
    else:
        return False
    if nb != ngenomes:
        if not prev:
            logger.error(("ERROR: family {} contains {} genomes in total instead of the {} "
                          "genomes in input.\n").format(num_fam, nb, ngenomes))
        return len_aln
    return True


def family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file,
                     num_fam, ngenomes, logger):
    """
    From a given family, align all its proteins with mafft, back-translate
    to nucleotides, and add missing genomes in this family.

    Parameters
    ----------
    prt_file : str
        path to file containing proteins extracted
    gen_file : str
        path to file containing genes extracted
    miss_file : str
        path to file containing list of genomes missing
    mafft_file : str
        path to file which will contain the protein alignment
    btr_file : str
        path to file which will contain the nucleotide alignment back-translated from protein
        alignment
    num_fam : int
        current family number
    ngenomes : int
        total number of genomes in dataset
    logger : logging.Logger
        logger with queueHandler to give logs to main logger

    Returns
    -------
    bool or str
        - False if problem with extractions or with alignment or with backtranslation
        - 'nb_seqs' = number of sequences aligned if everything went well (extractions and
          alignment ok, btr created without problem)
        - "OK" if extractions and alignments went well, and btr already exists and is ok
    """
    # Check number of proteins extracted
    nbfprt = check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger)
    nbfal = None
    # If problem with extractions, remove mafft and btr files if they exist, so that they will
    # be regenerated
    if not nbfprt:
        utils.remove(mafft_file)
        utils.remove(btr_file)
        return False
    # If mafft file already exists, check the number of proteins aligned corresponds to number of
    #  proteins extracted. If not, remove mafft and btr files.
    if os.path.isfile(mafft_file):
        message = ("fam {}: Will redo alignment, because found a different number of proteins "
                   "extracted in {} ({}) and proteins aligned in "
                   "existing {}").format(num_fam, prt_file, nbfprt, mafft_file)
        nbfal = check_nb_seqs(mafft_file, nbfprt, logger, message)
        if not nbfal:
            os.remove(mafft_file)
            utils.remove(btr_file)
    # If mafft file does not exist (removed because problem in its alignment, or just not generated
    # yet), remove btr (will be regenerated), and do alignment with mafft
    if not os.path.isfile(mafft_file):
        utils.remove(btr_file)
        nbfal = mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger)
    # If problem with alignment, return False
    if not nbfal:
        return False
    # If btr file already exists, means that it was already done before, and not removed because
    # extractions and mafft files are ok. So, return True, saying that btr file is done,
    # next step will be to check it, add missing genomes etc.
    if os.path.isfile(btr_file):
        message = ("fam {}: Will redo back-translation, because found a different number of "
                   "proteins aligned in {} ({}) and genes back-translated in "
                   "existing {}").format(num_fam, mafft_file, nbfal, btr_file)
        res = check_nb_seqs(btr_file, [nbfal, ngenomes], logger, message)
        if not res:
            utils.remove(btr_file)
        else:
            return "OK"
    # If btr file does not exist (removed because problem with mafft generated before,
    # or just not generated yet), do back-translation, and return True if it went well,
    # False otherwise
    return back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger)


def check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger):
    """
    Check that extractions went well for the given family:

    - check number of proteins and genes extracted compared to the
      number of genomes

    Parameters
    ----------
    num_fam : int
        current family number
    miss_file : str
        path to file containing the list of genomes missing for the current family
    prt_file : str
        path to file containing all proteins extracted
    gen_file : str
        path to file containing all genes extracted
    ngenomes : int
        total number of genomes in dataset
    logger : logging.Logger
        logger with queueHandler to give logs to main logger

    Returns
    -------
    bool or int
        False if any problem (nbmiss+prt != nbgenomes or nbmiss+gen != nbgenomes). If no
        problem, returns the number of proteins/genes extracted
    """
    logger.log(utils.detail_lvl(), "Checking extractions for family {}".format(num_fam))

    # Check that extractions went well
    nbmiss = utils.count(miss_file)
    nbfprt = utils.grep(prt_file, "^>", counts=True)
    nbfgen = utils.grep(gen_file, "^>", counts=True)
    if nbmiss + nbfprt != ngenomes:
        logger.error(("fam {}: wrong sum of missing genomes ({}) and prt extracted ({}) for {} "
                      "genomes in the dataset.").format(num_fam, nbmiss, nbfprt, ngenomes))
        return False
    if nbmiss + nbfgen != ngenomes:
        logger.error(("fam {}: wrong sum of missing genomes ({}) and gen extracted ({}) for {} "
                      "genomes in the dataset.").format(num_fam, nbmiss, nbfgen, ngenomes))
        return False
    return nbfprt


def mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger):
    """
    Align all proteins of the given family with mafft

    Parameters
    ----------
    num_fam : int
        current family number
    prt_file : str
        path to file containing all proteins extracted
    mafft_file : str
        path to file which will contain proteins alignment
    nbfprt : int
        number of proteins extracted in prt file
    logger : logging.Logger
        logger with queueHandler to give logs to main logger

    Returns
    -------
    bool
        True if no problem (alignment ok, same number of proteins extracted and aligned),
        False otherwise
    """
    logger.log(utils.detail_lvl(), "Aligning family {}".format(num_fam))
    cmd = "mafft --auto {}".format(prt_file)
    error = "Problem while trying to align fam {}".format(num_fam)
    stdout = open(mafft_file, "w")
    stderr = open(mafft_file + ".log", "w")
    logger.log(utils.detail_lvl(), cmd)
    ret = utils.run_cmd(cmd, error, stdout=stdout, stderr=stderr, logger=logger)
    stdout.close()
    if not isinstance(ret, int):
        ret = ret.returncode
    if ret != 0:
        os.remove(mafft_file)
        return False
    message = ("fam {}: different number of proteins extracted in {} ({}) and proteins "
               "aligned in {}").format(num_fam, prt_file, nbfprt, mafft_file)
    return check_nb_seqs(mafft_file, nbfprt, logger, message)


def back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger):
    """
    Backtranslate protein alignment to nucleotides

    Parameters
    ----------
    num_fam : int
        current family number. Used for log messages
    mafft_file : str
        path to file containing protein alignments by mafft
    gen_file : str
        path to file containing all sequences, not aligned, in nucleotides. It is used to
        convert the alignment in proteins into a nucleotide alignment
    btr_file : str
        path to the file that will contain the nucleotide alignment
    nbfal : int
        number of sequences aligned for the family by mafft
    logger : logging.Logger
        logger with queueHandler to give logs to main logger

    Returns
    -------
    bool
        - False if problem (back-translation, different number of families...)
        - number of sequences in btr file if everything went well
    """
    logger.log(utils.detail_lvl(), "Back-translating family {}".format(num_fam))
    curpath = os.path.dirname(os.path.abspath(__file__))
    awk_script = os.path.join(curpath, "prt2codon.awk")
    cmd = "awk -f {} {} {}".format(awk_script, mafft_file, gen_file)
    stdout = open(btr_file, "w")
    error = "Problem while trying to backtranslate {} to a nucleotide alignment".format(mafft_file)
    ret = utils.run_cmd(cmd, error, stdout=stdout, logger=logger)
    stdout.close()
    if not isinstance(ret, int):
        ret = ret.returncode
    if ret != 0:
        os.remove(btr_file)
        return False
    message = ("fam {}: different number of proteins aligned in {} ({}) and genes "
               "back-translated in {}").format(num_fam, mafft_file, nbfal, btr_file)
    return check_nb_seqs(mafft_file, nbfal, logger, message)


def check_nb_seqs(alnfile, nbfal, logger, message):
    """
    Check the number of sequences in the given alignment file

    Parameters
    ----------
    alnfile : str
        path to alignment file
    nbfal : int or []
        expected number of sequences, or list of expected numbers of sequences
    logger : logging.Logger
        logger with queueHandler to give logs to main logger
    message : str
        message to write when sequence numbers do not match

    Returns
    -------
    bool or int
        - False if not same number of sequences
        - nbseqs in align file if found among values in 'nbfal'
    """
    nbseqs = utils.grep(alnfile, "^>", counts=True)
    if isinstance(nbfal, int):
        nbfal = [nbfal]
    for num in nbfal:
        if nbseqs == num:
            return nbseqs

    logger.error(message + ' ({})'.format(nbseqs))
    return False


def check_lens(aln_file, num_fam, logger):
    """
    In the given alignment file, check that all sequences have the same length.
    If there is no problem, it returns the length of alignment, and the number
    of sequences aligned.

    Parameters
    ----------
    aln_file : str
        path to the alignment file to check
    num_fam : int
        current family number. Used for log message if problem
    logger : logging.Logger
        logger having a queue Handler to give logs to the main logger in the main process

    Returns
    -------
    bool or tuple
        False if there is a problem in the alignment (sequences do not all have the
        same length).
        If they all have the same length, returns this length and the number of sequences
    """
    nb_gen = 0
    all_sums = set()
    with open(aln_file, "r") as btrf:
        cur_sum = 0
        start = True
        for line in btrf:
            if line.startswith(">"):
                nb_gen += 1
                if not start:
                    all_sums.add(cur_sum)
                    cur_sum = 0
            else:
                start = False
                cur_sum += len(line.strip())
        all_sums.add(cur_sum)
    if len(all_sums) > 1:
        logger.error(("alignments for family {} do not all have the same "
                      "length. Lengths found are: {}\n").format(num_fam, all_sums))
        return False
    return list(all_sums)[0], nb_gen
