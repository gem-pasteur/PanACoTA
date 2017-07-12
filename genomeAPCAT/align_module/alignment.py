#!/usr/bin/env python3
# coding: utf-8

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

from genomeAPCAT import utils



def align_all_families(prefix, all_fams, ngenomes, quiet, threads):
    """
    For each family:
    - align all its proteins with mafft
    - back-translate to nucleotides
    - add missing genomes
    """
    main_logger = logging.getLogger("align.alignment")
    main_logger.info(("Starting alignment of all families: protein alignment, "
                      "back-translation to nucleotides, and add missing genomes in the family"))
    nbfam = len(all_fams)
    if not quiet:
        # Create progressbar
        widgets = ['Alignment: ', progressbar.Bar(marker='█', left='', right='', fill=' '),
                   ' ', progressbar.Counter(), "/{}".format(nbfam), ' (',
                   progressbar.Percentage(), ') - ', progressbar.Timer(), ' - '
                  ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbfam,
                                      term_width=100).start()
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
            while(True):
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
    return set(final) == set([True])


def handle_family(args):
    """
    For a given family:
    - align its proteins with mafft
    - back-translate to nucleotides
    - add missing genomes
    """
    prefix, num_fam, ngenomes, q = args
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('run_prokka')
    # Get file names
    prt_file = "{}-current.{}.prt".format(prefix, num_fam)
    gen_file = "{}-current.{}.gen".format(prefix, num_fam)
    miss_file = "{}-current.{}.miss.lst".format(prefix, num_fam)
    mafft_file = "{}-mafft-align.{}.aln".format(prefix, num_fam)
    btr_file = "{}-mafft-prt2nuc.{}.aln".format(prefix, num_fam)
    status1 = family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file,
                               num_fam, ngenomes, logger)
    # If it returned true, Add missing genomes
    status2 = False
    if status1:
        status2 = add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, logger)
    return status1 and status2


def add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, logger):
    """
    Once all family proteins are aligned, and back-translated to nucleotides,
    add missing genomes for the family to the alignment with '-'.
    """
    logger.log(utils.detail_lvl(), "Adding missing genomes for family {}".format(num_fam))
    len_aln, _ = check_lens(btr_file, num_fam, logger)
    with open(miss_file, "r") as missf, open(btr_file, "a") as btrf:
        for genome in missf:
            genome = genome.strip()
            toadd = ">" + genome + "\n" + "-" * len_aln + "\n"
            btrf.write(toadd)
    _, nb = check_lens(btr_file, num_fam, logger)
    if nb != ngenomes:
        logger.error(("ERROR: family {} contains {} in total instead of the {} "
                      "genomes in input.\n").format(num_fam, nb, ngenomes))
        return False
    return True


def family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file,
                     num_fam, ngenomes, logger):
    """
    From a given family, align all its proteins with mafft, back-translate
    to nucleotides, and add missing genomes in this family.

    Returns False if a problem occurred during alignment, checking, back-translation etc.
    Returns True if no problem found
    """
    nbfprt = check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger)
    if not nbfprt:
        return False
    nbfal = mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger)
    if not nbfal:
        return False
    return back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger)


def check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger):
    """
    Check that extractions went well for the given family:
    - check number of proteins and genes extracted compared to the
    number of genomes
    """
    logger.log(utils.detail_lvl(), "Checking extractions for family {}".format(num_fam))

    # Check that extractions went well
    nbmiss = utils.count(miss_file)
    nbfprt = utils.grep(prt_file, "^>", count=True)
    nbfgen = utils.grep(gen_file, "^>", count=True)
    if nbmiss + nbfprt != ngenomes:
        logger.error(("fam {}: wrong sum of missing genomes ({}) and prt "
                      "extracted ({}).").format(num_fam, nbmiss, nbfprt))
        return False
    if nbmiss + nbfgen != ngenomes:
        logger.error(("fam {}: wrong sum of missing genomes ({}) and gen "
                      "extracted ({}).").format(num_fam, nbmiss, nbfgen))
        return False
    return nbfprt


def mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger):
    """
    Align all proteins of the given family with mafft
    """
    logger.log(utils.detail_lvl(), "Aligning family {}".format(num_fam))
    cmd = "fftns --quiet {}".format(prt_file)
    error = "Problem while trying to align fam {}".format(num_fam)
    stdout = open(mafft_file, "w")
    ret = utils.run_cmd(cmd, error, stdout=stdout, logger=logger)
    if ret != 0:
        return False
    return check_mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger)


def check_mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger):
    """
    Check that mafft alignment went well: the number of proteins in the alignment
    is the same as the number of proteins extracted
    """
    nbfal = utils.grep(mafft_file, "^>", count=True)
    if nbfprt != nbfal:
        logger.error("fam {}: different number of proteins extracted in {} ({}) and proteins "
                     "aligned in {} ({}).".format(num_fam, prt_file, nbfprt, mafft_file, nbfal))
        return False
    return nbfal


def back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger):
    """
    Backtranslate protein alignment to nucleotides
    """
    logger.log(utils.detail_lvl(), "Back-translating family {}".format(num_fam))
    curpath = os.path.dirname(os.path.abspath(__file__))
    awk_script = os.path.join(curpath, "prt2codon.awk")
    cmd = "awk -f {} {} {}".format(awk_script, mafft_file, gen_file)
    stdout = open(btr_file, "w")
    error = "Problem while trying to backtranslate {} to a nucleotide alignment".format(mafft_file)
    ret = utils.run_cmd(cmd, error, stdout=stdout, logger=logger)
    if ret != 0:
        return False
    return check_backtranslate(num_fam, mafft_file, btr_file, nbfal, logger)


def check_backtranslate(num_fam, mafft_file, btr_file, nbfal, logger):
    """
    Check back-translation
    """
    nbfalnuc = utils.grep(btr_file, "^>", count=True)
    if nbfal != nbfalnuc:
        logger.error("fam {}: different number of proteins aligned in {} ({}) and genes "
                     "back-translated in {} ({})".format(num_fam, mafft_file, nbfal,
                                                         btr_file, nbfalnuc))
        return False
    return True


def check_lens(aln_file, num_fam, logger):
    """
    In the given alignment file, check that all sequences have the same length.
    If there is no problem, it returns the length of alignment, and the number
    of sequences aligned
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
    return list(all_sums)[0], nb_gen
