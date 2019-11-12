#!/usr/bin/env python3
# coding: utf-8

import sys
import os
import logging
import progressbar

from PanACoTA import utils

logger = logging.getLogger("align.extract")


def get_all_seqs(all_genomes, dname, dbpath, listdir, aldir, all_fams, quiet):
    """
    For all genomes, extract its proteins present in a persistent family to the file
    corresponding to this family.

    Parameters
    ----------
    all_genomes : []
        list of all genome names
    dname : str
        name of dataset
    dbpath : str
        path to folder containing 'Proteins' and 'Genes' folders
    listdir : str
        path to folder containing the lists of proteins/genes to extract
    aldir : str
        path to folder where extracted proteins/genes must be saved
    all_fams : []
        list of all family numbers
    quiet : bool
        True if nothin must be written to stdout/stderr, False otherwise
    """
    # Get list of files not already existing
    files_todo = check_existing_extract(all_fams, aldir, dname)
    if len(files_todo) == 0:
        logger.info(("All extraction files already existing (see detailed log for "
                     "more information)"))
        logger.warning("All prt and gene files for all families already exist. The program "
                       "will use them for the next step. If you want to re-extract a given "
                       "family, remove its prt and gen extraction files. If you want to "
                       "re-extract all families, use option -F (or --force).")
        return

    logger.info("Extracting proteins and genes from all genomes")
    nbgen = len(all_genomes)
    bar = None
    curnum = 1
    if not quiet:
        widgets = ['Gene and Protein extraction:',
                   progressbar.Bar(marker='█', left='', right='', fill=' '),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ") - ", progressbar.Timer(), ' ',
                   progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=150).start()
        curnum = 1
    for genome in all_genomes:
        ge_gen = os.path.join(listdir, dname + "-getEntry_gen_" + genome + ".txt")
        ge_prt = os.path.join(listdir, dname + "-getEntry_prt_" + genome + ".txt")
        logger.details("Extracting proteins and genes from {}".format(genome))
        prtdb = os.path.join(dbpath, "Proteins", genome + ".prt")
        gendb = os.path.join(dbpath, "Genes", genome + ".gen")
        get_genome_seqs(prtdb, ge_prt, files_todo)
        get_genome_seqs(gendb, ge_gen, files_todo)
        if not quiet:
            bar.update(curnum)
            curnum += 1
    if not quiet:
        bar.finish()


def check_existing_extract(all_fams, aldir, dname):
    """
    For each family, check if its prt and gen extraction file already exist.
    If both exist, no need to re-extract for those families.
    If only one or no one exists, put to list to extract.

    Parameters
    ----------
    all_fams : list
        list of all family numbers
    aldir : str
        path to directory where extraction files must be saved
    dname : str
        name of the dataset

    Returns
    -------
    []
        list of files that must be generated (prt and gen files)
    """
    extract_fams = []
    for fam in all_fams:
        genfile = os.path.join(aldir, "{}-current.{}.gen".format(dname, fam))
        prtfile = os.path.join(aldir, "{}-current.{}.prt".format(dname, fam))
        if not os.path.isfile(genfile) or not os.path.isfile(prtfile):
            # At least 1 file missing: re-extract all proteins and all genes
            utils.remove(genfile)
            utils.remove(prtfile)
            # As we re-extract proteins and genes, redo alignments
            mafft_file = os.path.join(aldir, "{}-mafft-align.{}.aln".format(dname, fam))
            btr_file = os.path.join(aldir, "{}-mafft-prt2nuc.{}.aln".format(dname, fam))
            utils.remove(mafft_file)
            utils.remove(btr_file)
            extract_fams.append(genfile)
            extract_fams.append(prtfile)
    # If we re-extract at least 1 family, redo the final files (concatenation and group by genome)
    if len(extract_fams) > 0:
        concat = os.path.join(aldir, "{}-complete.cat.aln".format(dname))
        outdir = os.path.split(aldir)[0]
        treedir = os.path.join(outdir, "Phylo-" + dname)
        grpfile = os.path.join(treedir, dname + ".grp.aln")
        utils.remove(concat)
        utils.remove(grpfile)
    return extract_fams


def get_genome_seqs(fasta, tabfile, files_todo, outfile=None):
    """
    From a fasta file, extract all sequences given in the tab file.
    The tab file can contain:

    - 1 sequence name per line -> all sequences will be extracted to the same file
    - 1 sequence + 1 filename per line -> each sequence will be extracted in the given file

    If outfile not given, the tab file must contain 2 columns (1 for the sequence name,
    1 for its output file). If an outfile is given, only the 1st column of tab file
    will be considered, and all sequences will be extracted to the given outfile.

    Parameters
    ----------
    fasta : str
        path to fasta file from which sequences must be extracted
    tabfile : str
        path to the tab file containing the names of sequences to extract
    files_todo : list
        list of files which must be generated (prt and gen files). Others
        already exist, so ignore them.
    outfile : str or None
        if None, the tab file must contain 2 columns (1 for the sequence name,
        1 for its output file). If an outfile is given (not None), only the 1st column of tab file
        will be considered, and all sequences will be extracted to the given outfile.
    """
    with open(tabfile, "r") as tabf:
        to_extract = get_names_to_extract(tabf, outfile)

    if outfile:
        if os.path.isfile(outfile):
            logger.warning("Sequences are already extracted in {}. This will "
                           "be used for next step. If you want "
                           "to re-extract all sequences, use option -F (or "
                           "--force)".format(outfile))
            return
        with open(fasta, "r") as fasf, open(outfile, "a") as outf:
            extract_sequences(to_extract, fasf, outf=outf)
    else:
        with open(fasta, "r") as fasf:
            extract_sequences(to_extract, fasf, files_todo=files_todo)


def get_names_to_extract(tabf, outfile):
    """
    From the tab file, get names of sequences to extract.

    Parameters
    ----------
    tabf : _io.TextIO
        open file containing names of sequences to extract
    outfile : str or None
        if None, the tab file must contain 2 columns (1 for the sequence name,
        1 for its output file). If an outfile is given (not None), only the 1st column of tab file
        will be considered, and all sequences will be given 'outfile' as output file

    Returns
    -------
    dict
        {sequence_to_extract: file_to_which_it_will_be_extracted}
    """
    to_extract = {}
    for line in tabf:
        if outfile:
            seq = line.split()[0].strip()
            out = outfile
        else:
            try:
                seq, out = line.split()[:2]
            except ValueError:
                logger.error(("Your file {} does not contain an output filename for {}. "
                              "Please give an output filename for each sequence to extract, "
                              "or give a general output filename where all sequences will "
                              "be extracted.").format(tabf.name, line.strip()))
                sys.exit(1)
        to_extract[seq] = out
    return to_extract


def extract_sequences(to_extract, fasf, files_todo=None, outf=None):
    """
    Extract sequences from an open fasta file 'fasf', and a list of sequences to
    extract

    Parameters
    ----------
    to_extract : dict or []
        {sequence_to_extract: file_to_which_it_will_be_extracted} or list of sequences to
        extract (if an output filename is also given in 'outf')
    fasf : _io.TextIO
        open file containing sequences in multi-fasta format
    files_todo : list or None
        list of files which must be generated (prt and gen files). Others
        already exist, so ignore them.
    outf : _io.TextIO or None
        If an outfile is given (not None), and 'to_extract' is a dict, only its keys will be
        considered, and all these sequences will be extracted to 'outfile' (if 'to_extract' is a
        list, will extract all sequences of this list). Otherwise, if None,
        each sequence will be extracted to its corresponding value in 'to_extract'.
    """
    if files_todo is None:
        files_todo = []
    out_given = (outf is not None)
    extract = False
    for line in fasf:
        if line.startswith(">"):
            seq = line.split(">")[1].split()[0]
            # Seq is part of sequences to extract
            if seq in to_extract:
                # if out_given, outf is already open. Otherwise, open it only if must be written
                if not out_given:
                    # get corresponding outfile
                    out = to_extract[seq]
                    if out in files_todo:
                        outf = open(out, "a")
                    else:
                        outf = None
                if outf:
                    outf.write(line)
                    extract = True
                else:
                    extract = False
            else:
                if not out_given and outf:
                    outf.close()
                extract = False
        else:
            if extract:
                outf.write(line)
