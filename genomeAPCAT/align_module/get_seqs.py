#!/usr/bin/env python3
# coding: utf-8

import sys
import os
import logging
import progressbar

logger = logging.getLogger("align.extract")


def get_all_seqs(all_genomes, dname, dbpath, listdir, force, quiet):
    """
    For all genomes, extract its proteins present in a persistent family to the file
    corresponding to this family.
    """
    logger.info("Extracting proteins and genes from all genomes")
    nbgen = len(all_genomes)
    if not quiet:
        widgets = ['Gene and Protein extraction:',
                   progressbar.Bar(marker='█', left='', right='', fill=' '),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ") - ", progressbar.Timer(), ' ',
                   progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=150).start()
        curnum = 1
    files_explored = []
    for genome in all_genomes:
        ge_gen = os.path.join(listdir, dname + "-getEntry_gen_" + genome + ".txt")
        ge_prt = os.path.join(listdir, dname + "-getEntry_prt_" + genome + ".txt")
        logger.details("Extracting proteins and genes from {}".format(genome))
        prtdb = os.path.join(dbpath, "Proteins", genome + ".prt")
        gendb = os.path.join(dbpath, "Genes", genome + ".gen")
        get_genome_seqs(prtdb, ge_prt, files_explored, force=force)
        get_genome_seqs(gendb, ge_gen, files_explored, force=force)
        if not quiet:
            bar.update(curnum)
            curnum += 1
    if not quiet:
        bar.finish()


def get_genome_seqs(fasta, tabfile, files_explored, outfile=None, force=False):
    """
    From a fasta file, extract all sequences given in the tab file.
    The tab file can contain:
    - 1 sequence name per line -> all sequences will be extracted to the same file
    - 1 sequence + 1 filename per line -> each sequence will be extracted in the given file

    If outfile not given, the tab file must contain 2 columns (1 for the sequence name,
    1 for its output file). If an outfile is given, only the 1st column of tab file
    will be considered, and all sequences will be extracted to the given outfile.
    """
    with open(tabfile, "r") as tabf:
        to_extract = get_names_to_extract(tabf, outfile)

    if outfile:
        if outfile in files_explored:
            mode = "a"
        elif force:
            mode = "w"
            files_explored.append(outfile)
        elif os.path.isfile(outfile):
            logger.warning("Sequences are already extracted in {}. This will "
                           "be used for next step. If you want "
                           "to re-extract all sequences, use option -F (or "
                           "--force)".format(outfile))
            return
        with open(fasta, "r") as fasf, open(outfile, mode) as outf:
            extract_sequences(to_extract, fasf, outf=outf)
    else:
        with open(fasta, "r") as fasf:
            extract_sequences(to_extract, fasf, force=force, explored_files=files_explored)


def get_names_to_extract(tabf, outfile):
    """
    From the tab file, get names of sequences to extract.
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
                              "be extracted.").format(tabfile, line.strip()))
                sys.exit(1)
        to_extract[seq] = out
    return to_extract


def extract_sequences(to_extract, fasf, force=False, explored_files=[], outf=None):
    """
    Extract sequences from an open fasta file 'fasf', and a list of sequences to
    extract

    my_files = files already seen and for which we must append new sequences, corresponding
    to the same kind of data as fasf (protein or gene)
    files_other = same as my_files, but if fasf is proteins, the file has been seen as gene,
    and if fasf is gene, the file has been seen as protein.
    """
    out_given = outf != None
    change_ext = {".prt": ".gen", ".gen": ".prt"}
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
                    out_path, ext = os.path.splitext(out)
                    out_base = os.path.basename(out_path)
                    if out_base in explored_files:
                        # file already seen before: add sequence to it
                        outf = open(out, "a")
                    elif (os.path.isfile(out) and os.path.isfile(out_path + change_ext[ext])
                          and not force):
                        # file never seen before, not force, and the 2 files already exist:
                        # keep them
                        logger.warning("Sequences are already extracted in {} and {}. This will "
                                       "be used for next step. If you want to re-extract only "
                                       "those sequences, erase the 2 files before. If you want "
                                       "to re-extract all sequences, use option -F (or "
                                       "--force)".format(out, out_path + change_ext[ext]))
                        outf = None
                    else:
                        # file not seen before. overwrite them (either force, or not force but
                        # at least one file does not exist)
                        open(out_path + change_ext[ext], "w").close()
                        outf = open(out, "w")
                        explored_files.append(out_base)
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


if __name__ == "__main__":
    if len(sys.argv) == 3:
        for _ in range(200):
            get_genome_seqs(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        get_seqs(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("please give fasta, tabfile, [outfile]")


