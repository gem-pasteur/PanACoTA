#!/usr/bin/env python3
# coding: utf-8

import sys
import logging

logger = logging.getLogger("pangnome.mmseqs")


def get_seqs(fasta, tabfile, outfile=None):
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
        with open(fasta, "r") as fasf, open(outfile, "w") as outf:
            extract_sequences(to_extract, fasf, outf)
    else:
        with open(fasta, "r") as fasf:
            extract_sequences(to_extract, fasf)


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


def extract_sequences(to_extract, fasf, outf=None):
    """
    Extract sequences from an open fasta file 'fasf', and a list of sequences to
    extract
    """
    out_given = outf != None
    extract = False
    for line in fasf:
        if line.startswith(">"):
            seq = line.split(">")[1].split()[0]
            if seq in to_extract:
                if not out_given:
                    out = to_extract[seq]
                    outf = open(out, "w")
                outf.write(line.strip() + '\n')
                extract = True
            else:
                if not out_given and outf:
                    outf.close()
                extract = False
        else:
            if extract:
                outf.write(line)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        get_seqs(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        get_seqs(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("please give fasta, tabfile, [outfile]")


