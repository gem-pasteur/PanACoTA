#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
Script to prepare the raw sequence files for the pipeline:
    - rename all contigs of a given multifasta file. Contigs are named with the genome
     name (= fasta filename) + the contig number
    - get the total number of contigs in each genome
    - get the total size (in nuc) of each genome

Example for contig names:
ELAN.0316.00001.0001 for the first contig of the first ELAN strain

Input:
    - database folder, containing all sequences in multifasta
    - a file with 2 columns: gembase name and original name of each genome

Output:
    - for each genome, a new multi-fasta in the database folder (with contig names changed)
    - a file LSTINFO-complete.lst with 4 columns: gembase name, original name, nb contigs,
    genome size

@author: GEM, Institut Pasteur

April 2016
"""

import glob
import os
import sys


def main(db_path, lstfile):
    """
    Main method, renaming fasta contigs of all multifasta files found in db_path,
    according to the genome names defined in lstfile.
    """
    corres = read_lstinfo(lstfile)
    # For each genome, read its sequence file, and:
    # get number of contigs, calculate total size, change contig names
    get_all_genomes_info(corres, db_path)
    # Write LSTINFO complete
    write_lstinfo(lstfile + "-complete.lst", corres)


def write_lstinfo(lstfile, corres):
    """
    From info (saved in 'corres') of all genomes, write it to the output file,
    called 'lstfile'.
    """
    corres_gembase = {gem: [ori, cont, size] for ori, [gem, cont, size] in corres.items()}
    with open(lstfile, "w") as lstf:
        lstf.write("gembase_name\tGenome_orig_name\tnb_contigs\tgenome_size\n")
        for gembase in sorted(corres_gembase,
                              key=lambda x: (x.split(".")[0], int(x.split(".")[-1]))):
            to_write = gembase
            for info in corres_gembase[gembase]:
                to_write += "\t" + str(info)
            lstf.write(to_write + "\n")
    print(" * LSTINFO completed")


def get_all_genomes_info(corres, db_path):
    """
    Method reading all genome sequences and getting info (nb contigs, size), and
    creating a new sequence file with the contig names in gembase format.
    """
    print(" * Reading and formatting input sequences...")
    for ori, new in corres.items():
        orifile = os.path.join(db_path, ori)
        add_genome_info(orifile, new)


def add_genome_info(orifile, new):
    """
    For a given genome (sequence in 'orifile', read its multifasta sequence,
    change contig names using 'new', and get its number of contigs and its total size.
    """
    if os.path.isfile(orifile):
        cleanfile = orifile + "-gembase.fna"
        nbcont = 0
        gensize = 0
        # print "   ->", os.path.basename(orifile)  #, "\033[K\r",
        # sys.stdout.flush()
        with open(orifile, "r") as orif:
            with open(cleanfile, "w") as ouf:
                for line in orif:
                    if line.startswith(">"):
                        nbcont += 1
                        res = ">" + new[0] + "." + str(nbcont).zfill(4) + "\n"
                        ouf.write(res)
                    else:
                        gensize += len(line.strip())
                        ouf.write(line)
        new += [nbcont, gensize]
    else:
        new += ['NA', 'NA']


def read_lstinfo(lstfile):
    """
    Read lstinfo file, to have the corresponding original and new names.
    """
    corres = {}
    with open(lstfile, "r") as lsf:
        for line in lsf:
            if "gembase" not in line and line != "\n":
                gembase = line.split()[0]
                ori = line.split()[1]
                corres[ori] = [gembase]
    return corres


def parse():
    """
    Method to create a parser for command-line options
    """
    import argparse
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-d", "--db", dest="dbpath",
                        help=("Path to the folder containing all multifasta files whose "
                              "contigs must be renamed."),
                        required=True)
    parser.add_argument("-l", dest="lstfile",
                        help=("File containing in a first column the new genome name (name which "
                              "will be given to the contigs), and in the second column the "
                              "original genome name (multifasta filename)."), required=True)
    return parser.parse_args()


if __name__ == '__main__':
    OPTIONS = parse()
    main(OPTIONS.dbpath, OPTIONS.lstfile)
