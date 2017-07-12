#!/usr/bin/env python3
# coding: utf-8

"""
Concatenate all alignment files of all families. Then,
group alignments by genome.

@author: GEM, Institut Pasteur
March 2017
"""

import os
import sys
import logging
from genomeAPCAT import utils

logger = logging.getLogger("align.post")


def post_alignment(fam_nums, all_genomes, prefix, outdir, dname, quiet):
    """
    After the alignment of all proteins by family:
    - concatenate all alignment files
    - group the alignment by genome
    """
    all_alns = concat_alignments(fam_nums, prefix, quiet)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(treedir, exist_ok=True)
    group_by_genome(all_genomes, all_alns, treedir, dname, quiet)


def concat_alignments(fam_nums, prefix, quiet):
    """
    Concatenate all family alignment files to a unique file
    """
    logger.info("Concatenating all alignment files")
    list_files = ["{}-mafft-prt2nuc.{}.aln".format(prefix, num_fam) for num_fam in fam_nums]
    output = "{}-complete.cat.aln".format(prefix)
    if quiet:
        utils.cat(list_files, output)
    else:
        utils.cat(list_files, output, title = "Concatenation")
    return output


def group_by_genome(all_genomes, all_alns, treedir, dname, quiet):
    """
    From the alignment file 'all_alns' containing all proteins, group the alignments of
    proteins by their genome (listed in 'all_genomes'), and save the result
    in 'treedir'
    """
    logger.info("Grouping alignments per genome")
    outfile = os.path.join(treedir, dname + ".grp.aln")
    sequences = read_alignments(all_alns, all_genomes)
    write_groups(outfile, sequences)


def read_alignments(all_alns, all_genomes):
    """
    Read alignment file, and assign each sequence to a genome
    """
    sequences = {}  # name: [list of sequences]
    genome = None
    seq = ""
    with open(all_alns, 'r') as alnf:
        for line in alnf:
            if line.startswith(">"):
                # If new header, write previous protein name/sequence to 'sequences'
                if genome and seq:
                    sequences[genome].append(seq)
                    seq = ""
                # Get new genome header
                genome = get_genome(line, all_genomes)
                if genome not in sequences:
                    sequences[genome] = []
            else:
                seq += line.strip()
    if genome and seq:
        sequences[genome].append(seq)
    per_genome = [len(seq) for seq in sequences.values()]
    if len(set(per_genome)) != 1:
        logger.error("Problems occurred while grouping alignments by genome: all genomes "
                     "do not have the same number of sequences. Check that each protein "
                     "name contains the name of the genome from which it comes.")
        sys.exit(1)
    logger.details("{} sequences found per genome".format(per_genome[0]))
    return sequences


def write_groups(outfile, sequences):
    """
    Writing alignments per genome to output file.
    """
    logger.details("Writing alignments per genome")
    with open(outfile, "w") as outf:
        for genome in sorted(sequences, key=utils.sort_genomes):
            # write header for genome
            outf.write(">" + genome + "\n")
            # Write all sequences
            outf.write("".join(sequences[genome]) + "\n")


def get_genome(header, all_genomes):
    """
    Find to which genome belongs 'header'
    """
    header = header.split(">")[1].split()[0]
    for genome in all_genomes:
        if genome in header:
            return genome
    logger.error("Protein {} does not correspond to any genome name given... {}".format(header))
    sys.exit(1)
