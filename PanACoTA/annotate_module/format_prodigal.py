#!/usr/bin/env python3
# coding: utf-8

"""
Functions to convert prodigal result files to gembase format.

    * Proteins: multifasta with all CDS in aa
    * Replicons: multifasta of genome
    * Genes: multifasta of all genes in nuc
    * gff3: gff files without sequence
    * LSTINFO: information on annotation. Columns are: "start end strand type locus
    gene_name | product | EC_number | inference2" and strain is C (complement) or D (direct). Locus is: `<genome_name>.<contig_num><i or b>_<protein_num>`

@author gem
April 2019
"""

def format_one_genome(gpath, name, prok_path, lst_dir, prot_dir, gene_dir,
                      rep_dir, gff_dir, logger):
    """
     Format the given genome, and create its corresponding files in the following folders:

    - Proteins
    - Genes
    - Replicons
    - LSTINFO

    Parameters
    ----------
    gpath : str
        path to the genome sequence which was given to prokka for annotation
    name : str
        gembase name of the genome
    prok_path : str
        directory where prokka folders are saved
    lst_dir : srt
        path to LSTINFO folder
    prot_dir : str
        path to Proteins folder
    gene_dir : str
        path to Genes folder
    rep_dir : str
        path to Replicons folder
    gff_dir : str
        path to gff3 folder
    logger : logging.Logger
        logger object to write log information

    Returns
    -------
    bool :
        True if genome was correctly formatted, False otherwise
    """

    logger.info("Format prodigal results")