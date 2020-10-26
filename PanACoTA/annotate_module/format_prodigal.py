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
Functions to convert prodigal result files to gembase format.

    * Proteins: multifasta with all CDS in aa
    * Replicons: (multi)fasta of genome sequences
    * Genes: multifasta of all genes in nuc
    * gff3: gff files without sequence
    * LSTINFO: information on annotation. Columns are: "start end strand type locus
    gene_name | product | EC_number | inference2" and strain is C (complement) or D (direct).
    Locus is: `<genome_name>.<contig_num><i or b>_<protein_num>`
    For prodigal: "start end strand type locus NA | NA | NA | NA", as there is no
    functional annotation.

@author gem
July 2019
"""

import os
import shutil
import glob
import logging

import PanACoTA.utils as utils
from PanACoTA.annotate_module import general_format_functions as gfunc

logger = logging.getLogger("annotate.prodigal_format")


def format_one_genome(gpath, name, prod_path, lst_dir, prot_dir, gene_dir,
                      rep_dir, gff_dir):
    """
    Format the given genome, and create its corresponding files in the following folders:

    - Proteins
    - Genes
    - Replicons
    - LSTINFO
    - gff

    Parameters
    ----------
    gpath : str
        path to the genome sequence which was given to prodigal for annotation
    name : str
        gembase name of the genome
    prod_path : str
        directory where all tmp_files for all sequences are saved (sequences cut at each set
        of 5N, prodigal results and logs)
    lst_dir : str
        path to LSTINFO folder
    prot_dir : str
        path to Proteins folder
    gene_dir : str
        path to Genes folder
    rep_dir : str
        path to Replicons folder
    gff_dir : str
        path to gff3 folder

    Returns
    -------
    bool :
        True if genome was correctly formatted, False otherwise
    """
    # Get directory where prodigal results for the current genome are saved
    prodigal_dir = os.path.join(prod_path, os.path.basename(gpath) + "-prodigalRes")

    # Get prodigal result files
    prot_file = glob.glob(os.path.join(prodigal_dir, "*.faa"))[0]
    gen_file = glob.glob(os.path.join(prodigal_dir, "*.ffn"))[0]
    gff_file = glob.glob(os.path.join(prodigal_dir, "*.gff"))[0]

    # Define names for generated gembase files
    res_prot_file = os.path.join(prot_dir, name + ".prt")
    res_gene_file = os.path.join(gene_dir, name + ".gen")
    res_lst_file = os.path.join(lst_dir, name + ".lst")
    res_rep_file = os.path.join(rep_dir, name + ".fna")
    res_gff_file = os.path.join(gff_dir, name + ".gff")

    # Generate replicon file (same as input sequence but with gembase formatted headers). From
    # this file, get contig names, to be used to generate gff file
    contigs, sizes = utils.get_genome_contigs_and_rename(name, gpath, res_rep_file, logger)
    if not contigs:
        try:
            os.remove(res_rep_file)
            os.remove(res_lst_file)
            os.remove(res_gff_file)
            os.remove(res_gene_file)
            os.remove(res_prot_file)
        except OSError:
            pass
        logger.error("Problems while generating Replicon file for {}".format(name))
        return False

    # First, create .gen and .lst files. If they could not be formatted,
    # remove those files, and return False with error message
    ok = create_gene_lst(contigs, gen_file, res_gene_file, res_lst_file, gpath, name)
    if not ok:
        try:
            os.remove(res_rep_file)
            os.remove(res_gene_file)
            os.remove(res_lst_file)
            os.remove(res_gff_file)
            os.remove(res_prot_file)
        except OSError:
            pass
        logger.error(f"Problems while generating .gen and .lst files for {name}")
        return False

    # Create gff files.
    ok = create_gff(gpath, gff_file, res_gff_file, res_lst_file, contigs, sizes)
    # If problem while formatting the genome (rep or gff file), remove all
    # already created files, and return False (genome not formatted) with error message.
    if not ok:
        try:
            os.remove(res_gene_file)
            os.remove(res_lst_file)
            os.remove(res_rep_file)
            os.remove(res_gff_file)
            os.remove(res_prot_file)
        except OSError:
            pass
        logger.error("Problems while generating .gff (gff3 folder) "
                     "file for {}".format(name))
        return False

    # Generate .prt files (in Proteins directory)
    ok = create_prt(prot_file, res_prot_file, res_lst_file)
    # If problem while formatting prt file, return False, delete all generated
    # formatted files, and write an error message to user.
    if not ok:
        try:
            os.remove(res_gene_file)
            os.remove(res_lst_file)
            os.remove(res_rep_file)
            os.remove(res_gff_file)
            os.remove(res_prot_file)
        except OSError: # pragma: no cover
            pass
        logger.error("Problems while generating .prt file (Proteins folder) "
                     "for {}".format(name))
        return False
    return ok


def create_gene_lst(contigs, gen_file, res_gen_file, res_lst_file, gpath, name):
    """
    Generate .gen file, from sequences contained in .ffn, but changing the
    headers to match with gembase format.
    At the same time, generate .lst file, from the information given in prodigal ffn headers

    Parameters
    ----------
    contigs : dict
        {original_contig_name: gembase_contig_name}
    gen_file : str
        .ffn file generated by prodigal
    res_gen_file : str
        generated .gen file, to write in Genes directory
    res_lst_file : str
        generated .lst file to write in LSTINFO directory
    gpath : str
        path to the genome given to prodigal. Only used for error message
    name : str
        gembase name of the genome to format
    logger : logging.Logger
        logger object to put information

    Returns
    -------
    bool :
        True if conversion went well, False otherwise
    """
    # Variable which will contain the current gene sequence
    seq = ""
    # number of the current gene (first gene is 1, 2nd is 2 etc. each number is unique: do not
    # re-start at 1 for each new contig)
    locus_num = 0
    # contig name of the last gene. To check if we are now in a new contig (-> loc = b) or not
    prev_cont_name = ""
    # Previous ontig number: contig number to use in gembase format
    prev_cont_num = 0
    contig_num = 0
    # Keep start, end, strand and informations (prodigal gives information on start_type,
    # gc_cont etc.) from the previous gene, before overwriting it with information
    # on the new one
    prev_start = ""
    prev_end = ""
    prev_strand = ""
    prev_info = ""
    # Update loc when contig changes ('b' if gene at the border of a contig, 'i' if it is inside)
    prev_loc = "b"
    # To start, the first gene is, by definition, at the border of the contig
    loc = "b"
    # Open files: .ffn prodigal to read, .gen and .lst gembase to create
    with open(gen_file, "r") as ffn, open(res_gen_file, "w") as r_gen,\
         open(res_lst_file, "w") as r_lst:
        # Read all lines in ffn file (sequences in nuc. for each gene)
        for lineffn in ffn:
            # If it is a sequence, save it and go to next line
            if not lineffn.startswith(">"):
                seq += lineffn
                continue
            # Otherwise:
            # - write header of previous sequence to .gen
            # - write previous sequence (in 'seq') to .gen
            # - write LSTINFO information to .lst
            # - update information (new start, end, contig number etc.) for next gene
            else:
                # Get information given for the new gene (by .ffn file from prodigal)
                (gname, start, end, strand, info) = lineffn.strip().split(">")[-1].split("#")
                # Get contig number from prodigal gene header: prodigal first part of header is:
                #  <original genome name contig name>_<protein number>
                contig_name = gname.strip().split("_")
                if len(contig_name) > 1:
                    contig_name = "_".join(contig_name[:-1])
                else:
                    contig_name = contig_name[0]
                # If new contig:
                # - previous gene was the last of its contig -> prev_loc = "b" ;
                # - the current gene is the first of its contig (loc = "b")
                # - we must increment the contig number
                if contig_name != prev_cont_name:
                    # Check that this contig name is in the list, and get its gembase contig number
                    if contig_name in contigs:
                        contig_num = contigs[contig_name].split(".")[-1]
                    # if not in the list, problem, return false
                    else:
                        logger.error(f"'{contig_name}' found in {gen_file} does not exist in "
                                     f"{gpath}.")
                        return False
                    prev_loc = 'b'
                    loc = 'b'
                # If not new contig. If prev_loc == 'b', previous gene is the first protein
                # of this contig.
                # Current gene will be inside the contig (except if new contig for the next gene,
                # meaning only 1 gene in the contig)
                else:
                    loc = 'i'

                # If it is not the first gene of the genome, write previous gene information
                if prev_start != "":
                    # Write line in LSTINFO file, + header and sequence to the gene file
                    lstline = gfunc.write_gene("CDS", locus_num, "NA", "NA",
                                               prev_loc, name, prev_cont_num, "NA", prev_info,
                                               "NA", prev_strand, prev_start, prev_end, r_lst)
                    gfunc.write_header(lstline, r_gen)
                    r_gen.write(seq)
                # -> get new information, save it for the next gene, and go to next line
                # Strands are 1/-1 in prodigal, while we use D,C -> convert, so that next time
                # we find a new gene, it writes this before updating for this new gene
                if int(strand) == 1:
                    strand = "D"
                else:
                    strand = "C"
                # Prepare variables for next gene
                locus_num += 1
                seq = ""
                prev_cont_num = contig_num
                prev_cont_name = contig_name
                prev_start = start
                prev_end = end
                prev_strand = strand
                prev_loc = loc
                prev_info = info
        # Write last gene of the genome (-> loc = 'b'),
        # Just check that there was at least 1 gene found (prev_start != "").
        # Otherwise, nothing to write
        if prev_start != "":
            prev_loc = "b"
            lstline = gfunc.write_gene("CDS", locus_num, "NA", "NA",
                                       prev_loc, name, prev_cont_num, "NA", prev_info, "NA",
                                       prev_strand, prev_start, prev_end, r_lst)
            gfunc.write_header(lstline, r_gen)
            r_gen.write(seq)
    return True


def create_gff(gpath, gff_file, res_gff_file, res_lst_file, contigs, sizes):
    """
    Create .gff3 file.

    Format:

    ##gff-version 3
    ##sequence-region contig1 start end
    ##sequence-region contig2 start end
    ...
    seqid(=contig) source   type start   end score  strand phase attributes

    All fields tab separated.
    Empty fields contain '.'

    For example:
    ESCO.1017.00200.00001    Prodigal:2.6    CDS start  end .   +   .   ID=ESCO.1017.00200.b0001_00001;locus_tag=ESCO.1017.00200.b0001_00001;product=hypothetical protein


    genome1_1   Prodigal_v2.6.3 CDS 213 1880    260.0   +   0   ID=1_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.534;conf=99.99;score=259.99;cscore=268.89;sscore=-8.89;rscore=-8.50;uscore=-4.34;tscore=3.95;

    Parameters
    ----------
        gpath : str
            path to the genome sequence given to prodigal. Only used for the error message
        gff_file : str
            path to gff file generated by prodigal
        res-gff_file : str
            path to the gff file that must be created in result database
        res-lst_file : str
            path to the lst file that was created in result database in the previous step
        contigs : dict
            dict of contig names with their size. ["original_name": "gembase_name"]
        sizes : dict
            dict of contig gembase names with their sizes {"gembase_name": size}

    Returns
    -------
    bool
        True if everything went well, False if any problem

    """
    # open gff generated by prodigal to read it
    # open file to write new gff file
    # open lst file to read all information saved from prodigal results
    with open(gff_file, 'r') as gf, open(res_gff_file, "w") as rgf, open(res_lst_file, "r") as rlf:
        # Write headers of gff3 file
        rgf.write("##gff-version  3\n")
        for ori_name, new_name in contigs.items():
            # Write the list of contigs, with their size
            end = sizes[new_name]
            rgf.write(f"##sequence-region\t{new_name}\t1\t{end}\n")

        # Now, convert each line of prodigal gff to gembase formatted gff line
        for linegff in gf:
            # Ignore gff lines starting by #. For the new gff, these are already all written at the
            # beginning of the file.
            if linegff.startswith("#"):
                continue

            # We do not write the sequences
            if linegff.startswith(">"):
                break

            # Get all information from prodigal gff line. Strip each of them as trailing whitespace
            # can be hidden (leading to information from gff considered as different from
            # information from lst)
            fields_g = linegff.split("\t")
            fields_g = [info.strip() for info in fields_g]
            (contig_name, source, type_g, start_g, end_g,
             score, strand_g, phase, attributes) = fields_g
            # Get information given to this same sequence from the lst file
            # (next lst line corresponds to next gff line without #), as, for each format,
            # there is 1 line per gene)
            linelst = rlf.readline()
            fields_l = linelst.split("\t")
            fields_l = [info.strip() for info in fields_l]
            start_l, end_l, strand_l, type_l, locus_l, _, _ = fields_l


            # Get gff and ffn filenames to give information to user if error message
            gff = os.path.basename(gff_file)
            ffn = ".".join(gff.split(".")[:-1]) + ".ffn"
            # Path where gff and ffn generated by prodigal are
            tmp = gpath + "-prodigalRes"
            # Get gene name given by prodigal to current gene
            gname = attributes.split("ID=")[1].split(";")[0]

            # Compare information from lst and information from prodigal gff (start,
            # end and type of feature). They should correspond
            for (elemg, eleml, label) in zip([start_g, end_g, type_g],
                                             [start_l, end_l, type_l],
                                             ["start", "end", "type"]):
                # If 1 element is different (start or end position, or type), print error
                # message and return False: this genome could not be converted to gff
                if elemg != eleml:
                    logger.error(f"Files {ffn} and {gff} (in prodigal tmp_files: {tmp}) do not have "
                                 f"the same {label} value for gene {gname} ({elemg} in gff, {eleml} in "
                                 "ffn))")
                    return False

            # Replace prodigal ID with the new gene name (in gembase format), found in the
            # corresponding lst line
            at_split = attributes.split(";")
            new = [atr if "ID" not in atr else f'ID={locus_l}' for atr in at_split]

            # Write new line of gff file
            # Write new contig name
            cname = contigs[contig_name]
            info = "\t".join([cname, source, type_g, start_g, end_g, score, strand_g,
                              phase, ";".join(new)])
            rgf.write(info + "\n")
    return True


def create_prt(prot_file, res_prot_file, res_lst_file):
    """
    Generate .prt file (gembase formatted gene names), from features contained in .lst file generated just before.

    Parameters
    ----------
    prot_file : str
        .faa file generated by prodigal
    res_prot_file : str
        output file, to write in Proteins directory
    res_lst_file : str
        .lst file to get all gene names in gembase format instead of re-generating them
    Returns
    -------
    bool :
        True if conversion went well, False otherwise
    """

    # Open:
    # - prot file to read gene sequences from prodigal results
    # - res_prot file to write sequences with gembase headers
    # - res_lst_file to get gene gembase names and other infos (strand, size...)

    with open(prot_file, "r") as faa, open(res_prot_file, "w") as r_prt,\
         open(res_lst_file, "r") as r_lst:
         # Read prt file generated by prodigal
        for lineprot in faa:
            # If protein sequence, write it
            if not lineprot.startswith(">"):
                r_prt.write(lineprot)
                continue
            # If header, replace by gembase header
            # For that, get next lst line (corresponding to next protein,
            # as there is 1 protein per line in .lst -> 1 protein per header in .prt)
            linelst = r_lst.readline().strip()
            # Try to get info from lstline.
            # If lstline empty, it means that the current protein
            # is missing from lst file. We already read the last protein of lst file.
            if linelst != '':
                # If not empty, check lst format, return False if not right format
                try:
                    # If ok, gembase name is in the fifth column of lst file
                    start, end, _, _, gem_name, _, _ = linelst.split("\t")
                except ValueError:
                    logger.error("Problem in format of lstline ({})".format(linelst))
                    return False
            else:
                logger.error("No more protein in lst file. We cannot get information on this "
                             "protein ({})! Check that you do not have more proteins than genes "
                             "in prodigal results".format(lineprot.strip()))
                return False
            # Write this gembase name as a new header
            # Size of protein sequence is the third of gene sequence. Check that it is an int.
            try :
                size_gen = (int(end) - int(start) + 1)
            except ValueError:
                logger.error("Start and/or end of protein {} position is not a number (start "
                             "= {}; end = {})".format(gem_name, start, end))
                return False
            # Find size of protein in number of aa
            # If number of nucleotides in gene cannot be divided by 3 to get number of amino-acids, there is a problem with this protein: return False to ignore this genome
            size_prot = size_gen/3
            if int(size_prot) != size_prot:
                logger.error("Gene {} has a number of nucleotides ({}) that is not divisible "
                             "by 3.".format(gem_name, size_gen))
                return False
            gfunc.write_header(linelst, r_prt)
            # new_header = "\t".join([gem_name, str(int(size_prot)), product, info])
            # r_prt.write(">" + new_header + "\n")
        # Check that there are no more proteins in lst than in this prt file
        linelst = r_lst.readline()
        if linelst.strip() != '':
            gem_name = linelst.strip().split("\t")[4]
            logger.error("Protein {} is in .lst file but its sequence is not in the protein "
                         "file generated by prodigal.".format(gem_name))
            return False
    return True
