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
Functions to:

- convert prokka tbl file to our tab file
- convert prokka ffn and faa headers to our format
- Create the database, with the following folders in the given "res_path":

    * Proteins: multifasta with all CDS in aa
    * Replicons: multifasta of genome
    * Genes: multifasta of all genes in nuc
    * gff3: gff files without sequence
    * LSTINFO: information on annotation. Columns are: "start end strand type locus\
    gene_name | product | EC_number | inference2" with the same types as prokka file,\
    and strain is C (complement) or D (direct). Locus is:\
    `<genome_name>.<contig_num><i or b>_<protein_num>`

@author gem
April 2019
"""

import os
import glob
import sys
import shutil
import logging
import PanACoTA.utils as utils
from PanACoTA.annotate_module import general_format_functions as general

logger = logging.getLogger("annotate.prokka_format")


def format_one_genome(gpath, name, prok_path, lst_dir, prot_dir, gene_dir,
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

    Returns
    -------
    bool :
        True if genome was correctly formatted, False otherwise
    """
    prokka_dir = os.path.join(prok_path, os.path.basename(gpath) + "-prokkaRes")
    # Get needed Prokka result files
    fna_file = glob.glob(os.path.join(prokka_dir, "*.fna"))[0]
    prokka_tbl_file = glob.glob(os.path.join(prokka_dir, "*.tbl"))[0]
    prokka_gff_file = glob.glob(os.path.join(prokka_dir, "*.gff"))[0]
    prokka_ffn_file = glob.glob(os.path.join(prokka_dir, "*.ffn"))[0]
    prokka_faa_file = glob.glob(os.path.join(prokka_dir, "*.faa"))[0]

    #  Define names for generated gembase files
    res_lst_file = os.path.join(lst_dir, name + ".lst")
    res_gff_file = os.path.join(gff_dir, name + ".gff")
    res_gene_file = os.path.join(gene_dir, name + ".gen")
    res_prt_file = os.path.join(prot_dir, name + ".prt")
    res_rep_file = os.path.join(rep_dir, name + ".fna")


    # Generate replicon file (same as input sequence but with gembase formatted headers). From
    # this file, get contig names, to be used to generate gff file
    contigs, sizes = utils.get_genome_contigs_and_rename(name, fna_file, res_rep_file, logger)
    if not contigs:
        try:
            os.remove(res_rep_file)
            os.remove(res_lst_file)
            os.remove(res_gff_file)
            os.remove(res_gene_file)
            os.remove(res_prt_file)
        except OSError:
            pass
        logger.error("Problems while generating Replicon file for {}".format(name))
        return False

    # Convert prokka tbl file to gembase .lst file format
    ok_tbl = tbl2lst(prokka_tbl_file, res_lst_file, contigs, name, fna_file)
    if not ok_tbl:
        try:
            os.remove(res_rep_file)
            os.remove(res_lst_file)
            os.remove(res_gff_file)
            os.remove(res_gene_file)
            os.remove(res_prt_file)
        except OSError:
            pass
        logger.error("Problems while generating LSTINFO file for {}".format(name))
        return False
    # Create gff3 file for annotations
    ok_gff = generate_gff(fna_file, prokka_gff_file, res_gff_file, res_lst_file, sizes, contigs)
    if not ok_gff:
        try:
            os.remove(res_rep_file)
            os.remove(res_lst_file)
            os.remove(res_gff_file)
            os.remove(res_gene_file)
            os.remove(res_prt_file)
        except OSError:
            pass
        logger.error("Problems while generating .gff file for {}".format(name))
        return False
    # create Genes file (and check no problem occurred with return code)
    ok_gene = create_gen(prokka_ffn_file, res_lst_file, res_gene_file)
    # If gene file not created because a problem occurred, return False:
    # format did not run for this genome
    if not ok_gene:
        try:
            os.remove(res_rep_file)
            os.remove(res_lst_file)
            os.remove(res_gff_file)
            os.remove(res_gene_file)
            os.remove(res_prt_file)
        except OSError:
            pass
        logger.error("Problems while generating .gen file for {}".format(name))
        return False


    # If gene file was created, create Proteins file
    ok_prt = create_prt(prokka_faa_file, res_lst_file, res_prt_file)
    # If protein file not created, return False: format did not run for this genome
    if not ok_prt:
        try:
            os.remove(res_lst_file)
            os.remove(res_gff_file)
            os.remove(res_gene_file)
            os.remove(res_prt_file)
            os.remove(res_rep_file)
            # Remove twice to be able to check that when there is a problem while removing files,
            # it generates the expected error
            os.remove(res_rep_file)
        except OSError:
            pass
        logger.error("Problems while generating .prt file for {}".format(name))
        return False
    return True


def tbl2lst(tblfile, lstfile, contigs, genome, gpath):
    """
    Read prokka tbl file, and convert it to the lst file.

    * prokka tbl file format::

        >Feature contig_name
        start end type
                [EC_number text]
                [gene text]
                inference ab initio prediction:Prodigal:2.6
                [inference text]
                locus_tag test
                product text

    where `type` can be CDS, tRNA, rRNA, etc ...
    lines between [] are not always present

    * lst file format::

        start end strand type locus gene_name | product | EC_number | inference 2 | db_xref

    with the same types as prokka file, and strain is C (complement) or D (direct)
    locus is: `<genome_name>.<contig_num><i or b>_<protein_num>`

    Parameters
    ----------
    tblfile : str
        name of prokka output tbl file to read
    lstfile : str
        name of lst file to generate
    contigs : dict
        {original_contig_name: gembase_contig_name}
    genome : str
        genome name (gembase format)
    gpath : str
        path to the genome given to prodigal. Only used for error message
    changed_name : bool
        True if contig names have been changed (cutn != 0) -> contig names end by '_num',
        False otherwise.

    Returns
    -------
    bool :
        True if genome name used in lstfile and prokka tblfile are the same, False otherwise
    """
    # Protein localisation in contig (b = border ; i = inside)
    cont_loc = "b"
    prev_cont_loc = "b"
    # Current contig number. Used to compare with new one, to know if protein is
    # inside or at the border of a contig
    prev_cont_num = -1
    prev_cont_name = ""
    tbl_cont_name = ""
    tbl_cont_num = 0  # Contig number in tbl file (<orig_name>_<cont_num>)
    lst_cont_num = 0  # Exp cont num: new feature means exp_cont_num += 1
    # prev_cont_num = -1  # Previous contig number seen
    # Information on current feature. At the beginning, everything empty, no information
    gene_name = "NA"
    locus_num = "NA"
    product = "NA"
    ecnum = "NA"
    inf2 = "NA"
    db_xref = "NA"
    start = "-1"
    end = "-1"
    strand = "D"
    # Feature type (CDS, tRNA...)
    feature_type = ""

    # Check that tblfile is not empty
    if os.stat(tblfile).st_size == 0:
        logger.error(f"{tblfile} is empty.")
        return False

    with open(tblfile, "r") as tblf, open(lstfile, "w") as lstf:
        for line in tblf:
            elems = line.strip().split("\t")
            if line.startswith(">Feature"):
                tbl_cont_name = line.split()[-1]
            else:
                if not tbl_cont_name:
                    logger.error(f"Wrong format for {tblfile}.")
                    return False
                # Get line type, and retrieve info according to it
                # If it is not the line with start, end, type, there are only 2 elements
                #  in the line:
                #  - information_type information_value
                if len(elems) == 2:
                    if "locus_tag" in elems[0]:
                        locus_num = elems[-1].split("_")[-1]
                    if "gene" in elems[0]:
                        gene_name = elems[1]
                    if "product" in elems[0]:
                        product = elems[1]
                    if "EC_number" in elems[0]:
                        ecnum = elems[1]
                    if "inference" in elems[0] and "ab initio prediction:Prodigal" not in line:
                        inf2 = elems[1]
                    if "db_xref" in elems[0]:
                        db_xref = elems[1]
                # Information on start, end, type
                else:
                    # new gene
                    # if new gene is not on the same contig as previously,
                    # get new contig and loc = 'b'
                    if tbl_cont_name != prev_cont_name:
                        # Check that this contig name is in the list, and get its gembase contig
                        # number
                        if tbl_cont_name in contigs:
                            contig_num = contigs[tbl_cont_name].split(".")[-1]
                        # if not in the list, problem, return false
                        else:
                            logger.error(f"'{tbl_cont_name}' found in {tblfile} does not exist in "
                                         f"{gpath}.")
                            return False
                        # Previous loc was 'i' (because we were in the same contig as the
                        # previous one). But now, we know the it was the last gene of
                        # its contig: we change loc to 'b'
                        prev_cont_loc = "b"
                        cont_loc = "b"
                    # Same contig as previously: this gene is inside the contig (cont_loc = "i").
                    # If, in fact, it was the last gene of this contig, it will be changed
                    # when discovering next gene.
                    else:
                        cont_loc = 'i'

                    # If not first gene of the contig, write the previous gene to .lst file
                    # The first gene will be written while reading the 2nd gene
                    if start != "-1" and end != "-1" and not crispr:
                        lstline = general.write_gene(feature_type, locus_num,
                                                     gene_name, product,
                                                     prev_cont_loc, genome,
                                                     prev_cont_num, ecnum, inf2,
                                                     db_xref, strand, start, end, lstf)

                    # Get new values for the next gene: start, end, strand and feature type
                    start, end, feature_type = elems
                    crispr = "CRISPR" in feature_type or "repeat_region" in feature_type
                    if crispr:
                        continue

                    # Get strain of gene
                    if int(end) < int(start):
                        start, end = end, start
                        strand = "C"
                    else:
                        strand = "D"
                    # Initialize variables for next feature
                    # (except start, end, strand and feature type that we just calculated).
                    # prev_cont_num = exp_cont_num
                    prev_cont_loc = cont_loc
                    prev_cont_num = contig_num
                    prev_cont_name = tbl_cont_name
                    locus_num = "NA"
                    gene_name = "NA"
                    product = "NA"
                    ecnum = "NA"
                    inf2 = "NA"
                    db_xref = "NA"

        # Write last feature
        if start != -1 and end != -1:
            prev_cont_loc = "b"
            general.write_gene(feature_type, locus_num, gene_name, product,
                               prev_cont_loc, genome, prev_cont_num,
                               ecnum, inf2, db_xref, strand, start, end, lstf)
    return True


def generate_gff(gpath, prokka_gff_file, res_gff_file, res_lst_file, sizes, contigs):
    """
    From the lstinfo file and contig names (retrieved from generation of Replicons files),
    generate a gff file.

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
        path to the genome sequence given to prokka. Only used for error message
    res-gff_file : str
        path to the gff file that must be created in result database
    res-lst_file : str
        path to the lst file that was created in result database in the previous step
    sizes : list
        dict of contig names with their size. {"gembase1": "size", "gembase2":"size2" ...]
    contigs : list
        dict of contig original and gembase names. {"contig1": "gembase1"...}

    Returns
    -------
    bool :
        True if conversion worked well, False otherwise
    """
    # open gff generated by prokka to read it
    # open file to write new gff file (in gff3 folder)
    # open lst file (from LSTINFO folder) to read all annotation information saved
    # from prokka results
    with open(prokka_gff_file, "r") as prokf, open(res_lst_file, "r") as lstf, \
            open(res_gff_file, "w") as gfff:
        # Write headers of gff3 file
        gfff.write("##gff-version  3\n")
        # Write all sequences with their size. Order by name in gembase format
        for old, new in sorted(contigs.items(), key=lambda items:items[1]):
            end = sizes[new]
            # Write the list of cobisntigs, with their size
            gfff.write(f"##sequence-region\t{new}\t{1}\t{end}\n")

        # Now, convert each line of prokka gff to gembase formatted gff line
        for linegff in prokf:
            # Ignore gff lines starting by #. For the new gff, these are already all written at the
            # beginning of the file.
            if linegff.startswith("#"):
                continue
            # We do not write sequences
            if linegff.startswith('>'):
                break
            # Try to extract information on new line
            try:
                # Get all information from prokka gff line. Strip each of them as trailing whitespace
                # can be hidden (leading to information from gff considered as different from
                # information from lst)
                fields_g = linegff.strip().split("\t")
                fields_g = [info.strip() for info in fields_g]
                # Ignore lines with sequence
                # if len(fields_g) != 9:
                #     continue
                (contig_name, source, type_g, start_g, end_g,
                 score, strand_g, phase, attributes) = fields_g
                # Ignore CRISPR
                if "CRISPR" in type_g or "repeat_region" in type_g:
                    continue
                # Get information given to this same sequence from the lst file
                # (next lst line corresponds to next gff line without #), as, for each format,
                # there is 1 line per gene)
                linelst = lstf.readline()
                fields_l = linelst.split("\t")
                fields_l = [info.strip() for info in fields_l]
                start_l, end_l, strand_l, type_l, locus_l, l_gene, l_info = fields_l
                # Get gff and ffn filenames to give information to user if error message
                gff = os.path.basename(prokka_gff_file)
                tbl = gff.replace(".gff", ".tbl")
                # Path where gff and ffn generated by prodigal are
                tmp = gpath + "-prokkaRes"
                # Get gene name given by prodigal to current gene
                gname = attributes.split("ID=")[1].split(";")[0]
                # Get locus_tag given by prokka to current feature (should be the same as ID)
                loc_name = attributes.split("locus_tag=")[1].split(";")[0]
                if loc_name != gname:
                    logger.error(f"Problem in {gff}: ID={gname} whereas locus_tag={loc_name}.")
                    return False
                # Compare information from lst and information from prodigal gff (start,
                # end and type of feature). They should correspond
                for (elemg, eleml, label) in zip([start_g, end_g, type_g],
                                                 [start_l, end_l, type_l],
                                                 ["start", "end", "type"]):
                    # If 1 element is different (start or end position, or type), print error
                    # message and return False: this genome could not be converted to gff
                    if elemg != eleml:
                        logger.error(f"Files {tbl} and {gff} (in prokka tmp_files: {tmp}) "
                                     f"do not have the same {label} value for gene {gname} ({elemg} "
                                     f"in gff, {eleml} in tbl)")
                        return False

                # Replace prokka ID and locus_tag with the new gene name (in gembase format),
                # found in the corresponding lst line
                at_split = attributes.split(";")
                newID = [atr if "ID" not in atr else f'ID={locus_l}' for atr in at_split]
                new = [atr if "locus_tag" not in atr else f'locus_tag={locus_l}' for atr in newID]

                # Write new line of gff file
                cname = contigs[contig_name]
                info = "\t".join([cname, source, type_g, start_g, end_g, score, strand_g,
                                  phase, ";".join(new)])
                gfff.write(info + "\n")
            except:
                logger.error(f"Wrong format for {prokka_gff_file}.")
                return False
        return True


def create_gen(ffnseq, lstfile, genseq):
    """
    Generate .gen file, from sequences contained in .ffn, but changing the
    headers using the information in .lst

    Parameters
    ----------
    ffnseq : str
        .ffn file generated by prokka
    lstfile : str
        lstfile converted from prokka tbl file
    genseq : str
        output file, to write in Genes directory
    logger : logging.Logger
        logger object to put information

    Returns
    -------
    bool :
        True if conversion went well, False otherwise
    """
    problem = False
    write = True  # Write next sequence
    with open(ffnseq) as ffn, open(lstfile) as lst, open(genseq, "w") as gen:
        for line_ffn in ffn:
            # Ignore gene that we do not want to write (should be a crispr)
            # If line of sequence, write it as is, and go to next line
            if not line_ffn.startswith(">"):
                # We just read a seq line. If we can write (write is True), do it and go
                # to next line
                # Otherwise, just go to next line
                if write:
                    gen.write(line_ffn)
                continue
            # Try to get gene ID. If does not work, ignore this gene (it may be a
            # CRISPR, and we ignore them
            test_gen_id = line_ffn.split()[0].split("_")[-1]
            if not test_gen_id.isdigit():
                # Maybe a CRISPR? Or wrong gene name? -> ignore
                logger.log(utils.detail_lvl(),
                           f"Unknown header format for {line_ffn.strip()}. "
                           "This gene will be ignored in .gen output file.")
                write = False
                continue
            # If ffn contains a gene header, find its information in lst file
            else:
                write = True
                lstline = lst.readline().strip()
                gen_id = int(test_gen_id)
                # genID exists, ffn header is for a gene. Check that it corresponds to
                # information in lst file.
                id_lst = lstline.split("\t")[4].split("_")[-1]
                gen_id_lst = int(id_lst)
                # in lst, find the same gene ID as in ffn (some gene IDs in lst can be absent
                # from ffn, if prokka do not give their sequence).
                # As they are ordered by increasing number, go to next lstline until
                # corresponding gene ID is found. However, if ffn ID > lst ID: ID does not
                # exist in .lst -> problem.
                while gen_id > gen_id_lst:
                    lstline = lst.readline().strip()
                    if not lstline:
                        gen_id_lst = "-1"
                        break
                    id_lst = lstline.split("\t")[4].split("_")[-1]
                    gen_id_lst = int(id_lst)
                # If it found the same gene ID, write info in gene file
                if gen_id == gen_id_lst:
                    general.write_header(lstline.strip(), gen)
                # If gene ID of ffn not found, write error message and stop
                else:
                    logger.error(f"Missing info for gene {line_ffn.strip()} "
                                 f"(from {ffnseq}) in {lstfile}. If it is actually present "
                                 "in the lst file, check that genes are ordered by increasing number in both lst and ffn files.")
                    return False
    return True


def create_prt(faaseq, lstfile, prtseq):
    """
    Generate .prt file, from sequences in .faa, but changing the headers
    using information in .lst

    **Note:** works if proteins are in increasing order (of number after "_" in their name)
    in faa and tbl (hence lst) files.

    If a header is not in the right format, or a protein exists in prt file but not in lstfile,
    conversion is stopped, an error message is output, and prt file is removed.

    Parameters
    ----------
    faaseq : str
        faa file output of prokka
    lstfile : str
        lstinfo converted from prokka tab file
    prtseq : str
        output file where converted proteins must be saved

    Returns
    -------
    bool :
        True if conversion went well, False otherwise
    """
    problem = False
    with open(faaseq) as faa, open(lstfile) as lst, open(prtseq, "w") as prt:
        for line in faa:
            # all header lines must start with PROKKA_<geneID>
            if line.startswith(">"):
                try:
                    # get gene ID
                    gen_id = int(line.split()[0].split("_")[-1])
                except ValueError as err:
                    logger.error(f"Unknown header format {line.strip()} in {faaseq}. "
                                 f"Gene ID is not a number.")
                    return False
                gen_id_lst = 0
                # get line of lst corresponding to the gene ID
                lstline = ""
                while gen_id > gen_id_lst:
                    lstline = lst.readline().strip()
                    id_lst = lstline.split("\t")[4].split("_")[-1]
                    # don't cast to int if info for a crispr
                    if id_lst.isdigit():
                        gen_id_lst = int(id_lst)
                # check that gen_id is the same as the lst line
                if gen_id == gen_id_lst:
                    general.write_header(lstline, prt)
                else:
                    logger.error(f"Missing info for protein {line.strip()} (from {faaseq}) "
                                 f"in {lstfile}. If it is actually present "
                                 "in the lst file, check that proteins are ordered by increasing "
                                 "number in both lst and faa files.")
                    return False
            # not header: inside sequence, copy it to the .prt file
            else:
                prt.write(line)
    return True


