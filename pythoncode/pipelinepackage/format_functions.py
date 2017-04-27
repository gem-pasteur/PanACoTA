#!/usr/bin/env python3
# coding: utf-8

"""
Functions to:
- convert prokka tbl file to our tab file
- convert prokka ffn and faa headers to our format
- Create the database, with the following folders in the given "res_path":
    - Proteins: multifasta with all CDS in aa
    - Replicons: multifasta of genome
    - Genes: multifasta of all genes in nuc
    - LSTINFO: information on annotation. Columns are: `start end strand type locus
    gene_name | product | EC_number | inference2` with the same types as prokka file,
    and strain is C (complement) or D (direct). Locus is:
    `<genome_name>.<i or b><contig_num>_<protein_num>`

@author gem
April 2017
"""

import os
import shutil
import logging
import progressbar
import glob
from pipelinepackage import genome_seq_functions as gfunc

logger = logging.getLogger()


def format_genomes(genomes, results, res_path, prok_path):
    """
    For all genomes which were annotated by prokka, reformat them
    in order to have, in 'res_path', the following folders:
    * LSTINFO: containing a .lst file for each genome, with all genes
    * Replicons: containing all multifasta files
    * Genes: containing 1 multi-fasta per genome, with all its genes in nuc
    * Proteins: containing 1 multi-fasta per genome, with all its proteins in aa

    - genomes = {genome: [name, gpath, size, nbcont, l90]}
    - results = {genome: bool} True if prokka ran well, False otherwise
    - res_path = path to folder where the 4 directories must be created
    - prok_path = path to folder named "<genome_name>-prokkaRes" where all prokka
    results are saved.
    """
    logger.info("Formatting all genomes")
    lst_dir = os.path.join(res_path, "LSTINFO")
    prot_dir = os.path.join(res_path, "Proteins")
    gene_dir = os.path.join(res_path, "Genes")
    rep_dir = os.path.join(res_path, "Replicons")
    os.makedirs(lst_dir, exist_ok=True)
    os.makedirs(prot_dir, exist_ok=True)
    os.makedirs(gene_dir, exist_ok=True)
    os.makedirs(rep_dir, exist_ok=True)

    nbgen = len(genomes)
    # Create progressbar
    widgets = ['Formatting genomes: ', progressbar.Bar(marker='█', left='', right='', fill=' '),
               ' ', progressbar.Percentage()]
    bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=100).start()
    skipped = []  # list of genomes skipped: no format step run
    skipped_format = []  # List of genomes for which forat step had problems
    for num, (genome, (name, gpath, _, _, _)) in enumerate(genomes.items()):
        # Ignore genomes with bad quality (not annotated)
        if genome not in results:
            bar.update(num + 1)
            continue
        # if prokka did not run well for a genome, don't format it
        if not results[genome]:
            skipped.append(genome)
            bar.update(num + 1)
            continue
        ok_format = format_one_genome(gpath, name, prok_path, lst_dir, prot_dir, gene_dir, rep_dir)
        if not ok_format:
            skipped_format.append(genome)
        bar.update(num + 1)
    bar.finish()
    return skipped, skipped_format


def format_one_genome(gpath, name, prok_path, lst_dir, prot_dir, gene_dir, rep_dir):
    """
    Format the given genome, and create its corresponding files in the following folders:
    - Proteins
    - Genes
    - Replicons
    - LSTINFO

    arguments:
    * gpath: path to the genome sequence which was given to prokka for annotation
    * name: gembase name of the genome
    * prok_path: directory where prokka folders are saved
    * lst_dir: path to LSTINFO folder
    * prot_dir: path to Proteins folder
    * gene_dir: path to Genes folder
    * rep_dir: path to Replicons folder

    Returns: True if format worked for this genome, False otherwise.
    """
    prokka_dir = os.path.join(prok_path, os.path.basename(gpath) + "-prokkaRes")
    # Get .tbl file
    tblgenome = glob.glob(os.path.join(prokka_dir, "*.tbl"))
    if len(tblgenome) == 0:
        logger.error(("No .tbl file in {}.").format(prokka_dir))
        return False
    elif len(tblgenome) > 1:
        logger.error("Several .tbl files in {}.".format(prokka_dir))
        return False
    else:
        tblgenome = tblgenome[0]
    lstgenome = os.path.join(lst_dir, name + ".lst")
    # convert tbl to lst, and check that genome name in tbl is the same as the given genome name
    # If the user ran prokka, and then changed the genome names (for example, by changing
    # L90 threshold), and wants to format again its genomes with the new names, but without
    # running prokka again, we must change here the contig names (containing genome name).
    same_name = tbl2lst(tblgenome, lstgenome, name)
    # create Genes file (and check no problem occurred with return code)
    ffngenome = glob.glob(os.path.join(prokka_dir, "*.ffn"))
    if len(ffngenome) == 0:
        logger.error(("No .ffn file in {}.").format(prokka_dir))
        os.remove(lstgenome)
        return False
    elif len(ffngenome) > 1:
        logger.error("Several .ffn files in {}.".format(prokka_dir))
        os.remove(lstgenome)
        return False
    else:
        ffngenome = ffngenome[0]
    gengenome = os.path.join(gene_dir, name + ".gen")
    ok_gene = create_gen(ffngenome, lstgenome, gengenome)
    # If gene file not created because a problem occurred, return False:
    # format did not run for this genome
    if not ok_gene:
        os.remove(lstgenome)
        return False
    # If gene file was created, create Proteins file
    faagenome = glob.glob(os.path.join(prokka_dir, "*.faa"))
    if len(faagenome) == 0:
        logger.error(("No .faa file in {}.").format(prokka_dir))
        os.remove(lstgenome)
        os.remove(gengenome)
        return False
    elif len(faagenome) > 1:
        logger.error("Several .faa files in {}.".format(prokka_dir))
        os.remove(lstgenome)
        os.remove(gengenome)
        return False
    else:
        faagenome = faagenome[0]
    prtgenome = os.path.join(prot_dir, name + ".prt")
    ok_prt = create_prt(faagenome, lstgenome, prtgenome)
    # If protein file not created, return False: format did not run for this genome
    if not ok_prt:
        os.remove(lstgenome)
        os.remove(gengenome)
        return False
    # Create Replicons file
    rep_file = os.path.join(rep_dir, name + ".fna")
    # If the genome name did not change after prokka run, copy genome sequence
    # given to prokka to Replicons folder.
    if same_name:
        shutil.copyfile(gpath, rep_file)
    # otherwise, change headers in the file.
    else:
        gfunc.rename_genome_contigs(name, gpath, rep_file)
    return True


def create_gen(ffnseq, lstfile, genseq):
    """
    Generate .gen file, from sequences contained in .ffn, but changing the
    headers using the information in .lst

    - ffnseq = .ffn file generated by prokka
    - lstfile = lstfile converted from prokka tbl file
    - genseq = output file, to write in Genes directory
    """
    problem = False
    crisprID = 1
    with open(ffnseq, "r") as ffn, open(lstfile, "r") as lst, open(genseq, "w") as gen:
        for line_ffn in ffn:
            # If line of sequence, write it as is, and go to next line
            if not line_ffn.startswith(">"):
                gen.write(line_ffn)
                continue
            # if lstline indicates a CRISPR, header in ffn file is the genome name. Hence,
            # it should not contain a '_' followed by a number.
            lstline = lst.readline().strip()
            if "CRISPR" in lstline:
                if '_' in line_ffn:
                    if line_ffn.strip().split("_")[-1].isdigit():
                        logger.error("According to lstinfo file, gene {} should be a CRISPR. "
                                     "However, its name has the same format as a gene name (not "
                                     "CRISPR). Format function will stop here, and gen file will "
                                     "not be created for {}.".format(line_ffn.strip(), ffnseq))
                        problem = True
                        break
                # check crispr ID is the same as in the current lst line
                crisprIDlst = int(lstline.split("\t")[4].split("_CRISPR")[1])
                if (crisprID == crisprIDlst):
                    write_header(lstline, gen)
                    crisprID += 1
                else:
                    logger.error(("Problem with CRISPR numbers in {}. CRISPR {} in ffn is "
                                  "CRISPR num {}, which is not found at this place in "
                                  "lstinfo file.").format(lstfile, line_ffn.strip(), crisprID))
                    problem = True
                    break
            else:
                # get geneID of this header, and the next line of the lst file
                try:
                    genID = int(line_ffn.split()[0].split("_")[-1])
                except Exception as err:
                    logger.error(("Unknown header format {} in {}. "
                               "Error: {}\nGen file will not be "
                               "created.").format(line_ffn.strip(), ffnseq, err))
                    problem = True
                    break
                try:
                    genIDlst = int(lstline.split("\t")[4].split("_")[-1])
                except Exception as err:
                    logger.error(("Unknown gene format {} in {}. "
                               "Error: {}\nGen file will not be "
                               "created.").format(lstline.strip(), lstfile, err))
                    problem = True
                    break
                # check that genID is the same as the lst line
                if (genID == genIDlst):
                    write_header(lstline.strip(), gen)
                else:
                    logger.error("Missing info for gene {} in {}. If it is actually present "
                                 "in the lst file, check that genes are ordered by increasing "
                                 "number in both lst and ffn files.\nGen file not created "
                                 "from {}.".format(line_ffn.strip(), lstfile, ffnseq))
                    problem = True
                    break
    if problem:
        os.remove(genseq)
    return not problem


def create_prt(faaseq, lstfile, prtseq):
    """
    Generate .prt file, from sequences in .faa, but changing the headers
    using information in .lst

    * faaseq: faa file output of prokka
    * lstfile: lstinfo converted from prokka tab file
    * prtseq: output file where converted proteins must be saved

    Note: works if proteins are in increasing order (of number after "_" in their name)
    in faa and tbl (hence lst) files.

    If a header is not in the right format, or a protein exists in prt file but not in lstfile,
    conversion is stopped, an error message is output, and prt file is removed.
    """
    problem = False
    with open(faaseq, "r") as faa, open(lstfile, "r") as lst, open(prtseq, "w") as prt:
        for line in faa:
            # all header lines must start with PROKKA_<geneID>
            if line.startswith(">"):
                try:
                    # get gene ID
                    genID = int(line.split()[0].split("_")[-1])
                except Exception as err:
                    logger.error(("Unknown header format {} in {}. "
                               "Error: {}\nPrt file not created "
                               "from {}.").format(line.strip(), faaseq, err, faaseq))
                    problem = True
                    break
                genIDlst = 0
                # get line of lst corresponding to the gene ID
                while (genID > genIDlst):
                    lstline = lst.readline().strip()
                    IDlst = lstline.split("\t")[4].split("_")[-1]
                    # don't cast to int if info for a crispr
                    if(IDlst.isdigit()):
                        genIDlst = int(IDlst)
                # check that genID is the same as the lst line
                if (genID == genIDlst):
                    write_header(lstline, prt)
                else:
                    logger.error(("Missing info for protein {} in {}. If it is actually present "
                                  "in the lst file, check that proteins are ordered by increasing "
                                  "number in both lst and faa files.\nPrt file not created "
                                  "from {}.").format(line.strip(), lstfile, faaseq))
                    problem = True
                    break
            # not header: inside sequence, copy it to the .prt file
            else:
                prt.write(line)
    if problem:
        os.remove(prtseq)
    return not problem


def write_header(lstline, outfile):
    """
    write heaader to output file. Header is generated from the lst line.
    """
    name = lstline.split("\t")[4]
    size = int(lstline.split("\t")[1]) - int(lstline.split("\t")[0]) + 1
    geneName = lstline.split("\t")[5]
    info = lstline.split("\t")[6]
    towrite = " ".join([name, str(size), geneName, info])
    outfile.write(">" + towrite + "\n")


def tbl2lst(tblfile, lstfile, genome):
    """
    Read prokka tbl file, and convert it to the lst file.

    genome = genome name (gembase format)

    * prokka tbl file format:
    ```
    >Feature contig_name
    start end type
            [EC_number text]
            [gene text]
            inference ab initio prediction:Prodigal:2.6
            [inference text]
            locus_tag test
            product text
    ```
    where `type` can be CDS, tRNA, rRNA, etc ...
    lines between [] are not always present

    * lst file format:
    ```
    start end strand type locus gene_name | product | EC_number | inference 2
    ```
    with the same types as prokka file, and strain is C (complement) or D (direct)
    locus is: `<genome_name>.<i or b><contig_num>_<protein_num>`

    """
    same_name = True
    crispr_num, cont_num = 1, 0
    start, end = -1, -1
    strand, gtype, locus_num = [""] * 3
    gene_name, product, ecnum, inf2 = ["NA"] * 4
    cont_loc = "b"
    with open(tblfile, "r") as tblf, open(lstfile, "w") as lstf:
        for line in tblf:
            elems = line.strip().split("\t")
            # New feature (= new contig)
            if line.startswith(">Feature"):
                # if there was something in the previous contig, print last gene of
                # previous contig (with contigLoc=b), and reinitiate for new contig/genes
                if start != -1 and end != -1:
                    cont_loc = "b"
                    crispr_num = write_gene(gtype, locus_num, gene_name, product,
                                            crispr_num, cont_loc, genome, cont_num,
                                            ecnum, inf2, strain, start, end, lstf)
                    # init for next feature
                    start, end = -1, -1
                    strand, gtype, locus_num = [""] * 3
                    gene_name, product, ecnum, inf2 = ["NA"] * 4
                # Read new feature
                if same_name:
                    contig = line.strip().split()[-1]
                    c_genome = ".".join(contig.split(".")[:-1])
                    if c_genome != genome:
                        same_name = False
                    else:
                        cont_num = int(contig.split(".")[-1])
                if not same_name:
                    cont_num += 1
            # Line indicating position of gene
            if len(elems) == 3:
                # if not first gene of the contig, write previous gene
                if start != -1 and end != -1:
                    crispr_num = write_gene(gtype, locus_num, gene_name, product,
                                            crispr_num, cont_loc, genome, cont_num,
                                            ecnum, inf2, strain, start, end, lstf)
                    gtype, strand, locus_num = [""] * 3
                    start, end = -1, -1
                    gene_name, product, ecnum, inf2 = ["NA"] * 4
                    cont_loc = "i"
                start, end, gtype = int(elems[0]), int(elems[1]), elems[2]
                # Get strain of gene
                if end < start:
                    start, end = end, start
                    strain = "C"
                else:
                    strain = "D"
            if "locus_tag" in elems[0] and len(elems) == 2:
                locus_num = elems[-1].split("_")[-1]
            if "gene" in elems[0] and len(elems) == 2:
                gene_name = elems[1]
            if "product" in elems[0] and len(elems) == 2:
                product = elems[1]
            if "EC_number" in elems[0] and len(elems) == 2:
                ecnum = elems[1]
            if ("inference" in elems[0] and "ab initio prediction:Prodigal" not in line
                and len(elems) == 2):
                inf2 = elems[1]
        # Write last gene:
        if start != -1 and end != -1:
            cont_loc = "b"
            write_gene(gtype, locus_num, gene_name, product, crispr_num, cont_loc,
                       genome, cont_num, ecnum, inf2, strain, start, end, lstf)
    return same_name


def write_gene(gtype, locus_num, gene_name, product, crispr_num, cont_loc,
               genome, cont_num, ecnum, inf2, strand, start, end, lstopenfile):
    """
    Write given gene to output file

    * gtype: type of feature (CDS, tRNA, etc.)
    * locus_num: number of locus given by prokka
    * gene_name: gene name found by prokka ("NA" if no gene name)
    * product: found by prokka, "NA" if no product
    * crispr_num: current crispr number. In prokka tbl, CRISPRs are not numbered, they all
    have the same name. We name them by adding a unique number to each CRISPR. If the current
    gene to add is a CRISPR, this number will be incremented and returned. If not, this same
    name will be returned.
    * cont_loc: 'i' if the gene is inside a contig, 'b' if its on the border (first or last gene
    of the contig)
    * genome: genome name (spegenus.date.strain_num)
    * cont_num: contig number
    * ecnum: EC number, found by prokka, or "NA" if no EC number
    * inf2: more information found by prokka, or "NA" if no more information
    * strand: C (complement) or D (direct)
    * start: start of gene in the contig
    * end: end of gene in the contig
    * lstopenfile: open file where lstinfo must be written
    """
    # if last gene was a crispr
    if gtype == "repeat_region":
        gtype = "CRISPR"
        locus_num = "CRISPR" + str(crispr_num)
        gene_name = "crispr"
        product = "crispr-array"
        crispr_num += 1
    locus_name = "{}.{}{}_{}".format(genome, cont_loc,
                                     str(cont_num).zfill(4),
                                     str(locus_num).zfill(5))
    # If '|' character found in those fields, replace by '_' to avoid problems while parsing
    more_info = "| {} | {} | {}".format(product.replace("|", "_"),
                                        ecnum.replace("|", "_"),
                                        inf2.replace("|", "_"))
    lst_line = "\t".join([str(start), str(end), strand, gtype,
                          locus_name, gene_name, more_info])
    lstopenfile.write(lst_line + "\n")
    return crispr_num
