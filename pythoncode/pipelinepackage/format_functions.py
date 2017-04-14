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

logger = logging.getLogger()


def format_genomes(genomes, results, res_path):
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
    skipped = []  # list of genomes skipped
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
        format_one_genome(gpath, name, lst_dir, prot_dir, gene_dir, rep_dir)
        bar.update(num + 1)
    bar.finish()
    return skipped


def format_one_genome(gpath, name, lst_dir, prot_dir, gene_dir, rep_dir):
    """
    Format the given genome, and create its corresponding files in the following folders:
    - Proteins
    - Genes
    - Replicons
    - LSTINFO

    arguments:
    * gpath: path to the genome sequence which was given to prokka for annotation
    * name: gembase name of the genome
    * lst_dir: path to LSTINFO folder
    * prot_dir: path to Proteins folder
    * gene_dir: path to Genes folder
    * rep_dir: path to Replicons folder
    """
    # convert tbl to lst
    tblgenome = os.path.join(gpath + "-prokkaRes", name + ".tbl")
    lstgenome = os.path.join(lst_dir, name + ".lst")
    tbl2lst(tblgenome, lstgenome)
    # copy fasta file to Replicons
    rep_file = os.path.join(rep_dir, name + ".fna")
    shutil.copyfile(gpath, rep_file)
    # create Genes file
    ffngenome = os.path.join(gpath + "-prokkaRes", name + ".ffn")
    gengenome = os.path.join(gene_dir, name + ".gen")
    create_gen(ffngenome, lstgenome, gengenome, name)
    # Create Proteins file
    faagenome = os.path.join(gpath + "-prokkaRes", name + ".faa")
    prtgenome = os.path.join(prot_dir, name + ".prt")
    create_prt(faagenome, lstgenome, prtgenome)



def create_gen(ffnseq, lstfile, genseq, genome_name):
    """
    Generate .gen file, from sequences contained in .ffn, but changing the
    headers using the information in .lst

    - ffnseq = .ffn file generated by prokka
    - lstfile = lstfile converted from prokka tbl file
    - genseq = output file, to write in Genes directory
    - genome_name = gembase name. With prokka, genes corresponding to crispr repeat are
    called with the genome name.
    """
    crisprID = 1
    with open(ffnseq, "r") as ffn, open(lstfile, "r") as lst, open(genseq, "w") as gen:
        for line in ffn:
            # header starting with genome_name: CRISPR
            if line.startswith(">" + genome_name):
                # get next line of lst file
                lstline = lst.readline()
                # check crispr ID is the same as in the current lst line
                crisprIDlst = int(lstline.split("\t")[4].split("_CRISPR")[1])
                if (crisprID == crisprIDlst):
                    write_header(lstline, gen)
                    crisprID += 1
                else:
                    logger.error("ERROR, missing info for gene: {} in {}".format(line, lstfile))
            # header starting with >something_<geneID>: not a crispr
            elif line.startswith(">"):
                # get geneID of this header, and the next line of the lst file
                try:
                    genID = line.split()[0].split("_")[-1]
                except Exception as err:
                    log.error(("Unknown header format {} in {}. "
                               "Error: {}").format(line, ffnseq, err))
                lstline = lst.readline()
                genIDlst = lstline.split("\t")[4].split("_")[1]
                # check that genID is the same as the lst line
                if (genID == genIDlst):
                    write_header(lstline, gen)
                else:
                    logger.error("ERROR, missing info for gene: {} in {}".format(line, lstfile))
            # not header: inside sequence, copy it to .gen file
            else:
                gen.write(line)


def create_prt(faaseq, lstfile, prtseq):
    """
    Generate .prt file, from sequences in .faa, but changing the headers
    using information in .lst

    * faaseq: faa file output of prokka
    * lstfile: lstinfo converted from prokka tab file
    * prtseq: output file where converted proteins must be saved
    """
    with open(faaseq, "r") as faa, open(lstfile, "r") as lst, open(prtseq, "w") as prt:
        for line in faa:
            # all header lines must start with PROKKA_<geneID>
            if line.startswith(">"):
                try:
                    # get gene ID
                    genID = int(line.split()[0].split("_")[1])
                except Exception as err:
                    log.error(("Unknown header format {} in {}. "
                               "Error: {}").format(line, faaseq, err))
                genIDlst = 0
                # get line of lst corresponding to the gene ID
                while (genID > genIDlst):
                    lstline = lst.readline()
                    IDlst = lstline.split("\t")[4].split("_")[1]
                    # don't cast to int if info for a crispr
                    if(IDlst.isdigit()):
                        genIDlst = int(IDlst)
                # check that genID is the same as the lst line
                if (genID == genIDlst):
                    write_header(lstline, prt)
                else:
                    logger.error("Missing info for protein {} in {}".format(line, lstfile))
            # not header: inside sequence, copy it to the .prt file
            else:
                prt.write(line)


def write_header(lstline, outfile):
    """
    write heaader to output file. Header is generated from the lst line.
    """
    name = lstline.split("\t")[4]
    size = int(lstline.split("\t")[1]) - int(lstline.split("\t")[0]) + 1
    geneName = lstline.split("\t")[5]
    info = lstline.split("\t")[6]
    towrite = " ".join([name, str(size), geneName, info])
    outfile.write(">" + towrite)


def tbl2lst(tblfile, lstfile):
    """
    Read prokka tbl file, and convert it to the lst file.

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
    crispr_num = 1
    start, end = -1, -1
    strand, gtype, genome, cont_num, locus_num = [""] * 5
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
                    strand, gtype, genome, cont_num, locus_num = [""] * 5
                    gene_name, product, ecnum, inf2 = ["NA"] * 4
                # Read new feature
                contig = line.strip().split()[-1]
                genome = ".".join(contig.split(".")[:-1])
                cont_num = contig.split(".")[-1]
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


def write_gene(gtype, locus_num, gene_name, product, crispr_num, cont_loc,
               genome, cont_num, ecnum, inf2, strain, start, end, lstopenfile):
    """
    Write given gene to output file
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
    more_info = "| {} | {} | {}".format(product, ecnum, inf2)
    lst_line = "\t".join([str(start), str(end), strain, gtype,
                          locus_name, gene_name, more_info])
    lstopenfile.write(lst_line + "\n")
    return crispr_num
