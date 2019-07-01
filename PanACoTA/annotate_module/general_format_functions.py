#!/usr/bin/env python3
# coding: utf-8

"""
Functions to convert prokka or prodigal results to gembase format:

    * Proteins: multifasta with all CDS in aa
    * Replicons: multifasta of genome
    * Genes: multifasta of all genes in nuc
    * gff3: gff files without sequence
    * LSTINFO:
        - if annotated by prokka: information on annotation. Columns are:
        "start end strand type locus gene_name | product | EC_number | inference2"
        with the same types as prokka file, and strain is C (complement) or D (direct).
        Locus is: `<genome_name>.<contig_num><i or b>_<protein_num>`
        - if annotated by prodigal

@author gem
May 2019
"""

import os
import logging
import logging.handlers
import progressbar
import multiprocessing
import threading
import PanACoTA.utils as utils
import PanACoTA.annotate_module.format_prokka as fprokka
import PanACoTA.annotate_module.format_prodigal as fprodigal


def format_genomes(genomes, results, res_path, annot_path, prodigal_only, threads=1, quiet=False):
    """
    For all genomes which were annotated (by prokka or prodigal), reformat them
    in order to have, in 'res_path', the following folders:

    * LSTINFO: containing a .lst file for each genome, with all genes
    * Replicons: containing all multifasta files
    * Genes: containing 1 multi-fasta per genome, with all its genes in nuc
    * Proteins: containing 1 multi-fasta per genome, with all its proteins in aa
    * gff: containing all gff files

    Parameters
    ----------
    genomes : dict
        {genome: [name, gpath, size, nbcont, l90]}
    results : list
        List of genome annotated by prokka/prodigal, that must be formatted
    res_path : str
        path to folder where the 4 directories must be created
    annot_path : str
        path to folder named "<genome_name>-[prokka, prodigal]Res" where all prokka/prodigal
        results are saved.
    prodigal_only: True if it was annotated by prodigal, False if annotated by prokka
    threads : int
        number of threads to use to while formatting genomes
    quiet : bool
        True if nothing must be sent to stderr/stdout, False otherwise

    Returns
    -------
    (skipped, skipped_format) : tuple

        * skipped : list of genomes skipped because they had a problem in annotation step
        * skipped_format : list of genomes skipped because they had a problem in format step
    """

    main_logger = logging.getLogger("qc_annotate.geneffunc")
    main_logger.info("Formatting all genomes")
    lst_dir = os.path.join(res_path, "LSTINFO")
    prot_dir = os.path.join(res_path, "Proteins")
    gene_dir = os.path.join(res_path, "Genes")
    rep_dir = os.path.join(res_path, "Replicons")
    gff_dir = os.path.join(res_path, "gff3")
    os.makedirs(lst_dir, exist_ok=True)
    os.makedirs(prot_dir, exist_ok=True)
    os.makedirs(gene_dir, exist_ok=True)
    os.makedirs(rep_dir, exist_ok=True)
    os.makedirs(gff_dir, exist_ok=True)

    # If this function goes until here, it means that there is at least 1 genome to annotate
    nbgen = len(genomes)
    bar = None
    if not quiet:
        # Create progressbar
        widgets = ['Formatting genomes: ',
                   progressbar.Bar(marker='█', left='', right=''),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ") - ", progressbar.Timer()]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=100).start()
    # Create a Queue to put logs from processes, and handle them after from a single thread
    m = multiprocessing.Manager()
    q = m.Queue()

    # if at least 1 genome ok, try to format it
    params = [(genome, name, gpath, annot_path, lst_dir, prot_dir, gene_dir,
               rep_dir, gff_dir, results, prodigal_only, q)
              for genome, (name, gpath, _, _, _) in genomes.items()]

    # Create pool and launch parallel formating steps
    pool = multiprocessing.Pool(threads)
    final = pool.map_async(handle_genome, params, chunksize=1)
    pool.close()
    lp = threading.Thread(target=utils.logger_thread, args=(q,))
    lp.start()
    if not quiet:
        while True:
            if final.ready():
                break
            remaining = final._number_left
            bar.update(nbgen - remaining)
        bar.finish()
    pool.join()
    q.put(None)
    lp.join()
    res = final.get()

    # List of genomes for which format step had problems: will be filled after
    skipped_format = []
    for output in res:
        if not output[0]:
            skipped_format.append(output[1])
    return skipped_format


def handle_genome(args):
    """
    For a given genome, check if it has been annotated (in results), if annotation
    (by prokka or prodigal) ran without problems (result = True). In that case,
    format the genome and get the output to see if everything went ok.

    Parameters
    ----------
    args : tuple
        (genome, name, gpath, prok_path, lst_dir, prot_dir,\
         gene_dir, rep_dir, gff_dir, results, q)\
         with:

         * genome : original genome name
         * name : gembase name of the genome
         * gpath : path to the genome sequence which was given to prokka/prodigal for annotation
         * annot_path : directory where prokka/prodigal folders are saved
         * lst_dir : path to 'LSTINFO' folder
         * prot_dir : path to 'Proteins' folder
         * gene_dit : path to 'Genes' folder
         * rep_dir : path to 'Replicons' folder
         * gff_dir : path to 'gff3' folder
         * results : {genome_name: <True if genome was correctly annotated, False otherwise>}
         * prodigal_only : True if annotated by prodigal, False if annotated by prokka
         * q : multiprocessing.managers.AutoProxy[Queue] queue to put logs during subprocess

    Returns
    -------
    (bool, str) :

        * True if genome was annotated as expected, False otherwise
        * genome name (used to get info from the pool.map_async)
    """
    (genome, name, gpath, annot_path, lst_dir, prot_dir,
     gene_dir, rep_dir, gff_dir, results, prodigal_only, q) = args

    # Define which formatting must be used, given the annotation software
    if prodigal_only:
        format_one_genome = fprodigal.format_one_genome
    else:
        format_one_genome = fprokka.format_one_genome

    # Set logger for this process
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('format.handle_genome')
    # Handle genome, if we have its informations
    if genome not in results:
        return "no_res", genome
    ok_format = format_one_genome(gpath, name, annot_path, lst_dir,
                                  prot_dir, gene_dir, rep_dir, gff_dir, logger)
    return ok_format, genome


def write_gene(gtype, locus_num, gene_name, product, crispr_num, cont_loc,
               genome, cont_num, ecnum, inf2, strand, start, end, lstopenfile):
    """
    Write given gene to output file

    Parameters
    ----------
    gtype : str
        type of feature (CDS, tRNA, etc.)
    locus_num : str
        number of the locus given by prokka/prodigal
    gene_name : str
        gene name found by prokka/prodigal ("NA" if no gene name)
    product : str
        found by prokka/Prodigal, "NA" if no product (always the case for prodigal)
    crispr_num : int
        current crispr number. In prokka tbl, CRISPRs are not numbered, they all
        have the same name. We name them by adding a unique number to each CRISPR. If the current
        gene to add is a CRISPR, this number will be incremented and returned. If not, this same
        name will be returned.
    cont_loc : str
        'i' if the gene is inside a contig, 'b' if its on the border (first or last gene
        of the contig)
    genome : str
        genome name (spegenus.date.strain_num)
    cont_num : int
        contig number
    ecnum : str
        EC number, found by prokka, or "NA" if no EC number
    inf2 : str
        more information found by prokka, or "NA" if no more information
    strand : str
        C (complement) or D (direct)
    start : int
        start of gene in the contig
    end : int
        end of gene in the contig
    lstopenfile : _io.TextIOWrapper
        open file where lstinfo must be written

    Returns
    -------
    int :
        Current crispr number
    """
    # if last gene was a crispr
    if gtype == "repeat_region":
        gtype = "CRISPR"
        locus_num = "CRISPR" + str(crispr_num)
        gene_name = "crispr"
        product = "crispr-array"
        crispr_num += 1
    locus_name = "{}.{}{}_{}".format(genome, str(cont_num).zfill(4), cont_loc,
                                     str(locus_num).zfill(5))
    # If '|' character found in those fields, replace by '_' to avoid problems while parsing
    more_info = "| {} | {} | {}".format(product.replace("|", "_"),
                                        ecnum.replace("|", "_"),
                                        inf2.replace("|", "_"))
    lst_line = "\t".join([str(start), str(end), strand, gtype,
                          locus_name, gene_name, more_info])
    lstopenfile.write(lst_line + "\n")
    return crispr_num, lst_line


def write_header(lstline, outfile):
    """
    write header to output file. Header is generated from the lst line.

    Parameters
    ----------
    lstline : str
        current line of lst file
    outfile : _io.TextIOWrapper
        open file where header must be written-

    """
    start, end, _, _, name, gene_name, info = lstline.split("\t")
    size = int(end) - int(start) + 1
    towrite = " ".join([name, str(size), gene_name, info])
    outfile.write(">" + towrite + "\n")


def get_contig_name(genome, cont_num):
    """
    From a genome name and a contig number, write the header corresponding to this contig
    (for 'Replicons' files)
    """
    return ">{}.{}".format(genome, str(cont_num).zfill(4))

