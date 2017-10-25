#!/usr/bin/env python3
# coding: utf-8

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
    `<genome_name>.<i or b><contig_num>_<protein_num>`

@author gem
April 2017
"""

import os
import shutil
import logging
import logging.handlers
import progressbar
import glob
import multiprocessing
import threading
import genomeAPCAT.utils as utils


def format_genomes(genomes, results, res_path, prok_path, threads=1, quiet=False):
    """
    For all genomes which were annotated by prokka, reformat them
    in order to have, in 'res_path', the following folders:

    * LSTINFO: containing a .lst file for each genome, with all genes
    * Replicons: containing all multifasta files
    * Genes: containing 1 multi-fasta per genome, with all its genes in nuc
    * Proteins: containing 1 multi-fasta per genome, with all its proteins in aa

    Parameters
    ----------
    genomes : dict
        {genome: [name, gpath, size, nbcont, l90]}
    results : dict
        {genome: bool} True if prokka ran well, False otherwise
    res_path : str
        path to folder where the 4 directories must be created
    prok_path : str
        path to folder named "<genome_name>-prokkaRes" where all prokka
        results are saved.
    threads : int
        number of threads to use to while formatting genomes
    quiet : bool
        True if nothing must be sent to stderr/stdout, False otherwise

    Returns
    -------
    (skipped, skipped_format) : tuple

        * skipped : list of genomes skipped because they had a problem in prokka step
        * skipped_format : list of genomes skipped because they had a problem in format step
    """
    main_logger = logging.getLogger("qc_annote.ffunc")
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
    skipped = []  # list of genomes skipped: no format step run
    skipped_format = []  # List of genomes for which forat step had problems
    params = [(genome, name, gpath, prok_path, lst_dir, prot_dir, gene_dir,
               rep_dir, gff_dir, results, q)
              for genome, (name, gpath, _, _, _) in genomes.items()]
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
    for output in res:
        if output[0] == "bad_prokka":
            skipped.append(output[1])
        elif not output[0]:
            skipped_format.append(output[1])
    return skipped, skipped_format


def handle_genome(args):
    """
    For a given genome, check if it has been annotated (in results), if prokka ran without
    problems (result = True). In that case, format the genome and get the output to
    see if everything went ok.

    Parameters
    ----------
    args : tuple
        (genome, name, gpath, prok_path, lst_dir, prot_dir,\
         gene_dir, rep_dir, gff_dir, results, q)\
         with:

         * genome : original genome name
         * name : gembase name of the genome
         * gpath : path to the genome sequence which was given to prokka for annotation
         * prok_path : directory where prokka folders are saved
         * lst_dir : path to 'LSTINFO' folder
         * prot_dir : path to 'Proteins' folder
         * gene_dit : path to 'Genes' folder
         * rep_dir : path to 'Replicons' folder
         * gff_dir : path to 'gff3' folder
         * results : {genome_name: <True if genome was formatted correctly by prokka,\
         False otherwise>}
         * q : multiprocessing.managers.AutoProxy[Queue] queue to put logs during subprocess

    Returns
    -------
    (bool, str) :

        * True if genome was annotated as expected, False otherwise
        * genome name (used to get info from the pool.map_async)
    """
    (genome, name, gpath, prok_path, lst_dir, prot_dir,
     gene_dir, rep_dir, gff_dir, results, q) = args
    # Set logger for this process
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('format.handle_genome')
    # Handle genome
    if genome not in results:
        return "no_res", genome
    # if prokka did not run well for a genome, don't format it
    if not results[genome]:
        return "bad_prokka", genome
    ok_format = format_one_genome(gpath, name, prok_path, lst_dir,
                                  prot_dir, gene_dir, rep_dir, gff_dir, logger)
    return ok_format, genome


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
    prokka_dir = os.path.join(prok_path, os.path.basename(gpath) + "-prokkaRes")
    # Get .tbl file
    tblgenome = glob.glob(os.path.join(prokka_dir, "*.tbl"))[0]
    lstgenome = os.path.join(lst_dir, name + ".lst")
    # convert tbl to lst, and check that genome name in tbl is the same as the given genome name
    # If the user ran prokka, and then changed the genome names (for example, by changing
    # L90 threshold), and wants to format again its genomes with the new names, but without
    # running prokka again, we must change here the contig names (containing genome name).
    same_name = tbl2lst(tblgenome, lstgenome, name)

    # Create gff3 file for annotations
    prokka_gff = glob.glob(os.path.join(prokka_dir, "*.gff"))[0]
    gffgenome = os.path.join(gff_dir, name + ".gff")
    ok_gff = generate_gff(prokka_gff, gffgenome, lstgenome, logger)
    if not ok_gff:
        os.remove(gffgenome)
        os.remove(lstgenome)
        return False

    # create Genes file (and check no problem occurred with return code)
    ffngenome = glob.glob(os.path.join(prokka_dir, "*.ffn"))[0]
    gengenome = os.path.join(gene_dir, name + ".gen")
    ok_gene = create_gen(ffngenome, lstgenome, gengenome, logger)
    # If gene file not created because a problem occurred, return False:
    # format did not run for this genome
    if not ok_gene:
        os.remove(gffgenome)
        os.remove(lstgenome)
        return False
    # If gene file was created, create Proteins file
    faagenome = glob.glob(os.path.join(prokka_dir, "*.faa"))[0]
    prtgenome = os.path.join(prot_dir, name + ".prt")
    ok_prt = create_prt(faagenome, lstgenome, prtgenome, logger)

    # If protein file not created, return False: format did not run for this genome
    if not ok_prt:
        os.remove(gffgenome)
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
        utils.rename_genome_contigs(name, gpath, rep_file)
    return True


def generate_gff(prokka_gff, gffgenome, lstgenome, logger):
    """
    From the lstinfo file, generate a gff file.

    Parameters
    ----------
    prokka_gff : str
        prokka output gff file
    gffgenome : str
        name of the gff file which will be created as output
    lstgenome : str
        name of lst file already created for this genome
    logger : logging.Logger
        logger object to add log information

    Returns
    -------
    bool :
        True if conversion worked well, False otherwise
    """
    with open(prokka_gff) as prokf, open(lstgenome) as lstf, \
            open(gffgenome, "w") as gfff:
        # write header line
        line_gff = prokf.readline()
        gfff.write(line_gff)
        while line_gff.startswith("#"):
            line_gff = prokf.readline()
        line_lst = lstf.readline()
        handle_line_gff(line_lst, line_gff, gfff)
        for line_lst, line_gff in zip(lstf, prokf):
            try:
                handle_line_gff(line_lst, line_gff, gfff)
            except IndexError:
                logger.error("Problem with your gff file. '{}' is not a gff entry "
                             "line, whereas it should correspond "
                             "to '{}'".format(line_gff.strip(), line_lst.strip()))
                return False
    return True


def handle_line_gff(line_lst, line_gff, gfff):
    """
    Read a line from prokka gff output and a line from lstinfo, and get needed information

    Parameters
    ----------
    line_lst : str
        a string with a line of generated lst file
    line_gff : str
        a string with a line of the prokka gff file
    gfff : _io.TextIOWrapper
        new gff file open

    """
    start, end, strand, tp, name, gene, other = line_lst.strip().split("\t")
    _, product, ecnum, pred = other.split("|")
    soft = line_gff.strip().split("\t")[1]
    cont_num = name.split(".")[-1].split("_")[0][1:]
    contig_name = ".".join(name.split(".")[:3]) + "." + cont_num
    strands = {"D": "+", "C": "-"}
    att = "ID={}".format(name)
    if ecnum.strip() != "NA":
        att += ";eC_number={}".format(ecnum.strip())
    if gene != "NA":
        att += ";Name={0};gene={0}".format(gene)
    if pred.strip() != "NA":
        att += ";inference={}".format(pred.strip())
    att += ";locus_tag={}".format(name)
    if product.strip() != "NA":
        att += ";product={}".format(product.strip())
    infos = [contig_name, soft, tp, start, end, ".", strands[strand], ".", att]
    gfff.write("\t".join(infos) + "\n")


def create_gen(ffnseq, lstfile, genseq, logger):
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
    crispr_id = 1
    with open(ffnseq) as ffn, open(lstfile) as lst, open(genseq, "w") as gen:
        for line_ffn in ffn:
            # If line of sequence, write it as is, and go to next line
            if not line_ffn.startswith(">"):
                gen.write(line_ffn)
                continue
            lstline = lst.readline().strip()
            # Try to get gene ID. If does not work, look if it is a CRISPR in lstinfo
            test_gen_id = line_ffn.split()[0].split("_")[-1]
            if not test_gen_id.isdigit():
                # If it is a CRISPR in lstline, and header of ffn does not have a gene format,
                # then ffn contains the CRISPR sequence
                if lstline.strip().split()[3] == "CRISPR":
                    crispr_id_lst = int(lstline.split("\t")[4].split("_CRISPR")[-1])
                    print(lstline, crispr_id_lst)
                    if crispr_id == crispr_id_lst:
                        write_header(lstline, gen)
                        crispr_id += 1
                    else:
                        logger.error(("Problem with CRISPR numbers in {}. CRISPR {} in ffn is "
                                      "CRISPR num {}, whereas it is annotated as CRISPR num {} "
                                      "in lst file.").format(lstfile, line_ffn.strip(), crispr_id,
                                                             crispr_id_lst))
                        problem = True
                        break
                # It is not a CRISPR in lstline, and header of ffn does not have a gene format:
                # problem
                else:
                    logger.error(("Unknown header format {} in {}."
                                  "\nGen file will not be "
                                  "created.").format(line_ffn.strip(), ffnseq))
                    problem = True
                    break
            # If ffn contains a gene header, find its information in lst file
            else:
                gen_id = int(test_gen_id)
                # genID exists, ffn header is for a gene. Check that it corresponds to
                # information in lst file.
                id_lst = lstline.split("\t")[4].split("_")[-1]
                # if line in lst corresponds to a gene -> get gene ID.
                # Otherwise, genID = 0 (CRISPR line in lst)
                if id_lst.isdigit():
                    gen_id_lst = int(id_lst)
                else:
                    gen_id_lst = 0
                # in lst, find the same gene ID as in ffn
                # as they are ordered by increasing number, stop if ffn ID is higher than lst ID
                while gen_id > gen_id_lst:
                    lstline = lst.readline().strip()
                    id_lst = lstline.split("\t")[4].split("_")[-1]
                    # don't cast to int if info for a crispr
                    if id_lst.isdigit():
                        gen_id_lst = int(id_lst)
                # If it found the same gene ID, write info in gene file
                if gen_id == gen_id_lst:
                    write_header(lstline.strip(), gen)
                # If gene ID of ffn not found, write error message and stop
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


def create_prt(faaseq, lstfile, prtseq, logger):
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
    logger : logging.Logger
        log object to add information

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
                except Exception as err:
                    logger.error(("Unknown header format {} in {}. "
                                  "Error: {}\nPrt file not created "
                                  "from {}.").format(line.strip(), faaseq, err, faaseq))
                    problem = True
                    break
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
    write header to output file. Header is generated from the lst line.

    Parameters
    ----------
    lstline : str
        current line of lst file
    outfile : _io.TextIOWrapper
        open file where header must be written-

    """
    name = lstline.split("\t")[4]
    size = int(lstline.split("\t")[1]) - int(lstline.split("\t")[0]) + 1
    gene_name = lstline.split("\t")[5]
    info = lstline.split("\t")[6]
    towrite = " ".join([name, str(size), gene_name, info])
    outfile.write(">" + towrite + "\n")


def tbl2lst(tblfile, lstfile, genome):
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

        start end strand type locus gene_name | product | EC_number | inference 2

    with the same types as prokka file, and strain is C (complement) or D (direct)
    locus is: `<genome_name>.<i or b><contig_num>_<protein_num>`

    Parameters
    ----------
    tblfile : str
        name of prokka output tbl file to read
    lstfile : str
        name of lst file to generate
    genome : str
        genome name (gembase format)

    Returns
    -------
    bool :
        True if genome name used in lstfile and prokka tblfile are the same, False otherwise
    """
    same_name = True
    crispr_num, cont_num = 1, 0
    start, end = -1, -1
    strand, gtype, locus_num = [""] * 3
    gene_name, product, ecnum, inf2 = ["NA"] * 4
    cont_loc = "b"
    strain = ""
    with open(tblfile) as tblf, open(lstfile, "w") as lstf:
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
            if "inference" in elems[0] and "ab initio prediction:Prodigal" not in line and\
               len(elems) == 2:
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

    Parameters
    ----------
    gtype : str
        type of feature (CDS, tRNA, etc.)
    locus_num : str
        number of locus given by prokka
    gene_name : str
        gene name found by prokka ("NA" if no gene name)
    product : str
        found by prokka, "NA" if no product
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
