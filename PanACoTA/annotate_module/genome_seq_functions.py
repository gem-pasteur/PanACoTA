#!/usr/bin/env python3
# coding: utf-8

"""
Functions to:

- analyse a genome (cut at stretch of N if asked, calc L90, nb cont, size...)
- if genome cut by stretch of N, write the new sequence to new file
- rename genome contigs and write new sequence to new file


@author gem
April 2017
"""
import os
import re
import numpy as np
import logging
import progressbar

from PanACoTA import utils

logger = logging.getLogger("qc_annotate.gseq")


def analyse_all_genomes(genomes, dbpath, tmp_path, nbn, prodigal_only, logger, quiet=False):
    """

    Parameters
    ----------
    genomes : dict
        {genome: spegenus.date}
    dbpath : str
        path to folder containing genomes
    tmp_path : str
        path to put out files
    nbn : int
        minimum number of 'N' required to cut into a new contig
    prodigal_only : bool
        True if we annotate with prodigal, False if we annotate with prokka
    logger : logging.Logger
        logger object to write log information
    quiet : bool
        True if nothing must be written to stdout/stderr, False otherwise

    Returns
    -------
    dict
        {genome: [spegenus.date, path, size, nbcont, l90]}

    """
    cut = nbn > 0
    pat = 'N' * nbn + "+"
    nbgen = len(genomes)
    bar = None
    curnum = None
    if cut:
        logger.info(("Cutting genomes at each stretch of at least {} 'N', "
                     "and then, calculating genome size, number of contigs and L90.").format(nbn))
    else:
        logger.info("Calculating genome size, number of contigs, L90")
    if not quiet:
        # Create progressbar
        widgets = ['Analysis: ', progressbar.Bar(marker='█', left='', right=''),
                   ' ', progressbar.Counter(), "/{}".format(nbgen), ' (',
                   progressbar.Percentage(), ') - ', progressbar.Timer(), ' - ',
                   progressbar.ETA()
                   ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=nbgen, term_width=100).start()
        curnum = 1
    toremove = []
    # Analyse genomes 1 by 1
    for genome, name in genomes.items():
        # If not quiet option, show progress bar
        if not quiet:
            bar.update(curnum)
            curnum += 1
        # analyse genome, and check everything went well
        res = analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes, prodigal_only, logger)
        # Problem while analysing genome -> genome ignored
        if not res:
            toremove.append(genome)
    # If there are some genomes to remove (analysis failed), remove them from genomes dict.
    if toremove:
        for gen in toremove:
            del genomes[gen]
    if not quiet:
        bar.finish()


def analyse_genome(genome, dbpath, tmp_path, cut, pat, genomes, prodigal_only, logger):
    """
    Analyse given genome:

    - if cut is asked:

        - cut its contigs at each time that 'pat' is seen
        - save cut genome in new file

    - calculate genome size, L90, nb contigs and save it into genomes

    Parameters
    ----------
    genome : str
        given genome to analyse
    dbpath : str
        path to the folder containing the given genome sequence
    tmp_path : str
        path to folder where output files must be saved.
    cut : bool
        True if contigs must be cut, False otherwise
    pat : str
        pattern on which contigs must be cut. ex: "NNNNN"
    genomes : dict
        {genome: [spegenus.date]} as input, and will be changed to\
         -> {genome: [spegenus.date, path, gsize, nbcont, L90]}
    prodigal_only : bool
        True if we annotate with prodigal, False if we annotate with prokka
    logger : logging.Logger
        logger object to write log information

    Returns
    -------
    bool
        True if genome analysis went well, False otherwise
    """
    gpath, grespath = get_output_dir(prodigal_only, dbpath, tmp_path, genome, cut, pat)

    # Open original sequence file
    with open(gpath) as genf:
        # If a new file must be created (sequences cut, and/or annotating with prokka), open it
        gresf = None
        if grespath:
            gresf = open(grespath, "w")
        # Initialize variables
        cur_contig_name = "" # header text (less than 20 characters if prokka used)
        contig_sizes = {}  # {header text: size}
        cur_seq = "" # sequence
        num = 1 # Used to get unique contig names

        # Read each line of original sequence
        for line in genf:
            #### NEW CONTIG
            # Line corresponding to a new contig
            if line.startswith(">"):
                # If not first contig, write info to  output file (of needed)
                if cur_seq != "":
                    num = format_contig(cut, pat, cur_seq, cur_contig_name, contig_sizes,
                                        gresf, num, logger)
                    # If problem while formatting contig, return False -> genome ignored
                    if num == -1:
                        return False
                # Get contig name for next turn, and reset sequence
                # If prodigal, contig name is as given by original sequence
                if prodigal_only:
                    cur_contig_name = line.strip()
                # If prokka, contig name is 1st word, 1st 15 characters
                else:
                    cur_contig_name = line.split()[0][:15]
                # Initialize for next contig
                cur_seq = ""
            # #### SEQUENCE LINE
            # If sequence line, keep it, all in upper case
            else:
                # Add this line without \n to sequence (1 sequence per line)
                cur_seq += line.strip().upper()

        # LAST CONTIG
        if cur_contig_name != "":
            num = format_contig(cut, pat, cur_seq, cur_contig_name, contig_sizes, gresf,
                                num, logger)
            if num == -1:
                return False

    # GLOBAL INFORMATION
    nbcont = len(contig_sizes)
    gsize = sum(contig_sizes.values())
    l90 = calc_l90(contig_sizes)
    if not l90:
        logger.error("Problem with genome {}. Not possible to get L90".format(genome))
        return False
    else:
        if grespath:
            genomes[genome] += [grespath, gsize, nbcont, l90]
        else:
            genomes[genome] += [gpath, gsize, nbcont, l90]
    # If we wrote a new sequence file, close it
    if grespath:
        gresf.close()
    return True


def get_output_dir(prodigal_only, dbpath, tmp_path, genome, cut, pat):
    """
    Get output file to put sequence cut and/or sequence with shorter contigs (prokka)

    Parameters
    ----------
    prodigal_only : bool
        True if we annotate with prodigal, False if we annotate with prokka
    dbpath : str
        path to the folder containing the given genome sequence
    tmp_path : str
        path to folder where output files must be saved.
    genomes : dict
        {genome: [spegenus.date]} as input, and will be changed to
        -> {genome: [spegenus.date, path, gsize, nbcont, L90]}
    cut : bool
        True if contigs must be cut, False otherwise
    pat : str
        pattern on which contigs must be cut. ex: "NNNNN"

    Return
    ------
    grespath : str
        path to ouput file. None if no need to create new sequence file
    """
    # Define soft name (to put in new sequence filename)
    if prodigal_only:
        soft = "prodigal"
    else:
        soft = "prokka"
    # Path to sequence to analyze
    gpath = os.path.join(dbpath, genome)
    # genome file is in dbpath except if it was in several files in dbpath,
    # in which case it has been concatenated to a file in tmp_path
    if not os.path.isfile(gpath):
        gpath = os.path.join(tmp_path, genome)
    # New file create if needed. If not (prodigal and not cut), empty filename
    grespath = None
    # If user asks to cut at each 'pat', need to create a new sequence file,
    # whatever the annotation soft used
    if cut:
        new_file = genome + "_{}-split{}N.fna".format(soft, len(pat) - 1)
        grespath = os.path.join(tmp_path, new_file)
    # If user does not want to cut, but annotates with prokka, need a new file with headers shorter
    elif not prodigal_only:
        new_file = genome + "_{}-shorter-contigs.fna".format(soft, len(pat) - 1)
        grespath = os.path.join(tmp_path, new_file)
    # If no cut and using prodigal, just keep original sequence, no need to create new file.
    # Just check that contigs have different names
    return gpath, grespath


def format_contig(cut, pat, cur_seq, cur_contig_name, contig_sizes, gresf, num, logger):
    """
    Format given contig, and save to output file if needed

    - if cut: cut it and write each subsequence
    - if prokka: write new contig (15 first characters of 1st word + contig_num)
    - if prodigal (and no cut), just check that contig names are different

    Parameters
    ----------
    cut : bool
        True if contigs must be cut, False otherwise
    pat : str
        pattern on which contigs must be cut. ex: "NNNNN"
    cur_seq : str
        current sequence (aa)
    cur_contig_name : str
        name of current contig
    cont_sizes : dict
        {contig_name : sequence length}
    gresf : io.TextIOWrappe
        open file to write new sequence
    num : int
        current contig number
    logger : logging.Logger
        logger object to write log information

    Returns
    -------
    bool
        True if contig has been written without problem, False if any problem
    """
    # "CUT" if cut: need to cut at each 'pat' -> write new header + seq in new file
    if cut:
        # Cut sequence and write header + sequence to res file
        num = split_contig(pat, cur_seq, cur_contig_name, contig_sizes, gresf, num)
    # PROKKA User does not want to cut, but will annotate with prokka, so we still
    # have to create a new sequence file
    elif gresf:
        new_contig_name = "{}_{}\n".format(cur_contig_name, num)
        gresf.write(new_contig_name)
        gresf.write(cur_seq + "\n")
        contig_sizes[new_contig_name] = len(cur_seq)
        num += 1
    # PRODIGAL No cut, and prodigal used -> no new file created, but check
    # contig unique names
    else:
        if cur_contig_name in contig_sizes.keys():
            logger.error("{} contig name is used for several contigs. Please put "
                         "different names for each contig. This genome will be "
                         "ignored.".format(cur_contig_name))
            return -1
        else:
            contig_sizes[cur_contig_name] = len(cur_seq)
    return num


def split_contig(pat, whole_seq, cur_contig_name, contig_sizes, gresf, num):
    """
    Save the contig read just before into dicts and write it to sequence file.
    Contig name must be at most 20 characters (required by prokka)

    Parameters
    ----------
    pat : str
        pattern to split a contig
    whole_seq : str
        sequence of current contig, to save once split according to pat
    cur_contig_name : str
        name of current contig to save once split according to pat
    contig_sizes : dict
        {name: size} save cur_contig once split according to pat
    gresf : _io.TextIOWrapper
        file open in w mode to save the split sequence
    num : int
        current contig number.

    Returns
    -------
    int
        new contig number, after giving number(s) to the current contig
    """
    # split contig each time a stretch of at least nbn 'N' is found (pattern pat)
    cont_parts = re.split(pat, whole_seq)
    # save contig parts
    for seq in cont_parts:
        # Only save non empty contigs (with some patterns, it could arrive that
        # we get empty contigs, if 2 occurrences of the pattern are side by side).
        if len(seq) == 0:
            continue
        new_contig_name = "{}_{}\n".format(cur_contig_name, num)
        contig_sizes[cur_name] = len(seq)
        gresf.write(new_contig_name)
        gresf.write(seq + "\n")
        num += 1
    return num


def calc_l90(contig_sizes):
    """
    Calc L90 of a given genome

    Parameters
    ----------
    contig_sizes : dict
        {name: size}

    Returns
    -------
    None or int
        if L90 found, returns L90. Otherwise, returns nothing
    """
    gsize = sum(contig_sizes.values())
    sizes = [contig_sizes[cont] for cont in contig_sizes]
    cum_sizes = np.cumsum(sorted(sizes, reverse=True))
    lim = 0.9 * gsize
    for num, val in enumerate(cum_sizes):
        if val >= lim:
            return num + 1


def rename_all_genomes(genomes):
    """
    Sort kept genomes by L90 and then nb contigs.
    For each genome, assign a strain number, and rename all its contigs.

    Parameters
    ----------
    genomes : dict
        {genome: [name, path, gsize, nbcont, L90]} as input, and will become\
        {genome: [gembase_name, path, gsize, nbcont, L90]} at the end

    """
    logger.info("Renaming kept genomes according to their quality.")
    last_name = ""
    last_strain = 0
    # "SAEN.1015.{}".format(str(last_strain).zfill(5))
    for genome, [name, _, _, _, _] in sorted(genomes.items(), key=sort_genomes):
        if last_name != name.split(".")[0]:
            last_strain = 1
            last_name = name.split(".")[0]
        else:
            last_strain += 1
        gembase_name = ".".join([name, str(last_strain).zfill(5)])
        genomes[genome][0] = gembase_name


def sort_genomes(x):
    """
    Sort all genomes with the following criteria:

    - sort by species (x[1][0] is species.date)
    - for each species, sort by l90
    - for same l90, sort by nb contigs

    Parameters
    ----------
    x : [[]]
        [genome_name, [species.date, path, gsize, nbcont, L90]]

    Returns
    -------
    tuple
        information on species, l90 and nb_contigs
    """
    return x[1][0].split(".")[0], x[1][-1], x[1][-2]


def plot_distributions(genomes, res_path, listfile_base, l90, nbconts):
    """
    Plot distributions of L90 and nbcontig values.

    Parameters
    ----------
    genomes : dict
        {genome: [name, path, size, nbcont, l90]}
    res_path : str
        path to put all output files
    listfile_base : str
        name of list file
    l90 : int
        L90 threshold
    nbconts : int
        nb contigs threshold

    Returns
    -------
    (l90_vals, nbcont_vals, dist1, dist2) :

    - l90_vals : list of l90 values for all genomes
    - nbcont_vals : list of nbcontigs for all genomes
    - dist1 : matplotlib figure of distribution of L90 values
    - dist2 : matplotlib figure of distribution of nb contigs values

    """
    """
    genomes:
    res_path:
    listfile_base:
    l90: max value of l90 allowed
    nbconts: max value of nb contigs allowed
    """
    logger.info("Generating distribution of L90 and #contigs graphs.")
    l90_vals = [val for _, (_, _, _, _, val) in genomes.items()]
    outl90 = os.path.join(res_path, "QC_L90-" + listfile_base + ".png")
    nbcont_vals = [val for _, (_, _, _, val, _) in genomes.items()]
    outnbcont = os.path.join(res_path, "QC_nb-contigs-" + listfile_base + ".png")
    dist1 = utils.plot_distr(l90_vals, l90, "L90 distribution for all genomes",
                             "max L90 =")
    dist2 = utils.plot_distr(nbcont_vals, nbconts,
                             "Distribution of number of contigs among all genomes",
                             "max #contigs =")
    dist1.savefig(outl90)
    dist2.savefig(outnbcont)
    return l90_vals, nbcont_vals, dist1, dist2
