#!/usr/bin/env python3

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
Functions helping for doing quality control on genomes in order to eliminate
bad quality sequences, and then run Mash loops in order to discard too close genomes.

@author gem
August 2019
"""

import os
import sys
import glob
import logging
import progressbar
import scipy.sparse
from scipy.sparse import dok_matrix

from PanACoTA import utils
from PanACoTA.annotate_module import genome_seq_functions as gfunc

logger = logging.getLogger("prepare.filter")


def check_quality(species_linked, db_path, tmp_dir, max_l90, max_cont, cutn):
    """
    Do a quality control of all genomes in db_path

    Parameters
    ----------
    outdir : str
        directory for all results (refseq downloads, database init etc)
    species_linked : str
        given NCBI species with '_' instead of spaces, or NCBI taxID if species
        name not given
    dbpath : str
        directory to 'Database_init' containing all .fna files
    tmp_dir : str
        directory where all tmp files must be saved (files cut at each stretch of 'x' N)
    max_l90 : int
        max L90 value tolerated to keep a genome
    max_cont : int
        Max number of contigs tolerated to keep a genome
    cutn : int
        cut at each stretch of this number of 'N'. Don't cut if equal to 0

    Returns
    -------
        genomes : {genome_file: [genome_name, orig_path, path_to_seq_to_annotate, size,
                                 nbcont, l90]}
        no need for small name, we won't annotate genomes. genome_name is the same as filename
        but without extension
    """
    # Check database folder exists
    if not os.path.isdir(db_path):
        logger.error(f"{db_path} does not exist.")
        sys.exit(1)
    if not os.path.isdir(tmp_dir):
        logger.error(f"{tmp_dir} does not exist.")
        sys.exit(1)
    # Get all genome filenames
    all_genomes = os.listdir(db_path)
    if len(all_genomes) == 0:
        logger.error(f"There is no genome in {db_path}.")
        sys.exit(1)
    # Get name of genomes without extension
    genomes = {g:[os.path.splitext(g)[0]] for g in all_genomes}
    logger.info("Total number of genomes for {}: {}".format(species_linked, len(all_genomes)))

    # cut at stretches of 'N' if asked, and get L90, nbcontig, size for all genomes
    # -> {genome_file: [genome_g, orig_path, to_annotate_path, size, nbcont, l90]}
    gfunc.analyse_all_genomes(genomes, db_path, tmp_dir, cutn, "prepare", logger, quiet=False)
    return genomes

def sort_genomes_minhash(genomes, max_l90, max_cont):
    """
    Sort genomes:
    - draft genomes, sorted by L90 and then nb_contigs

    Parameters
    ----------

    genomes : {genome_file: [genome_name, orig_name, path_to_seq_to_annotate, size,
                             nbcont, l90]}
    max_l90 : int
        max L90 value tolerated to keep a genome
    max_cont : int
        Max number of contigs tolerated to keep a genome

    Returns
    -------
    sorted_genomes: list of 'genome_file' for all genomes kept (L90 and nbcont ok),
    ordered by decreasing quality
    """
    logger.info("Sorting all {} genomes by quality".format(len(genomes)))
    sorted_genomes = []
    # Add all genomes, sorted by L90 and then nb_cont, if L90 <= 100 and nb_cont <= 999
    nbcomp = 0
    nb_disc = 0
    for ginfo in sorted(genomes.items(), key=utils.sort_genomes_l90_nbcont):
        gname, (_, _, _, _, nb_cont, L90) = ginfo
        if L90 > max_l90 or nb_cont > max_cont:
            nb_disc += 1
        else:
            sorted_genomes.append(gname)
    logger.info(f"{len(sorted_genomes)} genomes after quality control ({nb_disc} "
                 "discarded)")
    return sorted_genomes


def iterative_mash(sorted_genomes, genomes, outdir, species_linked, min_dist, max_dist,
                   threads, quiet):
    """
    Run mash all vs all, to get all pairwise distances.
    Then, take the first genome of the list, and remove those for which the distance to it
    is not between 1e-4 and 0.06. Restart with the next genome kept in the list, and so on
    until the last genome.

    Parameters
    ----------
    sorted_genomes: list
        list of 'genome_file' for all genomes kept (L90 and nbcont ok)
    genomes : dict
        {genome_file: [genome_name, orig_name, path_to_seq_to_annotate, size, nbcont, l90]}
    outdir : str
        path to directory where all results are saved
    species_linked : str
        species name if given, otherwise species taxID
    min_dist : float
        lower limit of distance between 2 genomes to keep them
    max_dist : float
        max limit of distance between 2 genomes to keep them
    threads :
        max number of threads to use
    quiet : bool
        True if nothing must be sent to stdout/stderr, False otherwise

    Returns
    -------

    genomes_removed : dict
        {genome_name: [ref_name, dist]} genome against which 'genome_name' is removed, and corresponding distance (justifying removal)
    """
    logger.info("Starting filtering steps according to distance between genomes.")
    # Run mash all vs all
    mash_dir = os.path.join(outdir, "mash_files")
    os.makedirs(mash_dir, exist_ok=True)
    # List of genome sequences to compare
    list_reps = os.path.join(mash_dir, f"list-to-sketch-{species_linked}.txt")
    # Mash logfile
    mash_log = os.path.join(mash_dir, f"mash-all-{species_linked}.log")
    # Binary file generated by minhash sketch (index of all sequences)
    out_msh = os.path.join(mash_dir, f"all-genomes-{species_linked}")
    # Matrix with pairwise distances between all genomes
    matrix = os.path.join(mash_dir, f"matrix-all-genomes-{species_linked}.txt")
    # Binary file to save matrix of pairwise distances
    sparse_mat = os.path.join(mash_dir, f"matrix-all-genomes-{species_linked}.npz")

    # Sketch genomes
    sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads)

    # Compute pairwise distances
    compare_all(out_msh, matrix, sparse_mat, mash_log, threads)
    # Iteratively discard genomes
    # List of genomes to compare to the next ones until a limit value is reached.
    # genomes ordered by decreasing L90/nbcont (used to pop elements in comparing step)
    to_try = sorted_genomes[::-1]
    # Put list of genomes removed by mash comparison, and why
    # (out of limits distance with which genome)
    genomes_removed = {}  # {genome: [compared_with, dist]}
    nbgen = len(to_try)

    # Get link between genome_file (genomes key) and place in sorted_genomes
    corresp_file = {genome: num for num, genome in enumerate(sorted_genomes)}
    # Read matrix (from npz file if existing, otherwise from txt file)
    if os.path.exists(sparse_mat):
        logger.info(f"Loading matrix contained in {sparse_mat}")
        # convert matrix returned by load_npz (coo format, as saved) to dok format
        mat_sp = scipy.sparse.load_npz(sparse_mat).todok()
    # Read matrix txt file generated by minhash, and save this python object matrix to a npz file.
    else:
        logger.info("Reading matrix from txt file generated by Mash.")
        mat_sp = read_matrix(genomes, sorted_genomes, matrix)
        # corresp_file = {genome_file : num of genome in sorted_genomes}
        logger.info("Saving matrix to npz file to be loaded quicker if needed later")
        # Convert dok_matrix to coo format, as dok format is not allowed by save_npz
        coo_mat = mat_sp.tocoo()
        scipy.sparse.save_npz(sparse_mat, coo_mat)

    # Iteratively discard genomes too close or too far
    logger.info("Starting iterative discarding steps")
    if not quiet:
        widgets = ['Genomes compared: ',
                   progressbar.Bar(marker='█', left='', right='', fill=' '), ' ',
                   progressbar.Counter(), "/{}".format(nbgen), ' ',
                   progressbar.Timer(), ' - '
                  ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=len(to_try), term_width=79).start()
        done = 0

    while len(to_try) > 1:
        mash_step(to_try, corresp_file, mat_sp, genomes_removed, min_dist, max_dist)
        if not quiet:
            done = nbgen - len(to_try)
            bar.update(done)
    if not quiet:
        bar.finish()
    logger.info("Final number of genomes in dataset: {}".format(nbgen - len(genomes_removed)))
    return genomes_removed


def sketch_all(genomes, sorted_genomes, outdir, list_reps, out_msh, mash_log, threads):
    """
    Sketch all genomes to a combined archive.

    Parameters
    ----------
    genomes : dict
        {genome_file: [genome_name, orig_name, path_to_seq_to_annotate, size, nbcont, l90]}
    sorted_genomes: list
        list of 'genome_file' for all genomes kept (L90 and nbcont ok), ordered by
        decreasing quality
    outdir : str
        path to directory where all results are saved
    list_reps : str
        file with list of genomes to sketch. File will be emptied if it contain something, and
        filled with the informations from 'genomes'.
    out_msh : str
        output of mash
    mash_log : str
        mash logfile
    threads :
        max number of threads to use

    Returns
    -------

    return value (0 if OK, 1 if error)

    """
    # If given outdir does not exist, close it
    if not os.path.isdir(outdir):
        logger.error(f"Your output directory '{outdir}' does not exist.")
        sys.exit(1)
    # Empty list_reps file
    open(list_reps, "w").close()
    # Complete paths to genomes to compare: 'path_to_seq_to_annotate' = genome_file[2]
    file_paths = [genomes[g][2] for g in sorted_genomes]
    # Write list of genomes to compare to a file
    utils.write_list(file_paths, list_reps)
    # Sketch all genome sequences if not already done
    if os.path.isfile(out_msh + ".msh"):
        logger.warning(f"Mash sketch file {out_msh}.msh already exists. PanACoTA will "
                        "use it for next step.")
        os.remove(list_reps)
        return 0
    logger.info("Sketching all genomes...")
    cmd_sketch = f"mash sketch -o {out_msh} -p {threads} -l {list_reps} -s 1e4"
    logger.details(cmd_sketch)
    error_sketch = (f"Error while trying to sketch {len(sorted_genomes)} genomes to combined "
                    "archive. Maybe some genome sequences in "
                    "'tmp_files' are missing! Check logfile: "
                    f"{mash_log}")

    outf = open(mash_log, "w")
    utils.run_cmd(cmd_sketch, error_sketch, eof=True, stdout=outf, stderr=outf, logger=logger)
    outf.close()
    return 0


def compare_all(out_msh, matrix, npz_matrix, mash_log, threads):
    """
    Comparing all pairwise genomes that are already been sketched in the given file.

    Parameters
    ----------
    out_msh : str
        output of mash
    matrix : str
        File to put generated matrix of pairwise distances between all genomes
    npz_matrix : str
        matrix of pairwise distances saved in a binary file
    mash_log : str
        mash logfile
    threads :
        max number of threads to use

    Returns
    -------

    return code
    """
    # txt matrix already exists
    if os.path.isfile(matrix):
        logger.warning("Matrix file {} already exists. The program will use this distance matrix "
                       "to filter all genomes according to their distances.".format(matrix))
        return 0
    # npz matrix already exists
    if os.path.isfile(npz_matrix):
        logger.warning("Matrix file {} already exists. The program will use this distance matrix "
                       "to filter all genomes according to their distances.".format(matrix))
        return 0
    logger.info("Computing pairwise distances between all genomes")
    cmd_dist = f"mash dist -p {threads} {out_msh}.msh {out_msh}.msh"
    logger.details(cmd_dist)
    # Open matfile to write matrix inside
    matfile = open(matrix, "w")
    # Open mash log to add log of 'mash dist' to log of 'mash sketch'
    outf = open(mash_log, "a")
    error_dist = ("Error while trying to estimate pairwise distances between all genomes. "
                  f"See {mash_log}.")
    utils.run_cmd(cmd_dist, error_dist, eof=True, stdout=matfile, stderr=outf)
    outf.close()
    matfile.close()
    return 0


def mash_step(to_try, corresp, mat_sp, genomes_removed, min_dist, max_dist):
    """
    Prepare a mash run, with a given genome as reference, and others to compare to.

    Parameters
    ----------
    to_try : list
        list of genome_file (keys of 'genomes') to compare, ordered by decreasing L90/nbcont
    corresp : dict
        {genome_file : num_of_genome in sorted_genomes}
    mat_sp : scipy.sparse.dok.dok_matrix
        triangle matrix containing pairwise distance comparisons
    genomes_removed : dict
        {genome_file: [ref_name, dist]} genome against which 'genome_name' is removed, and
        corresponding distance (justifying removal)
    min_dist : float
        lower limit of distance between 2 genomes to keep them
    max_dist : float
        max limit of distance between 2 genomes to keep them

    Returns
    -------

    to_try is updated (reference element and all genomes not compatible with it are removed)
    genomes_removed is updated
    return code (0 if no problem)

    """
    # Get last element (which is the 'best' genome), and remove it from the list
    ref_name = to_try.pop()
    # Line of genome in mat_sp
    ref_num = corresp[ref_name]
    # Genomes (ordered by increasing L90/nbcont) to compare to the selected element (ref_name)
    others = to_try[::-1]

    # For each genome, compare its distance to reference genome 'ref_name'
    for gname in others:
        # Column of genome in mat_sp
        other_num = corresp[gname]
        # Get distance between the 2 genomes
        if ref_num < other_num:
            dist = mat_sp[ref_num, other_num]
        else:
            logger.warning("Should never happen as mat_sp is a triangle matrix!")
            dist = mat_sp[other_num, ref_num]
        # If distance not in the limits, remove genome from to_try and add to genomes_removed list
        if not min_dist <= dist <= max_dist:
            to_try.remove(gname)
            genomes_removed[gname] = [ref_name, dist]
    return 0


def read_matrix(genomes, sorted_genomes, matrix):
    """
    Read the matrix of pairwise distances between all genomes, and save it to a sparse
    matrix (only upper triangle).

    Parameters
    ----------
    genomes : dict
        {genome_file: [genome_name, orig_name, path_to_seq_to_annotate, size, nbcont, l90]}
    sorted_genomes: list
        list of 'genome_file' for all genomes kept (L90 and nbcont ok)
    matrix : str
        File containing the matrix of pairwise distances between all genomes

    Returns
    -------

    mat_sp : str
        python dok_matrix object
    """
    if not os.path.isfile(matrix):
        logger.error(f"Matrix file {matrix} does not exist. We cannot read it "
                     "and do the next steps. Program ending.")
        sys.exit(1)

    nbgen = len(sorted_genomes)
    corresp_abs = {genomes[genome][2]: num for num, genome in enumerate(sorted_genomes)}
    # Create square matrix with nbgen cols/lines. dok format is a 'Dictionary Of Keys'
    # -> writes (0, 1) value
    mat_sp = dok_matrix((nbgen, nbgen), dtype=float)
    # Write matrix values
    with open(matrix, "r") as matf:
        for line in matf:
            path1, path2, dist = line.split()[:3]
            num1 = corresp_abs[path1]
            num2 = corresp_abs[path2]
            # only in lower triangle (no duplicate)
            if num1 < num2:
                mat_sp[num1, num2] = float(dist)
            else:
                mat_sp[num2, num1] = float(dist)
    return mat_sp


def write_outputfiles(genomes, sorted_genomes, genomes_removed, outdir, gspecies, min_dist, max_dist):
    """
    Write the list of genomes kept in a file, 1 genome per line -> will be the input file for
    annotation and next steps
    Write discarded genomes to another file, with, for each line:
    - genome name
    - problem when compared with which other genome
    - distance to this other genomes

    Parameters
    ----------
    genomes : dict
        {genome_file: [genome_name, orig_name, path_to_seq_to_annotate, size, nbcont, l90]}
    sorted_genomes: list
        list of 'genome_file' for all genomes kept (L90 and nbcont ok)
    genomes_removed : dict
        {genome_name: [ref_name, dist]} genome against which 'genome_name' is removed, and corresponding distance (justifying removal)
    outdir : str
        directory where those list files must be created
    gspecies : str
        species name if given, otherwise species taxID
    min_dist : float
        lower limit of distance between 2 genomes to keep them
    max_dist : float
        upper limit of distance between 2 genomes to keep them

    Returns
    -------
    return code
    """
    if not os.path.isdir(outdir):
        logger.error(f"The given output directory ({outdir}) does not exist. We cannot "
                      "create output files there")
        sys.exit(1)
    list_file = os.path.join(outdir, f"LSTINFO-{gspecies}-filtered-{min_dist}_{max_dist}.txt")
    kept_genomes = []
    discard_file = os.path.join(outdir, f"discarded-by-minhash-{gspecies}-{min_dist}_{max_dist}.txt")

    # Get list of kept genomes and write them in list_file
    for genome in sorted_genomes:
        if genome not in genomes_removed:
            kept_genomes.append(genome)

    # Write 4 columns in LSTINFO file (genomes kept) :
    # path to file analyzed, genome size, nbcont and L90
    with open(list_file, "w") as lf:
        lf.write('to_annotate\tgsize\tnb_conts\tL90\n')
        # For each genome in kept_genomes, find required information on this it, using 'genomes'
        for g in kept_genomes:
            _, _, analyzed, size, nbcont, l90 = genomes[g]
            towrite = utils.list_to_str([analyzed, size, nbcont, l90], sep="\t")
            lf.write(towrite)

    # Write list of discarded genomes and why they are discarded
    with open(discard_file, "w") as disf:
        disf.write("genome_name\tproblem_compared_with\tdist\n")
        for genome, info in genomes_removed.items():
            disf.write(utils.list_to_str([genome] + info))
    logger.info(f"Final list of genomes in the dataset: {list_file}")
    logger.info(f"List of genomes discarded by minhash steps: {discard_file}")
    return list_file
