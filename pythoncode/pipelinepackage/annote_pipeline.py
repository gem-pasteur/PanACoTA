#!/usr/bin/env python
# coding: utf-8

"""
Pipeline to annotate genomes. Steps are:
- optional: find stretches of at least 'n' N (default 5), and cut into a new contig at this stretch
- for each genome, calc L90 and number of contigs (after cut at N stretches if used)
- keep only genomes with:
    - L90 < x (default = 100)
    - #contig < y (default = 1000)
- rename those genomes and their contigs, with strain name increasing with quality (L90 and
 #contig)
- annotate kept genomes with prokka
- gembase format
- find essential genes and give a distribution graph/table so that the user can choose which
genomes he wants to remove on the 'essential genes' criteria

Input:
- list_file: list of genome filenames to annotate. Each genome is in a multi-fasta file. This file
contains 1 genome filename per line. It can contain a second column, with the species in 4 letters.
Genomes without this 2nd column will be renamed with the name given in species
- species: with 4 letters to rename genomes (except those whose species name is precised
in 2nd column of list file)
- dbpath: path to folder containing all multi-fasta sequences of genomes
- respath: path to folder where outputs must be saved (folders Genes, Replicons, Proteins,
LSTINFO and LSTINFO_dataset.lst file)
- threads: number of threads that can be used (default 1)

Output:
- In your given respath, you will find 4 folders: LSTINFO (information on each genome, with gene annotations), Genes (nuc. gene sequences), Proteins (aa proteins sequences), Replicons (input sequences but with formatted headers).
- In the database, folders with prokka results will be created for each input genome. If errors are generated during prokka step, you can look at the log file to see what was wrong.
- In your given respath, a file called log-<LSTINFO_file>-<current_date>.out will be generated. You can find there all logs: problems during annotation (hence no formatting step ran), and problems during formatting step. All steps start with a '*', and problems start with a ' - '.
- In your given respath, a file called err-<LSTINFO_file>-<current_date>.err will be generated, containing information on errors occured. If this file is empty, then annotation and formatting steps finished without any problem for all genomes.
- In your given respath, you will find a file called <LSTINFO_file>.lst with information on all
genomes: gembase_name, original_name, genome_size, L90, nb_contigs


@author gem
April 2017
"""

def main(list_file, db_path, res_path, name, l90, nbcont, cutn, threads):
    """

    """
    if cutn > 0:
        # run cut_n
        pass
    # calc L90, nbcontig and genome size for each genome
    # list of genomes kept, ordered by quality
    # rename contigs and generate LSTINFO for genomes kept



def parse():
    """
    Method to create a parser for command-line options
    """
    import argparse
    def gen_name(param):
        if len(gen_name) != 4:
            msg = ("The genome name must contain 4 characters. For example, this name can "
                   " correspond to the 2 first letters of genus, and 2 first letters of "
                   "species, e.g. ESCO for Escherichia Coli.")
            raise argparse.ArgumentTypeError(msg)
        return param

    parser = argparse.ArgumentParser(description=("Annotate all genomes"))
    # Create command-line parser for all options and arguments to give
    parser.add_argument(dest="list_file",
                        help=("File containing the list of genome filenames to annotate (1 genome"
                              " per line). Each genome is in multi-fasta format. You can "
                              "specify the species name (4 letters) you want to give to each "
                              "genome by adding it after the genome filename, separated by a "
                              "space. If not given, the species name will be the one given in "
                              "'species' argument. "))
    parser.add_argument("-d", dest="db_path", required=True,
                        help=("Path to folder containing all multifasta genome files"))
    parser.add_argument("-r", dest="res_path", required=True,
                        help=("Path to folder where output annotated genomes must be saved"))
    parser.add_argument("-s", dest="name", required=True, type=gen_name,
                        help=("Choose a name for your annotated genomes. This name should contain 4 letters. Generally, they correspond to the 2 first letters of genus, and 2 first letters of species, e.g. ESCO for Escherichia Coli."))
    parser.add_argument("--l90", dest="l90", type=int, default=100,
                        help=("Maximum value of L90 allowed to keep a genome. Default is 100."))
    parser.add_argument("--nbcont", dest="nbcont", type=int, default=999,
                        help=("Maximum number of contigs allowed to keep a genome. "
                              "Default is 999."))
    parser.add_argument("--cutN", dest="cutn", type=int, default=5,
                        help=("By default, each genome will be cut into new contigs at each "
                              "stretch of at least 5 'N' in its sequence. If you don't want to "
                              "cut genomes into new contigs when there are stretches of 'N', "
                              "put 0 to this option. If you want to cut from a different number "
                              "of 'N' stretches, put this value to this option."))
    parser.add_argument("--threads", dest="threads", type=int, default=1,
                        help=("Specify how many threads can be used (default=1)"))
    args = parser.parse_args()
    # if args.multi and args.mixed:
    #     parser.error("-M and -X options cannot be activated together. Choose if you want to:\n"
    #                  "- allow several members in any number of genomes of a family (-M)\n"
    #                  "- allow several members in only '1-tol'% of the genomes of a family "
    #                  "(other 'tol'% genomes must have exactely 1 member)")
    # if args.mixed and args.tol==1:
    #     parser.error("You are asking for mixed families, while asking for 100% of the genomes of "
    #                  "a family to have exactly one member, which is not compatible. Do you want "
    #                  "to \n- lower the percentage of genomes required to have exactly "
    #                  "1 member (-t tol)\n- not allow mixed families (remove -X option)")
    return args

if __name__ == '__main__':
    OPTIONS = parse()
    main(OPTIONS.list_file, OPTIONS.db_path, OPTIONS.res_path, OPTIONS.name, OPTIONS.l90,
         OPTIONS.nbcont, OPTIONS.cutn, OPTIONS.threads)
