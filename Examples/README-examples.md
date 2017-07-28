# Examples

## Genomes description

We provide 4 different fictive genomes:
- genome1: 
    * 1 file, with 1 contig (fasta file)
    * example of complete genome of a bacteria without plasmid. However, there is a stretch of 11 "N" in its sequence. It will then be split into 2 contigs with the default parameters.
- genome2: 
    * 1 file, with 3 contigs (multifasta)
    * example of a draft genome
- genome3:
    * 2 files
        - genome3-chromo.fst with 1 contig (fasta)
        - genome3-plasmid.fst with 1 contig (fasta)
    * example of a complete genome with 1 plasmid. Note that there is a stretch of 6 "N" in the plasmid sequence, so it will be split into 2 contigs with the default parameters.
- genome4:
    * 1 file, with 1 contigs (fasta)
    * example of a complete genome (no stretch of "N" in the sequence)

**Note:** the provided genomic sequences are taken from real genes, but then modified in order to have an example showing different situations. Hence, the examples results should not be interpreted biologically!


## Annotate step

To annotate genomes, you need to provide a list of genomes to annotate, in a text file. An example, corresponding to the genomes in `Examples/genomes` is provided in `Examples/input_files/list_genomes.lst`. You can see that each line contains the filename(s) corresponding to one genome. 

For genome2 and genome4, we specify a species name (GEN2 and GEN4 respectively). For genome1 and genome3, we do not specify a name: the default name given to the program (`-n <name>`) will be used to name them.  
For genome3 and genome4, we specify a date at which they were downloaded. For genome1 and genome2, we do not specify a date: if a default date is given to the program (`--date <date>`), it will be used for those 2 genomes. Otherwise, the current date will be used.

### Quality control

If you just want to do quality control on the dataset, type:

    genomeAPCAT annotate -d Examples/genomes -r my_results Examples/input_files/list_genomes.lst -Q  

This will create a folder `my_results`, containing:
- `QC_L90-list_genomes.png`: histogram of the L90 values of all genomes
- `QC_nb-contigs-list_genomes.png`: histogram of number of contigs in all genomes
- `discarded-list_genomes.lst`: should be empty. The default limits are L90 <= 100 and #contigs <= 999. In the png files, we can see that we are very far from those limits, so, no genome is discarded.
- `genomeAPCAT-annotate_list_genomes.log`: log file. See information on what happened during the run: traceback of stdout.
- `genomeAPCAT-annotate_list_genomes.log.err`: log file but only with Warnings and errors. If it is empty, everything went well!
- `genomeAPCAT-annotate_list_genomes.log.details`: same as `.log` file, but with more detailed information (for example, time of start/end of annotation of each individual genome). This file can be quite big if you have a lot of genomes.
- `info-genomes-list_genomes.lst`: file with information on each genome: size, number of contigs and L90.
- `tmp_files` folder: containing your genomic sequences, split at each stretch of at least 5 'N'. 
- In your `Examples/genomes` folder, you should now also have a file called `genome3-chromo.fst-all.fna`, containing the concatenation of the 2 files corresponding to `genome3`.

In the `QC_L90-list_genomes.png`, we can see that all genomes have a L90 <= 4.
Similarly, in `QC_nb-contigs-list_genomes.png`, we can see that all genomes have #contigs <= 4.

### Annotation 

Now that you have seen the distribution of L90 and #contig values in your genomes, and decided which limits you want to use (if you do not want to use the default ones), you can annotate the genomes which are under those limits:

    genomeAPCAT annotate -d Examples/genomes -r my_results Examples/input_files/list_genomes.lst -n GENO --l90 3 --nbcont 10

Here, we put the L90 limit to 7, which should lead to the removal of 1 genome (genome2, according to the `QC_L90-list_genomes.png` and `info-genomes-list_genomes.lst` files). The nbcont limit should not remove any genome. We put this limit just to show how the program works with your own limits, but they do not have any significance here, as a genome with L90 = 4 is not a bad quality genome!  
We also have to add an option, `-n <name>` to specify the default name to give to the genomes if it is not specified in the list file (here, for genome1 and genome3). We here choose `GENO`, as short for 'genome'... Choose more appropriate names!

In your `my_results` directory, you should now have:
- the png files as previously: the distribution of values is the same, but the new limit for L90 now appears as a red line, as it is in the same range as the values.
- `discarded-list_genomes.lst`: there is now one genome discarded: genome2. In this file, we can indeed see that its L90 is 4 (higher than the limit).
- `LSTINFO-list_genomes.lst`: get information on the 3 genomes which were annotated. You can check that 
    - genome1 was named using 'GENO' (default name), and the current date (MMYY)
    - genome3 was named using the default value 'GENO', and the date specified in list file (1216)
    - genome4 was named using 'GEN4' (specified in list file) and the date specified in list file (1111)
- log files as previously. Check in the `.log.err` file that no error occurred. Note that if you used the same output directory as for the previous step, and did not remove the log files, the new ones do not erase the existing ones: they now have a timestamp corresponding to the time/date when you launched this annotation step.
- in `tmp_files`, you still have the 'split5N' genomic sequence files, as well as prokka result folders.
- You have 4 new folders: `Replicons`, `LSTINFO`, `Genes`, `Proteins` each one containing 3 files (1 per genome) with your results.

## PanGenome step

To do a pangenome, you need to provide the list of genomes to consider, with 1 genome per line. Only the first column (genome name) will be considered, but you can use a file containing other columns...such as the one you already have in the result folder of annotation step: `LSTINFO-list_genomes.lst`! However, of course, if you want to do a pangenome of less genomes than the ones you annotated, you are free to create a new file with the genomes you want!  
Here, we are doing a pangenome of the 3 genomes annotated before. Here is the command line:

    genomeAPCAT pangenome -l my_results/LSTINFO-list_genomes.lst -n GENO3 -d my_results/Proteins -i 0.8 -o my_results/pangenome 

With:
- `-l`: the provided file (`input_for_pangenome.txt`) is the same as what you should have in `my_results/LSTINFO-list_genomes.lst`.
- `-n`: we name this dataset 'GENO3' (for 3 genomes of species 'GENO')
- `-d`: path to the folder containing all protein files
- `-i`: we want a pangenome with 80% identity
- `-o`: put all result and temporary files to this directory

In your `my_results/pangenome` folder, you should have your pangenome in `PanGenome-GENO3.All.prt-clust-0.8-model1_<date>.tsv.lst`. It contains 1 line per family. The first column is the family number, and others are all family members. You also have the qualitative (`.quali.txt`) and quantitative (`.quanti.txt`) matrix of this pangenome, as well as a summary file (`.summary.txt`).

You should also have other files that could help you if you need to investigate the results:
- `tmp_GENO3.All.prt<more_info>` folder: all temporary files used by MMseqs2 to cluster your proteins.
- `GENO3.All.prt-msDB*`: 5 files (alone, .index, .lookup, _h, _h.index) corresponding to the protein databank, in the format used by MMseqs2.
- `GENO3.All.prt-clust-0.8-mode1_<date>*`: 3 files (alone, .index, .tsv) generated by MMseqs2 corresponding to the clustering of your proteins
- `genomeAPCAT-pangenome_GENO3.log*`: the 3 log files as in the annotation step.
- `mmseq_GENO3.All.prt_0.8-mode1_<date>.log`: MMseqs2 log file.
- `PanGenome-GENO3.All.prt-clust-0.8-mode1-th10_2017-06-28_15-19-12.tsv.lst.bin`: binary file containing Python objects corresponding to the pangenome. Used by the program to do calculations faster the next time it needs this information (to generate Core or Persistent genome for example).

In your `my_results/Proteins` folder, you should have a new file, `GENO3.All.prt`, containing all proteins of the 3 genomes.

If you used the same dataset and parameters as in this README file, you should get a pangenome with 14 families.

## Core/Persistent Genome step

The core genome is inferred from the PanGenome. So, the only required file is your pangenome, obtained at last step. By default, it will generate a CoreGenome. Here is the command line to obtain the CoreGenome of our dataset:

    genomeAPCAT corepers -p my_results/pangenome/PanGenome-GENO3.All.prt-clust-0.8-model1_<date>.tsv.lst

**Replace `<date>` by your real filename**

The provided file (`Examples/input_files/input_for_corepers.txt`) is the same as the PanGenome you should get from the previous step (`my_results/pangenome/PanGenome-GENO3.All.prt-clust-0.8-model1_<date>.tsv.lst`), except for the date of genome1.fst, which corresponds to your current date.

You now have your Core or Persistent genome in `my_results/pangenome/PersGenome_<pangenome-filename>_1.lst`. With `_1` meaning that you asked for 100% of genomes present in each family.  
If you used the same dataset and parameters as in this README file, you should get a coregenome with 6 families.

If you want a Persistent Genome, specify the required options (minimum percentage of genomes in a family to be considered as persistent, allowing or not multi/mixed families...). You can see all options available by typing:
    
    genomeAPCAT corepers -h

## Alignment step

You can then do an alignment of all the proteins of each persistent family. For example, to align the 6 core families found in the previous step:

    genomeAPCAT align -c my_results/pangenome/PersGenome_<pangenome-filename>_1.lst -l my_results/LSTINFO-list_genomes.lst -n GENO3-1 -d my_results -o my_results/Phylogeny

**Replace `PersGenome_<pangenome-filename>_1.lst` by your real persistent genome filename**

with:
- `-c`: the core/persistent genome file generated in previous step
- `-l`: list of genomes in your dataset (generated by annotation step)
- `-n`: we name this dataset 'GENO3-1' (for 3 genomes of species 'GENO', and 100% of the genomes in each family)
- `-d`: path to the folder containing the directories 'Proteins' and 'Genes'
- `-o`: put all result and temporary files to this directory

In your output directory `my_results/Phylogeny`, you will find:
- `genomeAPCAT-align_GENO3-1.log*`: the 3 log files as in the other steps.
- a folder `List-GENO3-1`: contains, for each genome, the list of persistent proteins (that must be extracted to align them).
- a folder `Align-GENO3-1`: contains
    + for each family:
        + a file `GENO3-1-current.<fam>.gen` with all genes extracted
        + a file `GENO3-1-current.<fam>.prt` with all proteins extracted
        + a file `GENO3-1-current.<fam>.miss.lst` with the list of genomes not present in the family
    + `GENO3-1-complete.cat.aln` concatenation of all family alignments
- a folder `Phylo-GENO3-1`: contains `GENO3-1.grp.aln`, the alignment of all families then grouped by genome. This is the file you will need to infer a phylogenetic tree.

Again, see other options available by typing:
    
    genomeAPCAT align -h

## Tree step

You can infer a phylogenetic tree from the alignment of persistent families. By default, it uses FastTree to infer the phylogenetic tree, with a GTR DNA substitution model, and no bootstrap. To run this on the alignment generated by the previous step using 5 threads, use:

    genomeAPCAT tree -a my_results/Phylogeny/Phylo-GENO3-1/GENO3-1.grp.aln --threads 5

If you want to use FastME instead of FastTree, with all default parameters (TN93 substitution model, no bootstrap), using 5 threads, run:

    genomeAPCAT tree -a my_results/Phylogeny/Phylo-GENO3-1/GENO3-1.grp.aln -s fastme --threads 5

See `genomeAPCAT tree -h` to have an overview of all options available.

