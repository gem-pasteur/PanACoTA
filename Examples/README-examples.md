# Examples

## Genomes description

We provide 4 different fictive genomes:
- ACBA: 
    * 1 file, with 1 contig (fasta file)
    * example of complete genome of a bacteria without plasmid
- ESCO: 
    * 3 files
        - 1 chromosome, split into 8 contigs (multifasta)
        - 2 plasmids, with 1 contig each (2 fasta files)
    * example of a draft genome, with 2 plasmids (which are completely assembled)
- KLPN:
    * 1 file, with 5 contigs (multifasta)
    * examples of a draft genome (no plasmid, or plasmids in draft state, mixed with chromosome)
- SAEN:
    * 1 file, with 4 contigs (multifasta)
    * same as KLPN

## Annotate step

To annotate genomes, you need to provide a list of genomes to annotate, in a text file. An example, corresponding to the genomes in `Examples/genomes` is provided in `Examples/input_files/list_genomes.lst`. You can see that each line contains the filename(s) corresponding to one genome. 

For ACBA and ESCO, we specify a name (ACBA and ESCO). For KLPN and SAEN, we do not specify a name: the default name given to the program will be used to name them.  
For ESCO and SAEN, we specify a date at which they were downloaded. For ACBA and KLPN, we do not specify a name: if a default date is given to the program, it will be used for those 2 genomes. Otherwise, the current date will be used.

### Quality control

If you just want to do quality control on the dataset, type:

    genomeAPCAT annotate -d Examples/genomes -r my_results Examples/input_files/list_genomes.lst -Q  

This will create a folder `my_results`, containing:
- `QC_L90-list_genomes.png`: histogram of the L90 values of all genomes
- `QC_nb-contigs-list_genomes.png`: histogram of number of contigs in all genomes
- `discarded-list_genomes.lst`: should be empty. The default limits are L90 <= 100 and #contigs <= 999. In the png files, we can see that we are very far from those limits, so, no genome is discarded.
- `genomeAPCAT-annotate_list_genomes.log` and `genomeAPCAT-annotate_list_genomes.log.err`: log files. You can see what happened during the run.
- `tmp_files` folder: containing your genomic sequences, split at each stretch of at least 5 'N'. 
- In your `Examples/genomes` folder, you should now also have a file called `ESCO.chromo.fna-all.fna`, containing the concatenation of the 3 files corresponding to `ESCO`.

In the `QC_L90-list_genomes.png`, we can see that all genomes have a L90 <= 8.
Similarly, in `QC_nb-contigs-list_genomes.png`, we can see that all genomes have #contigs <= 10.

### Annotation 

Now that you have seen the distribution of L90 and #contig values in your genomes, and decided which limits you want to use (if you do not want to use the default ones), you can annotate the genomes which are under those limits:

    genomeAPCAT annotate -d Examples/genomes -r my_results Examples/input_files/list_genomes.lst -n DEFN --l90 7 --nbcont 10

Here, we put the L90 limit to 7, which should lead to the removal of 1 genome (according to the `QC_L90-list_genomes.png` file). The nbcont limit should not remove any genome. We put this limit just to show how the program works with your own limits, but they do not have any significance here, as a genome with L90 = 8 is not a bad quality genome!  
We also have to add an option, `-n <name>` to specify the default name to give to the genomes if it is not specified in the list file (here, for KLPN and SAEN). We here choose `DEFN`, as 'default name'.

In your `my_results` directory, you should now have:
- the png files as previously: the distribution of values is the same, but the new limits now appear as red lines, as they are in the same range as the values.
- `discarded-list_genomes.lst`: there is now one genome discarded: KLPN. In this file, we can indeed see that its L90 is 8 (higher than the limit).
- `LSTINFO-list_genomes.lst`: get information on the 3 genomes which were annotated. You can check that 
    - ACBA was named using 'ACBA' (specified in list file), and the current date (MMYY)
    - SAEN was named using the default value DEFN, and the date specified in list file (1116)
    - ESCO was named using 'ESCO' (specified in list file) and the date specified in list file (1116)
- log files as previously. Check in the `.log.err` file that no error occurred.
- in `tmp_files`, you now have the intermediate genomic sequence files, as well as prokka result folders.
- You have 4 new folders: `Replicons`, `LSTINFO`, `Genes`, `Proteins` each one containing 3 files (1 per genome) with your results.

