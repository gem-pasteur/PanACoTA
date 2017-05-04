# pipeline_annotation

[![build status](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/build.svg)](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/commits/master)
[![coverage report](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/coverage.svg)](http://aperrin.pages.pasteur.fr/pipeline_annotation/htmlcov)

Annotate genomes and format them to gembase format.  

## pipeline installation (final mode)

Dependencies:
 python3 (and pip)
 prokka

To install the package *pipelinepackage*, and all its dependances, from the root directory, just do::

    sudo pip3 install .

You will then be able to use the package from any directory in your computer,
just as any other software.

Warning: This must be done only if you downloaded a stable version of the package, and won't do any more changes on the scripts and modules. Indeed, by installing the package, the changes done after won't be taken into account while running the scripts.  
If you plan to work on the scripts, or to download a new version after,
choose the deployment installation (see below).


## pipeline installation

If you want to install the package while still working on modifying the scripts, or being able to download and run latest versions, just do::

    sudo pip3 install -e .

Your changes will then be taken into account. As you installed the package, you will be able to run the pipeline from any directory in your computer.


## Uninstalling pipeline

Whatever the way you installed the pipeline, you can uninstall it by running::

    sudo pip uninstall pipelinepackage


## Running Tests

If you want to work on pipeline scripts, you can use the tests provided with the software, used to check each of its functionalities. To run the tests, run, from the root of the project::

    PYTHONPATH+=. py.test

or, if you installed the package (final or development mode)::

    py.test

**Input:**
- File with list of genomes to format: 1 genome per line, 2 columns, no header. First column with gembase name of the genome, 2nd column with the current name of the genome (with file extension). The gembase name is: GGSS.mmyy.nnnnn with GGSS the 2 first letters of gender and 2 first letters of species, mmyy the month and year, nnnnn the strain number (with trainling 0s if needed).
- all genome sequences in multifasta (with current name) in a same folder `dbpath`.

**Output:**
- A folder `respath` with subfolders `Genes`, `Replicons`, `Proteins`, `LSTINFO` with all genomes in each.
- The list of genomes updated with 2 new columns: number of contigs and total length of sequence
- The log files .err and .out

## Prepare sequences
Prepares the raw genomic sequence files for the pipeline, with `prepare_sequences.py`:
- rename all contigs of a given multifasta file. Contigs are named with the genome name (= fasta filename) + the contig number + the size of the contig (bp)
- get the total number of contigs in each genome
- get the total size (in nuc) o each genome

Input: 
- database folder, containing all sequences in multifasta, with original names
- a file, "lstinfo-file" with 2 columns: gembase name and original name of each genome

Output: 
- for each genome, a new multi-fasta in the database folder (with contig names changed), called "original_name'-gembase.fna
- a file "lstinfo-file"-complete.lst with 4 columns: gembase name, original name, nb contigs, genome size

## Annotate genomes
used :
prokka 1.11 (with hmmer/3.1b1 aragorn/1.2.36 barrnap/0.4.2 minced/0.1.6 blast+/2.2.28 prodigal/2.60 infernal/1.1 ncbi_toolbox/20151127 signalp/4.0)

- Script `prokka_array.sh`
- Output: In the `Database` folder, for each genome, a subfolder `<genome_original_name>-prokka11Res` with prokka output files (see https://github.com/tseemann/prokka#output-files).

Prokka annotates each gene with 'PROKKA_xxxxx' with 'xxxxx' a unique ID for the gene in the whole genome (even if several contigs). This is the case for all CDS, tRNA, tmRNA, rRNA...   
But, for repeat_regions (CRISPR), there is no locus tag. The only information given in `.tbl`file is: rpt_family CRISPR and score (number of repetitions).  
In the `.ffn` file, gene headers are 'PROKKA_xxxxx + annotation' for CDS, tRNA etc., or only the 'contig header' for CRISPR regions (hence, if several CRISPR regions, several genes with the same header).  
If `--locustag` is given to prokka command, then PROKKA is replaced by the given locus tag in the gene IDs. But CRISPR regions are still named with the contig ID.

To check that Prokka ran well on all genomes: `check_prokka_run.sh` (called from `prokka_array.sh`). If problems occur, they are added to the error output file, called `err-<LSTINFO_file>-<current_date>.err`.


## Translate outputs to gembase format

This is done with the `generate_gembase.sh` script (called from `prokka_array.sh`):
- Parse .tbl prokka results and transform to .lst, with one line per gene, with the following fields: start, end, strand (C or D), type (CDS, RNA, CRISPR...), gembase_name, gene_name | product | EC_number | more information.
    + script: `tblToLst.awk` per genome.
    + output: `LSTINFO/<gembase_name>.lst`
- Copy `Database/<gembase_name>-gembase.fna` to `Replicons/<gembase_name.fna`
- Generate .prt (resp. .gen) fles, fromm .faa (resp .ffn) files generated by prokka in `<gembase_name>-prokka11Res`.
    + Script: `create_prt_gen.py` per genome
    + Output: `Genes/<gembase_name>.gen` and `Proteins/<gembase_name>.prt`.

Check that all went well (number of proteins/crispr/total genes corresponding between prokka results and generated files) with `check_gembase.sh`. If problems occur, they are added to the error output file, called `err-<LSTINFO_file>-<current_date>.err`.
