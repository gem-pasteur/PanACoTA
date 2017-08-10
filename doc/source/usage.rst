=============================
Running genomeAPCAT: tutorial
=============================

``genomeAPCAT`` contains 5 subcommands, for the different steps:
    - ``annotate`` (annotate all genomes of the dataset, after a quality control)
    - ``pangenome`` (generate pan-genome)
    - ``corepers`` (generate core-genome or persistent-genome)
    - ``align`` (align core/persistent families)
    - ``tree`` (infer phylogenetic tree from persistent genome)

You can run them by typing::

    genomeAPCAT <subcommand_name> <arguments_for_subcommand>

Each subcommand has its own options and inputs. To get the list of required arguments and other available options for the subcommand you want to run, type::

    genomeAPCAT <subcommand> -h


.. note:: In the example command lines, we put ``<>`` around the fields that you have to replace by the information corresponding to what you want. For example, if we write ``command -D <seqfile>`` and the sequence file you want to use is in your current directory and is called ``my_sequence.fa``, then you should write ``command -D my_sequence.fa``.

.. note:: In the example command lines, commands between ``[]`` are optional, meaning that you can run the command line without this part. For example, if we write ``command -D <seqfile> [-t <num> -i <percentage>]``, and your sequence file is the same as previously, and the default parameters are 10 for ``-t`` and 0.5 for ``-i``. You can run either ``command -D my_sequence`` (using default parameters for both options: 10 and 0.5), ``command -D my_sequence -t 8`` (specifying ``t=8`` option and default ``i=0.5``), ``command -D my_sequence -i 0.9`` (default ``t=10`` and specified ``i=0.9``) or ``command -D my_sequence -t 8 -i 0.9`` (specifying both options: ``t=8`` and ``i=0.9``) according to your needs (if default values are ok, you do not need to specify the option).



Here are the options shared by all subcommands:

    - ``-h`` or ``--help``: show help on subcommand, as described here above.
    - ``--quiet``: do not write anything on stdout nor stderr. Log files are still created, so you can check what is running.
    - ``-v`` or ``--verbose``: be more verbose:

        + ``-v`` will add warnings in stderr (by default, only errors are displayed in stderr, warning are just in log files),
        + ``-vv`` will do the same as ``-v``, and also add details to stdout (by default, only info is written to stdout)

We will now describe each subcommand, with its options.


``annotate`` subcommand
=======================

You can see all required arguments and available options with::

    genomeAPCAT annotate -h

The input for annotation is a set of genomes, in (multi-)fasta format. All files to annotate must be in a same directory, referred after by ``<db_path>``. However, this directory can also contain other files/sequences, not used in this study. The program will only use the files specified in the ``<list_file>``, which is the main file you have to provide for this step.

Input file formats
------------------

.. _lfile:

'list_file'
^^^^^^^^^^^

The ``list_file`` is a text file with the following format:

    - 1 genome per line. If a genome is contained in several (multi-)fasta files, give all filenames, separated by a space.
    - after the filename(s), you can specify more information on the genome. If you want to do so, add ``::`` to separate the genome filename(s) and the informations. Possible informations are:

        - the species name. Usually, we use the 2 first letters of genus and 2 first letters of species (e.g. ESCO for Escherichia coli). But you can choose any name, as long as it contains 4 alpha-numeric characters (letters or/and numbers). If the species name is not given in the genome line, the program will use the one given by the ``-n <name>`` option when running the command. Specifying the species name at a genome line in the ``list_file`` is useful when you want to annotate several genomes from different species. If all your dataset corresponds to the same species, just provide its name with the ``-n <name>`` option!

        - the date. Separate the species name and the date by a ``.``. If no species name given, just put this dot after the ``::`` separating filenames and information. This date allows you to specify when the genome was sequenced/retrieved, with 4 digits (MMYY). This can be useful if some genomes have not been sequenced at the same time as others, and you want to keep this information for later analyses. If not given, the program will use:

            + the date given with ``--date <date>`` option if given by user
            + today's date if not given

Example:

.. code-block:: text

    genome1.fasta
    genome2-chromo1.fasta genome2-pl.fst
    g3.fa :: ESCO
    gen4-contigs.fst :: ESCO.0217
    genome.fasta genome-plasmids.fasta :: .0217

We have here a dataset with 5 genomes:
    - the 1st genome's sequence is in the file called ``genome1.fasta`` (it can be either a fasta or multi-fasta, according to the assembly status - complete/draft - of the genome). Its species name and date will be the default ones given to the program
    - the 2nd genome's sequence is in 2 files: for example, its chromosome is in ``genome2-chromo1.fasta``, and its plasmid is in ``genome2-pl.fst``. Again, each of those files can contain complete or draft sequences. As the previous genome, its species name and date will be the default ones.
    - the 3rd genome's sequence is in ``g3.fa``. Its species name will be ``ESCO``, while its date will be the default one.
    - the 4th genome's sequence is in ``gen4-contigs.fst``. Its species name will be ``ESCO``, and its date ``0217`` (February 2017).
    - the 5th genome's sequence is in ``genome.fasta`` and ``genome-plasmids.fasta``. Its species name will be the default one, and the date will be ``0217``.


.. _seq:

sequence files
^^^^^^^^^^^^^^

Sequence files must be in fasta or multi-fasta format. A complete genome with only 1 chromosome will hence contain only 1 fasta entry. For example::

    >genome1
    ACCTTAGAGCGCTCTCGCGCATAG

If a genome contains several replicons (either chromosome and plasmids, either a draft genome with several contigs), it contains 1 fasta entry per replicon. For example::

    >genome1-chromo-contig1
    ACCGAAGCGCGCGAGAGTGTGTGGGA...
    >genome1-chromo-contig2
    ACCGAGAGCGCGCGCGGGAGAGAGAGAGC...
    >genome1-chromo-contig3
    ACACGAGCAATATACAGCAGACAGCAGACATATACTCAGACGACAG...
    >genome1-plasmid
    ACAGACGACATAAGAGACGACACAAAAAACACAGAGTTTATGA...

With some softwares, the different contigs of a draft genome are all concatenated in a same fasta entry, and their sequences are separated by stretches of ``N``. For example::

    >genome_seq
    AACACACGATCTCGGCAGCGCANNNNNNNNNNNNNACAGCATNNNNTCGCGCCGACGNNACTATAACAGCAGACNNNNNNNNNNCACACCGGGTATCAGCAGCAGACGACGACGAACGAANNNNNNNNNNACACAGCACTATACGNACAGCA...

This genome is a draft with 4 contigs. By default, ``genomeAPCAT`` will split the sequences each time there is stretch of at least 5 ``N``, in order to have 1 replicon per fasta entry. For example, with the previous file in input, it will create a new multi-fasta file with::

    >genome_seq_cont1
    AACACACGATCTCGGCAGCGCA
    >genome_seq_cont2
    ACAGCATNNNNTCGCGCCGACGNNACTATAACAGCAGAC
    >genome_seq_cont3
    CACACCGGGTATCAGCAGCAGACGACGACGAACGAA
    >genome_seq_cont4
    ACACAGCACTATACGNACAGCA

Stretches of less than 5 ``N`` are kept, while the longer ones are removed, and the 2 parts form 2 different entries.

If you want to deactivate this feature, or choose another minimal number of ``N`` to split, you can specify it with the option ``--cutN <number>`` (0 to deactivate) while running the program (see :ref:`options <option>`).

.. _outform:

Output file formats
-------------------

The annotation step will create 4 result folders. Here is a description of their content.

.. _lstinfof:

'LSTINFO_<list_file>.lst' file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains the list of all genomes annotated, sorted by species, and, in each species, by increasing L90 and number of contigs, with 5 columns:
    - new name of genome (called 'gembase_name'), with format ``<name>.<date>.<strain>`` with:

        - ``name`` given in ``-n <name>`` or line in list_file
        - ``date`` given in ``--date <date>``, line in list_file or current date
        - ``strain`` is a number with 5 digits, identifying the different genomes of a same species.
        - for example: ``ESCO.0217.00002`` for the 2nd strain of Escherichia coli.
    - original name of genome (as given in list_file)
    - genome size (number of bases)
    - number of contigs in genome
    - L90 of genome

Example:

.. code-block:: text

    gembase_name    orig_name                   gsize   nb_conts    L90
    ESCO.0817.00001 genome1.fst                 9808    2           2
    ESCO.1216.00002 genome3-chromo.fst-all.fna  8817    3           3
    GEN2.0817.00001 genome2.fst                 10711   4           4
    GEN4.1111.00001 genome4.fst                 7134    1           1


LSTINFO folder
^^^^^^^^^^^^^^

This folder contains 1 file per genome, called ``<genome_name>.lst``, containing 1 line per sequence annotated (gene, tRNA, rRNA etc.), with the following informations:
    - start position of sequence in the replicon
    - end position of sequence in the replicon
    - strand (D for direct, C for complement)
    - type of sequence (CDS, rRNA, CRISPR, etc.)
    - name of the sequence annotated. The name is ``<genome_name>.<place><contig>_<num>`` where:

        + ``<place>`` is ``i`` when the sequence is inside its replicon, or ``b`` when it is at the border of its replicon (first and last sequence of each replicon)
        + ``<contig>`` is the contig number, with 4 digits
        + ``<num>`` is the unique sequence number.
        + For example: ``ESCO.0217.00002.i0001_00005`` is a gene from the 2nd strain of E. coli, in contig 1 (not the first or last gene of this contig), and is the 5th sequence annotated in this genome.
    - gene name when applicable
    - more information on the sequence annotated (product, similar sequences in PFAM, etc.)

Example of a file which would be called ``ESCO.0417.00002.lst``:

.. code-block:: text

    34685   35866   C       CDS     ESCO.0417.00002.b0001_00001     thlA                | Acetyl-CoA acetyltransferase | 2.3.1.9 | similar to AA sequence:UniProtKB:P45359
    37546   40215   D       tRNA    ESCO.0417.00002.i0001_00002     NA                  | tRNA-Met(cat) | NA | COORDINATES:profile:Aragorn:1.2
    45121   47569   D       CDS     ESCO.0417.00002.i0001_00003     NA                  | Prophage CP4-57 regulatory protein (AlpA) | NA | protein motif:Pfam:PF05930.6
    50124   52465   D       CDS     ESCO.0417.00002.b0001_00004     P22 coat protein 5  | P22 coat protein - gene protein 5 | NA | protein motif:Pfam:PF11651.2
    1       2600    C       tRNA    ESCO.0417.00002.b0004_00005     NA                  | tRNA-Gly(ccc) | NA | COORDINATES:profile:Aragorn:1.2
    3500    5000    D       CDS     ESCO.0417.00002.i0004_00006     NA                  | hypothetical protein | NA | NA
    10000   10215   C       CRISPR  ESCO.0417.00002.b0004_CRISPR1   crispr              | crispr-array | NA | NA
    4568    5896    D       CDS     ESCO.0417.00002.b0006_00007     NA                  | hypothetical protein | NA | NA
    126     456     D       CDS     ESCO.0417.00002.b0007_00008     NA                  | hypothetical protein | NA | NA

Proteins folder
^^^^^^^^^^^^^^^

This folder contains 1 file per genome, called ``<genome_name>.prt``. This file is a multi-fasta file, and contains amino-acid sequences, corresponding to all CDS annotated.

Genes folder
^^^^^^^^^^^^

This folder contains 1 file per genome, called ``<genome_name>.gen``. This file, in multi-fasta format, contains nucleic sequences, corresponding to all sequences annotated (found in corresponding file in LSTINFO folder).

Replicons folder
^^^^^^^^^^^^^^^^

This folder contains 1 file per genome, called ``<genome_name>.fna``. It corresponds to the input file, containing all replicons of the genome, but with contigs renamed.

.. _qco:

Quality Control only
--------------------

Before annotating all genomes, we advise to run once the program with the ``-Q`` option, to do the quality control, but not the annotation. In that case, for each line of the list_file, it will:

    - concatenate sequences in 1 file if several are given
    - split concatenated contigs into different entries (see :ref:`sequences format <seq>`)
    - calculate the genome characteristics:

        + L90: minimum number of contigs needed to cover at least 90% of the sequence
        + number of contigs
        + sequence length

With this information, you will be able to see which genomes should be removed from the study, because of their bad quality. Then, you can annotate only the genomes you keep for the study.

You can run this quality control with (order of arguments does not matter)::

    genomeAPCAT annotate <list_file> -d <dbpath> -r <res_path> -Q

with:

    - ``<list_file>`` your list file as described in :ref:`input formats<lfile>`.
    - ``-d <dbpath>`` the path to the folder containing all your fasta files listed in list_file.
    - ``-r <res_path>`` path to the directory where you want to put the results (no need to create the directory before, the program will do it).
    - ``-Q`` specify that you only want the quality control

This will create a folder ``<res_path>``, with the following files inside:

    - ``QC_L90-<list_file>.png``: histogram of the L90 values of all genomes
    - ``QC_nb-contigs-<list_file>.png``: histogram of number of contigs in all genomes
    - ``discarded-<list_file>.lst``: list of genomes that would be discarded if you keep the default limits (L90 :math:`\leq` 100 and #contigs :math:`\leq` 999).
    - ``genomeAPCAT-annotate_<list_file>.log``: log file. See information on what happened during the run: traceback of stdout.
    - ``genomeAPCAT-annotate_<list_file>.log.err``: log file but only with Warnings and errors. If it is empty, everything went well!
    - ``genomeAPCAT-annotate_<list_file>.log.details``: same as ``.log`` file, but with more detailed information (for example, while running annotation, you can have the time of start/end of annotation of each individual genome). This file can be quite big if you have a lot of genomes.
    - ``info-genomes-<list_file>.lst``: file with information on each genome: size, number of contigs and L90.
    - ``tmp_files`` folder: containing your genomic sequences, split at each stretch of at least 5 ``N``.

.. _annot:

Annotation
----------

When you know the limits you want to use for the L90 and number of contigs, you can run the full annotation step, and not only the quality control. Use::

    genomeAPCAT annotate <list_file> -d <dbpath> -r <res_path> -n <name> [--l90 <num> --nbcont <num>]

with:
    - same arguments as before
    - ``-n <name>`` the default species name to use, for lines of the list_file which do not contain this information. This name must contain 4 alpha-numeric characters.
    - ``--l90 <num>``: *optional*. If the default value (max L90 = 100) does not fit your data, choose your own maximum limit.
    - ``--nbcont <num>``: *optional*. If the default value (max nb_contigs = 999) does not fit your data, choose your own maximum limit.

This command will run the same steps as described in quality control only, with additional steps:

    - annotation of each genome with prokka
    - formatting prokka results to the 4 output folders (see :ref:`output formats <outform>`)

This will also create a folder ``<res_path>``, with the following files inside:

    - same files as quality control only, except ``info-genomes-<list_file>.lst``.
    - ``LSTINFO_<list_file>.lst``: information on annotated genomes, as described :ref:`here<lstinfof>`
    - prokka result folders in your ``tmp_files`` directory
    - The 4 folders ``LSTINFO``, ``Replicons``, ``Genes`` and ``Proteins`` as described in :ref:`output file formats<outform>`.

.. _option:

Options
-------

Here is the list of options available when running ``genomeAPCAT annotate``:

    - ``-n <name>``: required when not running quality control only (see :ref:`annotation<annot>`)
    - ``-Q``: run quality control only (see :ref:`QC only<qco>`)
    - ``--l90 <l90>``: to specify the maximum L90 value accepted to keep a genome. Default is 100
    - ``--nbcont <number>``: to specify the maximum number of contigs allowed to keep a genome. Default is 999
    - ``--cutN <number>``: by default, each sequence is split at each stretch of at least 5 ``N`` (see :ref:`sequence format<seq>`). If you do not want to split sequences, put 0. If you want to change the condition, put the minimum number of ``N`` required to split the sequence.
    - ``--date <date>``: date used to name the genome (in gembase_format, see :ref:`first column of LSTINFO_file<lstinfof>`). If not given, and no information is given on a line in the list_file, the current date will be used.
    - ``--tmp <tmpdir>``: to specify where the temporary files must be saved. By default, they are saved in ``<res_path>/tmp_files``.
    - ``--prok <prok_dir>``: to specify where the prokka output folders must be saved. By default, they are saved in the same directory as ``<tmpdir>``. This can be useful if you want to run this step on a dataset for which some genomes are already annotated. For those genomes, it will use the already annotated results found in ``<prok_dir>`` to run the formatting steps, and it will only annotate the genomes not found.
    - ``-F`` or ``--force``: Force run: Add this option if you want to run prokka and formatting steps for all genomes even if their result folder (for prokka step) or files (for format step) already exist: override existing results. Without this option, if there already are results in the given result folder, the program stops. If there are no results, but prokka folder already exists, prokka won't run again, and the formating step will use the already existing folder if correct, or skip the genome if there are problems in prokka folder.
    - ``--threads <number>``: if you have several cores available, you can use them to run this step faster, by handling several genomes at the same time, in parallel. By default, only 1 core is used. You can specify how many cores you want to use, or put 0 to use all cores of your computer.


``pangenome`` subcommand
========================


dsfgrgr


``corepers`` subcommand
=======================

dgfdgdf

``align`` subcommand
====================


sdgdfgd

``tree`` subcommand
===================


fdgfdgdf

