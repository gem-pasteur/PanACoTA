# **PanACoTA**
 [![License (AGPL version 3)](https://img.shields.io/badge/license-GNU%20AGPL%20version%203-green.svg)](COPYING)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4724)
[![DOI:10.1101/2020.09.11.293472](https://zenodo.org/badge/DOI/10.1101/2020.09.11.293472.svg)](https://doi.org/10.1101/2020.09.11.293472)

[![GitHub release](https://img.shields.io/github/release/gem-pasteur/PanACoTA.svg)](https://github.com/gem-pasteur/PanACoTA/releases)
[![pipeline status](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/pipeline.svg)](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/-/commits)
[![coverage report](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/coverage.svg)](http://aperrin.pages.pasteur.fr/pipeline_annotation/htmlcov)

<!--
 To get right angles instead of rounded ones, 
add '?style=flat-square' after .svg -->


This README file provides some essential information to install/use PanACoTA. But it is better to read the [**full documentation**](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc), providing more details : [![](doc/source/images/manual.jpg) ](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc)

&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;

``PanACoTA``  (PANgenome with Annotations, COre identification, Tree and corresponding Alignments) is a software providing tools for large scale bacterial comparative genomics. You can download all refseq genomes for a given species, or use your set of complete and/or draft genomes, to:

- Do a quality control of your strains, to eliminate poor quality genomes, which would not give any information for the comparative study
- Uniformly annotate all genomes (with functional annotation, or only syntactic annotation, according to your needs)
- Do a Pan-genome
- Do a Core or Persistent genome
- Align all Core/Persistent families
- Infer a phylogenetic tree from the Core/Persistent families

&nbsp;
&nbsp;
&nbsp;

If you use PanACoTA, please cite:

Amandine PERRIN, Eduardo P.C. ROCHA (2020).  PanACoTA: A modular tool for massive microbial comparative genomics. BioRxiv. 
[![DOI:10.1101/2020.09.11.293472](https://zenodo.org/badge/DOI/10.1101/2020.09.11.293472.svg)](https://doi.org/10.1101/2020.09.11.293472)


# Installation

## Dependencies

PanACoTA is written in **python3**. So, you need python3 (and pip3 for installation) to run it.

PanACoTA has several dependencies. Install only the one(s) you need, according to the module you want to use: 
- For prepare module: [**mash**](https://mash.readthedocs.io/en/latest/) (to filter genomes)
- For annotate module: [**prokka**](https://github.com/tseemann/prokka) and/or [**prodigal**](https://github.com/hyattpd/Prodigal) (to uniformly annotate your genomes) 
- For pangenome module: [**mmseqs**](https://github.com/soedinglab/MMseqs2) (to generate pangenomes)
- For align module: [**mafft**](http://mafft.cbrc.jp/alignment/software/) (to align persistent genome)
- For tree module: At least one of those softwares, to infer a phylogenetic tree:
    - [**IQ Tree**](http://www.iqtree.org/)
    - [**FastTreeMP**](http://www.microbesonline.org/fasttree/#Install): We advise to download C code, and compile as described here above.
    - [**FastME**](http://www.atgc-montpellier.fr/fastme/binaries.php)
    - [**Quicktree**](https://github.com/tseemann/quicktree/releases)

To be able to install the dependencies, make sure you already have:

- `tar`
- `git`
- `wget`
- bioperl, java and some other base packages required for prokka: see [Prokka README](https://github.com/tseemann/prokka) for more information.

For FastTree, we advise to download C code from [here](http://www.microbesonline.org/fasttree/#Install), and compile it using:

    gcc -DOPENMP -fopenmp -DUSE_DOUBLE -Wall -O3 -finline-functions -funroll-loops -o FastTreeMP FastTree-2.1.9.c -lm

You can then add the output `FastTreeMP` to your `$PATH` to be able to run it from everywhere.


## Downloading and updating `PanACoTA`

You can download `PanACoTA` source code by downloading a [zip file](https://github.com/gem-pasteur/PanACoTA/archive/master.zip), or by cloning its github repository (which is more recommended). By cloning the github repository, you will then be able to update the code to new versions very easily and quickly. Here is how to clone the repository:

    git clone https://github.com/gem-pasteur/PanACoTA.git

This will create a repository called `PanACoTA`. Go inside this repository to install `PanACoTA`, as described hereafter.

If a new version of `PanACoTA` is released, and you want to use it, type the following command to update the source code:

    git pull

Then, you will be able to upgrade to the new version (see bellow).

## <a name="install"></a> Installing `PanACoTA` (final mode)


To install `PanACoTA`, and all its dependencies, from the root directory, type:

    ./make

or

    ./make install

You will then be able to use the package from any directory in your computer,
just as any other software.

If you have permission issues, you can either use 'sudo' before the previous command lines to install it as root, or add the `--user` option to install it locally.

**Warning:** If you plan to work on the scripts, choose the development installation (see [Development section](#develop)).

## <a name="uninstall"></a> Uninstalling `PanACoTA`

If you don't want `PanACoTA` anymore, or if you want to install a newer version, uninstall it by typing:

    ./make uninstall

## Upgrade to new version

If you want to install a new version of `PanACoTA` (and you got it by cloning the gitlab repository):
- update source code to the new version (`git pull`)
- upgrade installation to the new version (`./make upgrade`)

If you installed it by downloading a zip file, [Uninstall it](#uninstall), and [install](#install) the new version (by cloning gitlab repository, or downloading the new zip file).


# Running `PanACoTA`

## Quick run

`PanACoTA` contains 6 different subcommands:

- `prepare` (download genomes from refseq if you want to, or give your input database, to run a filtering quality control)
- `annotate` (annotate all genomes of the dataset, after a quality control)
- `pangenome` (generate pan-genome)
- `corepers` (generate core-genome or persistent-genome)
- `align` (align core/persistent families)
- `tree` (infer phylogenetic tree from persistent genome)

You can run them by typing:

    PanACoTA <subcommand_name> <arguments_for_subcommand>

Each subcommand has its own options and inputs. To get the list of required arguments and other available options for the subcommand you want to run, type:

    PanACoTA <subcommand> -h

## Examples

We provide a folder, `Examples`, containing genomic sequences (in `Examples/genomes`) and examples of input files (in `Examples/input_files`) for the software.
In the [example part of documentation](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc/examples.html), you will find information explaining you how to run the different modules of `PanACoTA` with this dataset, so that you can try the software. We also describe the results that should be created by each command line.

**Note:** the provided genomic sequences are taken from real genomes, but then modified and drastically shortened in order to have an example showing different situations, but running very fast. Hence, the examples results should not be interpreted biologically!

## Documentation

You can find more information in [PanACoTA documentation](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc)!


## On cluster

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4724)

We provide a singularity image, to help running PanACoTA on a cluster.

First, download the singularity image:

    singularity pull --name panacota.img shub://gem-pasteur/PanACoTA[:<version>]

If you want a specific version, like version 1.0, specify `shub://gem-pasteur/PanACoTA:1.0`. If you want the latest version, use `shub://gem-pasteur/PanACoTA`or `shub://gem-pasteur/PanACoTA:latest`

Then, you can run PanACoTA as previously:

    ./panacota.img -h  # to get help on the whole PanACoTA program
    ./panacota.img <subcommand_name> <arguments_for_subcommand>  # to run a module of PanACoTA on your data.


# <a name="develop"></a>  Development

This part is for people who want to work on developing `PanACoTA` package. In the documentation, there is a part dedicated to [developers](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc/develop.html).

PanACoTA is also hosted in gitlab, where all CI is done. Here is the link: https://gitlab.pasteur.fr/aperrin/pipeline_annotation

## Installing `PanACoTA` (development mode)

If you want to install `PanACoTA` while still working on modifying the scripts, type:

    ./make develop

Your changes will then be taken into account. As you installed the package, you will be able to run it from any directory in your computer.

If you don't want to install the software, you can still test it, and contribute to the tests and documentation by installing the libraries needed for the software, and those 
needed for development by running:

    pip3 install -r requirements.txt  # dependencies used by PanACoTA
    pip3 install -r requirements-dev.txt  # libraries used to run tests, generate documentation etc.

**Note:** biopython is only used for 'tree' subcommand, with option ``--soft fastme`` or ``--soft quicktree``. If you do not plan to use this, you do not need to install biopython. You can comment the ``biopython>=1.60`` line in ``requirements.txt`` (add a ``#`` at the beginning of the line). 

## Running Tests

If you want to work on the scripts, you can use the tests provided with the software, used to check each of its functionalities. To run the tests, run, from the root of the project:

    PYTHONPATH+=. py.test test/test_unit
    PYTHONPATH+=. py.test test/test_functional

or, if you installed the package (final or development mode):

    py.test test/test_unit
    py.test test/test_functional

Add ``-sv`` for more details on each test.

## Contributing to documentation

The full documentation, found [here](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc) is generated with
 [sphinx](http://www.sphinx-doc.org/en/stable/).
You can add your contribution to it. To generate the html documentation locally, go to ``doc/sources`` directory, and run:

    make html
    
Then, open ``doc/build/html/index.html`` on your browser.

The online version will be automatically updated when pushed on master branch.
