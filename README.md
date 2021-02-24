# **PanACoTA**
[![GitHub release](https://img.shields.io/github/release/gem-pasteur/PanACoTA.svg)](https://github.com/gem-pasteur/PanACoTA/releases)
[![PyPI version](https://badge.fury.io/py/PanACoTA.svg)](https://badge.fury.io/py/PanACoTA)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/panacota/badges/version.svg)](https://anaconda.org/bioconda/panacota)

<!--
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4724)
-->

[![DOI:10.1093/nargab/lqaa106](https://zenodo.org/badge/DOI/10.1093/nargab/lqaa106.svg)](https://academic.oup.com/nargab/article/3/1/lqaa106/6090162)

 [![License (AGPL version 3)](https://img.shields.io/badge/license-GNU%20AGPL%20version%203-green.svg)](COPYING)
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

Amandine Perrin, Eduardo P. C. Rocha, PanACoTA: a modular tool for massive microbial comparative genomics, *NAR Genomics and Bioinformatics*, Volume 3, Issue 1, March 2021 
[![DOI:10.1093/nargab/lqaa106](https://zenodo.org/badge/DOI/10.1093/nargab/lqaa106.svg)](https://academic.oup.com/nargab/article/3/1/lqaa106/6090162)


**Content of this README:**

Installation
- [Dependences](#dep)
- [pip](#pypi)
- [cloning github repository](#clone)
- [singularity](#singularity)
- [conda](#conda)

Running
- [Quick run](#run)
- [Examples](#example)
- [Documentation](#doc)

[Development](#develop)

# Installation

## <a name="dep"></a> Dependencies

PanACoTA is written in **python3**. So, you need python3 (and pip3 for installation) to run it.

Then, PanACoTA has several external dependencies. If you use [`singularity`](#singularity) installation (for ex. to run on a cluster), you do not need to install any dependency. Otherwise, install only the one(s) you need, according to the module(s) you want to use: 
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

## <a name="install"></a> Installing `PanACoTA` and update

You have different possibilities to install `PanACoTa`.

**Warning:** If you plan to work on the scripts, choose the development installation (see [Development section](#develop)).

### <a name="pypi"></a> From pip

[![PyPI version](https://badge.fury.io/py/PanACoTA.svg)](https://badge.fury.io/py/PanACoTA)

A very simple way to install the last stable version. This will install files in your python site-packages folder.

    pip install panacota

And to get new version:

    pip install panacota --upgrade

If you have permission issues, you can either use 'sudo' before the previous command lines to install it as root, or add the `--user` option to install it locally.

### <a name="clone"></a> From github repository

This allows you to get the very last version, and be able to test the last enhancements before they are uploaded to the other platforms. For that, go to where you want to install it `(<your_dir>)`, and type: 

    git clone https://github.com/gem-pasteur/PanACoTA.git

This will create a repository called `PanACoTA`, containing the content of this github repository. To install PanACoTA:

    cd PanACoTA 
    ./make

If you have permission issues, you can either use 'sudo' before the previous command lines to install it as root, or add the `--user` option to install it locally.

To update to new version, go back to your repository:

    cd <your_dir>/PanACoTA
    git pull
    ./make upgrade


### <a name="singularity"></a> From singularity image

Very useful if you do not have permission rights on the computer, such as, for example, on a cluster. The other advantage is that you do not need to install any dependence (except singularity itself of course). Singularity image includes all of them. You just have to download 1 file, and nothing will be installed anywhere on your computer.

First, download the singularity image:

    singularity pull --name panacota.img docker://gempasteur/panacota[:<version>] 

If you want a specific version, like version 1.0, specify `docker://gempasteur/panacota:1.0`. 

To get latest version:

    singularity pull --name panacota.img docker://gempasteur/panacota

(This is the same as `singularity pull --name panacota.img docker://gempasteur/panacota:latest`)

It will replace your file panacota.img by a new one corresponding to the latest version.

### <a name="conda"></a> From conda

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/panacota/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/panacota/badges/version.svg)](https://anaconda.org/bioconda/panacota)

Be careful while using conda, especially if you are not familiar with it. We advise to install PanACoTA in a dedicated conda environment, in order to avoid unwanted interactions with other softwares (like needed versions of dependencies automatically installed by conda). To install the package, use `conda install -c bioconda panacota`. But, as described in [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#id6), we advise to install it with:

    # Create an environment: This creates the 'myenv' environment in '/envs/'. No packages will be installed in this environment.
    conda create --name myenv
    # Activate the environment
    conda activate myenv
    # Install PanACoTA
    conda install -c bioconda panacota
    # When you have finished using PanACoTA, deactivate environment
    conda deactivate

To update to new version:

    conda update panacota


### From zip version

For people wanting to download source code of a specific version, we provide releases. You can download last one here: [![GitHub release](https://img.shields.io/github/release/gem-pasteur/PanACoTA.svg)](https://github.com/gem-pasteur/PanACoTA/releases)

## <a name="uninstall"></a> Uninstalling `PanACoTA`

If you don't want `PanACoTA` anymore uninstall it by typing:

    pip unintall panacota  # If you installed from pip
    ./make uninstall  # If you installed from github repository

Or, if you used singularity, just remove the downloaded image: `rm -r panacota.img`

# <a name="run"></a> Running `PanACoTA`

## Quick run

`PanACoTA` contains 6 different subcommands:

- `prepare` (download genomes from refseq if you want to, or give your input database, to run a filtering quality control). To help you find NCBI species taxid you need, you can use their [taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi).
- `annotate` (annotate all genomes of the dataset, after a quality control)
- `pangenome` (generate pan-genome)
- `corepers` (generate core-genome or persistent-genome)
- `align` (align core/persistent families)
- `tree` (infer phylogenetic tree from persistent genome)

You can run them by typing:

    PanACoTA <subcommand_name> <arguments_for_subcommand>

Each subcommand has its own options and inputs. To get the list of required arguments and other available options for the subcommand you want to run, type:

    PanACoTA <subcommand> -h


When using singularity, just replace `PanACoTA` by `./panacota.img`:

    ./panacota.img <subcommand_name> <arguments_for_subcommand>  
    ./panacota.img -h 

It also provides a subcommand `PanACoTA all` to run all modules in a row.

## <a name="example"></a> Examples

We provide a folder, `Examples`, containing genomic sequences (in `Examples/genomes`) and examples of input files (in `Examples/input_files`) for the software.
In the [example part of documentation](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc/examples.html), you will find information explaining you how to run the different modules of `PanACoTA` with this dataset, so that you can try the software. We also describe the results that should be created by each command line.

**Note:** the provided genomic sequences are taken from real genomes, but then modified and drastically shortened in order to have an example showing different situations, but running very fast. Hence, the examples results should not be interpreted biologically!

## <a name="doc"></a> Documentation

You can find more information in [PanACoTA documentation](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc)!

# <a name="develop"></a>  Development

This part is for people who want to work on developing `PanACoTA` package. In the documentation, there is a part dedicated to [developers](http://aperrin.pages.pasteur.fr/pipeline_annotation/html-doc/develop.html).

PanACoTA is also hosted in gitlab, where all CI is done. Here is the link: https://gitlab.pasteur.fr/aperrin/pipeline_annotation

## Installing `PanACoTA` (development mode)

If you want to install `PanACoTA` while still working on modifying the scripts, use `./make develop` instead of `./make install` once you have cloned the repository.

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
