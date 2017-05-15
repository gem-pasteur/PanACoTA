# **-- genomeAPCAT --**

[![build status](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/build.svg)](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/commits/master)
[![coverage report](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/coverage.svg)](http://aperrin.pages.pasteur.fr/pipeline_annotation/htmlcov)

genome APCAT is a software providing tools for large scale comparative genomics:
- annotation of genomes
- pan-genome
- persistent genome
- phylogenetic tree from persistent genome 


# Installation

## Dependencies

genomeAPCAT is written in **python3**. So, you need python3 (and pip3 for installation) to run it.

Its external dependencies are:
- prokka (to annotate the genomes). 

You can either install the external dependencies yourself, with the version you want, or use the installation script `make.py`, which will install the dependencies.

To be able to install the dependencies (by yourself, or with the installation script), make sure you have: `tar`, `git` and `wget`.

Then, for prokka installation, you need to install some system packages, as well as bioperl and java, if not already done (see [Prokka README](https://github.com/tseemann/prokka) for more information).

## Downloading `genomeAPCAT`

You can download `genomeAPCAT` source code by cloning its gitlab repository. Then, go to the directory created to install `genomeAPCAT` as described bellow. Here is an example to download the code and move to the created repository (called 'root directory' in the next sections):

    git clone https://aperrin@gitlab.pasteur.fr/aperrin/pipeline_annotation.git
    cd pipeline_annotation

## Installing `genomeAPCAT` (final mode)

To install `genomeAPCAT`, and all its dependencies, from the root directory, type:

    ./make.py 

or 

    ./make.py install

You will then be able to use the package from any directory in your computer,
just as any other software.

**Warning:** If you plan to work on the scripts, or to download a new version after, choose the development installation (see below).


## Installing `genomeAPCAT` (development mode)

If you want to install `genomeAPCAT` while still working on modifying the scripts, type:

    ./make.py develop

Your changes will then be taken into account. As you installed the package, you will be able to run it from any directory in your computer.

## Uninstalling `genomeAPCAT`

If you don't want `genomeAPCAT` anymore, or if you want to install a newer version, uninstall it by typing:

    ./make.py uninstall

## Update to new version

If you want to install a new version of `genomeAPCAT`:
- uninstall the previous version (`./make.py uninstall`)
- download the new version
- install the new version (`./make.py`)

## Cleaning dependencies

If you installed the dependencies (such as prokka) via our installation script, but now want to install your own version, you can remove all dependencies downloaded and installed by `make.py` by doing:

    ./make.py clean

# Running `genomeAPCAT`

## Quick run

`genomeAPCAT` contains 4 different modules: `annotate`, `pan-genome`, `pers-genome` and `tree`. You can run then by typing:

    genomeAPCAT <module_name> <arguments_for_module>

Each module has its own options and inputs. Type: `genomeAPCAT help <module_name>` or `genomeAPCAT <module_name> -h` to get the list of required arguments and other available options for the module you want to run.

## Examples

Provide an example dataset, with command lines for all modules... ?

# Development

This part is for people who want to work on developing `genomeAPCAT` package.

## Running Tests

If you want to work on the scripts, you can use the tests provided with the software, used to check each of its functionalities. To run the tests, run, from the root of the project:

    PYTHONPATH+=. py.test test/test_unit
    PYTHONPATH+=. py.test test/test_functional

or, if you installed the package (final or development mode)::

    py.test test/test_unit
    py.test test/test_functional
