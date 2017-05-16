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

## Downloading and updating `genomeAPCAT`

You can download `genomeAPCAT` source code by downloading an archive (zip, tar.gz), or by cloning its gitlab repository. By cloning the gitlab repository, you will then be able to update the code to new versions very easily and quickly. Here is how to clone the repository:

    git clone https://gitlab.pasteur.fr/aperrin/pipeline_annotation

Give your gitlab login, and password.

This will create a repository called `pipeline_annotation`. Go inside this repository to install `genomeAPCAT`, as described hereafter.

If a new version of `genomeAPCAT` is released, and you want to use it, type the following command to update the source code:

    git pull

Then, you will be able to install the new version (see bellow).

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
- update to the new version (`git pull`)
- install the new version (`./make.py`)

## Cleaning dependencies

If you installed the dependencies (such as prokka) via our installation script, but now want to install your own version, you can remove all dependencies downloaded and installed by `make.py` by doing:

    ./make.py clean

# Running `genomeAPCAT`

## Quick run

`genomeAPCAT` contains 4 different modules: 
- `annotate` (annotate all genomes of the dataset, after a quality control)
- `pan-genome` (generate pan-genome)
- `pers-genome` (generate persistent-genome)
- `tree` (infer phylogenetic tree from persistent genome)

You can run then by typing:

    genomeAPCAT <module_name> <arguments_for_module>

Each module has its own options and inputs. To get the list of required arguments and other available options for the module you want to run, type: 

    genomeAPCAT help <module_name>
    # or 
    genomeAPCAT <module_name> -h

## Examples

We provide a folder, `Examples`, containing genomic sequences (in `Examples/genomes`) and examples of input files (in `Examples/input_files`) for the software.  
In this folder, you will also find a README file, explaining you how to run the different modules of `genomeAPCAT` with this dataset, so that you can try the software. We also describe the results that should be created by each command line.

**Note:** the provided genomic sequences are taken from real genomes, but then modified and shortened in order to have an example showing different situations. Hence, the examples results should not be interpreted biologically!

## Documentation

Not done yet, there will be a documentation describing each module and its arguments/options...

# Development

This part is for people who want to work on developing `genomeAPCAT` package.

## Running Tests

If you want to work on the scripts, you can use the tests provided with the software, used to check each of its functionalities. To run the tests, run, from the root of the project:

    PYTHONPATH+=. py.test test/test_unit
    PYTHONPATH+=. py.test test/test_functional

or, if you installed the package (final or development mode)::

    py.test test/test_unit
    py.test test/test_functional
