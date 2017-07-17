Using genomeAPCAT
*****************


``genomeAPCAT`` is a Python package, developed in Python 3.

Installation: '``./make``' and its options
========================================================

Dependencies
------------

``genomeAPCAT`` is written in **python3**. So, you need python3 (and pip3 for installation) to run it.

Its external dependencies are:

- `prokka <https://github.com/tseemann/prokka>`_  (to annotate the genomes). If you do not already have it, this software will be automatically installed by `make` script, along with `genomeAPCAT` (see :ref:`Install section <installing>`)
- `mmseqs <https://github.com/soedinglab/MMseqs2>`_  (to generate pangenomes)
- `fftns <http://mafft.cbrc.jp/alignment/software/>`_ (from mafft, to align persistent genome)
- `FastTreeMP <http://www.microbesonline.org/fasttree/#Install>`_ (to infer a phylogenetic tree). We advise to download C code, and compile as described here above.

You can either install the external dependencies yourself, with the version you want, or use the installation script ``make``, for the dependencies that it installs.

To be able to install the dependencies (by yourself, or with the installation script), make sure you already have:

- ``tar``
- ``git``
- ``wget``
- bioperl, java and some other base packages required for prokka: see `Prokka README <https://github.com/tseemann/prokka>`_ for more information.

For FastTree, we advise to download C code from `here <http://www.microbesonline.org/fasttree/#Install>`_, and compile it using:

.. code-block:: bash

    gcc -DOPENMP -fopenmp -DUSE_DOUBLE -Wall -O3 -finline-functions -funroll-loops -o FastTreeMP -lm FastTree-2.1.9.c

You can then add the output ``FastTreeMP`` to your ``$PATH`` to be able to run it from everywhere.


Downloading and updating ``genomeAPCAT``
----------------------------------------

You can download ``genomeAPCAT`` source code by downloading an archive (zip, tar.gz), or by cloning its gitlab repository. By cloning the gitlab repository, you will then be able to update the code to new versions very easily and quickly. Here is how to clone the repository:

.. code-block:: bash

    git clone https://gitlab.pasteur.fr/aperrin/pipeline_annotation

Give your gitlab login, and password.

This will create a repository called ``pipeline_annotation``. Go inside this repository to install ``genomeAPCAT``, as described hereafter.

If a new version of ``genomeAPCAT`` is released, and you want to use it, type the following command to update the source code:

.. code-block:: bash

    git pull

Then, you will be able to upgrade to the new version (see bellow).

.. _installing:

Installing ``genomeAPCAT`` (final mode)
---------------------------------------


To install ``genomeAPCAT``, and all its dependencies, from the root directory, type:

.. code-block:: bash

    ./make

or

.. code-block:: bash

    ./make install

You will then be able to use the package from any directory in your computer,
just as any other software.

If you have permission issues, you can either use ``sudo`` before the previous command lines to install it as root, or, if you do not have root access, use ``./make --user`` to install it locally.

.. warning:: If you plan to work on the scripts, choose the development installation (see below).

.. note:: Dependencies installed by ``make``: prokka


Installing ``genomeAPCAT`` (development mode)
---------------------------------------------

If you want to install ``genomeAPCAT`` while still working on modifying the scripts, type:

.. code-block:: bash

    ./make develop

Your changes will then be taken into account. As you installed the package, you will be able to run it from any directory in your computer.

Uninstalling ``genomeAPCAT``
----------------------------

If you don't want ``genomeAPCAT`` anymore, uninstall it by typing:

.. code-block:: bash

    ./make uninstall

Upgrade to new version
----------------------

If you want to install a new version of ``genomeAPCAT``:

.. code-block:: bash

    git pull         # update source code to the new version
    ./make upgrade   # upgrade to the new version


Cleaning dependencies
---------------------

If you installed the dependencies (such as prokka) via our installation script, but now want to install your own version, you can remove all dependencies downloaded and installed by ``make`` by doing:

.. code-block:: bash

    ./make clean

Running ``genomeAPCAT``
=======================

## Quick run

`genomeAPCAT` contains 5 different subcommands:
- `annotate` (annotate all genomes of the dataset, after a quality control)
- `pangenome` (generate pan-genome)
- `corepers` (generate core-genome or persistent-genome)
- `align` (align core/persistent families)
- `tree` (infer phylogenetic tree from persistent genome)

You can run them by typing:

    genomeAPCAT <subcommand_name> <arguments_for_subcommand>

Each subcommand has its own options and inputs. To get the list of required arguments and other available options for the subcommand you want to run, type:

    genomeAPCAT <subcommand> -h

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
