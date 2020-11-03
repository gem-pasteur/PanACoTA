======================
Starting with PanACoTA
======================


``PanACoTA`` is a Python package, developed in Python 3.6.


Installation
============

Dependencies
------------

``PanACoTA`` is written in **python3**. So, you need python3 (and pip3 for installation) to run it.

Then, ``PanACoTA`` has several external dependencies. If you use :ref:`Singularity <singularity>` installation (for ex. to run on a cluster), you do not need to install any dependency. Otherwise, install only the one(s) you need, according to the module(s) you want to use:

- For prepare module: `mash <https://mash.readthedocs.io/en/latest/>`_ (to filter genomes)
- For annotate module: `prokka <https://github.com/tseemann/prokka>`_  and/or `prodigal <https://github.com/hyattpd/Prodigal>`_  (to uniformly annotate your genomes)
- For pangenome module: `mmseqs <https://github.com/soedinglab/MMseqs2>`_  (to generate pangenomes)
- For align module: `mafft <http://mafft.cbrc.jp/alignment/software/>`_ (to align persistent genome)
- For tree module: At least one of those softwares:

    - `IQtree <http://www.iqtree.org/>`_
    - `FastTreeMP <http://www.microbesonline.org/fasttree/#Install>`_: We advise to download C code, and compile as described :ref:`here above <fasttree>`.
    - `FastME <http://www.atgc-montpellier.fr/fastme/binaries.php>`_
    - `Quicktree <https://github.com/tseemann/quicktree/releases>`_

To be able to install the dependencies, make sure you already have:

- ``tar``
- ``git``
- ``wget``
- bioperl, java and some other base packages required for prokka: see `Prokka README <https://github.com/tseemann/prokka>`_ for more information.

.. _fasttree:

For FastTree, we advise to download C code from `here <http://www.microbesonline.org/fasttree/#Install>`_, and compile it using:

.. code-block:: bash

    gcc -DOPENMP -fopenmp -DUSE_DOUBLE -Wall -O3 -finline-functions -funroll-loops -o FastTreeMP -lm FastTree.c

You can then add the output ``FastTreeMP`` to your ``$PATH`` to be able to run it from everywhere.

Installation and update:
------------------------

You have different possibilities to install ``PanACoTa``.

.. warning:: If you plan to work on the scripts, choose the development installation (see :doc:`Developer documentation <develop>`).


From pip
********

|pip|

.. |pip| image:: https://badge.fury.io/py/PanACoTA.svg
    :target: https://badge.fury.io/py/PanACoTA

A very simple way to install the last stable version. This will install files in your python site-packages folder.

.. code-block:: bash

    pip install panacota

And to get new version

.. code-block:: bash

    pip install panacota --upgrade

If you have permission issues, you can either use 'sudo' before the previous command lines to install it as root, or add the ``--user`` option to install it locally.


From Github repository
**********************

.. _clone:

This allows you to get the very last version, and be able to test the last enhancements before they are uploaded to the other platforms (pip, conda, singularity...). For that, go to where you want to install it ``(<your_dir>)``, and type:

.. code-block:: bash

    git clone https://github.com/gem-pasteur/PanACoTA.git

This will create a repository called ``PanACoTA``, containing the content of this Github repository. To install PanACoTA, and be able to launch it from anywhere:

.. code-block:: bash

    cd PanACoTA
    ./make


If you have permission issues, you can either use 'sudo' before the previous command lines to install it as root, or add the ``--user`` option to install it locally.

To upload to new version, go back to your repository:

.. code-block:: bash

    cd <your_dir>/PanACoTA
    git pull
    ./make upgrade

From singularity image
**********************

.. _singularity:

|singularity|

.. |singularity| image:: https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg
   :target: https://singularity-hub.org/collections/4724

Very useful if you do not have permission rights on the computer, such as, for example, on a cluster. The other advantage is that you do not need to install any dependence (except singularity itself of course). Singularity image includes all of them. You just have to download 1 file, and nothing will be installed anywhere on your computer.

First, download the singularity image:

.. code-block:: bash

    singularity pull --name panacota.img shub://gem-pasteur/PanACoTA[:<version>]

If you want a specific version, like version 1.0, specify ``shub://gem-pasteur/PanACoTA:1.0``.

To get latest version:

.. code-block:: bash

    singularity pull --name panacota.img shub://gem-pasteur/PanACoTA

(This is the same as ``singularity pull --name panacota.img shub://gem-pasteur/PanACoTA:latest``)

It will replace your file ``panacota.img`` by a new one corresponding to the latest version.


From zip version
****************

For people wanting to download source code of a specific version, we provide releases. You can download last one here:

|zip|

.. |zip| image:: https://img.shields.io/github/release/gem-pasteur/PanACoTA.svg
    :target: https://github.com/gem-pasteur/PanACoTA/releases


.. _installing:

Uninstalling ``PanACoTA``
-------------------------

.. _uninstall:

If you don't want ``PanACoTA`` anymore uninstall it by typing:

.. code-block:: bash

    pip uninstall panacota   # If you installed from pip
    ./make uninstall         # If you installed from github repository

Or, if you used singularity, just remove the downloaded image: ``rm -r panacota.img``.


Quick run
=========

``PanACoTA`` contains 6 different subcommands:

- ``prepare`` (download assemblies from refseq if you want to, or give your input database, to run a filtering quality control). To help you find NCBI species taxid you need, you can use their `taxonomy browser <https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi>`_
- ``annotate`` (annotate all genomes of the dataset, after a quality control)
- ``pangenome`` (generate pan-genome)
- ``corepers`` (generate core-genome or persistent-genome)
- ``align`` (align core/persistent families)
- ``tree`` (infer phylogenetic tree from persistent genome)

You can run them by typing:

.. code-block:: bash

    PanACoTA <subcommand_name> <arguments_for_subcommand>

Each subcommand has its own options and inputs. To get the list of required arguments and other available options for the subcommand you want to run, type:

.. code-block:: bash

    PanACoTA <subcommand> -h

When using singularity, just replace ``PanACoTA`` by ``./panacota.img``:

.. code-block:: bash

    ./panacota.img <subcommand_name> <arguments_for_subcommand>
    ./panacota.img -h

It also provides a subcommand ``PanACoTA all`` to run all modules in a row.
