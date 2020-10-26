======================
Starting with PanACoTA
======================


``PanACoTA`` is a Python package, developed in Python 3.6.


Dependencies
============

``PanACoTA`` is written in **python3**. So, you need python3 (and pip3 for installation) to run it.

There are several external dependencies. Install only the one(s) you need, according to the module(s) you want to use:

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

Installation:
=============

Downloading and updating
------------------------


You can get ``PanACoTA`` source code by downloading an archive of a given release (zip, tar.gz), or by cloning its github repository. By cloning the github repository, you will then be able to update the code to new versions very easily and quickly. Here is how to clone the repository:

.. code-block:: bash

    git clone https://github.com/gem-pasteur/PanACoTA.git
ve your github login, and password.

This will create a repository called ``PanACoTA``. Go inside this repository (``cd PanACoTA``) to install the software, as described hereafter.

If a new version of ``PanACoTA`` is released, and you want to use it, type the following command to update the source code:

.. code-block:: bash

    git pull

Then, you will be able to upgrade to the new version (see :ref:`Upgrade section <upgrade>`).


.. _installing:

Installing ``PanACoTA``
--------------------------

To install ``PanACoTA``, from its directory, type:

.. code-block:: bash

    ./make

or

.. code-block:: bash

    ./make install

You will then be able to use the package from any directory in your computer,
just as any other software.

.. note:: If you have permission issues, you can either use ``sudo`` before the previous command lines to install it as root, or, if you do not have root access (or prefer a local installation), use ``./make --user`` to install it locally.

.. warning:: If you plan to work on the scripts, choose the development installation (see :doc:`Developer documentation <develop>`).

.. _uninstall:

Uninstalling ``PanACoTA``
----------------------------

If you don't want ``PanACoTA`` anymore, uninstall it by typing:

.. code-block:: bash

    ./make uninstall

.. note:: If you have permission issues, and installed the package as root, use ``sudo`` before the previous command line to uninstall it.

.. _upgrade:

Upgrade to new version
----------------------

If you want to install a new version of ``PanACoTA``:

.. code-block:: bash

    git pull         # update source code to the new version
    ./make upgrade   # upgrade to the new version

.. note:: If you have permission issues, and installed the package as root, use ``sudo`` before the second command line (``sudo ./make upgrade``) to upgrade. Or, if you installed the package locally, use ``./make upgrade --user`` to upgrade this local version.

If you installed it by downloading a zip file, :ref:`Uninstall it <uninstall>`, and install the new version (by cloning gitlab repository, or downloading the new zip file).


Quick run
=========

``PanACoTA`` contains 6 different subcommands:

- ``prepare`` (download genomes from refseq if you want to, or give your input database, to run a filtering quality control)
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

Running with singularity image
==============================

We provide a singularity image, to help running PanACoTA on a cluster.

First, download the singularity image::

    singularity pull --name panacota.img shub://gem-pasteur/PanACoTA[:version]

If you want a specific version, like version 1.0, specify ``shub://gem-pasteur/PanACoTA:1.0``. If you want the latest version, use ``shub://gem-pasteur/PanACoTA`` or ``shub://gem-pasteur/PanACoTA:latest``.

Then, you can run PanACoTA in the same way as previously, using:

.. code-block:: bash

    ./panacota.img -h  # to get help on the whole PanACoTA program
    ./panacota.img <subcommand_name> <arguments_for_subcommand>  # to run a module of PanACoTA on your data.
