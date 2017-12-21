========================
Work on genomeAPCAT code
========================

This part is for people who want to work on developing `genomeAPCAT` package: adding new features, correcting bugs etc.


Installing ``genomeAPCAT`` in development mode
==============================================

If you want to install ``genomeAPCAT`` while still working on modifying the scripts, type:

.. code-block:: bash

    ./make develop

Your changes will then be taken into account. As you installed the package, you will be able to run it from any directory in your computer.

If you don't want to install the software, you can still work on it, test it, and contribute to the tests and documentation by only installing the libraries needed for the software, and those
needed for development by running:

.. code-block:: bash

    pip3 install -r requirements.txt  # dependencies used by genomeAPCAT
    pip3 install -r requirements-dev.txt  # libraries used to run tests, generate documentation etc.

.. note:: biopython is only used for 'tree' subcommand, with option ``--soft fastme`` or ``--soft quicktree``. If you do not plan to use this, you do not need to install biopython. You can comment  the ``biopython>=1.60`` line in ``requirements.txt`` (add a ``#`` at the beginning of the line).


Running Tests
=============

If you want to work on the scripts, you can use the tests provided with the software, used to check each of its functionalities. Tests are done with ``pytest`` framework.

To run the tests, run, from the root of the project::

    py.test test/test_unit
    py.test test/test_functional

or, if you did not install the package::

    PYTHONPATH+=. py.test test/test_unit
    PYTHONPATH+=. py.test test/test_functional

Add ``-v`` to get more detailed information on each test run.

This will also generate the coverage report. Open ``htmlcov/index.html`` on your browser if you want to check code coverage of your new function/module. The on-line version can be found `here <http://aperrin.pages.pasteur.fr/pipeline_annotation/htmlcov/>`_, and is automatically updated at each push on master branch.

.. warning:: If you add new features, or modify existing scripts please complete/update the tests!




Contributing to documentation
=============================

This documentation is generated with `sphinx <http://www.sphinx-doc.org/en/stable/>`_. You can add your contribution to it. To generate the html documentation locally, go to ``doc/sources`` directory, and run:

.. code-block:: bash

    make html

Then, open ``doc/build/html/index.html`` on your browser.

The on-line version will be automatically updated when modifications on rst files are pushed on master branch.
