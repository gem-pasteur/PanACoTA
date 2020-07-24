=====================
Work on PanACoTA code
=====================

This part is for people who want to work on developing `PanACoTA` package: adding new features, correcting bugs etc.


Installing ``PanACoTA`` in development mode
===========================================

If you want to install ``PanACoTA`` while still working on modifying the scripts, type:

.. code-block:: bash

    ./make develop

Your changes will then be taken into account. As you installed the package, you will be able to run it from any directory in your computer.

If you don't want to install the software, you can still work on it, test it, and contribute to the tests and documentation by only installing the libraries needed for the software, and those
needed for development by running:

.. code-block:: bash

    pip3 install -r requirements.txt  # dependencies used by PanACoTA
    pip3 install -r requirements-dev.txt  # libraries used to run tests, generate documentation etc.

.. note:: biopython is only used for 'tree' subcommand, with option ``--soft fastme`` or ``--soft quicktree``. If you do not plan to use this, you do not need to install biopython. You can comment  the ``biopython>=1.60`` line in ``requirements.txt`` (add a ``#`` at the beginning of the line).


Running Tests
=============

If you want to work on the scripts, you can use the tests provided with the software, used to check each of its functionalities. Tests are done with `pytest <https://docs.pytest.org/en/latest/>`_ framework.

To run the tests, run, from the root of the project::

    py.test test/test_unit
    py.test test/test_functional

or, if you did not install the package::

    PYTHONPATH+=. py.test test/test_unit
    PYTHONPATH+=. py.test test/test_functional

Add ``-v`` to get more detailed information on each test run.

If you want to run only a specific test file, run::

    py.test test/test_<unit or functional>/<test_file.py>

If you want to run only a specific test, run::

    py.test test/test_<unit or functional>/<test_file.py>::<test_name>

When you run tests (all of them or individual ones), it will also always generate the coverage report. Open ``htmlcov/index.html`` on your browser if you want to check code coverage of your new function/module. The online version can be found `here <http://aperrin.pages.pasteur.fr/pipeline_annotation/htmlcov/>`_, and is automatically updated at each push on master branch.

.. warning:: If you add new features, or modify existing scripts please complete/update the tests!

We created one test file per module. If you create a new module, please create the corresponding test file in ``test/test_unit`` for unit tests, and ``test/test_functional`` for functional tests. The only condition is that your test filename must start with ``test_``, and each test name must start with ``test_``. This is required so that they are automatically run with ``py.test``, and is useful to differentiate tests and helper functions. If you modified a module, please modify/update the corresponding tests.


Contributing to documentation
=============================

This documentation is generated with `sphinx <http://www.sphinx-doc.org/en/stable/>`_. You can add your contribution to it. To generate the html documentation locally, go to ``doc/sources`` directory, and run:

.. code-block:: bash

    make html

Then, open ``doc/build/html/index.html`` on your browser.

The online version will be automatically updated when modifications on rst files are pushed on master branch.
