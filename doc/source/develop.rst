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


Running Tests
=============

If you want to work on the scripts, you can use the tests provided with the software, used to check each of its functionalities. Tests are done with ``pytest`` package.

To run the tests, run, from the root of the project::

    PYTHONPATH+=. py.test test/test_unit
    PYTHONPATH+=. py.test test/test_functional

or, if you installed the package (final or development mode)::

    py.test test/test_unit
    py.test test/test_functional

If you add new features, please complete the tests!

