#!/usr/bin/env python3
# coding: utf-8

"""
Setup script
"""

import genomeAPCAT
try:
    from setuptools import setup
    from setuptools.command.test import test as TestCommand
except ImportError:
    from distutils.core import setup
    from distutils.core import Command as TestCommand


class PyTest(TestCommand):
    def initialize_options(self):
        pass

    def finalize_options(self):
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        import sys
        errno = pytest.main(self.test_args)
        sys.exit(errno)


def parse_requirements(requirements):
    """
    Returns the list of requirements found in the requirements file
    """
    with open(requirements, "r") as req_file:
        return [l.strip('\n') for l in req_file if l.strip('\n')
                and not l.startswith('#')]

packages = ['genomeAPCAT', 'genomeAPCAT.qcAnnote', 'genomeAPCAT.subcommands']
requires = parse_requirements("requirements.txt")
scripts = ['genomeAPCAT/subcommands/qc_and_annote.py']

classifiers = [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: ???",
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

with open('README.md') as f:
    long_description = f.read()

setup(
    name='genomeAPCAT',
    packages=packages,
    version=genomeAPCAT.__version__,
    description="Large scale comparative genomics tools: annotate genomes, do pangenome, "
                "core/persistent genome, align core/persistent families, infer phylogenetic tree.",
    long_description=long_description,
    author='Amandine Perrin',
    author_email='amandine.perrin@pasteur.fr',
    license='???',
    platforms='OS Independent',
    package_data={'': ['LICENSE']},
    download_url='??',
    url='??',
    scripts=scripts,
    include_package_data=True,
    install_requires=requires,
    tests_require=['pytest'],
    cmdclass={'test': PyTest},
    classifiers=classifiers
)
