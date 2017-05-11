# pipeline_annotation

[![build status](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/build.svg)](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/commits/master)
[![coverage report](https://gitlab.pasteur.fr/aperrin/pipeline_annotation/badges/master/coverage.svg)](http://aperrin.pages.pasteur.fr/pipeline_annotation/htmlcov)

Annotate genomes and format them to gembase format.  

# Installation

totomain is written in python3. So, you need python3 (and pip3 for installation) to run it.

Its external dependencies are:
- prokka (to annotate the genomes). 

You can either install the external dependencies with the version you want, or use the installation script `make.py`, which will install the dependencies.

To be able to install the dependencies (by yourself, or with the installation script), make sure you have:
- tar
- git
- wget

Then, for prokka installation, you need to install some system packages, as well as bioperl and java, if not already done. Here are command lines to install the required packages. For more information, see [Prokka README](https://github.com/tseemann/prokka):

**Centos/Fedora/RHEL (RPM)**
```
sudo yum install perl-Time-Piece perl-XML-Simple perl-Digest-MD5 git java perl-CPAN perl-Module-Build
sudo cpan -i Bio::Perl  # if you don't have Bioperl installed (it will be tedious)
```

**Ubuntu/Debian/Mint (APT)**
```
sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
```

**Mac OS X**
```
sudo cpan Time::Piece XML::Simple Digest::MD5 Bio::Perl
```

## totomain installation (final mode)

To install the package *totomain*, and all its dependencies, from the root directory, just do:

    ./make.py 

or 

    ./make.py install

You will then be able to use the package from any directory in your computer,
just as any other software.

Warning: If you plan to work on the scripts, or to download a new version after, choose the development installation (see below).


## totomain installation (development mode)

If you want to install the package while still working on modifying the scripts, or being able to download and run latest versions, just do:

    ./make.py develop

Your changes will then be taken into account. As you installed the package, you will be able to run the pipeline from any directory in your computer.

## Uninstalling totomain

If you don't want totomain anymore, or if you want to install a newer version, uninstall it by typing:

    ./make.py uninstall

## Update to new version

If you want to install a new version of totomain:
- uninstall the previous version (`./make.py uninstall`)
- download the new version
- install the new version (`./make.py`)

## Cleaning dependencies

If you installed the dependencies (such as prokka) via our installation script, but now want to install your own version, you can remove all dependencies downloaded and installed by `make.py` by doing:

    ./make.py clean

# Running totomain

## Input/Output

## Fast run

## Examples


# Development

This part is for people who want to work on developing totomain package.

## Running Tests

If you want to work on pipeline scripts, you can use the tests provided with the software, used to check each of its functionalities. To run the tests, run, from the root of the project:

    PYTHONPATH+=. py.test test/test_unit
    PYTHONPATH+=. py.test test/test_functional

or, if you installed the package (final or development mode)::

    py.test test/test_unit
    py.test test/test_functional
