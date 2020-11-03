#!/usr/bin/env python3
# coding: utf-8

# ###############################################################################
# This file is part of PanACOTA.                                                #
#                                                                               #
# Authors: Amandine Perrin                                                      #
# Copyright © 2018-2020 Institut Pasteur (Paris).                               #
# See the COPYRIGHT file for details.                                           #
#                                                                               #
# PanACOTA is a software providing tools for large scale bacterial comparative  #
# genomics. From a set of complete and/or draft genomes, you can:               #
#    -  Do a quality control of your strains, to eliminate poor quality         #
# genomes, which would not give any information for the comparative study       #
#    -  Uniformly annotate all genomes                                          #
#    -  Do a Pan-genome                                                         #
#    -  Do a Core or Persistent genome                                          #
#    -  Align all Core/Persistent families                                      #
#    -  Infer a phylogenetic tree from the Core/Persistent families             #
#                                                                               #
# PanACOTA is free software: you can redistribute it and/or modify it under the #
# terms of the Affero GNU General Public License as published by the Free       #
# Software Foundation, either version 3 of the License, or (at your option)     #
# any later version.                                                            #
#                                                                               #
# PanACOTA is distributed in the hope that it will be useful, but WITHOUT ANY   #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS     #
# FOR A PARTICULAR PURPOSE. See the Affero GNU General Public License           #
# for more details.                                                             #
#                                                                               #
# You should have received a copy of the Affero GNU General Public License      #
# along with PanACOTA (COPYING file).                                           #
# If not, see <https://www.gnu.org/licenses/>.                                  #
# ###############################################################################

"""
corepers is a subcommand of PanACoTA

Generate a core genome (families containing 1 member in all genomes of the dataset)
or a persistent genome (families with a given % of genomes having exactly 1 member).
You can also allow:

- mixed families: exactly 1 member in the given percentage of genomes, but the other genomes
  can contain 0 or several members
- multi families: allow several members in any genome.


@author gem
June 2017
"""

import os

import sys


def main_from_parse(args):
    """
    Call main function from the arguments given by parser

    Parameters
    ----------
    args : argparse.Namespace
        result of argparse parsing of all arguments in command line
    """
    cmd = "PanACoTA " + ' '.join(args.argv)
    main(cmd, args.pangenome, args.tol, args.multi, args.mixed, args.outputdir,
         floor=args.floor, verbose=args.verbose, quiet=args.quiet)


def main(cmd, pangenome, tol, multi, mixed, outputdir, floor=False, verbose=0, quiet=False):
    """
    Read pangenome and deduce Persistent genome according to the user criteria

    Parameters
    ----------
    pangenome : str
        file containing pangenome
    tol : float
        min % of genomes present in a family to consider it as persistent (between 0 and 1)
    multi : bool
        True if multigenic families are allowed, False otherwise
    mixed : bool
        True if mixed families are allowed, False otherwise
    outputdir : str or None
        Specific directory for the generated persistent genome. If not given, pangenome directory is used.
    floor : bool
        Require at least floor(nb_genomes*tol) genomes if True, ceil(nb_genomes*tol) if False
    verbose : int
        verbosity:
        - defaut 0 : stdout contains INFO, stderr contains ERROR.
        - 1: stdout contains INFO, stderr contains WARNING and ERROR
        - 2: stdout contains (DEBUG), DETAIL and INFO, stderr contains WARNING and ERROR
        - >=15: Add DEBUG in stdout
    quiet : bool
        True if nothing must be sent to stdout/stderr, False otherwise
    """
    # import needed packages
    import logging
    from PanACoTA import utils
    from PanACoTA import utils_pangenome as utilsp
    import PanACoTA.corepers_module.persistent_functions as pers
    from PanACoTA import __version__ as version

    # get pangenome name info
    _, base_pan = os.path.split(pangenome)
    # Define output filename
    output_name = "PersGenome_" + base_pan + "_"
    if floor:
        output_name += "F"
    output_name += str(tol)
    if multi:
        output_name += "-multi.lst"
    elif mixed:
        output_name += "-mixed.lst"
    else:
        output_name += ".lst"
    # Define output directory and filename path
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    outputfile = os.path.join(outputdir, output_name)
    logfile_base = os.path.join(outputdir, "PanACoTA-corepers")
    # level is the minimum level that will be considered.
    # for verbose = 0 or 1, ignore details and debug, start from info
    if verbose <= 1:
        level = logging.INFO
    # for verbose = 2, ignore only debug
    if verbose >= 2 and verbose < 15:
        level = 15 # int corresponding to detail level
    # for verbose >= 15, write everything
    if verbose >= 15:
        level = logging.DEBUG
    utils.init_logger(logfile_base, level, 'corepers', verbose=verbose, quiet=quiet)
    logger = logging.getLogger("corepers")
    logger.info(f'PanACoTA version {version}')
    logger.info("Command used\n \t > " + cmd)

    logger.info(get_info(tol, multi, mixed, floor))

    # Read pangenome
    fams_by_strain, families, all_strains = utilsp.read_pangenome(pangenome, logger)
    # Generate persistent genome
    fams = pers.get_pers(fams_by_strain, families, len(all_strains), tol, multi, mixed, floor)
    # Write persistent genome to file
    pers.write_persistent(fams, outputfile)
    logger.info("Persistent genome step done.")
    return outputfile


def get_info(tol, multi, mixed, floor):
    """
    Get a string corresponding to the information that will be given to logger.

    Parameters
    ----------
    tol : float
        min % of genomes present in a family to consider it as persistent (between 0 and 1)
    multi : bool
        True if multigenic families are allowed, False otherwise
    mixed : bool
        True if mixed families are allowed, False otherwise
    floor : bool
        Require at least floor(nb_genomes*tol) genomes if True, ceil(nb_genomes*tol) if False

    Returns
    -------
    str
        Information to give to logger
    """
    if tol == 1:
        return "Will generate a CoreGenome."
    else:
        if floor:
            floorstr = "floor"
        else:
            floorstr = "ceil"
        toprint = (f"Will generate a Persistent genome with member(s) in at least {100*tol}"
                   f"% of all genomes in each family.\n")
        if multi:
            toprint += ("Multigenic families are allowed (several members in "
                        "any genome of a family).")
        elif mixed:
            toprint += ("Mixed families are allowed. To be considered as persistent, "
                        f"a family must have exactly 1 member in {tol*100}% of the genomes, "
                        f"but in the remaining {round((1-tol)*100,3)}% genomes, there can be 0, 1 or "
                        "several members.")
        else:
            toprint += ("To be considered as persistent, a family must contain exactly 1 member "
                        f"in at least {tol*100}% of all genomes. The other genomes are absent from the "
                        "family.")
        return toprint


def build_parser(parser):
    """
    Method to create a parser for command-line options

    Parameters
    ----------
    parser : argparse.ArgumentParser
        parser to configure in order to extract command-line arguments
    """
    import argparse
    from PanACoTA import utils_argparse

    # Create command-line parser for all options and arguments to give
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-p", dest="pangenome", required=True,
                          help="PanGenome file (1 line per family, first column is fam number)")
    required.add_argument("-o", dest="outputdir", required=True,
                          help=("Specify the output directory for your core/persistent genome."),
                          default=".")
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-t", "--tol", dest="tol", default=1, type=utils_argparse.percentage,
                          help=("min %% of genomes having at least 1 member in a family to "
                                "consider the family as persistent (between 0 and 1, "
                                "default is 1 = 100%% of genomes = Core genome)."
                                "By default, the minimum number of genomes will be "
                                "ceil('tol'*N) (N being the total number of genomes). If "
                                "you want to use floor('tol'*N) instead, add the '-F' option."))
    optional.add_argument("-M", dest="multi", action='store_true',
                          help=("Add this option if you allow several members in any genome "
                                "of a family. By default, only 1 (or 0 if tol<1) member "
                                "per genome are allowed in all genomes. If you want to allow "
                                "exactly 1 member in 'tol'%% of the genomes, and 0, 1 "
                                "or several members in the '1-tol'%% left, use the option -X "
                                "instead of this one: -M and -X options are not compatible."))
    optional.add_argument("-X", dest="mixed", action='store_true',
                          help="Add this option if you want to allow families having several "
                               "members only in '1-tol'%% of the genomes. In the other genomes, "
                               "only 1 member exactly is allowed. This option is not compatible "
                               "with -M (which is allowing multigenic families: having several "
                               "members in any number of genomes).")
    optional.add_argument("-F", dest="floor", action="store_true",
                          help="When you specify the '-tol' option, with a number lower "
                               "than 1, you can add this option to use floor('tol'*N) "
                               "as a minimum number of genomes instead of ceil('tol'*N) "
                               "which is the default behavior.")

    helper = parser.add_argument_group('Others')
    helper.add_argument("-v", "--verbose", dest="verbose", action="count", default=0,
                        help="Increase verbosity in stdout/stderr.")
    helper.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                        help=("Do not display anything to stdout/stderr. log files will "
                              "still be created."))
    helper.add_argument("-h", "--help", dest="help", action="help",
                        help="show this help message and exit")


def check_args(parser, args):
    """
    Check that arguments given to parser are as expected.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser used to parse command-line
    args : argparse.Namespace
        Parsed arguments

    Returns
    -------
    argparse.Namespace or None
        The arguments parsed, updated according to some rules. Exit program
        with error message if error occurs with arguments given.
    """
    if args.multi and args.mixed:
        parser.error("-M and -X options cannot be activated together. Choose if you want to:\n"
                     "- allow several members in any number of genomes of a family (-M)\n"
                     "- allow several members in only '1-tol'% of the genomes of a family "
                     "(other 'tol'% genomes must have exactly 1 member) (-X)")
    if args.mixed and args.tol == 1:
        parser.error("You are asking for mixed families, while asking for 100% of the genomes of "
                     "a family to have exactly one member, which is not compatible. Do you want "
                     "to \n- lower the percentage of genomes required to have exactly "
                     "1 member (-t tol)\n- not allow mixed families (remove -X option)")
    if args.floor and args.tol == 1:
        parser.error("You are asking to use floor('tol'*N) as a minimum number of genomes "
                     "present in a family, but with 'tol'=1: the minimum number of genomes "
                     "will always be equal to N, using floor or the default ceil! Either "
                     "use a 'tol' lower than 1, or remove the '-F' option.")
    return args


def parse(parser, argu):
    """
    Parse arguments given to parser

    Parameters
    ----------
    parser : argparse.ArgumentParser
        the parser used
    argu : [str]
        command-line given by user, to parse using parser

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    args = parser.parse_args(argu)
    return check_args(parser, args)


if __name__ == '__main__':
    import argparse
    myparser = argparse.ArgumentParser(description="Compute core or persistent genome",
                                       add_help=False)
    build_parser(myparser)
    OPTIONS = parse(myparser, sys.argv[1:])
    main_from_parse(OPTIONS)
