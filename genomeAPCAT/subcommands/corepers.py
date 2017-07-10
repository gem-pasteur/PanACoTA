#!/usr/bin/env python3
# coding: utf-8

"""
corepers is a subcommand of genomeAPCAT


@author gem
June 2017
"""

import os


def main_from_parse(args):
    """
    Call main function from the arguments given by parser
    """
    main(args.pangenome, args.tol, args.multi, args.mixed, outputfile=args.outfile)


def main(pangenome, tol, multi, mixed, outputfile=None):
    """ Read pangenome and deduce Persistent genome according to the user criteria

    pangenome: file containing pangenome
    tol: min number of genomes present in a family to consider it as persistent
    multi: True if multigenic families are allowed, False otherwise
    multi: True if mixed families are allowed, False otherwise
    outputfile: name of persistent genome file
    """
    # import needed packages
    import logging
    from genomeAPCAT import utils
    from genomeAPCAT import utils_pangenome as utilsp
    import genomeAPCAT.corepers_module.persistent_functions as pers

    # name logfile, add timestamp if already existing
    path_pan, base_pan = os.path.split(pangenome)
    logfile_base = os.path.join(path_pan, "genomeAPCAT-corepers")
    level = logging.DEBUG
    utils.init_logger(logfile_base, level, '') #, verbose=verbose, quiet=quiet)
    logger = logging.getLogger()

    # Define output filename
    if not outputfile:
        outputfile = os.path.join(path_pan, "PersGenome_" + base_pan + "_" + str(tol))
        if multi:
            outputfile += "-multi.lst"
        elif mixed:
            outputfile += "-mixed.lst"
        else:
            outputfile += ".lst"
    else:
        outputfile = os.path.join(path_pan, outputfile)
    logger.info(get_info(tol, multi, mixed))

    # Read pangenome
    fams_by_strain, families, all_strains = utilsp.read_pangenome(pangenome)
    # Generate persistent genome
    fams = pers.get_pers(fams_by_strain, families, len(all_strains), tol, multi, mixed)
    pers.write_persistent(fams, outputfile)


def get_info(tol, multi, mixed):
    """
    Get a string corresponding to the information that will be given to logger.
    """
    if tol == 1:
        return "Will generate a CoreGenome."
    else:
        toprint = ("Will generate a Persistent genome with member(s) in at least {}"
                   "% of all genomes in each family.\n").format(100*tol)
        if multi:
            toprint += ("Multigenic families are allowed (several members in "
                        "any genome of a family).")
        elif mixed:
            toprint += ("Mixed families are allowed. To be considered as persistent, "
                        "a family must have exactly 1 member in {}% of the genomes, "
                        "but in the remaining {}% genomes, there can be 0, 1 or "
                        "several members.".format(tol*100, round((1-tol)*100), 3))
        else:
            toprint += ("To be considered as persistent, a family must contain exactly 1 member "
                        "in at least {}% of all genomes.").format(tol*100)
        return toprint


def build_parser(parser):
    """
    Method to create a parser for command-line options
    """
    import argparse

    def percentage(param):
        try:
            param = float(param)
        except Exception:
            msg = "argument -t tol: invalid float value: {}".format(param)
            raise argparse.ArgumentTypeError(msg)
        if param < 0 or param > 1:
            msg = ("The minimum %% of genomes required in a family to be persistent must "
                   "be in [0, 1]. Invalid value: {}".format(param))
            raise argparse.ArgumentTypeError(msg)
        return param

    # Create command-line parser for all options and arguments to give
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-p", dest="pangenome", required=True,
                        help=("PanGenome file (1 line per family, first column is fam number)"))

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("-t", dest="tol", default=1, type=percentage,
                          help=("min %% of genomes having at least 1 member in a family to "
                                "consider the family as persistent (between 0 and 1, "
                                "default is 1 = 100%% of genomes = Core genome)"))
    optional.add_argument("-M", dest="multi", action='store_true',
                          help=("Add this option if you allow several members in any genome "
                                "of a family. By default, only 1 (or 0 if tol<1) member "
                                "per genome are allowed in all genomes. If you want to allow "
                                "exactly 1 member in 'tol'%% of the genomes, and 0, 1 "
                                "or several members in the '1-tol'%% left, use the option -X "
                                "instead of this one: -M and -X options are not compatible."))
    optional.add_argument("-X", dest="mixed", action='store_true',
                          help=("Add this option if you want to allow families having several "
                                "members only in '1-tol'%% of the genomes. In the other genomes, "
                                "only 1 member exactly is allowed. This option is not compatible "
                                "with -M (which is allowing multigenic families: having several "
                                "members in any number of genomes)."))
    optional.add_argument("-o", dest="outfile",
                          help=("You can specify an output filename (without path)"
                                " for the Persistent genome "
                                "deduced from Pan. If not given, it will be saved as "
                                "PersGenome_<pangenome>_tol[-multi][-mixed].lst"))

    helper = parser.add_argument_group('Others')
    helper.add_argument("-v", "--verbose", dest="verbose", action="count", default=0,
                        help=("Increase verbosity in stdout/stderr."))
    helper.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                        help=("Do not display anything to stdout/stderr. log files will "
                              "still be created."))
    helper.add_argument("-h", "--help", dest="help", action="help",
                        help="show this help message and exit")


def check_args(parser, args):
    """
    Check that arguments given to parser are as expected.
    """
    if args.multi and args.mixed:
        parser.error("-M and -X options cannot be activated together. Choose if you want to:\n"
                     "- allow several members in any number of genomes of a family (-M)\n"
                     "- allow several members in only '1-tol'% of the genomes of a family "
                     "(other 'tol'% genomes must have exactly 1 member) (-X)")
    if args.mixed and args.tol==1:
        parser.error("You are asking for mixed families, while asking for 100% of the genomes of "
                     "a family to have exactly one member, which is not compatible. Do you want "
                     "to \n- lower the percentage of genomes required to have exactly "
                     "1 member (-t tol)\n- not allow mixed families (remove -X option)")
    return args


def parse(parser, argu):
    """
    Parse arguments given to parser
    """
    args = parser.parse_args(argu)
    return check_args(parser, args)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=("Compute core or persistent genome"),
                                     add_help=False)
    build_parser(parser)
    OPTIONS = parse(parser, sys.argv[1:])
    main_from_parse(OPTIONS)