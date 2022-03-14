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


import sys
from textwrap import dedent

from PanACoTA import __version__ as version

from PanACoTA.subcommands import all_modules
from PanACoTA.subcommands import prepare
from PanACoTA.subcommands import annotate
from PanACoTA.subcommands import pangenome
from PanACoTA.subcommands import corepers
from PanACoTA.subcommands import align
from PanACoTA.subcommands import tree

def main():
    """
    Start program according to arguments given by user.
    """
    action, args = parse_arguments(sys.argv[1:])
    action(args)


def parse_arguments(argv):
    """
    Extract command-line arguments for different actions.
    """
    import argparse

    header = '''



 ___                 _____  ___         _____  _____
(  _`\              (  _  )(  _`\      (_   _)(  _  )
| |_) )  _ _   ___  | (_) || ( (_)   _   | |  | (_) |
| ,__/'/'_` )/' _ `\|  _  || |  _  /'_`\ | |  |  _  |
| |   ( (_| || ( ) || | | || (_( )( (_) )| |  | | | |
(_)   `\__,_)(_) (_)(_) (_)(____/'`\___/'(_)  (_) (_)


       Large scale comparative genomics tools

     -------------------------------------------
     '''

    footer = "For more details, see PanACoTA documentation."

    # Create main parser
    # TITLE with ascii art PanACoTA
    # footer for doc
    parser = argparse.ArgumentParser(epilog=footer,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=dedent(header))


    parser.add_argument('-V', '--version', action='version',
                        version='PanACoTA - v. ' + str(version),
                        help="Print the version number and exit")

    # Create subparsers, for all submodules
    subparsers = parser.add_subparsers(dest='subparser_called')
    # dest: to be able to get the subparser called with args.subparser_called
    actions = {}  # to add the action to do according to the subparser called
    checks = {}  # to add the function to call to check the subparser arguments

    # Running all modules at once. Start with ASCII art title, + small description of subcommand
    parser_all = subparsers.add_parser('all',
                                        formatter_class=argparse.RawDescriptionHelpFormatter,
                                        description=(dedent(header) +
                                        "\n=> Run all PanACoTA modules"),
                                        epilog=footer,
                                        help="Run all PanACoTA modules",
                                        add_help=False)
    all_modules.build_parser(parser_all)
    actions["all"] = all_modules.main_from_parse
    checks["all"] = all_modules.check_args
    # checks["all_modules"] = all_modules.check_args

    # Preparation part. Start with ASCII art title, + small description of subcommand
    parser_prepare = subparsers.add_parser('prepare',
                                            formatter_class=argparse.RawDescriptionHelpFormatter,
                                            description=(dedent(header) +
                                            "\n=> Prepare draft dataset"),
                                            epilog=footer,
                                            help="Prepare draft dataset",
                                            add_help=False)
    prepare.build_parser(parser_prepare)
    actions["prepare"] = prepare.main_from_parse
    checks["prepare"] = prepare.check_args

    # QC and annotation part. Start with ASCII art title, + small description of subcommand
    parser_annotate = subparsers.add_parser('annotate',
                                            formatter_class=argparse.RawDescriptionHelpFormatter,
                                            description=(dedent(header) +
                                            "\n=> Quality control and annotation of genomes"),
                                            epilog=footer,
                                            help="Quality control and annotation of genomes",
                                            add_help=False)
    annotate.build_parser(parser_annotate)
    actions["annotate"] = annotate.main_from_parse
    checks["annotate"] = annotate.check_args

    # Pan genome part. Start with ASCII art title, + small description of subcommand
    parser_pan = subparsers.add_parser('pangenome',
                                       formatter_class=argparse.RawDescriptionHelpFormatter,
                                       description=(dedent(header) +
                                       "\n=> Generate a pan-genome of your dataset"),
                                       epilog=footer, help="Generate a pan-genome of your dataset",
                                       add_help=False)
    pangenome.build_parser(parser_pan)
    actions["pangenome"] = pangenome.main_from_parse

    # Persistent genome part. Start with ASCII art title, + small description of subcommand
    parser_corepers = subparsers.add_parser('corepers',
                                            formatter_class=argparse.RawDescriptionHelpFormatter,
                                            description=(dedent(header) +
                                            "\n=> Compute a Core or Persistent genome of your "
                                            "dataset"),
                                            epilog=footer,
                                            help="Compute a Core or Persistent genome of your "
                                                 "dataset",
                                            add_help=False)
    corepers.build_parser(parser_corepers)
    actions["corepers"] = corepers.main_from_parse
    checks["corepers"] = corepers.check_args

    # Alignment part. Start with ASCII art title, + small description of subcommand
    parser_align = subparsers.add_parser('align',
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description=(dedent(header) +
                                         "\n=> Align Core/Persistent families"),
                                         epilog=footer,
                                         help="Align Core/Persistent families",
                                         add_help=False)
    align.build_parser(parser_align)
    actions["align"] = align.main_from_parse

    # tree part. Start with ASCII art title, + small description of subcommand
    parser_tree = subparsers.add_parser('tree',
                                        formatter_class=argparse.RawDescriptionHelpFormatter,
                                        description=(dedent(header) +
                                        "\n=> Infer phylogenetic tree based on "
                                        "core/persistent genome"),
                                        epilog=footer,
                                        help="Infer phylogenetic tree based on "
                                             "core/persistent genome",
                                        add_help=False)
    tree.build_parser(parser_tree)
    actions["tree"] = tree.main_from_parse
    checks["tree"] = tree.check_args

    # Parse arguments and execute corresponding action
    arguments = parser.parse_args(argv)
    arguments.argv = argv
    action_called = arguments.subparser_called
    # If checks are needed, do it (if some arguments are not compatible etc.)
    if action_called in checks:
        checks[action_called](parser, arguments)

    # If subparser called does not exist, error
    if action_called not in actions:
        parser.error("too few arguments. Use '-h' to get help.")
    return actions[action_called], arguments


if __name__ == '__main__':
    main()
