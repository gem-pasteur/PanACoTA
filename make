#!/usr/bin/env python3
# coding: utf-8

"""
Installation script for PanACoTA
Targets available are:

- `./make.py` = `./make.py install`: install all dependencies if not already present,
and install PanACoTA
- `./make.py upgrade`: upgrade PanACoTA package.
- `./make.py develop`: same as install, but install PanACoTA in development mode, to be
able to change the script and run it.
- `./make.py uninstall`: uninstall PanACoTA

@author gem
May 2017
update February 2019
"""

import subprocess
import os
import logging
from logging.handlers import RotatingFileHandler
import sys
import shlex
import shutil


def check_path(install_dir):
    """
    Check that given install dir is in $PATH. If not, close and ask user to provide
    an installation directory which is in the path, or to add /usr/local/bin to
    its path.
    """
    elems = os.environ["PATH"].split(os.pathsep)
    if install_dir not in elems:
        logger.error("Your installation directory ({}) is not in your $PATH. Please, provide "
                     "an installation directory which is in your PATH with '--prefix "
                     "<install_dir>', or add the default installation directory "
                     "(/usr/local/bin) to your PATH.".format(install_dir))
        sys.exit(1)


def uninstall():
    """
    Uninstall PanACoTA python package
    """
    logger.info("Uninstalling PanACoTA...")
    cmd = "pip3 uninstall -y PanACoTA"
    error = ("A problem occurred while trying to uninstall PanACoTA. If you have "
             "permission errors, try to add 'sudo' before your command line.")
    run_cmd(cmd, error)
    link_dest = os.path.join(os.sep + "usr", "local", "bin", "PanACoTA")
    if os.path.exists(link_dest):
        os.remove(link_dest)


def upgrade(user=False):
    """
    Upgrade PanACoTA module to the latest version. User must update the sources
    before running upgrade!
    """
    logger.info("Upgrading PanACoTA...")
    if user:
        cmd = "pip3 install --upgrade --no-deps --user ."
    else:
        cmd = "pip3 install --upgrade --no-deps ."
    error = ("An error occurred while trying to update PanACoTA. If you have "
             "permission errors, try to add 'sudo' before your command line.")
    run_cmd(cmd, error, eof=True)


def install_all(install_dir, target, dev=False, user=False):
    """
    Install PanACoTA.
    check dependencies and show missing ones
    dev: install PanACoTA in development mode if true. Otherwise, install in final mode
    user: install in user mode if True
    """
    to_install_user = check_dependencies(target)
    if user:
        logger.info("Installing PanACoTA in user mode...")
        opt = "--user"
    else:
        logger.info("Installing PanACoTA...")
        opt = ""
    if dev:
        logger.info("Installing developer packages needed for PanACoTA")
        cmd = "pip3 install " + opt + " -e ."
        cmd2 = "pip3 install -r requirements-dev.txt"
        error2 = ("Problem while trying to install developer tools. If you have "
                  "permission errors, try to add 'sudo' before your command line. If "
                  "you do not have root access, install with the '--user' option")
        run_cmd(cmd2, error2, eof=True)
    else:
        cmd = "pip3 install " + opt + " ."
    error = ("A problem occurred while trying to install PanACoTA. If you have "
             "permission errors, try to add 'sudo' before your command line. If "
             "you do not have root access, install with the '--user' option")
    run_cmd(cmd, error, eof=True)
    if user:
        gapcat_bin = os.path.join(os.getcwd(), "PanACoTA", "bin", "run_panacota.py")
        os.symlink(gapcat_bin, os.path.join(install_dir, "PanACoTA"))
    if to_install_user:
        msg = ("Some dependencies needed for some subcommands of PanACoTA are not installed. "
               "Here is the list of missing dependencies, and for what they are used. If you plan "
               "to use the subcommands hereafter, first install required dependencies:\n")

        if "mash" in to_install_user:
            msg += ("\t- mash (for prepare subcommand, to filter genomes) \n")
        if "prodigal" in to_install_user:
            msg += ("\t- prodigal : for annotate subcommand, you at least need prodigal (for "
                    "syntaxic annotation only). If you even need functional annotation, also "
                    "install prokka\n")
        if "prokka" in to_install_user:
            msg += ("\t- prokka (for annotate subcommand, with syntaxic + functional annotation). "
                    "If you only need syntaxic annotation, prodigal is enough.\n")
        if "barrnap" in to_install_user:
            msg += ("\t- barrnap. If you use Prokka for functional annotation, it will not predict"
                    " RNA.\n")
        if "mmseqs" in to_install_user:
            msg += "\t- mmseqs (for pangenome subcommand)\n"
        if "mafft" in to_install_user:
            msg += ("\t- mafft (to align persistent genomes in order to infer a phylogenetic tree "
                    "after)\n")
        if "trees" in to_install_user:
            msg += ("\t- One of the 4 following softwares, used to infer a phylogenetic tree:\n"
                    "\t\t* FastTree (see README or documentation for more information on how to "
                    "install it)\n"
                    "\t\t* FastME\n"
                    "\t\t* Quicktree\n"
                    "\t\t* IQtree (or IQtree2)")
        msg += ("See more information on how to download/install those softwares in README or in "
                "documentation.")
        logger.warning(msg)


def check_dependencies(target):
    """
    Check that all required tools are available.
    - for each step (download-prepare, annotate-format, pangenome, core-pers, align, tree), check if the tools are installed. If not, return a message for user, so that he can install them if needed.
    - check that pip3 is available
    """
    to_install = []
    to_install_user = []
    msg = None
    if target == "install" or target == "develop":
        if not cmd_exists("pip3"):
            logger.error("You need pip3 to install PanACoTA.")
            sys.exit(1)
        if not cmd_exists("mash"):
            to_install_user.append("mash")
        if not cmd_exists("barrnap"):
            to_install_user.append("barrnap")
        if not cmd_exists("prodigal") and not cmd_exists("prokka"):
            to_install_user.append("prodigal")
        if not cmd_exists("prokka"):
            to_install_user.append("prokka")
        if not cmd_exists("mmseqs"):
            to_install_user.append("mmseqs")
        if not cmd_exists("mafft"):
            to_install_user.append("mafft")
        if not cmd_exists("FastTreeMP") and not cmd_exists("fastme") and not cmd_exists(
                "quicktree") and not cmd_exists("iqtree") and not cmd_exists("iqtree2"):
            to_install_user.append("trees")
    return to_install_user


def run_cmd(cmd, error, eof=False):
    """
    Run the given command line. If the return code is not 0, print error message
    if eof (exit on fail) is True, exit program if error code is not 0.
    """
    retcode = subprocess.call(shlex.split(cmd))
    if retcode != 0:
        logger.error(error)
        if eof:
            sys.exit(retcode)
    return retcode


def cmd_exists(cmd):
    """
    Check if the command is in $PATH and can then be executed
    """
    torun = "which " + cmd
    trying = subprocess.Popen(torun.split(), stdout=subprocess.PIPE)
    out, _ = trying.communicate()
    if trying.returncode == 0:
        if os.path.isfile(out.strip()):
            return True
    return False


def rmdir(mydir):
    """
    Remove directory if existing
    """
    if os.path.isdir(mydir):
        shutil.rmtree(mydir)


def parse():
    """
    Method to create a parser for command-line options
    """
    import argparse
    parser = argparse.ArgumentParser(description=("Script to install, clean or uninstall "
                                                  "PanACoTA"))
    # Create command-line parser for all options and arguments to give
    targets = ['install', 'develop', 'upgrade', 'uninstall']
    parser.add_argument("target", default='install', choices=targets, nargs='?',
                        help=("Choose what you want to do:\n"
                              " - install: install PanACoTA and its dependencies. If not "
                              "already installed by user, dependencies packages "
                              "will be downloaded "
                              "and built in 'dependencies' folder, and their binary files will be "
                              "put to 'binaries' folder.\n"
                              " - upgrade: upgrade PanACoTA module. \n"
                              " - develop: same as install, but PanACoTA will be installed "
                              "in development mode, so that you can modify the script and "
                              "take the changes into account while running.\n"
                              " - uninstall: uninstall PanACoTA.\n"
                              "Default is %(default)s."))
    parser.add_argument("--prefix", dest="install_dir",
                        help=("By default, all scripts will be installed in /usr/local/bin. "
                              "If you want to install them in another directory, give it "
                              "with this --prefix option."))
    parser.add_argument("--user", dest="user", action="store_true",
                        help="Install package in user mode.")
    args = parser.parse_args()
    if args.user and args.target not in ["install", "upgrade"]:
        parser.error("--user option can only be used when installing (or upgrading) the package.")
    return args


if __name__ == '__main__':
    # Init logger
    logger = logging.getLogger()
    level = logging.INFO
    logger.setLevel(level)
    # create formatter for log messages: "timestamp :: level :: message"
    formatterFile = logging.Formatter('[%(asctime)s] :: %(levelname)s :: %(message)s',
                                      '%Y-%m-%d %H:%M:%S')
    formatterStream = logging.Formatter(' * [%(asctime)s] :: %(message)s', '%Y-%m-%d %H:%M:%S')
    # Create handler 1: writing to 'logfile'
    logfile = "install.log"
    open(logfile, "w").close()  # empty logfile if existing
    logfile_handler = RotatingFileHandler(logfile, 'w', 1000000, 5)
    # set level to the same as the logger level
    logfile_handler.setLevel(level)
    logfile_handler.setFormatter(formatterFile)  # add formatter
    logger.addHandler(logfile_handler)  # add handler to logger
    # Create handler 2: write to stderr
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)  # write any message
    stream_handler.setFormatter(formatterStream)
    logger.addHandler(stream_handler)  # add handler to logger

    # Get user arguments
    OPTIONS = parse()
    user_mode = OPTIONS.user
    if not OPTIONS.install_dir:
        my_install_dir = os.path.join(os.sep + "usr", "local", "bin")
    else:
        my_install_dir = OPTIONS.install_dir
    my_target = OPTIONS.target
    # Execute target
    if my_target == "install":
        # Check that installation directory is in $PATH
        check_path(my_install_dir)
        install_all(my_install_dir, my_target, user=user_mode)
    elif my_target == "develop":
        # Check that installation directory is in $PATH
        check_path(my_install_dir)
        install_all(my_install_dir, my_target, dev=True)
    elif my_target == "upgrade":
        upgrade(user=user_mode)
    elif my_target == "uninstall":
        uninstall()

    # Print end of process (to get time)
    logger.info("DONE")
