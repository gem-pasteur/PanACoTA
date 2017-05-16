#!/usr/bin/env python3
# coding: utf-8

"""
Installation script for genomeAPCAT
Targets available are:

- `./make.py` = `./make.py install` : install all dependencies if not already present,
and install genomeAPCAT
- `./make.py develop` : same as install, but install genomeAPCAT in development mode, to be
able to change the script and run it.
- `./make.py clean` : clean all dependencies that were installed by this script
(uninstall them, remove source and bin folders)
- `./make.py uninstall` : clean dependencies + uninstall genomeAPCAT

By default, the dependencies are installed in /usr/local/bin. If the user wants
to install them in another directory (which is in the path), he can specify it with
`--prefix <directory>`.

@author gem
May 2017
"""

import subprocess
import os
import logging
from logging.handlers import RotatingFileHandler
import sys
import shlex
import shutil
import glob


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


def clean_dependencies(install_dir):
    """
    Remove links to files in 'binaries', and remove 'binaries' and 'dependencies' folders with
    all their content.
    """
    logger.info("Cleaning dependencies...")
    binpath = os.path.join(os.getcwd(), "binaries")
    srcpath = os.path.join(os.getcwd(), "dependencies")
    # Remove barrnap downloaded archive if not already done
    if os.path.isfile("0.8.tar.gz"):
        os.remove("0.8.tar.gz")
    # Remove prokka downloaded folder if it was not moved to 'dependencies'
    shutil.rmtree("prokka", ignore_errors=True)
    # Remove barrnap initial folder if it was not moved to 'dependencies'
    shutil.rmtree("barrnap-0.8", ignore_errors=True)
    # If there are binaries in 'binaries' folder, remove their link to install_dir
    if os.path.isdir(binpath):
        for binf in glob.glob(os.path.join(binpath, "*")):
            linked = os.path.join(install_dir, os.path.basename(binf))
            if os.path.isfile(linked):
                os.remove(linked)
    # Remove 'dependencies' and 'binaries' folders
    shutil.rmtree(binpath, ignore_errors=True)
    shutil.rmtree(srcpath, ignore_errors=True)


def uninstall():
    """
    Uninstall genomeAPCAT python package
    """
    logger.info("Uninstalling genomeAPCAT...")
    cmd = "pip3 uninstall -y genomeAPCAT"
    error = "A problem occurred while trying to uninstall genomeAPCAT."
    run_cmd(cmd, error)


def install_all(install_dir, dev=False):
    """
    Install all needed software(s).
    First, check if dependencies are installed.
    If at least one dependency is not installed, install it.
    Then, install genomeAPCAT.

    dev: install genomeAPCAT in development mode if true. Otherwise, install in final mode
    """
    to_install = check_dependencies()
    if to_install != []:
        binpath = os.path.join(os.getcwd(), "binaries")
        os.makedirs("dependencies", exist_ok=True)
        os.makedirs(binpath, exist_ok=True)
        if "barrnap" in to_install:
            ret = install_barrnap()
            if ret != 0:
                logger.warning("Barrnap was not installed (see above). Prokka will "
                               "not predict RNA.")
        if "prokka" in to_install:
            install_prokka()
        logger.info("Finalizing dependencies installation...")
        for binf in glob.glob(os.path.join(binpath, "*")):
            os.symlink(binf, os.path.join(install_dir, os.path.basename(binf)))
    logger.info("Installing genomeAPCAT...")
    if dev:
        cmd = "pip3 install -e ."
    else:
        cmd = "pip3 install ."
    error = "A problem occurred while trying to install genomeAPCAT."
    run_cmd(cmd, error)


def check_dependencies():
    """
    Check that all required tools are available.

    - if prokka is not installed:
        - if barrnap is not installed: check that wget and tar are available
        - check that git is available
    - check that pip3 is available
    """
    to_install = []
    if target == "install":
        if not cmd_exists("prokka"):
            if not cmd_exists("barrnap"):
                if not cmd_exists("wget"):
                    logger.error(("You need wget to install barrnap, the RNA predictor"
                                  " used by prokka."))
                    sys.exit(1)
                if not cmd_exists("tar"):
                    logger.error(("You need tar to install barrnap, the RNA predictor"
                                  " used by prokka."))
                    sys.exit(1)
                to_install.append("barrnap")
            if not cmd_exists("git"):
                logger.error("You need git to install prokka.")
                sys.exit(1)
            to_install.append("prokka")
        if not cmd_exists("pip3"):
            logger.error("You need pip3 to install genomeAPCAT.")
            sys.exit(1)
    return to_install


def install_barrnap():
    """
    Install barrnap, the RNA predictor used by prokka
    """
    logger.info("Installing barrnap...")
    cmd = "wget https://github.com/tseemann/barrnap/archive/0.8.tar.gz"
    error = "A problem occurred while trying to download barrnap. See log above."
    ret = run_cmd(cmd, error)
    if ret != 0:
        return ret
    cmd = "tar -xf 0.8.tar.gz"
    error = "A problem occurred while untaring barrnap. See log above."
    ret = run_cmd(cmd, error)
    if ret != 0:
        return ret
    cmd = "rm 0.8.tar.gz"
    error = "A problem occurred while removing barrnap archive. See log above."
    ret = run_cmd(cmd, error)
    binpath = os.path.join(os.getcwd(), "binaries")
    srcpath = os.path.join(os.getcwd(), "dependencies")
    cmd = "mv barrnap-0.8 " + srcpath
    error = ("A problem occurred while moving barrnap package to "
             "dependencies folder. See log above.")
    ret = run_cmd(cmd, error)
    if ret != 0:
        return ret
    os.symlink(os.path.join(srcpath, "barrnap-0.8", "bin", "barrnap"),
               os.path.join(binpath, "barrnap"))
    return 0


def install_prokka():
    """
    Install prokka
    """
    logger.info("Installing prokka...")
    cmd = "git clone https://github.com/tseemann/prokka.git"
    error = "A problem occurred while trying to download prokka. See log above."
    ret = run_cmd(cmd, error, eof=True)
    cmd = "mv prokka dependencies"
    error = "A problem occurred while moving prokka to 'dependencies'. See log above."
    ret = run_cmd(cmd, error, eof=True)
    binpath = os.path.join(os.getcwd(), "binaries")
    srcpath = os.path.join(os.getcwd(), "dependencies")
    cmd = os.path.join(srcpath, "prokka", "bin", "prokka") +  " --setupdb"
    error = "A problem occurred while initializing prokka db. See log above."
    ret = run_cmd(cmd, error, eof=True)
    os.symlink(os.path.join(srcpath, "prokka", "bin", "prokka"),
               os.path.join(binpath, "prokka"))


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


def parse():
    """
    Method to create a parser for command-line options
    """
    import argparse
    parser = argparse.ArgumentParser(description=("Script to install, clean or uninstall "
                                                  "genomeAPCAT"))
    # Create command-line parser for all options and arguments to give
    targets = ['install', 'develop', 'clean', 'uninstall']
    parser.add_argument("target", default='install', choices=targets, nargs='?',
                        help=("Choose what you want to do:\n"
                              " - install: install genomeAPCAT and its dependencies. If not "
                              "already installed by user, dependencies packages "
                              "will be downloaded "
                              "and built in 'dependencies' folder, and their binary files will be "
                              "put to 'binaries' folder.\n"
                              " - develop: same as install, but genomeAPCAT will be installed "
                              "in development mode, so that you can modify the script and "
                              "take the changes into account while running.\n"
                              " - clean: clean dependencies: for dependencies "
                              "which were installed "
                              "via this script, uninstall them, and remove their downloaded "
                              "sources from 'dependencies' folder. Can be used if the user wants "
                              "to install another version of the dependencies.\n"
                              " - uninstall: uninstall genomeAPCAT, as well as the dependencies "
                              "which were installed for it (in 'dependencies' folder).\n"
                              "Default is %(default)s."))
    parser.add_argument("--prefix", dest="install_dir",
                        help=("By default, all scripts will be installed in /usr/local/bin. "
                              "If you want to install them in another directory, give it "
                              "with this --prefix option."))
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # Init logger
    logger = logging.getLogger()
    level = logging.INFO
    # logging.basicConfig(level=level,
                        # datefmt='%a, %d %b %Y %H:%M:%S')
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
    if not OPTIONS.install_dir:
        install_dir = os.path.join(os.sep + "usr", "local", "bin")
    else:
        install_dir = OPTIONS.install_dir
    target = OPTIONS.target
    # Execute target
    if target == "install":
        # Check that installation directory is in $PATH
        check_path(install_dir)
        install_all(install_dir)
    elif target == "develop":
        # Check that installation directory is in $PATH
        check_path(install_dir)
        install_all(install_dir, dev=True)
    elif target == "clean":
        clean_dependencies(install_dir)
    elif target == "uninstall":
        clean_dependencies(install_dir)
        uninstall()

