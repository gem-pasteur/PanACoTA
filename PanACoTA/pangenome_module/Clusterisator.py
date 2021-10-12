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
Abstract class for launching different clusterisators.

@author RedSnail
October 2021
"""
import logging
import os
import progressbar
import threading
import time

from PanACoTA import utils
from PanACoTA import utils_pangenome as utils_pan
from abc import ABC, abstractmethod, abstractproperty

def infinite_progressbar(func):
    def wrapper(self):
        try:
            stop_bar = False
            if self.quiet:
                widgets = []
            # If not quiet, start a progress bar while clustering proteins. We cannot guess
            # how many time it will take, so we start an "infinite" bar, and send it a signal
            # when it has to stop. If quiet, we start a thread that will immediatly stop
            else:
                widgets = [progressbar.BouncingBar(marker=progressbar.RotatingMarker(markers="◐◓◑◒")),
                            "  -  ", progressbar.Timer()]
            x = threading.Thread(target=utils.thread_progressbar, args=(widgets, lambda: stop_bar,))
            x.start()
            res = func(self)
            # except KeyboardInterrupt: # pragma: no cover
        except:  # pragma: no cover
            stop_bar = True
            x.join()
            sys.exit(1)

        stop_bar = True
        x.join()
        return res

    return wrapper

class Clusterisator(ABC):
    """
    The Clusterisator object is created for each launch of clusterisator (mmseq or proteinortho) and stores parameters
    of the launch

    Parameters
    ----------
    threads : int
        number of threads which can be used.
    outdir : str
        directory where output cluster file must be saved.
    prt_path : str
        path to file containing all proteins to cluster.
    panfile : str or None
        name for output pangenome file. Otherwise, will use default name
    quiet : bool
        True if nothing must be written on stdout, False otherwise.

    Attributes
    ----------
        threads : int
            Where threads parameter is stored.
        outdir : str
            Where outdir parameter is stored.
        prt_path : str
            Where prt_path parameter is stored.
        quiet : bool
            Where quiet parameter is stored.
        logger : logging.Logger
            A logger object which is used by methods for logging.
        log_path : str
            Path to write logs.
        tmpdir : str
            Path to store temporary files.
        panfile : str
            Name for output pangenome file.
    """
    def __init__(self, threads, outdir, prt_path, panfile, quiet):
        self.threads = threads
        self.outdir = outdir
        self.prt_path = prt_path
        self.quiet = quiet
        self.logger = logging.getLogger(f"pangenome.{self.method_name}")

        prt_bank = os.path.basename(self.prt_path)
        self.log_path = os.path.join(outdir, f"{self.method_name}_" + prt_bank + "_" + self.infoname + ".log")
        self.tmpdir = os.path.join(self.outdir, "tmp_" + prt_bank + "_" + self.infoname)

        if not panfile:
            self.panfile = f"PanGenome-{prt_bank}-clust-{self.infoname}.lst"

        self.panfile = os.path.join(self.outdir, self.panfile)

    @property
    @abstractmethod
    def method_name(self):
        """
        Getter for the name of the method that is used for clusterisation.

        Returns
        -------
        method_name : str
            Name of the command that is used for clustering and pangenome building.
        """
        pass

    def check_installed(self):
        """
        Check if required software is installed.

        Returns
        -------
        installed : bool
            Is it installed.
        """
        return utils.check_installed(self.method_name)

    @property
    @abstractmethod
    def info_string(self):
        """
        A string with method name and method-specific parameters information that is printed to logs.

        Returns
        -------
        info_string : str
            Method-specific info string for logging
        """
        pass

    @property
    def info_with_threads(self):
        """
        A string with method name, parameters and number of threads that is printed to logs.
        Returns
        -------
        info_string : str
            An info string that is printed to logs
        """
        return self.info_string if self.threads == 1 else self.info_string + f"\n\t- {self.threads} threads"
    
    @property
    @abstractmethod
    def infoname(self):
        """
        A string with information about parameters that is used for temporary files and directories naming

        Returns
        -------
        infoname : str
            An info string that is used for temporary files' and dirs' naming
        """
        pass

    @property
    @abstractmethod
    def expected_files(self):
        """
        Paths of the files that are required for launch of clustering.

        Returns
        -------
        expected : List[str]
            list of paths
        """
        pass

    @property
    @abstractmethod
    def tmp_files_cmd(self):
        """
        Tuple of command and error message that are used for tmp files creation.

        Returns
        -------
        cmd : str
            bash command
        msg : str
            error message
        """
        pass

    @infinite_progressbar
    def create_tmp_files(self):
        """
        Creates temporary files and directories that are required for clustering method launch.

        Returns
        -------
        status : bool
            If files were re-done.
        """
        files_existing = list(filter(os.path.isfile, self.expected_files))

        if len(files_existing) == len(self.expected_files):
            self.logger.warning(f"{self.method_name} temporary files already exists. The program will "
                                "use them.")
            return False

        if len(files_existing) > 0:
            self.logger.warning(f"Some, but not all {self.method_name} temporary files exist, the program will remove"
                                f"remaining.")
            for file in files_existing:
                os.remove(file)
                logger.details(f"Removing '{file}'.")


        self.logger.debug(f"Existing files: {len(files_existing)}")
        self.logger.debug(f"Expected files: {len(self.expected_files)}")
        cmd, msg = self.tmp_files_cmd

        self.logger.details(f"{self.method_name} command: {cmd}")
        with open(self.log_path, "w") as logf:
            utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
        return True

    @property
    @abstractmethod
    def clustering_files(self):
        """
        Files that are created by clustering algorithm

        Returns
        -------
        downstream_files : List[str]
        """
        pass

    def do_pangenome(self, status):
        """
        Run clustering command and parse the result.

        Parameters
        ----------
        status : bool
            If the files were re-done.

        Returns
        -------
        families : {fam_num: [all members]}
            List of families
        """
        present_downstream_files = list(filter(os.path.isfile, self.clustering_files + [self.panfile]))

        if len(present_downstream_files) != len(self.clustering_files) + 1:
            status = True

        if status and len(present_downstream_files) > 0:
            self.logger.details("Removing existing clustering and/or pangenome files.")
            map(utils.remove, present_downstream_files)

        if status:
            self.logger.info("Clustering proteins...")
            self.run_clust()
        else:
            self.logger.warning((f"mmseqs clustering {self.mmseqclust} already exists. The program will now convert "
                                 "it to a pangenome file."))

        families = self.parse_to_pangenome() # here should edit parsing
        return families

    @infinite_progressbar
    def run_clust(self):
        cmd, msg = self.clust_cmd
        with open(self.log_path, "a") as logm:
            utils.run_cmd(cmd, msg, eof=False, stdout=logm, stderr=logm)

    @property
    @abstractmethod
    def clust_cmd(self):
        """
        Command to perform clustering and error message in case of failure

        Returns
        -------
        cmd : str
            command that is used for clustering
        msg : str
            an error message that is displayed in case of failure
        """
        pass

    @abstractmethod
    def parse_to_pangenome(self):
        """
        Parse clustering results to pangenome.

        Returns
        -------
        families : List[List[str]]
            List of families
        """
        pass

    def write_panfile(self, families):
        """
        Write panfile.

        Parameters
        ----------
        families

        Returns
        -------
        families : {fam_num: [all members]}
            List of families
        """
        with open(self.panfile, "w") as panf:
            for i, fam in enumerate(families):
                panf.write(" ".join([str(i)] + fam))

    def run(self):
        """
        Run all steps to build a pangenome:

        - create temporary files required for launch from protein bank
        - cluster proteins
        - convert to pangenome
        Returns
        -------
        families : {fam_num: [all members]}
            List of families in format {fam_num: [all members]}
        outfile : str
            Pangenome file path.
        """
        start = time.strftime('%Y-%m-%d_%H-%M-%S')
        with open(self.log_path, "a") as logm:
            logm.write(f"Start: {start}")

        # Get general information and file/directory names
        self.logger.info(self.info_with_threads)

        # If pangenome file already exists, read it to get families
        if os.path.isfile(self.panfile):
            self.logger.warning(f"Pangenome file {self.panfile} already exists. PanACoTA will read it to get families.")
            _, families, _ = utils_pan.read_pan_file(self.panfile, self.logger)
        else:
            os.makedirs(self.tmpdir, exist_ok=True)
            self.logger.info("Creating temporary files")
            status = self.create_tmp_files()

            families = self.do_pangenome(status)

        self.write_panfile(families)

        end = time.strftime('%Y-%m-%d_%H-%M-%S')
        with open(self.log_path, "a") as logm:
            logm.write(f"End: {end}")

        return dict(zip(range(1, len(families) + 1), families)), self.panfile

