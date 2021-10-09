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
Functions to use proteinortho to create a pangenome

@author RedSnail
October 2021
"""
import logging
import os

from PanACoTA import utils
from PanACoTA import utils_pangenome as utils_pan
from abc import ABC, abstractmethod, abstractproperty

class Clusterisator(ABC):
    def __init__(self, threads, outdir, prt_path, panfile, quiet):
        self.threads = threads
        self.outdir = outdir
        self.prt_path = prt_path
        self.quiet = quiet
        self.logger = logging.getLogger(f"pangenome.{self.method_name}")

        print(prt_path)
        self.prt_bank = os.path.basename(self.prt_path)
        self.infoname = self.get_info()
        self.log_path = os.path.join(outdir, f"{self.method_name}_" + self.prt_bank + "_" + self.infoname + ".log")
        self.tmpdir = os.path.join(self.outdir, "tmp_" + self.prt_bank + "_" + self.infoname)

        if not panfile:
            self.panfile = f"PanGenome-{self.prt_bank}-clust-{self.infoname}.lst"

        self.panfile = os.path.join(self.outdir, self.panfile)

    @property
    @abstractmethod
    def method_name(self) -> str: pass

    def check_installed(self):
        return utils.check_installed(self.method_name)

    @abstractmethod
    def info_string(self) -> str: pass

    @abstractmethod
    def get_info(self) -> str: pass

    @abstractmethod
    def arrange_files(self, tmp_dir) -> bool: pass

    @abstractmethod
    def do_pangenome(self, status) -> tuple: pass

    def run(self):
        # Get general information and file/directory names
        information = self.info_string()
        if self.threads > 1:
            information += f"\n\t- {self.threads} threads"
        self.logger.info(information)

        # If pangenome file already exists, read it to get families
        if os.path.isfile(self.panfile):
            logger.warning(f"Pangenome file {self.panfile} already exists. PanACoTA will read it to get families.")
            _, families, _ = utils_pan.read_pan_file(self.panfile, self.logger)
        else:
            os.makedirs(self.tmpdir, exist_ok=True)
            status = self.arrange_files(self.tmpdir)

            families = self.do_pangenome(status)

        return families, self.panfile


