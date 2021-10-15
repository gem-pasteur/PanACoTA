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
A class for launching mmseqs clusterisation

@author RedSnail
October 2021
"""

import progressbar
import threading
import os
import sys
import copy
import time

from PanACoTA.pangenome_module.Clusterisator import Clusterisator, infinite_progressbar
from PanACoTA import utils

class MMseq(Clusterisator):
    def __init__(self, min_id, clust_mode, outdir, prt_path, threads, panfile, quiet):
        self.min_id = min_id
        self.clust_mode = clust_mode
        super(MMseq, self).__init__(threads, outdir, prt_path, panfile, quiet)

        prt_bank = os.path.basename(self.prt_path)

        if panfile is None:
            self.panfile_tmp = f"PanGenome-{prt_bank}-clust-{self.infoname}.lst"
        else:
            self.panfile_tmp = panfile

        self.panfile_tmp = os.path.join(outdir, self.panfile_tmp)

        self.mmseqdb = os.path.join(self.tmpdir, prt_bank + "-msDB")
        self.mmseqclust = os.path.join(self.tmpdir, prt_bank + "-clust-" + self.infoname)
        self.mmseqstsv = self.mmseqclust + ".tsv"

    @property
    def panfile(self):
        return self.panfile_tmp

    @property
    def method_name(self) -> str:
        return "mmseqs"

    @property
    def info_string(self) -> str:
        return ("Will run MMseqs2 with:\n"
                   f"\t- minimum sequence identity = {self.min_id*100}%\n"
                   f"\t- cluster mode {self.clust_mode}")

    @property
    def infoname(self) -> str:
        if self.threads != 1:
            threadinfo = "-th" + str(self.threads)
        else:
            threadinfo = ""
        infoname = str(self.min_id) + "-mode" + str(self.clust_mode) + threadinfo
        return infoname

    @property
    def expected_files(self):
        return list(map(lambda ext: self.mmseqdb + ext,
                        ["", ".index", ".dbtype", ".lookup", "_h", "_h.index", "_h.dbtype"]))

    @property
    def tmp_files_cmds(self):
        return [(f"mmseqs createdb {self.prt_path} {self.mmseqdb}",
                f"Problem while trying to convert database {self.prt_path} to mmseqs "
               "database format.")]

    @property
    def clustering_files(self):
        return list(map(lambda ext: self.mmseqclust + ext, [".0", ".1", ".dbtype", ".tsv"]))

    @property
    def clust_cmds(self):
        clust_cmd = (f"mmseqs cluster {self.mmseqdb} {self.mmseqclust} {self.tmpdir} "
                    f"--min-seq-id {self.min_id} --threads {self.threads} --cluster-mode "
                    f"{self.clust_mode}")
        tsv_cmd = f"mmseqs createtsv {self.mmseqdb} {self.mmseqdb} {self.mmseqclust} {self.mmseqstsv}"
        return [(clust_cmd, f"En error occured while clustering was done. See logs in {self.log_path}"),
                (tsv_cmd, f"En error occured while converting to tsv. See logs in {self.log_path}")]

    def parse_to_pangenome(self):
        clusters = {}  # {representative: [all members]}
        with open(self.mmseqstsv) as mmsf:
            for line in mmsf:
                repres, other = line.strip().split()
                if repres in clusters:
                    clusters[repres].append(other)
                else:
                    clusters[repres] = [repres]

        families = list(map(lambda mems: sorted(mems, key=utils.sort_proteins), clusters.values()))
        return dict(zip(range(1, len(families) + 1), families))


