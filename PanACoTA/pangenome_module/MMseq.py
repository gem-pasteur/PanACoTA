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

import progressbar
import threading
import os
import sys
import copy
import time

from PanACoTA.pangenome_module.Clusterisator import Clusterisator
from PanACoTA import utils

class MMseq(Clusterisator):
    def __init__(self, min_id, clust_mode, outdir, prt_path, threads, panfile, quiet):
        self.min_id = min_id
        self.clust_mode = clust_mode
        super(MMseq, self).__init__(threads, outdir, prt_path, panfile, quiet)

        self.mmseqdb = os.path.join(self.tmpdir, self.prt_bank + "-msDB")
        self.mmseqclust = os.path.join(self.tmpdir, self.prt_bank + "-clust-" + self.infoname)
        self.mmseqstsv = self.mmseqclust + ".tsv"

    @property
    def method_name(self) -> str:
        return "mmseq"

    def info_string(self) -> str:
        return ("Will run MMseqs2 with:\n"
                   f"\t- minimum sequence identity = {self.min_id*100}%\n"
                   f"\t- cluster mode {self.clust_mode}")

    def get_info(self) -> str:
        if self.threads != 1:
            threadinfo = "-th" + str(self.threads)
        else:
            threadinfo = ""
        infoname = str(self.min_id) + "-mode" + str(self.clust_mode) + threadinfo
        return infoname

    def arrange_files(self, tmp_dir) -> int:
        self.logger.info("Creating database")
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
            res = self.create_mmseqs_db()
        # except KeyboardInterrupt: # pragma: no cover
        except:  # pragma: no cover
            stop_bar = True
            x.join()
            sys.exit(1)
        # Clustering done, stop bar and join (if quiet, it was already finished, so we just join it)
        stop_bar = True
        x.join()
        return res

    def create_mmseqs_db(self):
        outext = ["", ".index", ".dbtype", ".lookup", "_h", "_h.index", "_h.dbtype"]
        files_existing = []
        if os.path.isfile(self.mmseqdb):
            for file in [self.mmseqdb + ext for ext in outext]:
                if not os.path.isfile(file):
                    continue
                files_existing.append(file)
            if len(files_existing) != len(outext):
                self.logger.warning(f"mmseqs database {self.mmseqdb} already exists, but at least 1 associated "
                               "file (.dbtype, .index etc). is missing. The program will "
                               "remove existing files and recreate the database.")
                files_remaining = copy.deepcopy(files_existing)
                for file in files_existing:
                    os.remove(file)  # Delete file
                    files_remaining.remove(file)  # Remove file from list of existing files
                    logger.details(f"Removing '{file}'.")
                files_existing = copy.deepcopy(files_remaining)
            else:
                self.logger.warning(f"mmseqs database {mself.mseqdb} already exists. The program will "
                               "use it.")
                return False
        self.logger.debug("Existing files: {}".format(len(files_existing)))
        self.logger.debug("Expected extensions: {}".format(len(outext)))
        cmd = f"mmseqs createdb {self.prt_path} {self.mmseqdb}"
        msg = (f"Problem while trying to convert database {self.prt_path} to mmseqs "
               "database format.")
        self.logger.details(f"MMseqs command: {cmd}")
        with open(self.log_path, "w") as logf:
            utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
        return True

    def do_pangenome(self, status) -> tuple:
        # If we just made the database, we must redo all next steps
        # -> if existing, remove
        # mmseqsclust (created by run_mmseqs_clust)
        # mmseqstsv (created by mmseqs_to_pangenome)
        # pangenome file
        if status and os.path.isfile(self.mmseqclust) or os.path.isfile(self.mmseqstsv) or os.path.isfile(self.panfile):
            self.logger.details("Removing existing clustering and/or pangenome files.")
            utils.remove(self.mmseqclust)
            utils.remove(self.mmseqstsv)
            utils.remove(self.panfile)
        bar = None
        self.logger.debug(self.mmseqclust)
        if os.path.isfile(self.mmseqclust):
            self.logger.warning((f"mmseqs clustering {self.mmseqclust} already exists. The program will now convert "
                            "it to a pangenome file."))
        else:
            self.logger.info("Clustering proteins...")
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
                self.run_mmseqs_clust()
            # except KeyboardInterrupt: # pragma: no cover
            except:  # pragma: no cover
                stop_bar = True
                x.join()
                sys.exit(1)
            # Clustering done, stop bar and join (if quiet, it was already finished, so we just join it)
            stop_bar = True
            x.join()
        # Convert output to tsv file (one line per comparison done)
        #  # Convert output to tsv file (one line per comparison done)
        # -> returns (families, outfile)
        families = self.mmseqs_to_pangenome()
        return families

    def run_mmseqs_clust(self):
        cmd = (
            f"mmseqs cluster {self.mmseqdb} {self.mmseqclust} {self.tmpdir} " 
            f"--min-seq-id {self.min_id} --threads {self.threads} --cluster-mode "
            f"{self.clust_mode}")
        self.logger.details(f"MMseqs command: {cmd}")
        msg = f"Problem while clustering proteins with mmseqs. See log in {self.log_path}"
        with open(self.log_path, "a") as logm:
            utils.run_cmd(cmd, msg, eof=False, stdout=logm, stderr=logm)

    def mmseqs_to_pangenome(self):
        cmd = f"mmseqs createtsv {self.mmseqdb} {self.mmseqdb} {self.mmseqclust} {self.mmseqstsv}"
        msg = "Problem while trying to convert mmseq result file to tsv file"
        self.logger.details(f"MMseqs command: {cmd}")
        with open(self.log_path, "a") as logf:
            utils.run_cmd(cmd, msg, eof=True, stdout=logf, stderr=logf)
        # Convert the tsv file to a 'pangenome' file: one line per family
        families = self.mmseqs_tsv_to_pangenome()
        return families

    def mmseqs_tsv_to_pangenome(self):
        self.logger.info("Converting mmseqs results to pangenome file")
        clusters = self.mmseq_tsv_to_clusters()
        families = self.clusters_to_file(clusters)
        end = time.strftime('%Y-%m-%d_%H-%M-%S')
        with open(self.log_path, "a") as logm:
            logm.write(f"End: {end}")
        return families

    def mmseq_tsv_to_clusters(self):
        clusters = {}  # {representative: [all members]}
        with open(self.mmseqstsv) as mmsf:
            for line in mmsf:
                repres, other = line.strip().split()
                if repres in clusters:
                    clusters[repres].append(other)
                else:
                    clusters[repres] = [repres]
        self.logger.info(f"Pangenome has {len(clusters)} families.")
        return clusters

    def clusters_to_file(self, clust):
        families = {}  # {famnum: [members]}
        with open(self.panfile, "w") as fout:
            num = 1
            for _, fam in clust.items():
                families[num] = []
                fout.write(str(num))
                for mem in sorted(fam, key=utils.sort_proteins):
                    families[num].append(mem)
                    fout.write(" " + mem)
                fout.write("\n")
                num += 1
        return families

