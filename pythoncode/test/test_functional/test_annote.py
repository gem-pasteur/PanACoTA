#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for annote_pipeline.py
"""

import pytest
import os
import subprocess
import shutil
import time


def test_annote_allDefault():
    """
    Test that when we call the pipeline with all default parameters, all expected output files
    are created. Check the content of result files (LSTINFO, Genes, Proteins, Replicons),
    lstinfo, discarded and log files.
    """
    date = time.strftime("%m%y")
    fulldate = time.strftime("%Y-%m-%d")
    list_file = os.path.join("test", "data", "test_files", "list_genomes-func-test-default.txt")
    dbpath = os.path.join("test", "data", "genomes")
    respath = os.path.join("test", "data", "res_test_funcDefault")
    name = "GENO"
    cmd = "annote_pipeline.py {} -d {} -r {} -n {}".format(list_file, dbpath, respath, name)
    ret = subprocess.call(cmd.split())
    assert ret == 0
    # Get output files
    log_files = [os.path.join(respath, "annote-genomes-list_genomes-func-test-default" + ext)
                 for ext in [".log", ".log.err"]]
    lstfile = [os.path.join(respath, "LSTINFO-list_genomes-func-test-default.lst")]
    discfile = [os.path.join(respath, "discarded-list_genomes-func-test-default.lst")]
    QCfiles = [os.path.join(respath, qc + "list_genomes-func-test-default.png")
               for qc in ["QC_L90-", "QC_nb-contigs-"]]
    genomes = {"B2_A3_5.fasta-changeName.fna": "GENO.{}.00003".format(date),
               "H299_H561.fasta": "GENO.{}.00002".format(date),
               "A_H738.fasta": "GENO.{}.00001".format(date)}
    split5Nfiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna") for g in genomes]
    gembasefiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna")
                    for g in genomes]
    proklogfiles = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna-prokka.log")
                    for g in genomes]
    prokka_files = [os.path.join(respath, "tmp_files", g + "-split5N.fna-gembase.fna-prokkaRes",
                                 genomes[g])
                    for g in genomes]
    lstinffiles = [os.path.join(respath, "LSTINFO", name + ".lst")
                   for name in list(genomes.values())]
    prtfiles = [os.path.join(respath, "Proteins", name + ".prt")
                for name in list(genomes.values())]
    genfiles = [os.path.join(respath, "Genes", name + ".gen")
                for name in list(genomes.values())]
    repfiles = [os.path.join(respath, "Replicons", name + ".fna")
                for name in list(genomes.values())]
    # Check all output files exist
    for f in (log_files + lstfile + discfile + QCfiles + split5Nfiles + gembasefiles +
              proklogfiles + lstinffiles + prtfiles + genfiles + repfiles):
        assert os.path.isfile(f)
    for f in prokka_files:
        assert os.path.isfile(f + ".tbl")
        assert os.path.isfile(f + ".faa")
        assert os.path.isfile(f + ".ffn")
    # Check content of result database files
    exp_dir = os.path.join("test", "data", "exp_files", "results_test_func-default")
    for f in (lstinffiles + prtfiles + genfiles + repfiles):
        exp_file = os.sep.join(f.split(os.sep)[-2:])
        exp_path = os.path.join(exp_dir, exp_file)
        with open(exp_path, "r") as expf, open(f, "r") as outf:
            for line_exp, line_out in zip(expf, outf):
                assert line_exp == line_out
    # Check that err file is empty
    with open(log_files[1], "r") as errf:
        lines = errf.readlines()
        assert lines == []
    # Check that log file contains expected information (info level)
    with open(log_files[0], "r") as logf:
        infos = []
        for line in logf:
            assert line.startswith("[" + fulldate)
            assert line.count("::") == 2
            infos.append(line.strip().split("::")[-1].strip())
        assert "Reading genomes" in infos
        assert ("Cutting genomes at each stretch of at least 5 'N', and then, calculating "
                "genome size, number of contigs and L90.") in infos
        assert ("Annotating all genomes with prokka") in infos
        assert ("Formatting all genomes") in infos
    # Check lstinfo file content
    explst = os.path.join(exp_dir, "LSTINFO-list_genomes-func-test-default.lst")
    with open(lstfile[0], "r") as expf, open(explst, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    shutil.rmtree(respath)