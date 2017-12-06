#!/usr/bin/env python3

"""
Functional tests for genomeAPCAT pangenome
"""
import os

import shutil
import subprocess

import pytest

import genomeAPCAT.subcommands.align as al


# Define common variables
ALDIR = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALDIR, "exp_files")
TESTPATH = os.path.join(ALDIR, "test_files")
LOGFILE_BASE = "logs_test_postalign"


def setup_module():
    """
    create logger at start of this test module
    """
    utils.init_logger(LOGFILE_BASE, 0, '', verbose=1)
    print("Createc logger")


def teardown_module():
    """
    Remove log files at the end of this test module
    """
    os.remove(LOGFILE_BASE + ".log")
    os.remove(LOGFILE_BASE + ".log.details")
    os.remove(LOGFILE_BASE + ".log.err")
    print("Remove log files")


def test_main(caplog):
    """
    Test that when giving a database, a persistent genome and a list of genomes, it extracts
    expected proteins by family, aligns each family, back-translates them, concatenates all
    families into one file and groups them by genome.
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main"
    threads = 1
    force = False
    al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check creation of the 3 subdirectories
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    assert os.path.isdir(aldir)
    assert os.path.isdir(listdir)
    assert os.path.isdir(treedir)
    # Check content of listdir
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen)))
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen)))
    # Check content of aldir
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.gen').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.prt').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam))
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    exp_concat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    same_sequences(out_concat, exp_concat)
    # Check content of treedir
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    same_sequences(out_grp, exp_grp)
    # Check presence of log files, and log.err is empty
    base_log = os.path.join(outdir, "genomeAPCAT-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    with open(base_log + ".err", "r") as bf:
        assert bf.readlines() == []
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    assert "Extracting proteins and genes from all genomes" in caplog.text
    for gen in genomes:
        assert "Extracting proteins and genes from {}".format(gen) in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert "Aligning family {}".format(fam) in caplog.text
        assert "Back-translating family {}".format(fam) in caplog.text
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_main_exist_ok(caplog):
    """
    Test main all files exist and are ok, no force -> end without error, with warnings on re-use
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4exists"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_allexist"
    threads = 1
    force = False
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outgen = os.path.join(aldir, '{}-current.{}.gen').format(dname, fam)
        refgen = os.path.join(ex_aldir, "current.{}.gen").format(fam)
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(aldir, '{}-current.{}.prt').format(dname, fam)
        refprt = os.path.join(ex_aldir, "current.{}.prt").format(fam)
        shutil.copyfile(refprt, outprt)
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    outcat = os.path.join(aldir, dname + "-complete.cat.aln")
    refcat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    shutil.copyfile(refcat, outcat)
    # Create content of treedir
    outgrp = os.path.join(treedir, dname + ".grp.aln")
    refgrp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    shutil.copyfile(refgrp, outgrp)
    # Run align module
    al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    for gen in genomes:
        assert ("For genome {0}, "
                "test_align_main_allexist/List-TEST4exists/TEST4exists-getEntry_prt_{0}.txt and "
                "test_align_main_allexist/List-TEST4exists/TEST4exists-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen)) in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in caplog.text
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in caplog.text
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in caplog.text
    assert "Alignments already concatenated" in caplog.text
    assert ("Alignments already concatenated in "
            "test_align_main_allexist/Align-TEST4exists/TEST4exists-complete.cat.aln. Program will "
            "use it for next steps. If you want to redo it, remove it before "
            "running.")in caplog.text
    assert "Alignments already grouped by genome" in caplog.text
    assert ("Alignments already grouped by genome in "
            "test_align_main_allexist/Phylo-TEST4exists/TEST4exists.grp.aln. Program will "
            "end.") in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_main_exist_emptygrp(caplog):
    """
    test main all files exist but empty grp -> does nothing, grp is not checked if everything
    before was ok
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4empty-grp"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_allexist_emptygrp"
    threads = 1
    force = False
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outgen = os.path.join(aldir, '{}-current.{}.gen').format(dname, fam)
        refgen = os.path.join(ex_aldir, "current.{}.gen").format(fam)
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(aldir, '{}-current.{}.prt').format(dname, fam)
        refprt = os.path.join(ex_aldir, "current.{}.prt").format(fam)
        shutil.copyfile(refprt, outprt)
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    outcat = os.path.join(aldir, dname + "-complete.cat.aln")
    refcat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    shutil.copyfile(refcat, outcat)
    # Create content of treedir
    outgrp = os.path.join(treedir, dname + ".grp.aln")
    refgrp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    shutil.copyfile(refgrp, outgrp)
    # Run align module
    al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in caplog.text
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in caplog.text
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in caplog.text
    assert "Alignments already concatenated" in caplog.text
    assert ("Alignments already concatenated in {1}/Align-{0}/{0}-complete.cat.aln. Program will "
            "use it for next steps. If you want to redo it, remove it before "
            "running.".format(dname, outdir))in caplog.text
    assert "Alignments already grouped by genome" in caplog.text
    assert ("Alignments already grouped by genome in "
            "{1}/Phylo-{0}/{0}.grp.aln. Program will "
            "end.".format(dname, outdir)) in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_main_exist_emptycat(caplog):
    """
    test main all files exist but "1 sentence" concat -> do nothing, as we do not check concat file
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4empty-cat"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_allexist_emptycat"
    threads = 1
    force = False
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outgen = os.path.join(aldir, '{}-current.{}.gen').format(dname, fam)
        refgen = os.path.join(ex_aldir, "current.{}.gen").format(fam)
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(aldir, '{}-current.{}.prt').format(dname, fam)
        refprt = os.path.join(ex_aldir, "current.{}.prt").format(fam)
        shutil.copyfile(refprt, outprt)
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    outcat = os.path.join(aldir, dname + "-complete.cat.aln")
    with open(outcat, "w") as outf:
        outf.write("Hello !!")
    # Create content of treedir
    outgrp = os.path.join(treedir, dname + ".grp.aln")
    refgrp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    shutil.copyfile(refgrp, outgrp)
    # Run align module
    al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in caplog.text
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in caplog.text
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in caplog.text
    assert "Alignments already concatenated" in caplog.text
    assert ("Alignments already concatenated in {1}/Align-{0}/{0}-complete.cat.aln. Program will "
            "use it for next steps. If you want to redo it, remove it before "
            "running.".format(dname, outdir)) in caplog.text
    assert "Alignments already grouped by genome" in caplog.text
    assert ("Alignments already grouped by genome in "
            "{1}/Phylo-{0}/{0}.grp.aln. Program will "
            "end.".format(dname, outdir)) in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_main_exist_emptyp2n(caplog):
    """
    test main all files exist but empty prt2nuc for 1 fam -> redo this prt2nuc and concat and grp
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4empty-p2n"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_allexist_emptyp2n"
    threads = 1
    force = False
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outgen = os.path.join(aldir, '{}-current.{}.gen').format(dname, fam)
        refgen = os.path.join(ex_aldir, "current.{}.gen").format(fam)
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(aldir, '{}-current.{}.prt').format(dname, fam)
        refprt = os.path.join(ex_aldir, "current.{}.prt").format(fam)
        shutil.copyfile(refprt, outprt)
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    outbtr1 = os.path.join(aldir, '{}-mafft-prt2nuc.1.aln').format(dname)
    with open(outbtr1, "w") as outf:
        outf.write("Hello !!")
    outcat = os.path.join(aldir, dname + "-complete.cat.aln")
    refcat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    shutil.copyfile(refcat, outcat)
    # Create content of treedir
    outgrp = os.path.join(treedir, dname + ".grp.aln")
    refgrp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    shutil.copyfile(refgrp, outgrp)
    # Run align module
    al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check btr file for fam 1
    ref_btr1 = os.path.join(ex_aldir, "mafft-prt2nuc.1.aln")
    same_files(outbtr1, ref_btr1)
    # Check concat and group files
    same_files(refcat, outcat)
    same_sequences(outgrp, refgrp)
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in caplog.text
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in caplog.text
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    assert "Checking extractions for family 1" in caplog.text
    assert ("fam 1: Will redo back-translation, because found a different number of proteins "
            "aligned in {0}/Align-{1}/{1}-mafft-align.1.aln (4) and genes back-translated in "
            "existing {0}/Align-{1}/"
            "{1}-mafft-prt2nuc.1.aln".format(outdir, dname)) in caplog.text
    assert "Back-translating family 1" in caplog.text
    fams2 = fams[1:]
    for fam in fams2:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in caplog.text
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_main_emptymafft1(caplog):
    """
    Test main all files exist but mafft empty for 1 fam -> redo mafft and prt2nuc for this fam,
    + concat + grp
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-emptm1"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_emptym1"
    threads = 1
    force = False
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outgen = os.path.join(aldir, '{}-current.{}.gen').format(dname, fam)
        refgen = os.path.join(ex_aldir, "current.{}.gen").format(fam)
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(aldir, '{}-current.{}.prt').format(dname, fam)
        refprt = os.path.join(ex_aldir, "current.{}.prt").format(fam)
        shutil.copyfile(refprt, outprt)
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    outmafft1 = os.path.join(aldir, '{}-mafft-align.1.aln').format(dname)
    with open(outmafft1, "w") as outf:
        outf.write("Hello !!")
    al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check content of mafft and prt2nuc files for fam 1
    ref_mafft1 = os.path.join(ex_aldir, "mafft-align.1.aln")
    same_files(outmafft1, ref_mafft1)
    out_p2n1 = os.path.join(aldir, '{}-mafft-prt2nuc.1.aln').format(dname)
    ref_p2n1 = os.path.join(ex_aldir, "mafft-prt2nuc.1.aln")
    same_files(out_p2n1, ref_p2n1)
    # Check content of concatenated and grouped files
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    exp_concat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    same_sequences(out_concat, exp_concat)
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    same_sequences(out_grp, exp_grp)
    # Check presence of log files, and log.err is empty
    base_log = os.path.join(outdir, "genomeAPCAT-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in caplog.text
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in caplog.text
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    assert "Checking extractions for family 1" in caplog.text
    assert ("fam 1: Will redo alignment, because found a different number of proteins "
            "extracted in {0}/Align-{1}/{1}-current.1.prt (4) and proteins aligned in "
            "existing {0}/Align-{1}/{1}-mafft-align.1.aln".format(outdir, dname)) in caplog.text
    assert "Aligning family 1" in caplog.text
    assert "Back-translating family 1" in caplog.text
    fams2 = fams[1:]
    for fam in fams2:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in caplog.text
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text
    assert "END" in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_main_emptyprtgen(caplog):
    """
    Test main all files exist but empty prt extraction for 1 fam, empty gen extract for 1 other fam
    -> exits with error message because wrong extractions for family 1 and 4 (and not redone as
    not force). Should have removed align and prt files for fams 1 and 4, and concat and group files
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-empt-prtgen"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_empty-prtgen"
    threads = 1
    force = False
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outgen = os.path.join(aldir, '{}-current.{}.gen').format(dname, fam)
        refgen = os.path.join(ex_aldir, "current.{}.gen").format(fam)
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(aldir, '{}-current.{}.prt').format(dname, fam)
        refprt = os.path.join(ex_aldir, "current.{}.prt").format(fam)
        shutil.copyfile(refprt, outprt)
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    # Empty prt file for fam 1, empty gen file for fam 4
    outprt1 = os.path.join(aldir, '{}-current.1.prt').format(dname)
    outgen4 = os.path.join(aldir, '{}-current.4.gen').format(dname)
    open(outprt1, "w").close()
    open(outgen4, "w").close()
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    # Empty concat and group files
    open(out_concat, "w").close()
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    open(out_grp, "w").close()
    with pytest.raises(SystemExit):
        al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check that align and prt2nuc files are removed for fam 1 and 4
    for fam in [1, 4]:
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        assert not os.path.isfile(outaln)
        assert not os.path.isfile(outbtr)
    # Check concatenated and grouped files are removed
    assert not os.path.isfile(out_concat)
    assert not os.path.isfile(out_grp)
    # Check presence of log files
    base_log = os.path.join(outdir, "genomeAPCAT-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in caplog.text
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in caplog.text
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    assert "Checking extractions for family 1" in caplog.text
    assert ("fam 1: wrong sum of missing genomes (0) and prt extracted (0) for 4 genomes in "
            "the dataset.") in caplog.text
    assert "Checking extractions for family 4" in caplog.text
    assert ("fam 4: wrong sum of missing genomes (1) and gen extracted (0) for 4 genomes in "
            "the dataset.") in caplog.text
    fams2 = fams[2:]
    for fam in fams2:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in caplog.text
    assert ("At least one alignment did not run well. See detailed log file for more information. "
            "Program will stop here, alignments won't be grouped by genome.") in caplog.text
    shutil.rmtree(outdir)


def test_main_emptylist(caplog):
    """
    Test main all files exist but wrong list of proteins for 1 genome (and no files of extracted
    prt nor gen )-> exits with error, and deleted all gen, prt, align, concat and group files
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-wronglist"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_wronglist"
    threads = 1
    force = False
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    prefix = os.path.join(aldir, dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes[1:]:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        replace_path(outgen, prefix)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
        replace_path(outprt, prefix)
    # for genome 1, empty lists
    outgen = os.path.join(listdir, "{}-getEntry_gen_GEN2.1017.00001.txt".format(dname))
    outprt = os.path.join(listdir, "{}-getEntry_prt_GEN2.1017.00001.txt".format(dname))
    open(outgen, "w").close()
    open(outprt, "w").close()
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    # Empty concat and group files
    open(out_concat, "w").close()
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    open(out_grp, "w").close()
    with pytest.raises(SystemExit):
        al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # # Check that prt, gen, align and prt2nuc files are removed for all fams
    for fam in fams:
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        assert not os.path.isfile(outaln)
        assert not os.path.isfile(outbtr)
    # Check concatenated and grouped files are removed
    assert not os.path.isfile(out_concat)
    assert not os.path.isfile(out_grp)
    # Check presence of log files
    base_log = os.path.join(outdir, "genomeAPCAT-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in caplog.text
    assert "Extracting proteins and genes from all genomes" in caplog.text
    for gen in genomes:
        assert "Extracting proteins and genes from {}".format(gen) in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    assert "Checking extractions for family 4" in caplog.text
    assert ("fam 4: wrong sum of missing genomes (1) and prt extracted (2) for 4 genomes in "
            "the dataset.") in caplog.text
    assert "Checking extractions for family 8" in caplog.text
    assert ("fam 8: wrong sum of missing genomes (1) and prt extracted (2) for 4 genomes in "
            "the dataset.") in caplog.text
    for fam in [1, 6, 10, 11, 13, 14]:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert ("fam {}: wrong sum of missing genomes (0) and prt extracted (3) for 4 genomes in "
                "the dataset.".format(fam)) in caplog.text
    assert ("At least one alignment did not run well. See detailed log file for more information. "
            "Program will stop here, alignments won't be grouped by genome.") in caplog.text
    shutil.rmtree(outdir)


def test_main_force(caplog):
    """
    Test main all files exist but force -> redo everything
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-wronglist"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_wronglist"
    threads = 1
    force = True
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    prefix = os.path.join(aldir, dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes[1:]:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(ex_listdir, "getEntry_gen_{}".format(gen))
        shutil.copyfile(refgen, outgen)
        replace_path(outgen, prefix)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(ex_listdir, "getEntry_prt_{}".format(gen))
        shutil.copyfile(refprt, outprt)
        replace_path(outprt, prefix)
    # for genome 1, empty lists
    outgen = os.path.join(listdir, "{}-getEntry_gen_GEN2.1017.00001.txt".format(dname))
    outprt = os.path.join(listdir, "{}-getEntry_prt_GEN2.1017.00001.txt".format(dname))
    open(outgen, "w").close()
    open(outprt, "w").close()
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outmiss = os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam)
        refmiss = os.path.join(ex_aldir, "current.{}.miss.lst").format(fam)
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        refaln = os.path.join(ex_aldir, "mafft-align.{}.aln").format(fam)
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        refbtr = os.path.join(ex_aldir, "mafft-prt2nuc.{}.aln").format(fam)
        shutil.copyfile(refbtr, outbtr)
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    # Empty concat and group files
    open(out_concat, "w").close()
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    open(out_grp, "w").close()
    # Run align
    al.main(corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check content of listdir
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen)))
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen)))
    # Check content of aldir
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.gen').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.prt').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam))
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    exp_concat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    same_sequences(out_concat, exp_concat)
    # Check content of treedir
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    same_sequences(out_grp, exp_grp)
    # Check presence of log files, and log.err is empty
    base_log = os.path.join(outdir, "genomeAPCAT-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    assert "Extracting proteins and genes from all genomes" in caplog.text
    for gen in genomes:
        assert "Extracting proteins and genes from {}".format(gen) in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert "Aligning family {}".format(fam) in caplog.text
        assert "Back-translating family {}".format(fam) in caplog.text
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_main_from_parse(caplog):
    """
    Test main when we give the output of the parser
    """
    import argparse
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-fromparse"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_fromparse"
    threads = 1
    force = True
    args = argparse.Namespace()
    args.corepers = corepers
    args.list_genomes = list_genomes
    args.dataset_name = dname
    args.dbpath = dbpath
    args.outdir = outdir
    args.threads = threads
    args.force = force
    args.verbose = 1
    args.quiet = False
    al.main_from_parse(args)
    # Check creation of the 3 subdirectories
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    assert os.path.isdir(aldir)
    assert os.path.isdir(listdir)
    assert os.path.isdir(treedir)
    # Check content of listdir
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen)))
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen)))
    # Check content of aldir
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.gen').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.prt').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam))
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    exp_concat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    same_sequences(out_concat, exp_concat)
    # Check content of treedir
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    same_sequences(out_grp, exp_grp)
    # Check presence of log files, and log.err is empty
    base_log = os.path.join(outdir, "genomeAPCAT-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    with open(base_log + ".err", "r") as bf:
        assert bf.readlines() == []
    # Check logs
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in caplog.text
    assert "Extracting proteins and genes from all genomes" in caplog.text
    for gen in genomes:
        assert "Extracting proteins and genes from {}".format(gen) in caplog.text
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in caplog.text
        assert "Aligning family {}".format(fam) in caplog.text
        assert "Back-translating family {}".format(fam) in caplog.text
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text
    assert "END" in caplog.text
    shutil.rmtree(outdir)


def test_align_all():
    """
    Test when calling align from command line, it runs and gives expected output files
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-fromparse"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = "test_align_main_fromparse"

    cmd = "genomeAPCAT align -c {} -l {} -n {} -d {} -o {}".format(corepers, list_genomes, dname,
                                                                   dbpath, outdir)
    ret = subprocess.call(cmd.split())
    assert ret == 0
    # Check creation of the 3 subdirectories
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    assert os.path.isdir(aldir)
    assert os.path.isdir(listdir)
    assert os.path.isdir(treedir)
    # Check content of listdir
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen)))
        assert os.path.isfile(os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen)))
    # Check content of aldir
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.gen').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.prt').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-current.{}.miss.lst').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam))
        assert os.path.isfile(os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam))
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    exp_concat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    same_sequences(out_concat, exp_concat)
    # Check content of treedir
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    same_sequences(out_grp, exp_grp)
    # Check presence of log files, and log.err is empty
    base_log = os.path.join(outdir, "genomeAPCAT-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    with open(base_log + ".err", "r") as bf:
        assert bf.readlines() == []
    shutil.rmtree(outdir)


def replace_path(filein, prefix):
    """
    From a file containing list of genes/proteins, and path to the file from which it must be
    extracted, replace the current path by the given prefix

    Parameters
    ----------
    filein : str
        path to file containing gen/prt list
    prefix : str
        what should be written instead of current path
    """
    shutil.copyfile(filein, filein + ".save")
    with open(filein + ".save", "r") as inf, open(filein, "w") as outf:
        for line in inf:
            outf.write(line.replace("test_align_main/Align-TEST4/TEST4", prefix))
    os.remove(filein + ".save")


def same_sequences(file_out, file_exp):
    """
    Check that the 2 files have the same content.

    Parameters
    ----------
    file_out : str
        file generated by the test
    file_exp : str
        file containing what should be generated
    """
    seq_out = get_seqs(file_out)
    seq_exp = get_seqs(file_exp)
    assert seq_out.keys() == seq_exp.keys()
    for name, seq in seq_out.items():
        assert len(seq) == len(seq_exp[name])
        assert set(seq) == set(seq_exp[name])
    assert seq_out == seq_exp


def get_seqs(filein):
    """
    Return dict with headers as keys and corresponding sequence as value

    Parameters
    ----------
    filein : multi fasta file to read

    Returns
    -------
    dict
        {header: "sequence"}
    """
    seq_out = {}
    with open(filein, "r") as fo:
        cur_seq = None
        for line in fo:
            if line.startswith(">"):
                cur_seq = line.strip()
                seq_out[cur_seq] = ""
            else:
                seq_out[cur_seq] += line.strip()
    return seq_out


def same_files(file_out, file_exp):
    """
    Check that the 2 files have the same content.

    Parameters
    ----------
    file_out : str
        file generated by the test
    file_exp : str
        file containing what should be generated
    """
    with open(file_out, "r") as fo, open(file_exp, "r") as fe:
        lines_out = fo.readlines()
        lines_exp = fe.readlines()
        assert len(lines_exp) == len(lines_out)
        for linout, linexp in zip(lines_out, lines_exp):
            assert linout == linexp
