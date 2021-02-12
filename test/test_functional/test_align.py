#!/usr/bin/env python3

"""
Functional tests for 'align' subcommand
"""
import os
import shutil
import subprocess
import pytest
import logging

import PanACoTA.subcommands.align as al
import test.test_unit.utilities_for_tests as tutil
# from PanACoTA import utils


# Define common variables
ALDIR = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALDIR, "exp_files")
TESTPATH = os.path.join(ALDIR, "test_files")
LOGFILE_BASE = "logs_test_postalign"
GENEPATH = os.path.join(ALDIR, "generated_by_func_tests")

@pytest.fixture(autouse=True)
def setup_teardown_module():
    """
    Remove log files at the end of this test module

    Before each test:
    - init logger
    - create directory to put generated files

    After:
    - remove all log files
    - remove directory with generated results
    """
    if not os.path.isdir(GENEPATH):
        os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")


def test_main():
    """
    Test that when giving a database, a persistent genome and a list of genomes, it extracts
    expected proteins by family, aligns each family, back-translates them, concatenates all
    families into one file and groups them by genome.
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = GENEPATH
    threads = 1
    force = False
    cmd = "cmd"
    al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force)
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
        assert os.path.isfile(os.path.join(listdir, f"{dname}-getEntry_gen_{gen}.txt"))
        assert os.path.isfile(os.path.join(listdir, f"{dname}-getEntry_prt_{gen}.txt"))
    # Check content of aldir
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        assert os.path.isfile(os.path.join(aldir, f'{dname}-current.{fam}.gen'))
        assert os.path.isfile(os.path.join(aldir, f'{dname}-current.{fam}.prt'))
        assert os.path.isfile(os.path.join(aldir, f'{dname}-current.{fam}.miss.lst'))
        assert os.path.isfile(os.path.join(aldir, f'{dname}-mafft-align.{fam}.aln'))
        assert os.path.isfile(os.path.join(aldir, f'{dname}-mafft-prt2nuc.{fam}.aln'))
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    exp_concat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    assert tutil.compare_order_content(out_concat, exp_concat)
    # Check content of treedir
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    assert tutil.compare_order_content(out_grp, exp_grp)
    # Check presence of log files, and log.err is empty
    base_log = os.path.join(outdir, "PanACoTA-align_" + dname + ".log")
    assert os.path.isfile(base_log)
    assert os.path.isfile(base_log + ".details")
    assert os.path.isfile(base_log + ".err")
    with open(base_log + ".err", "r") as bf:
        assert bf.readlines() == []
    # Check logs
    with open(base_log + ".details", "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in " ".join(log_content)
    assert "Extracting proteins and genes from all genomes" in " ".join(log_content)
    for gen in genomes:
        assert "Extracting proteins and genes from {}".format(gen) in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in " ".join(log_content)
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in " ".join(log_content)
        assert "Aligning family {}".format(fam) in " ".join(log_content)
        assert "Back-translating family {}".format(fam) in " ".join(log_content)
    assert "Concatenating all alignment files" in " ".join(log_content)
    assert "Grouping alignments per genome" in " ".join(log_content)
    assert "END" in " ".join(log_content)


def test_main_exist_ok():
    """
    Test main all files exist and are ok, no force -> end without error, with warnings on re-use
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4exists"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_exist_ok")
    threads = 1
    force = False
    cmd = "cmd test_main_exist_ok"
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
        outgen = os.path.join(listdir, f"{dname}-getEntry_gen_{gen}.txt")
        refgen = os.path.join(ex_listdir, f"getEntry_gen_{gen}")
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, f"{dname}-getEntry_prt_{gen}.txt")
        refprt = os.path.join(ex_listdir, f"getEntry_prt_{gen}")
        shutil.copyfile(refprt, outprt)
    # Create content of aldir
    ex_aldir = os.path.join(EXPPATH, "exp_aldir-pers")
    fams = [1, 4, 6, 8, 10, 11, 13, 14]
    for fam in fams:
        outgen = os.path.join(aldir, f'{dname}-current.{fam}.gen')
        refgen = os.path.join(ex_aldir, f"current.{fam}.gen")
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(aldir, f'{dname}-current.{fam}.prt')
        refprt = os.path.join(ex_aldir, f"current.{fam}.prt")
        shutil.copyfile(refprt, outprt)
        outmiss = os.path.join(aldir, f'{dname}-current.{fam}.miss.lst')
        refmiss = os.path.join(ex_aldir, f"current.{fam}.miss.lst")
        shutil.copyfile(refmiss, outmiss)
        outaln = os.path.join(aldir, f'{dname}-mafft-align.{fam}.aln')
        refaln = os.path.join(ex_aldir, f"mafft-align.{fam}.aln")
        shutil.copyfile(refaln, outaln)
        outbtr = os.path.join(aldir, f'{dname}-mafft-prt2nuc.{fam}.aln')
        refbtr = os.path.join(ex_aldir, f"mafft-prt2nuc.{fam}.aln")
        shutil.copyfile(refbtr, outbtr)
    outcat = os.path.join(aldir, dname + "-complete.cat.aln")
    refcat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    shutil.copyfile(refcat, outcat)
    # Create content of treedir
    outgrp = os.path.join(treedir, dname + ".grp.aln")
    refgrp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    shutil.copyfile(refgrp, outgrp)
    # Run align module
    al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check logs
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4exists.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in " ".join(log_content)
    for gen in genomes:
        assert (f"For genome {gen}, test/data/align/generated_by_func_tests/test_main_exist_ok/"
                f"List-TEST4exists/TEST4exists-getEntry_prt_{gen}.txt and test/data/align/"
                "generated_by_func_tests/test_main_exist_ok/List-TEST4exists/"
                f"TEST4exists-getEntry_gen_{gen}.txt already exist. The program "
                "will use them to extract proteins and genes. If you prefer to rewrite "
                "them, use option -F ") in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in " ".join(log_content)
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in " ".join(log_content)
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in " ".join(log_content)
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in " ".join(log_content)
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in " ".join(log_content)
    assert ("Alignments already concatenated in "
            "test/data/align/generated_by_func_tests/test_main_exist_ok/Align-TEST4exists/"
            "TEST4exists-complete.cat.aln. Program will "
            "use it for next steps. If you want to redo it, remove it before "
            "running.") in " ".join(log_content)
    assert ("Alignments already grouped by genome in "
            "test/data/align/generated_by_func_tests/test_main_exist_ok/Phylo-TEST4exists/"
            "TEST4exists.grp.aln. Program will "
            "end.") in " ".join(log_content)
    assert "END" in " ".join(log_content)


def test_main_exist_emptygrp(capsys):
    """
    test main all files exist but empty grp -> does nothing, grp is not checked if everything
    before was ok. grp must still be empty
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4empty-grp"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_exist_ok")
    threads = 1
    force = False
    cmd = "cmd test_main_exist_emptygrp"
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
    # Create content of treedir, with empty grp
    outgrp = os.path.join(treedir, dname + ".grp.aln")
    open(outgrp, "w").close()
    # Run align module
    al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force, verbose=2)
    # Check grp still empty
    assert os.stat(outgrp).st_size == 0
    # Check logs
    out, err = capsys.readouterr()
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4empty-grp.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in out
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in out
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in out
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in " ".join(log_content)
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in " ".join(log_content)
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in " ".join(log_content)
    assert "Alignments already concatenated" in out
    assert ("Alignments already concatenated in {1}/Align-{0}/{0}-complete.cat.aln. Program will "
            "use it for next steps. If you want to redo it, remove it before "
            "running.").format(dname, outdir) in " ".join(log_content)
    assert "Alignments already grouped by genome" in " ".join(log_content)
    assert ("Alignments already grouped by genome in "
            "test/data/align/generated_by_func_tests/test_main_exist_ok/"
            "Phylo-{0}/{0}.grp.aln. Program will end.").format(dname) in " ".join(log_content)
    assert "END" in out


def test_main_exist_emptycat(capsys):
    """
    test main all files exist but "1 sentence" concat -> do nothing, as we do not check concat file
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4empty-cat"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_exist_emptycat")
    threads = 1
    force = False
    cmd = "cmd test_main_exist_emptycat"
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
    with open(outgrp, "w") as outf:
        outf.write("It's me !!")
    # Run align module
    al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force, verbose=16)
    # Check concat and grp did not change
    with open(outcat, "r") as of:
        assert of.readlines() == ["Hello !!"]
    with open(outgrp, "r") as of:
        assert of.readlines() == ["It's me !!"]
    # Check logs
    out, err = capsys.readouterr()
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4empty-cat.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in out
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in out
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in out
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in " ".join(log_content)
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in " ".join(log_content)
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in " ".join(log_content)
    assert "Alignments already concatenated" in out
    assert ("Alignments already concatenated in {1}/Align-{0}/{0}-complete.cat.aln. Program will "
            "use it for next steps. If you want to redo it, remove it before "
            "running.".format(dname, outdir)) in " ".join(log_content)
    assert "Alignments already grouped by genome" in out
    assert ("Alignments already grouped by genome in "
            "{1}/Phylo-{0}/{0}.grp.aln. Program will "
            "end.".format(dname, outdir)) in " ".join(log_content)
    assert "END" in out


def test_main_exist_emptyp2n(capsys):
    """
    test main all files exist but empty prt2nuc for 1 fam -> redo this prt2nuc and concat and grp
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4empty-p2n"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_exist_emptyp2n")
    threads = 1
    force = False
    cmd = "cmd test_main_exist_emptycat"
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
    with open(outcat, "w") as oc:
        oc.write("complete cat aln")
    # Create content of treedir
    outgrp = os.path.join(treedir, dname + ".grp.aln")
    with open(outgrp, "w") as og:
        og.write("grouped grp aln")
    # Run align module
    al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check btr file for fam 1
    ref_btr1 = os.path.join(ex_aldir, "mafft-prt2nuc.1.aln")
    assert tutil.compare_order_content(ref_btr1, outbtr1)
    # Check concat and group files
    refcat = os.path.join(ex_aldir, "complete.cat.aln")
    assert tutil.compare_order_content(refcat, outcat)
    refgrp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    assert tutil.compare_order_content(refgrp, outgrp)
    # Check logs
    out, err = capsys.readouterr()
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4empty-p2n.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in out
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in " ".join(log_content)
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in out
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in out
    assert "Checking extractions for family 1" in " ".join(log_content)
    assert ("fam 1: Will redo back-translation, because found a different number of proteins "
            "aligned in {0}/Align-{1}/{1}-mafft-align.1.aln (4) and genes back-translated in "
            "existing {0}/Align-{1}/"
            "{1}-mafft-prt2nuc.1.aln".format(outdir, dname)) in " ".join(log_content)
    assert "Back-translating family 1" in " ".join(log_content)
    fams2 = fams[1:]
    for fam in fams2:
        assert "Checking extractions for family {}".format(fam) in " ".join(log_content)
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in " ".join(log_content)
    assert "Concatenating all alignment files" in out
    assert "Grouping alignments per genome" in out
    assert "END" in out


def test_main_emptymafft1(capsys):
    """
    Test main all files exist but mafft empty for 1 fam -> redo mafft and prt2nuc for this fam,
    + concat + grp
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-emptym1"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_exist_emptym1")
    threads = 1
    force = False
    cmd = "cmd test_main_exist_empty-m1"
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
    al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check content of mafft and prt2nuc files for fam 1
    ref_mafft1 = os.path.join(ex_aldir, "mafft-align.1.aln")
    assert tutil.compare_order_content(outmafft1, ref_mafft1)
    out_p2n1 = os.path.join(aldir, '{}-mafft-prt2nuc.1.aln').format(dname)
    ref_p2n1 = os.path.join(ex_aldir, "mafft-prt2nuc.1.aln")
    assert tutil.compare_order_content(out_p2n1, ref_p2n1)
    # Check content of concatenated and grouped files
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    exp_concat = os.path.join(EXPPATH, "exp_pers4genome-complete.cat.aln")
    assert tutil.compare_order_content(out_concat, exp_concat)
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    assert tutil.compare_order_content(out_grp, exp_grp)
    # Check presence of log files, and log.err is empty
    out, err = capsys.readouterr()
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4-emptym1.log")
    assert os.path.isfile(logfile)
    assert os.path.isfile(logfile + ".details")
    assert os.path.isfile(logfile + ".err")
    # Check logs
    with open(logfile + ".details", "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in out
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in " ".join(log_content)
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in out
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in out
    assert "Checking extractions for family 1" in " ".join(log_content)
    assert ("fam 1: Will redo alignment, because found a different number of proteins "
            "extracted in {0}/Align-{1}/{1}-current.1.prt (4) and proteins aligned in "
            "existing {0}/Align-{1}/"
            "{1}-mafft-align.1.aln".format(outdir, dname)) in " ".join(log_content)
    assert "Aligning family 1" in " ".join(log_content)
    assert "Back-translating family 1" in " ".join(log_content)
    fams2 = fams[1:]
    for fam in fams2:
        assert "Checking extractions for family {}".format(fam) in " ".join(log_content)
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in " ".join(log_content)
    assert "Concatenating all alignment files" in out
    assert "Grouping alignments per genome" in out
    assert "END" in out


def test_main_emptyprtgen(capsys):
    """
    Test main all files exist but empty prt extraction for 1 fam, empty gen extract for 1 other fam
    -> exits with error message because wrong extractions for family 1 and 4 (and not redone as
    not force). Should have removed align and prt2nuc files for fams 1 and 4, and concat and group files
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-empty-prtgen"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_exist_empty-prtgen")
    threads = 1
    force = False
    cmd = "cmd test_main_exist_empty-m1"
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
    # Empty concat and group files
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    open(out_concat, "w").close()
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    open(out_grp, "w").close()
    # Run alignment, but problem
    with pytest.raises(SystemExit):
        al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check that align and prt2nuc files are removed for fam 1 and 4
    for fam in [1, 4]:
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        assert not os.path.isfile(outaln)
        assert not os.path.isfile(outbtr)
    # Check concatenated and grouped files are removed
    assert not os.path.isfile(out_concat)
    assert not os.path.isfile(out_grp)
    # Check logs
    out, err = capsys.readouterr()
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4-empty-prtgen.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in out
    for gen in genomes:
        assert ("For genome {0}, "
                "{2}/List-{1}/{1}-getEntry_prt_{0}.txt and {2}/List-{1}/{1}-getEntry_gen_{0}.txt "
                "already exist. The program will use them to extract proteins and genes. If you "
                "prefer to rewrite them, use option -F "
                "(or --force)".format(gen, dname, outdir)) in " ".join(log_content)
    assert ("All extraction files already existing (see detailed log for more "
            "information)") in out
    assert ("All prt and gene files for all families already exist. The program will use them "
            "for the next step. If you want to re-extract a given family, remove its prt and "
            "gen extraction files. If you want to re-extract all families, use option -F "
            "(or --force).") in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in out
    assert "Checking extractions for family 1" in " ".join(log_content)
    assert ("fam 1: wrong sum of missing genomes (0) and prt extracted (0) for 4 genomes in "
            "the dataset.") in " ".join(log_content)
    assert "Checking extractions for family 4" in " ".join(log_content)
    assert ("fam 4: wrong sum of missing genomes (1) and gen extracted (0) for 4 genomes in "
            "the dataset.") in " ".join(log_content)
    fams2 = fams[2:]
    for fam in fams2:
        assert "Checking extractions for family {}".format(fam) in " ".join(log_content)
        assert ("Alignment already done for family {}. The program will use it for next "
                "steps").format(fam) in " ".join(log_content)
    assert ("At least one alignment did not run well. See detailed log file for more information. "
            "Program will stop here, alignments won't be grouped by genome.") in err


def test_main_emptylist(capsys):
    """
    Test main all files exist but wrong list of proteins for 1 genome (and no files of extracted
    prt nor gen )-> exits with error, and deleted all gen, prt, align, concat and group files
    """
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-wronglist"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_exist_wronglist")
    threads = 1
    force = False
    cmd = "cmd test_main_exist_empty-m1"
    # Create output directories and files
    aldir = os.path.join(outdir, "Align-" + dname)
    listdir = os.path.join(outdir, "List-" + dname)
    treedir = os.path.join(outdir, "Phylo-" + dname)
    os.makedirs(aldir)
    os.makedirs(listdir)
    os.makedirs(treedir)
    # Create content of listdir
    ex_listdir = os.path.join(EXPPATH, "exp_listdir-pers")
    genomes = ["GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        outgen = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        refgen = os.path.join(EXPPATH, f"exp_wronglist-getEntry_gen_{gen}.txt")
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(EXPPATH, f"exp_wronglist-getEntry_prt_{gen}.txt")
        shutil.copyfile(refprt, outprt)
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
        al.main(cmd, corepers, list_genomes, dname, dbpath, outdir, threads, force)
    # Check that prt, gen, align and prt2nuc files are removed for all fams
    for fam in fams:
        outaln = os.path.join(aldir, '{}-mafft-align.{}.aln').format(dname, fam)
        outbtr = os.path.join(aldir, '{}-mafft-prt2nuc.{}.aln').format(dname, fam)
        assert not os.path.isfile(outaln)
        assert not os.path.isfile(outbtr)
    # Check logs
    out, err = capsys.readouterr()
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4-wronglist.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in out
    for gen in genomes:
        assert (f"For genome {gen}, test/data/align/generated_by_func_tests/"
                "test_main_exist_wronglist/List-TEST4-wronglist/"
                f"TEST4-wronglist-getEntry_prt_{gen}.txt and test/data/align/"
                "generated_by_func_tests/test_main_exist_wronglist/List-TEST4-wronglist/"
                f"TEST4-wronglist-getEntry_gen_{gen}.txt already exist. "
                "The program will use them to extract proteins and genes. If you prefer to "
                "rewrite them, use option -F (or --force)") in " ".join(log_content)
    assert "Extracting proteins and genes from all genomes" in out
    for gen in genomes:
        assert "Extracting proteins and genes from {}".format(gen) in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in out
    assert "Checking extractions for family 4" in " ".join(log_content)
    assert ("fam 4: wrong sum of missing genomes (1) and prt extracted (2) for 4 genomes in "
            "the dataset.") in " ".join(log_content)
    assert "Checking extractions for family 8" in " ".join(log_content)
    assert ("fam 8: wrong sum of missing genomes (1) and prt extracted (2) for 4 genomes in "
            "the dataset.") in " ".join(log_content)
    for fam in [1, 6, 10, 11, 13, 14]:
        assert "Checking extractions for family {}".format(fam) in " ".join(log_content)
        assert ("fam {}: wrong sum of missing genomes (0) and prt extracted (3) for 4 genomes in "
                "the dataset.".format(fam)) in " ".join(log_content)
    assert ("At least one alignment did not run well. See detailed log file for more information. "
            "Program will stop here, alignments won't be grouped by genome.") in err


def test_main_from_parse_force(capsys):
    """
    Test main all files exist but force -> redo everything
    """
    import argparse
    corepers = os.path.join(TESTPATH, "test_pers0.99FX.lst")
    list_genomes = os.path.join("test", "data", "pangenome", "test_files", "list_to_pan.txt")
    dname = "TEST4-fromparse"
    dbpath = os.path.join("test", "data", "pangenome", "test_files", "example_db")
    outdir = os.path.join(GENEPATH, "test_main_from_parse_force")
    threads = 1
    force = True
    # Create parser
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
    args.argv = "cmd test_main_exist_empty-m1"
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
    for gen in genomes[1:]:
        outgen = os.path.join(listdir, f"{dname}-getEntry_gen_{gen}.txt")
        refgen = os.path.join(EXPPATH, f"exp_wronglist-getEntry_gen_{gen}.txt")
        shutil.copyfile(refgen, outgen)
        outprt = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        refprt = os.path.join(EXPPATH, f"exp_wronglist-getEntry_prt_{gen}.txt")
        shutil.copyfile(refprt, outprt)
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
    al.main_from_parse(args)
    # Check content of listdir
    genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    for gen in genomes:
        assert os.path.isfile(os.path.join(listdir, f"{dname}-getEntry_gen_{gen}.txt"))
        assert os.path.isfile(os.path.join(listdir, f"{dname}-getEntry_prt_{gen}.txt"))
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
    assert tutil.compare_order_content(out_concat, exp_concat)
    # Check content of treedir
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    assert tutil.compare_order_content(out_grp, exp_grp)
    # Check logs
    out, err = capsys.readouterr()
    logfile = os.path.join(outdir, "PanACoTA-align_TEST4-fromparse.log.details")
    with open(logfile, "r") as lc:
        log_content = lc.readlines()
    assert ("Reading PersGenome and constructing lists of missing genomes in "
            "each family") in out
    assert "Extracting proteins and genes from all genomes" in out
    for gen in genomes:
        assert "Extracting proteins and genes from {}".format(gen) in " ".join(log_content)
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in out
    for fam in fams:
        assert "Checking extractions for family {}".format(fam) in " ".join(log_content)
        assert "Aligning family {}".format(fam) in " ".join(log_content)
        assert "Back-translating family {}".format(fam) in " ".join(log_content)
    assert "Concatenating all alignment files" in out
    assert "Grouping alignments per genome" in out
    assert "END" in out
