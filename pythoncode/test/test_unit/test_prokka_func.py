#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import pytest
import os
import shutil
import glob
import pipelinepackage.prokka_functions as pfunc


def test_count_tbl():
    """
    Count the different features found in the tbl file, and return
    nbcont, nbCDS, nbGene, nbCRISPR
    """
    tblfile = os.path.join("test", "data", "test_files", "prokka_out_for_test.tbl")
    ncont, ncds, ngene, ncris = pfunc.count_tbl(tblfile)
    assert ncont == 7
    assert ncds == 13
    assert ngene == 15
    assert ncris == 2


def test_count_headers():
    """
    Count how many sequences there are in the given multi-fasta file
    """
    seqfile = os.path.join("test", "data", "genomes", "genome4.fasta")
    nb = pfunc.count_headers(seqfile)
    assert nb == 5


def test_check_prokka_no_outdir(capsys):
    """
    Test that prokka returns the right error message when output directory does not exist
    """
    outdir = "toto"
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "Prokka could not run properly. Look at prokka.log for more information.\n"


def test_check_prokka_notbl(capsys):
    """
    Check that check_prokka returns false when a tbl file is missing, and an error message
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-misstbl"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-misstbl original_name.fna: no .tbl file\n"
    os.remove(os.path.join(ori_dir, name + ".faa"))
    os.remove(os.path.join(ori_dir, name + ".ffn"))


def test_check_prokka_nofaa(capsys):
    """
    Check that check_prokka returns false when a faa file is missing, and an error message
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-missfaa"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".tbl"), os.path.join(ori_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missfaa original_name.fna: no .faa file\n"
    os.remove(os.path.join(ori_dir, name + ".tbl"))
    os.remove(os.path.join(ori_dir, name + ".ffn"))


def test_check_prokka_noffn(capsys):
    """
    Check that check_prokka returns false when a ffn file is missing, and an error message
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-missffn"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".tbl"), os.path.join(ori_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missffn original_name.fna: no .ffn file\n"
    os.remove(os.path.join(ori_dir, name + ".tbl"))
    os.remove(os.path.join(ori_dir, name + ".faa"))


def test_check_prokka_wrong_cont(capsys):
    """
    Check that check_prokka returns an error message when the number of contigs in tbl
    file is not as expected
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 10
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == ("prokka_out_for_test original_name.fna: no matching number of contigs; "
                   "nbcontig=10; in tbl =7\n")


def test_check_prokka_wrong_tblCDS(capsys):
    """
    Check that check_prokka returns an error message when the number of CDS in tbl
    file is different from the number of headers in faa file
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-wrongCDS"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err.split("\n")[0] == ("prokka_out_for_test-wrongCDS original_name.fna: "
                                  "no matching number of proteins between tbl and faa; "
                                  "faa=13; in tbl =12")
    assert err.split("\n")[1] == ("prokka_out_for_test-wrongCDS original_name.fna: "
                                  "no matching number of genes between tbl and ffn; "
                                  "ffn=17; in tbl =14genes 2CRISPR")
    os.remove(os.path.join(ori_dir, name + ".ffn"))
    os.remove(os.path.join(ori_dir, name + ".faa"))


def test_check_prokka_wrong_tblCRISPR(capsys):
    """
    Check that check_prokka returns an error message when the number of headers in ffn
    file is different from the number of CDS + CRISPR in tbl file (1CRISPR in tbl, 2 in ffn)
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test-wrongtblCRISP"
    ori_dir = os.path.join("test", "data", "test_files")
    ori_file = "prokka_out_for_test"
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".ffn"), os.path.join(ori_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_dir, ori_file + ".faa"), os.path.join(ori_dir, name + ".faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == ("prokka_out_for_test-wrongtblCRISP original_name.fna: "
                   "no matching number of genes between tbl and ffn; "
                   "ffn=17; in tbl =15genes 1CRISPR\n")
    os.remove(os.path.join(ori_dir, name + ".ffn"))
    os.remove(os.path.join(ori_dir, name + ".faa"))


def test_check_prokka_ok():
    """
    Check that everything is ok with prokka results (tbl, faa and ffn files exist,
    and number of CDS, CRISPR and genes correspond between them)
    """
    outdir = os.path.join("test", "data", "test_files")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert pfunc.check_prokka(outdir, logf, name, gpath, nbcont)


def test_run_prokka_out_exists_ok(capsys):
    """
    Test that when the output directory already exists, and files inside are OK,
    run_prokka returns True, with a warning message indicating that prokka did not rerun.
    """
    gpath = "path/to/nogenome/original_name.fna"
    outdir = os.path.join("test", "data", "test_files")
    cores_prokka = 1
    name = "prokka_out_for_test"
    force = False
    nbcont = 7
    arguments = (gpath, outdir, outdir, cores_prokka, name, force, nbcont)
    assert pfunc.run_prokka(arguments)
    _, err = capsys.readouterr()
    assert err.endswith("Prokka results folder already exists. Prokka did not run again, "
                        "formatting step used already generated results of Prokka in "
                        "test/data/test_files. If you want to re-run prokka, first remove this "
                        "result folder, or use '-F' or '--force' option if you want to rerun "
                        "prokka for all genomes.\n")


def test_run_prokka_out_exists_error(capsys):
    """
    Test that when the output directory already exists, and some file is missing,
    run_prokka returns False, and writes the warning message saying that prokka did not
    rerun, + the warning message for the missing file(s).
    """
    gpath = "path/to/nogenome/original_name.fna"
    outdir = os.path.join("test", "data", "test_files")
    cores_prokka = 1
    name = "prokka_out_for_test-wrongCDS"
    force = False
    nbcont = 7
    arguments = (gpath, outdir, outdir, cores_prokka, name, force, nbcont)
    assert not pfunc.run_prokka(arguments)
    _, err = capsys.readouterr()
    assert ("Prokka results folder already exists. Prokka did not run again, "
            "formatting step used already generated results of Prokka in "
            "test/data/test_files. If you want to re-run prokka, first remove this "
            "result folder, or use '-F' or '--force' option if you want to rerun "
            "prokka for all genomes.") in err
    assert "prokka_out_for_test-wrongCDS original_name.fna: no .ffn file"
    assert "prokka_out_for_test-wrongCDS original_name.fna: no .faa file"


def test_run_prokka_out_exists_force():
    """
    Test that when the output directory already exists with wrong files, but force is on,
    prokka is rerun and outputs the right files
    """
    gpath = os.path.join("test", "data", "genomes", "H299_H561.fasta")
    outdir = os.path.join("test", "data")
    res_dir = outdir
    name = "test_runprokka_H299"
    # Put empty tbl, faa, ffn files in prokka output dir, to check that they are overridden
    open(os.path.join(outdir, name + ".tbl"), "w").close()
    open(os.path.join(outdir, name + ".faa"), "w").close()
    open(os.path.join(outdir, name + ".ffn"), "w").close()
    cores_prokka = 5
    force = "--force"
    nbcont = 3
    arguments = (gpath, outdir, res_dir, cores_prokka, name, force, nbcont)
    assert pfunc.run_prokka(arguments)
    for filename in glob.glob(os.path.join(outdir, name + ".*")):
        os.remove(filename)
    os.remove(os.path.join(outdir, "H299_H561.fasta-prokka.log"))


def test_run_prokka_out_doesnt_exist():
    """
    Test that when the output directory does not exist, it creates it, and runs prokka
    with all expected outfiles
    """
    gpath = os.path.join("test", "data", "genomes", "H299_H561.fasta")
    outdir = os.path.join("test", "data", "prokkaRes")
    res_dir = os.path.join("test", "data")
    cores_prokka = 5
    name = "test_runprokka_H299"
    force = False
    nbcont = 3
    arguments = (gpath, outdir, res_dir, cores_prokka, name, force, nbcont)
    assert pfunc.run_prokka(arguments)
    shutil.rmtree(outdir)
    os.remove(os.path.join(res_dir, "H299_H561.fasta-prokka.log"))


def test_run_prokka_out_problem_running(capsys):
    """
    Check that when a problem occurs while trying to run prokka, run_prokka returns False,
    and the error message indicating to read in the log why it couldn't run
    """
    gpath = os.path.join("test", "data", "genomes", "H299 H561.fasta")
    outdir = os.path.join("test", "data", "prokkaRes")
    res_dir = os.path.join("test", "data")
    cores_prokka = 5
    name = "test_runprokka_H299"
    force = False
    nbcont = 3
    arguments = (gpath, outdir, res_dir, cores_prokka, name, force, nbcont)
    assert not pfunc.run_prokka(arguments)
    _, err = capsys.readouterr()
    assert ("Prokka could not run properly. Look at test/data/H299 H561.fasta-prokka.log "
            "for more information.") in err
    os.remove(os.path.join(res_dir, "H299 H561.fasta-prokka.log"))


def test_run_all_1by1():
    """
    Check that when running with 3 threads (not parallel), prokka runs as expected,
    and returns True for each genome
    """
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join("test", "data", "genomes", genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join("test", "data", "genomes", genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, 456464645, 1, 465]}
    threads = 3
    force = False
    prok_folder = os.path.join("test", "data")
    res_dir = os.path.join("test", "data")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder, res_dir)
    assert final[genome1]
    assert final[genome2]
    shutil.rmtree(os.path.join(prok_folder, genome1 + "-prokkaRes"))
    shutil.rmtree(os.path.join(prok_folder, genome2 + "-prokkaRes"))
    os.remove(os.path.join(res_dir, genome1 + "-prokka.log"))
    os.remove(os.path.join(res_dir, genome2 + "-prokka.log"))


def test_run_all_parallel_more_threads():
    """
    Check that there is no problem when running with more threads than genomes (each genome
    uses nb_threads/nb_genome threads)
    """
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join("test", "data", "genomes", genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join("test", "data", "genomes", genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    prok_folder = os.path.join("test", "data")
    res_dir = os.path.join("test", "data")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder, res_dir)
    assert final[genome1]
    assert final[genome2]
    shutil.rmtree(os.path.join(prok_folder, genome1 + "-prokkaRes"))
    shutil.rmtree(os.path.join(prok_folder, genome2 + "-prokkaRes"))
    os.remove(os.path.join(res_dir, genome1 + "-prokka.log"))
    os.remove(os.path.join(res_dir, genome2 + "-prokka.log"))


def test_run_all_parallel_less_threads():
    """
    Check that there is no problem when running with less threads than genomes (each genomes
    uses 2 threads)
    All genomes should run well, except for "genome3", where there are 3 contigs, and
    we announce 1 contig, so check_prokka should return false.
    """
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    gnames = ["H299_H561.fasta", "A_H738.fasta", "genome1.fasta", "genome2.fasta", "genome3.fasta"]
    gpaths = [os.path.join("test", "data", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["test_runall_1by1_1", gpaths[0], 12656, 3, 1],
               gnames[1]: ["test_runall_1by1_2", gpaths[1], 456464645, 1, 1],
               gnames[2]: ["test_runall_1by1_2", gpaths[2], 456464645, 4, 1],
               gnames[3]: ["test_runall_1by1_2", gpaths[3], 456464645, 3, 1],
               gnames[4]: ["test_runall_1by1_2", gpaths[4], 456464645, 1, 1]
              }
    threads = 4
    force = False
    prok_folder = os.path.join("test", "data")
    res_dir = os.path.join("test", "data")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder, res_dir)
    assert final[gnames[0]]
    assert final[gnames[1]]
    assert final[gnames[2]]
    assert final[gnames[3]]
    assert not final[gnames[4]]
    for name in gnames:
        shutil.rmtree(os.path.join(prok_folder, name + "-prokkaRes"))
        os.remove(os.path.join(res_dir, name + "-prokka.log"))
