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
    tblfile = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes",
                           "prokka_out_for_test.tbl")
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
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "out_test_notbl")
    name = "prokka_out_for_test-misstbl"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    logf = "prokka.log"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-misstbl original_name-error.fna: no .tbl file\n"
    shutil.rmtree(out_dir)


def test_check_prokka_sevtbl(capsys):
    """
    Check that check_prokka returns false when there is more than 1 tbl file,
    and an error message
    """
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "out_test_notbl")
    name = "prokka_out_for_test-misstbl"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + "2.tbl"))
    logf = "prokka.log"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-misstbl original_name-error.fna: several .tbl files\n"
    shutil.rmtree(out_dir)


def test_check_prokka_nofaa(capsys):
    """
    Check that check_prokka returns false when a faa file is missing, and an error message
    """
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "out_test_nofaa")
    name = "prokka_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missfaa original_name.fna: no .faa file\n"
    shutil.rmtree(out_dir)


def test_check_prokka_sevfaa(capsys):
    """
    Check that check_prokka returns false when there is more than 1 faa file,
    and an error message
    """
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "out_test_nofaa")
    name = "prokka_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + "2.faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missfaa original_name.fna: several .faa files\n"
    shutil.rmtree(out_dir)


def test_check_prokka_noffn(capsys):
    """
    Check that check_prokka returns false when a ffn file is missing, and an error message
    """
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "out_test_noffn")
    name = "prokka_out_for_test-missffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missffn original_name.fna: no .ffn file\n"
    shutil.rmtree(out_dir)


def test_check_prokka_sevffn(capsys):
    """
    Check that check_prokka returns false when there is more than 1 ffn file,
    and an error message
    """
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "out_test_noffn")
    name = "prokka_out_for_test-missffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + "2.ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == "prokka_out_for_test-missffn original_name.fna: several .ffn files\n"
    shutil.rmtree(out_dir)


def test_check_prokka_wrong_cont(capsys):
    """
    Check that check_prokka returns an error message when the number of contigs in tbl
    file is not as expected
    """
    outdir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
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
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "res_checkProkkaWrongTbl")
    os.makedirs(out_dir)
    name = "prokka_out_for_test-wrongCDS"
    tblfile = os.path.join("test", "data", "test_files", name + ".tbl")
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(tblfile, os.path.join(out_dir, name + ".tbl"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err.split("\n")[0] == ("prokka_out_for_test-wrongCDS original_name.fna: "
                                  "no matching number of proteins between tbl and faa; "
                                  "faa=13; in tbl =12")
    assert err.split("\n")[1] == ("prokka_out_for_test-wrongCDS original_name.fna: "
                                  "no matching number of genes between tbl and ffn; "
                                  "ffn=17; in tbl =14genes 2CRISPR")
    shutil.rmtree(out_dir)


def test_check_prokka_wrong_tblCRISPR(capsys):
    """
    Check that check_prokka returns an error message when the number of headers in ffn
    file is different from the number of CDS + CRISPR in tbl file (1CRISPR in tbl, 2 in ffn)
    """
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "res_checlProkkaWrongCRISPR")
    os.makedirs(out_dir)
    name = "prokka_out_for_test-wrongtblCRISP"
    tblfile = os.path.join("test", "data", "test_files", name + ".tbl")
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(tblfile, os.path.join(out_dir, name + ".tbl"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont)
    _, err = capsys.readouterr()
    assert err == ("prokka_out_for_test-wrongtblCRISP original_name.fna: "
                   "no matching number of genes between tbl and ffn; "
                   "ffn=17; in tbl =15genes 1CRISPR\n")
    shutil.rmtree(out_dir)


def test_check_prokka_ok():
    """
    Check that everything is ok with prokka results (tbl, faa and ffn files exist,
    and number of CDS, CRISPR and genes correspond between them)
    """
    outdir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
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
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont)
    assert pfunc.run_prokka(arguments)
    _, err = capsys.readouterr()
    assert err.endswith("Prokka results folder already exists. Prokka did not run again, "
                        "formatting step used already generated results of Prokka in "
                        "test/data/test_files/original_name.fna-prokkaRes. If you want "
                        "to re-run prokka, first remove this "
                        "result folder, or use '-F' or '--force' option if you want to rerun "
                        "prokka for all genomes.\n")


def test_run_prokka_out_exists_error(capsys):
    """
    Test that when the output directory already exists, and some file is missing,
    run_prokka returns False, and writes the warning message saying that prokka did not
    rerun, + the warning message for the missing file(s).
    """
    ori_prok_dir = os.path.join("test", "data", "test_files", "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    new_prok_dir = os.path.join("test", "data", "test_files", "original_name-error-prokkaRes")
    name = "prokka_out_for_test-wrongCDS"
    os.makedirs(new_prok_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(new_prok_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(new_prok_dir, name + ".faa"))
    gpath = "path/to/nogenome/original_name-error"
    outdir = os.path.join("test", "data", "test_files")
    cores_prokka = 1
    force = False
    nbcont = 7
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont)
    assert not pfunc.run_prokka(arguments)
    _, err = capsys.readouterr()
    assert ("Prokka results folder already exists. Prokka did not run again, "
            "formatting step used already generated results of Prokka in "
            "test/data/test_files/original_name-error-prokkaRes. If you want to re-run prokka, first remove this "
            "result folder, or use '-F' or '--force' option if you want to rerun "
            "prokka for all genomes.") in err
    assert "prokka_out_for_test-wrongCDS original_name.fna: no .ffn file"
    assert "prokka_out_for_test-wrongCDS original_name.fna: no .faa file"
    shutil.rmtree(new_prok_dir)


def test_run_prokka_out_exists_force():
    """
    Test that when the output directory already exists with wrong files, but force is on,
    prokka is rerun and outputs the right files
    """
    gpath = os.path.join("test", "data", "genomes", "H299_H561.fasta")
    outdir = os.path.join("test", "data")
    out_prokdir = os.path.join(outdir, "H299_H561.fasta-prokkaRes")
    name = "test_runprokka_H299"
    # Put empty tbl, faa, ffn files in prokka output dir, to check that they are overridden
    os.makedirs(out_prokdir)
    open(os.path.join(out_prokdir, name + ".tbl"), "w").close()
    open(os.path.join(out_prokdir, name + ".faa"), "w").close()
    open(os.path.join(out_prokdir, name + ".ffn"), "w").close()
    cores_prokka = 5
    force = True
    nbcont = 3
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont)
    assert pfunc.run_prokka(arguments)
    # Check content of tbl, ffn and faa files
    exp_dir = os.path.join("test", "data", "exp_files", "H299_H561.fasta-gembase.fna-prokkaRes",
                           "test_runprokka_H299")
    out_tbl = os.path.join(out_prokdir, name + ".tbl")
    out_faa = os.path.join(out_prokdir, name + ".faa")
    out_ffn = os.path.join(out_prokdir, name + ".ffn")
    assert os.path.isfile(out_tbl)
    with open(exp_dir + ".tbl", "r") as expf, open(out_tbl, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    assert os.path.isfile(out_faa)
    with open(exp_dir + ".faa", "r") as expf, open(out_faa, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    assert os.path.isfile(out_ffn)
    with open(exp_dir + ".ffn", "r") as expf, open(out_ffn, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    shutil.rmtree(out_prokdir)
    os.remove(os.path.join(outdir, "H299_H561.fasta-prokka.log"))


def test_run_prokka_out_doesnt_exist():
    """
    Test that when the output directory does not exist, it creates it, and runs prokka
    with all expected outfiles
    """
    gpath = os.path.join("test", "data", "genomes", "H299_H561.fasta")
    prok_dir = os.path.join("test", "data")
    out_dir = os.path.join(prok_dir, "H299_H561.fasta-prokkaRes")
    cores_prokka = 5
    name = "test_runprokka_H299"
    force = False
    nbcont = 3
    arguments = (gpath, prok_dir, cores_prokka, name, force, nbcont)
    assert pfunc.run_prokka(arguments)
    shutil.rmtree(out_dir)
    os.remove(os.path.join(prok_dir, "H299_H561.fasta-prokka.log"))


def test_run_prokka_out_problem_running(capsys):
    """
    Check that when a problem occurs while trying to run prokka, run_prokka returns False,
    and the error message indicating to read in the log why it couldn't run
    """
    gpath = os.path.join("test", "data", "genomes", "H299 H561.fasta")
    outdir = os.path.join("test", "data")
    cores_prokka = 5
    name = "test_runprokka_H299"
    force = False
    nbcont = 3
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont)
    assert not pfunc.run_prokka(arguments)
    _, err = capsys.readouterr()
    assert ("Prokka could not run properly. Look at test/data/H299 H561.fasta-prokka.log "
            "for more information.") in err
    os.remove(os.path.join(outdir, "H299 H561.fasta-prokka.log"))


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
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
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
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
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
    Genomes H299 and A_H738 should run well, but genomes genome* have problems (no CDS found),
    so check_prokka should return false.
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
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
    assert final[gnames[0]]
    assert final[gnames[1]]
    assert not final[gnames[2]]
    assert not final[gnames[3]]
    assert not final[gnames[4]]
    for name in gnames:
        shutil.rmtree(os.path.join(prok_folder, name + "-prokkaRes"))
        os.remove(os.path.join(res_dir, name + "-prokka.log"))
