#!/usr/bin/env python3

"""
Functional tests for genomeAPCAT pangenome
"""
import glob
import os
import shutil
import subprocess

import genomeAPCAT.subcommands.pangenome as pan
from genomeAPCAT import utils


LOGFILE_BASE = "func_test_pangenome"
TEST_FILES = os.path.join("test", "data", "pangenome", "test_files")
DBPATH = os.path.join(TEST_FILES, "example_db", "Proteins")
PATH_EXP_FILES = os.path.join("test", "data", "pangenome", "exp_files")


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
    Test that from empty directory, it creates all expected files and
    returns correct logs
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "testPAN4"
    min_id = 0.8
    outdir = "test_main_pangenome_dir"
    clust_mode = 1
    spe_dir = None
    threads = 1

    pan.main(lstinfo, name, DBPATH, min_id, outdir, clust_mode, spe_dir, threads)
    prtbank = os.path.join(DBPATH, "testPAN4.All.prt")
    assert os.path.isfile(prtbank)
    os.remove(prtbank)
    assert os.path.isdir(outdir)
    # Check presence of mmseq DB files
    msdb = os.path.join(outdir, "testPAN4.All.prt-msDB")
    assert os.path.isfile(msdb)
    assert os.path.isfile(msdb + ".index")
    assert os.path.isfile(msdb + ".lookup")
    assert os.path.isfile(msdb + "_h")
    assert os.path.isfile(msdb + "_h.index")
    # Check presence of mmseq cluster files
    cluster = os.path.join(outdir, "testPAN4.All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 3
    clust_base = [cl for cl in clust_files if not cl.endswith(".index")
                  and not cl.endswith(".tsv")][0]
    assert os.path.isfile(clust_base + ".index")
    assert os.path.isfile(clust_base + ".tsv")
    assert os.path.isfile(clust_base)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_testPAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_base = os.path.join(outdir, "PanGenome-{}.tsv.lst".format(os.path.basename(clust_base)))
    assert os.path.isfile(pan_base)
    assert os.path.isfile(pan_base + ".quali.txt")
    assert os.path.isfile(pan_base + ".quanti.txt")
    assert os.path.isfile(pan_base + ".summary.txt")
    # Check content of pangenome
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pan_base, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check log content
    assert "Building bank with all proteins to testPAN4.All.prt" in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text
    # Remove output directory
    shutil.rmtree(outdir)


def test_main_prt_exist(caplog):
    """
    Test that when the prt bank already exists, it writes it to the logs, and
    continues using it.
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "test2PAN4"
    min_id = 0.8
    outdir = "test_prtExist_main_pangenome_dir"
    clust_mode = 1
    spe_dir = None
    threads = 1

    # Copy prt bank to database folder, so that it is not created again
    src_prt_bank = os.path.join(PATH_EXP_FILES, "exp_EXEM.All.prt")
    dest_prt_bank = os.path.join(DBPATH, "test2PAN4.All.prt")
    shutil.copyfile(src_prt_bank, dest_prt_bank)

    pan.main(lstinfo, name, DBPATH, min_id, outdir, clust_mode, spe_dir, threads)
    os.remove(dest_prt_bank)

    assert os.path.isdir(outdir)
    # Check presence of mmseq DB files
    msdb = os.path.join(outdir, "test2PAN4.All.prt-msDB")
    assert os.path.isfile(msdb)
    assert os.path.isfile(msdb + ".index")
    assert os.path.isfile(msdb + ".lookup")
    assert os.path.isfile(msdb + "_h")
    assert os.path.isfile(msdb + "_h.index")
    # Check presence of mmseq cluster files
    cluster = os.path.join(outdir, "test2PAN4.All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 3
    clust_base = [cl for cl in clust_files if not cl.endswith(".index")
                  and not cl.endswith(".tsv")][0]
    assert os.path.isfile(clust_base + ".index")
    assert os.path.isfile(clust_base + ".tsv")
    assert os.path.isfile(clust_base)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_test2PAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_base = os.path.join(outdir, "PanGenome-{}.tsv.lst".format(os.path.basename(clust_base)))
    assert os.path.isfile(pan_base)
    assert os.path.isfile(pan_base + ".quali.txt")
    assert os.path.isfile(pan_base + ".quanti.txt")
    assert os.path.isfile(pan_base + ".summary.txt")
    # Check content of pangenome
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pan_base, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check log content
    assert ("Protein bank test/data/pangenome/test_files/example_db/"
            "Proteins/test2PAN4.All.prt already exists. It will "
            "be used by mmseqs.") in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text
    # Remove output directory
    shutil.rmtree(outdir)


def test_main_spedir(caplog):
    """
    Test that from empty directory, it creates all expected files and
    returns correct logs
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "test3PAN4"
    min_id = 0.8
    outdir = "test_spedir_main_pangenome_dir"
    clust_mode = 1
    spe_dir = outdir + "-spe"
    threads = 1

    pan.main(lstinfo, name, DBPATH, min_id, outdir, clust_mode, spe_dir, threads)
    prtbank = os.path.join(spe_dir, "test3PAN4.All.prt")
    assert os.path.isfile(prtbank)
    shutil.rmtree(spe_dir)
    assert os.path.isdir(outdir)
    # Check presence of mmseq DB files
    msdb = os.path.join(outdir, "test3PAN4.All.prt-msDB")
    assert os.path.isfile(msdb)
    assert os.path.isfile(msdb + ".index")
    assert os.path.isfile(msdb + ".lookup")
    assert os.path.isfile(msdb + "_h")
    assert os.path.isfile(msdb + "_h.index")
    # Check presence of mmseq cluster files
    cluster = os.path.join(outdir, "test3PAN4.All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 3
    clust_base = [cl for cl in clust_files if not cl.endswith(".index")
                  and not cl.endswith(".tsv")][0]
    assert os.path.isfile(clust_base + ".index")
    assert os.path.isfile(clust_base + ".tsv")
    assert os.path.isfile(clust_base)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_test3PAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_base = os.path.join(outdir, "PanGenome-{}.tsv.lst".format(os.path.basename(clust_base)))
    assert os.path.isfile(pan_base)
    assert os.path.isfile(pan_base + ".quali.txt")
    assert os.path.isfile(pan_base + ".quanti.txt")
    assert os.path.isfile(pan_base + ".summary.txt")
    # Check content of pangenome
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pan_base, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check log content
    assert "Building bank with all proteins to test3PAN4.All.prt" in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text
    # Remove output directory
    shutil.rmtree(outdir)


def test_main_outfile(caplog):
    """
    Test that from empty directory, it creates all expected files and
    returns correct logs
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "test4PAN4"
    min_id = 0.8
    outdir = "test_outfile_main_pangenome_dir"
    clust_mode = 1
    spe_dir = outdir + "-spe"
    threads = 1
    outfile = "mypangenome.txt"

    pan.main(lstinfo, name, DBPATH, min_id, outdir, clust_mode, spe_dir, threads, outfile=outfile)
    prtbank = os.path.join(spe_dir, "test4PAN4.All.prt")
    assert os.path.isfile(prtbank)
    shutil.rmtree(spe_dir)
    assert os.path.isdir(outdir)
    # Check presence of mmseq DB files
    msdb = os.path.join(outdir, "test4PAN4.All.prt-msDB")
    assert os.path.isfile(msdb)
    assert os.path.isfile(msdb + ".index")
    assert os.path.isfile(msdb + ".lookup")
    assert os.path.isfile(msdb + "_h")
    assert os.path.isfile(msdb + "_h.index")
    # Check presence of mmseq cluster files
    cluster = os.path.join(outdir, "test4PAN4.All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 3
    clust_base = [cl for cl in clust_files if not cl.endswith(".index")
                  and not cl.endswith(".tsv")][0]
    assert os.path.isfile(clust_base + ".index")
    assert os.path.isfile(clust_base + ".tsv")
    assert os.path.isfile(clust_base)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_test4PAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    outf = os.path.join(outdir, outfile)
    assert os.path.isfile(outf)
    assert os.path.isfile(outf + ".quali.txt")
    assert os.path.isfile(outf + ".quanti.txt")
    assert os.path.isfile(outf + ".summary.txt")
    # Check content of pangenome
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(outf, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check log content
    assert "Building bank with all proteins to test4PAN4.All.prt" in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text
    # Remove output directory
    shutil.rmtree(outdir)


def test_main_from_parse():
    """
    Test main when we give the output of the parser
    """
    import argparse
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "testFromParsePAN4"
    min_id = 0.8
    outdir = "test_all_main-from-parse_pangenome_dir"
    args = argparse.Namespace()
    args.lstinfo_file = lstinfo
    args.dataset_name = name
    args.dbpath = DBPATH
    args.min_id = min_id
    args.outdir = outdir
    args.clust_mode = 1
    args.spedir = None
    args.threads = 1
    args.outfile = None
    args.verbose = 0
    args.quiet = False
    pan.main_from_parse(args)

    prtbank = os.path.join(DBPATH, name + ".All.prt")
    assert os.path.isfile(prtbank)
    os.remove(prtbank)
    assert os.path.isdir(outdir)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_testFromParsePAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    cluster = os.path.join(outdir, name + ".All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 3
    clust_base = [cl for cl in clust_files if not cl.endswith(".index")
                  and not cl.endswith(".tsv")][0]
    pan_base = os.path.join(outdir, "PanGenome-{}.tsv.lst".format(os.path.basename(clust_base)))
    assert os.path.isfile(pan_base)
    assert os.path.isfile(pan_base + ".quali.txt")
    assert os.path.isfile(pan_base + ".quanti.txt")
    assert os.path.isfile(pan_base + ".summary.txt")
    # Check content of pangenome
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pan_base, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Remove output directory
    shutil.rmtree(outdir)


def test_pangenome_all():
    """
    Test when calling pangenome from command line, it runs and gives expected output files
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "testAllPAN4"
    min_id = 0.8
    outdir = "test_all_main_pangenome_dir"

    cmd = "genomeAPCAT pangenome -l {} -n {} -d {} -i {} -o {}".format(lstinfo, name, DBPATH,
                                                                       min_id, outdir)
    ret = subprocess.call(cmd.split())
    assert ret == 0

    prtbank = os.path.join(DBPATH, "testAllPAN4.All.prt")
    assert os.path.isfile(prtbank)
    os.remove(prtbank)
    assert os.path.isdir(outdir)
    # Check presence of mmseq DB files
    msdb = os.path.join(outdir, "testAllPAN4.All.prt-msDB")
    assert os.path.isfile(msdb)
    assert os.path.isfile(msdb + ".index")
    assert os.path.isfile(msdb + ".lookup")
    assert os.path.isfile(msdb + "_h")
    assert os.path.isfile(msdb + "_h.index")
    # Check presence of mmseq cluster files
    cluster = os.path.join(outdir, "testAllPAN4.All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 3
    clust_base = [cl for cl in clust_files if not cl.endswith(".index")
                  and not cl.endswith(".tsv")][0]
    assert os.path.isfile(clust_base + ".index")
    assert os.path.isfile(clust_base + ".tsv")
    assert os.path.isfile(clust_base)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_testAllPAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_base = os.path.join(outdir, "PanGenome-{}.tsv.lst".format(os.path.basename(clust_base)))
    assert os.path.isfile(pan_base)
    assert os.path.isfile(pan_base + ".quali.txt")
    assert os.path.isfile(pan_base + ".quanti.txt")
    assert os.path.isfile(pan_base + ".summary.txt")
    # Check content of pangenome
    exp_pan = os.path.join(PATH_EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pan_base, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check presence of log files, and that .err is empty
    log_base = os.path.join(outdir, "genomeAPCAT-pangenome_testAllPAN4.log")
    assert os.path.isfile(log_base)
    assert os.path.isfile(log_base + ".details")
    assert os.path.isfile(log_base + ".err")
    with open(log_base + ".err") as errf:
        lines = errf.readlines()
    assert lines == []
    # Remove output directory
    shutil.rmtree(outdir)
