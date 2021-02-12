#!/usr/bin/env python3

"""
Functional tests for pangenome subcommand
"""
import glob
import os
import shutil
import subprocess
import pytest

import PanACoTA.subcommands.pangenome as pan

LOGFILE_BASE = "func_test_pangenome"
PAN_DIR = os.path.join("test", "data", "pangenome")
TEST_FILES = os.path.join(PAN_DIR, "test_files")
DBPATH = os.path.join(TEST_FILES, "example_db", "Proteins")
EXP_FILES = os.path.join(PAN_DIR, "exp_files")
GENEPATH = os.path.join(PAN_DIR, "generated_by_func-tests")


@pytest.fixture(autouse=True)
def setup_teardown_module():
    """
    Remove log files at the end of this test module
    """
    # Init logger to level detail (15)
    # utils.init_logger(LOGFILE_BASE, 0, 'test_func_pan', verbose=1)
    if not os.path.isdir(GENEPATH):
        os.mkdir(GENEPATH)
    print("setup")

    yield
    # os.remove(LOGFILE_BASE + ".log")
    # os.remove(LOGFILE_BASE + ".log.details")
    # os.remove(LOGFILE_BASE + ".log.err")
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("Removed log files")


def test_main_from_parse():
    """
    Test main when we give the output of the parser
    """
    import argparse
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "testFromParsePAN4"
    min_id = 0.8
    args = argparse.Namespace()
    args.lstinfo_file = lstinfo
    args.dataset_name = name
    args.dbpath = os.path.join(GENEPATH, "database")
    # copy db_path folder to output folder, as it will modify it
    shutil.copytree(DBPATH, args.dbpath)
    args.min_id = min_id
    args.outdir = GENEPATH
    args.clust_mode = 1
    args.spedir = None
    args.threads = 1
    args.outfile = None
    args.verbose = 0
    args.quiet = False
    args.argv = ["pangenome", "pan.py", "test_main_from_parse"]
    # Run main_from_parse
    pan.main_from_parse(args)
    # Check rprt bank was created, and in expected location
    prtbank = os.path.join(args.dbpath, name + ".All.prt")
    assert os.path.isfile(prtbank)

    # Check presence of tmp folder
    tmp_base = os.path.join(GENEPATH, "tmp_testFromParsePAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # check presence of mmseq cluster files
    cluster = os.path.join(GENEPATH, name + ".All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 4
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_files = glob.glob(os.path.join(GENEPATH, "PanGenome-testFromParsePAN4*"))
    to_check = [".tsv.lst", ".tsv.lst.quali.txt", ".tsv.lst.quanti.txt", ".tsv.lst.summary.txt"]
    found = []
    pangenome_file = ""
    for f in pan_files:
        for c in to_check:
            if f.endswith(c):
                found.append(c)
            if f.endswith(".tsv.lst"):
                pangenome_file = f
    assert set(found) == set(to_check)
    # Check content of pangenome
    exp_pan = os.path.join(EXP_FILES, "exp_pangenome-4genomes.lst")
    # Check that all families are as expected. Compare lines without the family number
    with open(exp_pan, "r") as ep, open(pangenome_file, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)


def test_main(caplog):
    """
    Test that from empty directory, it creates all expected files and
    returns correct logs
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "testPAN4"
    min_id = 0.8
    outdir = GENEPATH
    clust_mode = 1
    spe_dir = None
    threads = 1
    cmd = "cmd"
    used_dbpath = os.path.join(GENEPATH, "database")
    # copy db_path folder to output folder, as it will modify it
    shutil.copytree(DBPATH, used_dbpath)
    out_panfile = os.path.join(outdir, "PanGenome-testPAN4.All.prt-clust-0.8-mode1_")
    assert pan.main(cmd, lstinfo, name, used_dbpath, min_id, outdir, clust_mode,
                    spe_dir, threads, verbose=2).startswith(out_panfile)
    # Checl creation of prt bank
    prtbank = os.path.join(used_dbpath, "testPAN4.All.prt")
    assert os.path.isfile(prtbank)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_testPAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_files = glob.glob(os.path.join(GENEPATH, "PanGenome-testPAN4*"))
    to_check = [".tsv.lst", ".tsv.lst.quali.txt", ".tsv.lst.quanti.txt", ".tsv.lst.summary.txt"]
    found = []
    pangenome_file = ""
    for f in pan_files:
        for c in to_check:
            if f.endswith(c):
                found.append(c)
            if f.endswith(".tsv.lst"):
                pangenome_file = f
    assert set(found) == set(to_check)
    # Check content of pangenome
    exp_pan = os.path.join(EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pangenome_file, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check log content
    assert ("Building bank with all proteins to test/data/pangenome/"
            "generated_by_func-tests/database/testPAN4.All.prt") in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text


def test_main_prt_exist(caplog):
    """
    Test that when the prt bank already exists, it writes it to the logs, and
    continues using it.
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "test2PAN4"
    min_id = 0.8
    outdir = GENEPATH
    clust_mode = 1
    spe_dir = None
    threads = 1
    cmd = "cmd"
    used_dbpath = os.path.join(GENEPATH, "database")
    # copy db_path folder to output folder, as it will modify it
    shutil.copytree(DBPATH, used_dbpath)

    # Copy prt bank to database folder, so that it is not created again
    src_prt_bank = os.path.join(EXP_FILES, "exp_EXEM.All.prt")
    dest_prt_bank = os.path.join(used_dbpath, "test2PAN4.All.prt")
    shutil.copyfile(src_prt_bank, dest_prt_bank)

    out_panfile = os.path.join(outdir, "PanGenome-test2PAN4.All.prt-clust-0.8-mode1_")
    assert pan.main(cmd, lstinfo, name, used_dbpath, min_id, outdir, clust_mode, spe_dir,
                    threads, verbose=15).startswith(out_panfile)

    # Check presence of mmseq DB files
    msdb = os.path.join(GENEPATH, "test2PAN4.All.prt-msDB")
    assert os.path.isfile(msdb)
    assert os.path.isfile(msdb + ".index")
    assert os.path.isfile(msdb + ".lookup")
    assert os.path.isfile(msdb + "_h")
    assert os.path.isfile(msdb + "_h.index")
    # Check presence of mmseq cluster files
    cluster = os.path.join(outdir, "test2PAN4.All.prt-clust-0.8-mode1_*")
    clust_files = glob.glob(cluster)
    assert len(clust_files) == 4
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_files = glob.glob(os.path.join(GENEPATH, "PanGenome-test2PAN4*"))
    to_check = [".tsv.lst", ".tsv.lst.quali.txt", ".tsv.lst.quanti.txt", ".tsv.lst.summary.txt"]
    found = []
    pangenome_file = ""
    for f in pan_files:
        for c in to_check:
            if f.endswith(c):
                found.append(c)
            if f.endswith(".tsv.lst"):
                pangenome_file = f
    assert set(found) == set(to_check)
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_test2PAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check content of pangenome
    exp_pan = os.path.join(EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pangenome_file, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check log content
    assert ("Protein bank test/data/pangenome/generated_by_func-tests/"
            "database/test2PAN4.All.prt already exists. It will "
            "be used by mmseqs.") in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text


def test_main_spedir(caplog):
    """
    Test that from empty directory, it creates all expected files and
    returns correct logs
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "test3PAN4"
    min_id = 0.8
    outdir = GENEPATH
    clust_mode = 1
    spe_dir = os.path.join(GENEPATH, "spedir")
    threads = 1
    cmd = "cmd"
    used_dbpath = os.path.join(GENEPATH, "database")
    # copy db_path folder to output folder, as it will modify it
    shutil.copytree(DBPATH, used_dbpath)

    out_panfile = os.path.join(outdir, "PanGenome-test3PAN4.All.prt-clust-0.8-mode1_")
    assert pan.main(cmd, lstinfo, name, used_dbpath, min_id, outdir, clust_mode, spe_dir,
                    threads, verbose=15).startswith(out_panfile)
    # Checl creation of prt bank
    prtbank = os.path.join(spe_dir, "test3PAN4.All.prt")
    assert os.path.isfile(prtbank)
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
    assert len(clust_files) == 4
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_test3PAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_files = glob.glob(os.path.join(GENEPATH, "PanGenome-test3PAN4*"))
    to_check = [".tsv.lst", ".tsv.lst.quali.txt", ".tsv.lst.quanti.txt", ".tsv.lst.summary.txt"]
    found = []
    pangenome_file = ""
    for f in pan_files:
        for c in to_check:
            if f.endswith(c):
                found.append(c)
            if f.endswith(".tsv.lst"):
                pangenome_file = f
    assert set(found) == set(to_check)
    # Check content of pangenome
    exp_pan = os.path.join(EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pangenome_file, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check log content
    assert ("Building bank with all proteins to test/data/pangenome/"
            "generated_by_func-tests/spedir/test3PAN4.All.prt") in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text


def test_main_outfile(caplog):
    """
    Test that when giving a name for pangenome file, it creates expected files
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "test4PAN4"
    min_id = 0.8
    outdir = GENEPATH
    clust_mode = 1
    spe_dir = None
    threads = 1
    cmd = "cmd"
    outfile = "my_pangenome"
    used_dbpath = os.path.join(GENEPATH, "database")
    # copy db_path folder to output folder, as it will modify it
    shutil.copytree(DBPATH, used_dbpath)

    assert pan.main(cmd, lstinfo, name, used_dbpath, min_id, outdir, clust_mode, spe_dir,
                    threads, outfile=outfile) == os.path.join(outdir, outfile)

    prtbank = os.path.join(used_dbpath, "test4PAN4.All.prt")
    assert os.path.isfile(prtbank)
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
    assert len(clust_files) == 4
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
    exp_pan = os.path.join(EXP_FILES, "exp_pangenome-4genomes.lst")
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
    assert ("Building bank with all proteins to test/data/pangenome/"
            "generated_by_func-tests/database/test4PAN4.All.prt") in caplog.text
    assert "Creating database" in caplog.text
    assert "Clustering proteins..." in caplog.text
    assert "Converting mmseqs results to pangenome file" in caplog.text
    assert "Pangenome has 16 families" in caplog.text
    assert "Retrieving information from pan families" in caplog.text
    assert "Saving all information to a binary file for later use" in caplog.text
    assert "Generating qualitative and quantitative matrix, and summary file" in caplog.text


def test_pangenome_all():
    """
    Test when calling pangenome from command line, it runs and gives expected output files
    """
    lstinfo = os.path.join(TEST_FILES, "list_to_pan.txt")
    name = "testAllPAN4"
    min_id = 0.8
    outdir = GENEPATH
    clust_mode = 1
    spe_dir = None
    threads = 1
    used_dbpath = os.path.join(GENEPATH, "database")
    # copy db_path folder to output folder, as it will modify it
    shutil.copytree(DBPATH, used_dbpath)

    cmd = f"PanACoTA pangenome -l {lstinfo} -n {name} -d {used_dbpath} -o {outdir} -vv"
    ret = subprocess.call(cmd.split())
    assert ret == 0

    prtbank = os.path.join(used_dbpath, "testAllPAN4.All.prt")
    assert os.path.isfile(prtbank)
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
    assert len(clust_files) == 4
    # Check presence of tmp folder
    tmp_base = os.path.join(outdir, "tmp_testAllPAN4.All.prt_0.8-mode1_*")
    assert len(glob.glob(tmp_base)) == 1
    # Check presence of pangenome files (pangenome, matrices, summary)
    pan_files = glob.glob(os.path.join(GENEPATH, "PanGenome-testAllPAN4*"))
    to_check = [".tsv.lst", ".tsv.lst.quali.txt", ".tsv.lst.quanti.txt", ".tsv.lst.summary.txt"]
    found = []
    pangenome_file = ""
    for f in pan_files:
        for c in to_check:
            if f.endswith(c):
                found.append(c)
            if f.endswith(".tsv.lst"):
                pangenome_file = f
    assert set(found) == set(to_check)
    # Check content of pangenome
    exp_pan = os.path.join(EXP_FILES, "exp_pangenome-4genomes.lst")
    with open(exp_pan, "r") as ep, open(pangenome_file, "r") as panf:
        lines_exp = []
        lines_out = []
        for line_exp in ep:
            lines_exp.append(" ".join(line_exp.split()[1:]))
        for line_out in panf:
            lines_out.append(" ".join(line_out.split()[1:]))
    assert len(lines_exp) == len(lines_out)
    assert set(lines_exp) == set(lines_out)
    # Check presence of log files, and that .err is empty
    log_base = os.path.join(outdir, "PanACoTA-pangenome_testAllPAN4.log")
    assert os.path.isfile(log_base)
    assert os.path.isfile(log_base + ".details")
    assert os.path.isfile(log_base + ".err")
    with open(log_base + ".err") as errf:
        lines = errf.readlines()
    assert lines == []
