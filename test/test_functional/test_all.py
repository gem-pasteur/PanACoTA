#!/usr/bin/env python3

"""
Functional tests for PanACoTA 'all' module
"""
import pytest
import os
import shutil
import glob
import argparse

import PanACoTA.subcommands.all_modules as allm

DATADIR = os.path.join("test", "data", "all")
GENEPATH = os.path.join(DATADIR, "generated_by_func-tests")

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
    os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    print("teardown")


def test_main_default_qc_only():
    """
    Test downloading 104099 genomes and analysis until QC, with verbose==2
    """
    cmd = "cmd"
    outdir = os.path.join(GENEPATH, "main_all_default_qc")

    # Common arguments: outdir, threads, verbose, quiet
    args_all = (outdir, 1, 2, False)
    # args for prepare:
    # NCBI_species_taxid (int), NCBI_species (str), levels (str), tmp_dir (str),
    # norefseq (bool), db_dir (str), only_mash (bool), info_file (str), l90 (int),
    # nbcont (int), cutn (int), min_dist (float), max_dist (float)
    args_prepare = ("104099", "", "", "all", "refseq", "", False, "", False, "", 100, 999, 5, 1e-4, 0.06)
    # args for annotate:
    # name (str), qc_only (bool), date (str), prodigal_only (bool), small (bool)
    args_annot = ("TEST", True, "2101", False, False)
    # args for pangenome:
    # min_id (float), clust_mode (int), spe_dir (str), outfile (str)
    args_pan = (0.8, 1, "", "")
    # args for corepers
    # tol (float), mixed (bool), multi (bool), floor (bool)
    args_corepers = (1, False, False, False)
    #  args for tree module
    #  soft (str), model (str), boot (bool), write_boot (bool), memory (str), fast (bool)
    args_tree = ("iqtree2", "GTR", False, False, "", True)

    # Run 'all' module
    out = allm.main(cmd, args_all, args_prepare, args_annot, args_pan, args_corepers, args_tree)
    assert out == "QC_only done"
    # Check that there are 3 log files (log, err and details)
    log_files = glob.glob(os.path.join(outdir, "*log*"))
    assert len(log_files) == 3
    # Check that there are 5 files/folders : 3 logs + 2 result folders
    assert len(os.listdir(outdir)) == 5
    # Check result folder names
    prep_dir = os.path.join(outdir, "1-prepare_module")
    annot_dir = os.path.join(outdir, "2-annotate_module")
    assert os.path.isdir(prep_dir)
    assert os.path.isdir(annot_dir)
    # Check presence of some key files for prepare module
    assert(len(glob.glob(os.path.join(prep_dir, "*log*")))) == 3
    lst1 = os.path.join(prep_dir, "LSTINFO-104099-filtered-0.0001_0.06.txt")
    ass1 = os.path.join(prep_dir, "assembly_summary-104099.txt") 
    assert os.path.isfile(lst1)
    assert os.path.isfile(ass1)
    assert os.path.isdir(os.path.join(prep_dir, "refseq"))
    # Check presence of key files in prepare module
    assert(len(glob.glob(os.path.join(annot_dir, "*log*")))) == 3
    lst2 = os.path.join(annot_dir, "ALL-GENOMES-info-LSTINFO-104099-filtered-0.0001_0.06.lst")
    png1 = os.path.join(annot_dir, "QC_L90-LSTINFO-104099-filtered-0.0001_0.06.png")
    png2 = os.path.join(annot_dir, "QC_nb-contigs-LSTINFO-104099-filtered-0.0001_0.06.png")
    disc = os.path.join(annot_dir, "discarded-LSTINFO-104099-filtered-0.0001_0.06.lst")
    for f in (lst2, png1, png2, disc):
        assert os.path.isfile(f)


def test_main_norefseq():
    """
    Test with norefseq (4 genomes given, giving 13 families).
    """
    cmd = "cmd"
    outdir = os.path.join(GENEPATH, "main_all_default_qc")

    # Common arguments: outdir, threads, verbose, quiet
    args_all = (outdir, 1, 0, False)
    # args for prepare:
    # NCBI_species_taxid (int), NCBI_species (str), levels (str), tmp_dir (str),
    # norefseq (bool), db_dir (str), only_mash (bool), info_file (str), l90 (int),
    # nbcont (int), cutn (int), min_dist (float), max_dist (float)
    # db_dir = "test/data/pangenome/test_files/example_db/Replicons"
    # db_dir = "104099/Database_init"
    db_dir = os.path.join(DATADIR, "genomes")
    args_prepare = ("104099", "", "", "all", "refseq", "", True, db_dir, False, "", 100, 999, 5, 1e-4, 1)
    # args for annotate:
    # name (str), qc_only (bool), date (str), prodigal_only (bool), small (bool)
    args_annot = ("TEST", False, "2101", True, False)
    # args for pangenome:
    # min_id (float), clust_mode (int), spe_dir (str), outfile (str)
    args_pan = (0.8, 1, "", "")
    # args for args_corepers
    # tol (float), mixed (bool), multi (bool), floor (bool)
    args_corepers = (1, False, False, False)
    #  args for tree module
    #  soft (str), model (str), boot (bool), write_boot (bool), memory (str), fast (bool)
    args_tree = ("iqtree2", "GTR", False, False, "", True)

    # Run 'all' module
    out = allm.main(cmd, args_all, args_prepare, args_annot, args_pan, args_corepers, args_tree)
    assert out == 0
    # Check that there are 2 log files (log, err)
    log_files = glob.glob(os.path.join(outdir, "*log*"))
    assert len(log_files) == 2
    # Check that there are 5 files/folders : 2 logs + 4 result folders
    assert len(os.listdir(outdir)) == 8
    # Check result folder names
    prep_dir = os.path.join(outdir, "1-prepare_module")
    annot_dir = os.path.join(outdir, "2-annotate_module")
    pan_dir = os.path.join(outdir, "3-pangenome_module")
    core_dir = os.path.join(outdir, "4-corepers_module")
    ali_dir = os.path.join(outdir, "5-align_module")
    tree_dir = os.path.join(outdir, "6-tree_module")
    for d in (prep_dir, annot_dir, pan_dir, core_dir, ali_dir, tree_dir):
        assert os.path.isdir(d)
    # CHECK PREPARE
    # Check presence of some key files for prepare module
    assert(len(glob.glob(os.path.join(prep_dir, "*log*")))) == 3
    lst1 = os.path.join(prep_dir, "LSTINFO-104099-filtered-0.0001_1.txt")
    assert os.path.isfile(lst1)
    # CHECK ANNOTATE
    # Check presence of key files in prepare module
    assert(len(glob.glob(os.path.join(annot_dir, "*log*")))) == 3
    lst2 = os.path.join(annot_dir, "LSTINFO-LSTINFO-104099-filtered-0.0001_1.lst")
    png1 = os.path.join(annot_dir, "QC_L90-LSTINFO-104099-filtered-0.0001_1.png")
    png2 = os.path.join(annot_dir, "QC_nb-contigs-LSTINFO-104099-filtered-0.0001_1.png")
    for f in (lst2, png1, png2):
        assert os.path.isfile(f)
    # Check presence of proteins, genes, gff3, replicons and lstinfo folders, and 4 genomes in each
    prot_fold = os.path.join(annot_dir, "Proteins")
    gen_fold = os.path.join(annot_dir, "Genes")
    rep_fold = os.path.join(annot_dir, "Replicons")
    gff_fold = os.path.join(annot_dir, "gff3")
    lstinfo_fold = os.path.join(annot_dir, "LSTINFO")
    for d in (gen_fold, rep_fold, gff_fold, lstinfo_fold):
        assert os.path.isdir(d)
        files = os.listdir(d)
        assert len(files) == 4
    assert os.path.isdir(prot_fold)
    assert len(os.listdir(prot_fold)) == 5  # 4 genomes + concatenated DB
    # CHECK PAN
    assert(len(glob.glob(os.path.join(pan_dir, "PanACoTA*log*")))) == 3
    assert(len(glob.glob(os.path.join(pan_dir, "PanGenome-TEST*")))) == 5
    # CHECK COREPERS
    assert(len(glob.glob(os.path.join(core_dir, "PanACoTA*log*")))) == 2
    core_file = glob.glob(os.path.join(core_dir, "PersGenome*"))
    assert len(core_file) == 1
    # CHECK PAN
    assert(len(glob.glob(os.path.join(ali_dir, "PanACoTA*log*")))) == 3
    assert os.path.isdir(os.path.join(ali_dir, "Align-TEST_4"))
    assert os.path.isdir(os.path.join(ali_dir, "List-TEST_4"))
    assert os.path.isdir(os.path.join(ali_dir, "Phylo-TEST_4"))
    assert os.path.isfile(os.path.join(ali_dir, "Phylo-TEST_4", "TEST_4.grp.aln"))
    # CHECK TREE
    assert(len(glob.glob(os.path.join(tree_dir, "PanACoTA*log*")))) == 3
    assert os.path.isfile(os.path.join(tree_dir, "TEST_4.grp.aln.iqtree_tree.iqtree"))
    assert os.path.isfile(os.path.join(tree_dir, "TEST_4.grp.aln.iqtree_tree.treefile"))


def test_main_from_parse():
    """
    same test (from norefseq), but from 'main_from_parse' function, and with verbose == debug
    """
    args = argparse.Namespace()
    args.argv = ["all", "test_all_modules"]
    # common params
    outdir = os.path.join(GENEPATH, "from_parse")
    args.outdir = outdir
    args.threads = 1
    args.verbose = 15
    args.quiet = False
    # prepare params
    args.ncbi_species_name = ""
    args.ncbi_species_taxid = "104099"
    args.ncbi_taxid = ""
    args.levels = ""
    args.ncbi_section = "refseq"
    args.tmp_dir = ""
    args.norefseq = True
    args.db_dir = os.path.join(DATADIR, "genomes")
    args.only_mash = False
    args.info_file = ""
    args.l90 = 100
    args.nbcont = 999
    args.cutn = 5
    args.min_dist = 1e-4
    args.max_dist = 1
    # annotate params
    args.qc_only = False
    args.date = "2101"
    args.prodigal_only = False
    args.small = False
    args.name = "TEST"
    # pangenome params
    args.min_id = 0.8
    args.clust_mode = 1
    args.spedir = ""
    args.outfile = ""
    # corepers params
    args.tol = 1
    args.mixed = False
    args.multi = False
    args.floor = False
    # params tree
    args.soft = "iqtree2"
    args.model = "GTR"
    args.boot = False
    args.write_boot = False 
    args.memory = ""
    args.fast = True

    allm.main_from_parse(args)

    # Check that there are 2 log files (log, err)
    log_files = glob.glob(os.path.join(outdir, "*log*"))
    assert len(log_files) == 4
    # Check that there are 5 files/folders : 4 logs + 6 result folders
    assert len(os.listdir(outdir)) == 10
    # Check result folder names
    prep_dir = os.path.join(outdir, "1-prepare_module")
    annot_dir = os.path.join(outdir, "2-annotate_module")
    pan_dir = os.path.join(outdir, "3-pangenome_module")
    core_dir = os.path.join(outdir, "4-corepers_module")
    ali_dir = os.path.join(outdir, "5-align_module")
    tree_dir = os.path.join(outdir, "6-tree_module")
    for d in (prep_dir, annot_dir, pan_dir, core_dir, ali_dir, tree_dir):
        assert os.path.isdir(d)
    # CHECK PREPARE
    # Check presence of some key files for prepare module
    assert(len(glob.glob(os.path.join(prep_dir, "*log*")))) == 4
    lst1 = os.path.join(prep_dir, "LSTINFO-104099-filtered-0.0001_1.txt")
    assert os.path.isfile(lst1)
    # CHECK ANNOTATE
    # Check presence of key files in prepare module
    assert(len(glob.glob(os.path.join(annot_dir, "*log*")))) == 4
    lst2 = os.path.join(annot_dir, "LSTINFO-LSTINFO-104099-filtered-0.0001_1.lst")
    png1 = os.path.join(annot_dir, "QC_L90-LSTINFO-104099-filtered-0.0001_1.png")
    png2 = os.path.join(annot_dir, "QC_nb-contigs-LSTINFO-104099-filtered-0.0001_1.png")
    for f in (lst2, png1, png2):
        assert os.path.isfile(f)
    # Check presence of proteins, genes, gff3, replicons and lstinfo folders, and 4 genomes in each
    prot_fold = os.path.join(annot_dir, "Proteins")
    gen_fold = os.path.join(annot_dir, "Genes")
    rep_fold = os.path.join(annot_dir, "Replicons")
    gff_fold = os.path.join(annot_dir, "gff3")
    lstinfo_fold = os.path.join(annot_dir, "LSTINFO")
    for d in (gen_fold, rep_fold, gff_fold, lstinfo_fold):
        assert os.path.isdir(d)
        files = os.listdir(d)
        assert len(files) == 4
    assert os.path.isdir(prot_fold)
    assert len(os.listdir(prot_fold)) == 5  # 4 genomes + concatenated DB
    # CHECK PAN
    assert(len(glob.glob(os.path.join(pan_dir, "PanACoTA*log*")))) == 4
    assert(len(glob.glob(os.path.join(pan_dir, "PanGenome-TEST*")))) == 5
    # CHECK COREPERS
    assert(len(glob.glob(os.path.join(core_dir, "PanACoTA*log*")))) == 4
    core_file = glob.glob(os.path.join(core_dir, "PersGenome*"))
    assert len(core_file) == 1
    # CHECK PAN
    assert(len(glob.glob(os.path.join(ali_dir, "PanACoTA*log*")))) == 4
    assert os.path.isdir(os.path.join(ali_dir, "Align-TEST_4"))
    assert os.path.isdir(os.path.join(ali_dir, "List-TEST_4"))
    assert os.path.isdir(os.path.join(ali_dir, "Phylo-TEST_4"))
    assert os.path.isfile(os.path.join(ali_dir, "Phylo-TEST_4", "TEST_4.grp.aln"))
    # CHECK TREE
    assert(len(glob.glob(os.path.join(tree_dir, "PanACoTA*log*")))) == 4
    assert os.path.isfile(os.path.join(tree_dir, "TEST_4.grp.aln.iqtree_tree.iqtree"))
    assert os.path.isfile(os.path.join(tree_dir, "TEST_4.grp.aln.iqtree_tree.treefile"))