#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the post_align submodule in align module
"""
import os
import pytest
import shutil
import logging

import PanACoTA.align_module.post_align as pal
import test.test_unit.utilities_for_tests as tutil
from PanACoTA import utils

# Define common variables
ALPATH = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALPATH, "exp_files")
TESTPATH = os.path.join(ALPATH, "test_files")
GENEPATH = os.path.join(ALPATH, "generated_by_unit-tests")
LOGFILE_BASE = "logs_test_postalign"
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]

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
    utils.init_logger(LOGFILE_BASE, 0, 'test_postalign', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    for f in LOGFILES:
        if os.path.exists(f):
            os.remove(f)
    print("teardown")


def test_get_genome():
    """
    Given a header and a list of genomes, check that it returns the expected genome
    """
    header = ">TOTO.0215.00002.i006_00065"
    genomes = ["TOTO.0315.00001", "ESCO.0215.00002", "ESCO.0215.00001", "TOTO.0215.00002"]
    assert pal.get_genome(header, genomes) == "TOTO.0215.00002"


def test_get_genome_not_start():
    """
    Given a header and a list of genomes, check that it returns the expected genome. The genome
    name is not at the beginning of the protein name
    """
    header = ">mongenome,TOTO.0215.00002.i006_00065"
    genomes = ["TOTO.0315.00001", "ESCO.0215.00002", "ESCO.0215.00001", "TOTO.0215.00002"]
    assert pal.get_genome(header, genomes) == "TOTO.0215.00002"


def test_get_genome_notfound(caplog):
    """
    Test that when header given, but no genome of the list corresponds to it, it exits with an
    error message.
    """
    header = ">TOTO.0215.00002.i006_00065"
    genomes = ["TOTO.0315.00001", "ESCO.0215.00002", "ESCO.0215.00001"]
    assert pal.get_genome(header, genomes) is None
    assert ("Protein TOTO.0215.00002.i006_00065 does not correspond to any genome name given... "
            "['TOTO.0315.00001', 'ESCO.0215.00002', 'ESCO.0215.00001']") in caplog.text


def test_write_groups(caplog):
    """
    Check that giving a list of genomes and corresponding sequences, it writes all of them in
    fasta format in the given output file.
    """
    caplog.set_level(logging.DEBUG)
    outfile = os.path.join(GENEPATH, "test_write_groups.fa")
    sequences = {"ESCO.0216.00001": ["AAAAA-", "TTTTT", "CCCCC", "GGGGG"],
                 "ESCO.0216.00002": ["AAAAT-", "TTTTA", "-CCCG", "GGGGC"],
                 "ESCO.0216.00003": ["AAAAAA", "TT-TT", "CCC-C", "GCGG-"]}
    pal.write_groups(outfile, sequences)
    assert "Writing alignments per genome" in caplog.text
    with open(outfile, "r") as outf:
        lines = outf.readlines()
    assert len(lines) == 6
    exp_lines = [">ESCO.0216.00001\n", "AAAAA-TTTTTCCCCCGGGGG\n",
                 ">ESCO.0216.00002\n", "AAAAT-TTTTA-CCCGGGGGC\n",
                 ">ESCO.0216.00003\n", "AAAAAATT-TTCCC-CGCGG-\n"]
    assert lines == exp_lines


def test_read_alignment(caplog):
    """
    Giving a file with proteins aligned, and a list of genomes, returns, for each genome,
    the list of sequences found.
    """
    caplog.set_level(logging.DEBUG)
    alnfile = os.path.join(TESTPATH, "complete.cat.fictive4genomes.aln")
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    seqs = pal.read_alignments(alnfile, all_genomes)
    exp_seq = {"GEN2.1017.00001": ["AAATCTCGGAGATCTCGCGCGGATATATCGCGCTATAGCGCT",
                                   "---TTCGCGCGGAGAGGCGCTCTCGGCG",
                                   "----------------------------------------------------"],
               "GEN4.1111.00001": ["AAATCTCGGAGTTCTCGCGCGGGTA---CGCGCTATAGCGCT",
                                   "---TTCGCGCGGAGAGGCGCTCTCG---",
                                   "CGGTATAGGGCGCTAC----CCGGATATAGGGCTCTCGGAATTCGCAAC-CC"],
               "GENO.1017.00001": ["-------GGAGATCTCGCGCGGATATATCGCGCTATAGCGCT",
                                   "AAATTCGCGCGGAGAGGCGCTCTCGAAA",
                                   "CGGTATAGGGCGCTACTTTTCCGGATATAGGGCTCTCGGAATTCG--ACACC"],
               "GENO.1216.00002": ["AAATCTCGGAGATCTCGCGCGGATATATAAAACTATAGCGCT",
                                   "----------------------------",
                                   "CGGTATGGGGCGCTAC----CCGGATATAGGGCTCTCGGAATTCGCAAC-CG"]
               }
    assert seqs == exp_seq
    assert "3 sequences found per genome" in caplog.text


def test_read_alignment_diffnbseq(caplog):
    """
    Giving a file with proteins aligned, and a list of genomes, with 3 proteins for each genome,
    except for 1 genome which has 4 proteins. Exits with an error message.
    """
    caplog.set_level(logging.DEBUG)
    alnfile = os.path.join(TESTPATH, "complete.cat.fictive4genomes-diffnbseq.aln")
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    assert pal.read_alignments(alnfile, all_genomes) is None
    assert ("Problems occurred while grouping alignments by genome: all genomes do not have the "
            "same number of sequences. Check that each protein name contains the name of the "
            "genome from which it comes.") in caplog.text


def test_read_alignment_nogenome(caplog):
    """
    Giving a file with proteins aligned, and a list of genomes, but with one genome missing,
    check that it returns None, and an error message specifying the protein not having a
    corresponding genome.
    """
    caplog.set_level(logging.DEBUG)
    alnfile = os.path.join(TESTPATH, "complete.cat.fictive4genomes-diffnbseq.aln")
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001"]
    assert pal.read_alignments(alnfile, all_genomes) is None
    assert ("Protein GENO.1216.00002.i0001_00003 does not correspond to any genome name given... "
            "['GEN2.1017.00001', 'GEN4.1111.00001', 'GENO.1017.00001']") in caplog.text


def test_group_by_genome(caplog):
    """
    Test that giving a file with all proteins aligned, a list of genomes, and an output
    filename, it writes in output the alignment grouped by genome and returns True
    """
    caplog.set_level(logging.DEBUG)
    alnfile = os.path.join(TESTPATH, "complete.cat.fictive4genomes.aln")
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    outgrp = os.path.join(GENEPATH, "test_group_by_genome")
    args = (all_genomes, alnfile, outgrp)
    assert pal.group_by_genome(args)
    exp_grp = os.path.join(EXPPATH, "exp_fictive.grp.aln")
    assert tutil.compare_order_content(outgrp, exp_grp)
    assert "3 sequences found per genome" in caplog.text
    assert "Writing alignments per genome" in caplog.text


def test_group_by_genome_nogenome(caplog):
    """
    Test that giving a file with all proteins aligned, a list of genomes but with 1 genome missing,
    it returns False, and an error message with the protein not having a corresponding genome
    """
    caplog.set_level(logging.DEBUG)
    alnfile = os.path.join(TESTPATH, "complete.cat.fictive4genomes.aln")
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001"]
    outgrp = os.path.join(GENEPATH, "test_group_by_genome")
    args = (all_genomes, alnfile, outgrp)
    assert pal.group_by_genome(args) is False
    assert not os.path.isfile(outgrp)
    assert ("Protein GENO.1216.00002.i0001_00003 does not correspond to any genome name given... "
            "['GEN2.1017.00001', 'GEN4.1111.00001', 'GENO.1017.00001']") in caplog.text


def test_launch_gbg(caplog):
    """
    Giving an alignment file, and a list of genomes, and a tree directory, check that it
    generates the expected filename, with the expected content, and returns True.
    """
    caplog.set_level(logging.DEBUG)
    all_genomes = ["GENO.1216.00002", "GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001"]
    alnfile = os.path.join(TESTPATH, "complete.cat.pers4genomes.aln")
    status = "Done"
    treedir = os.path.join(GENEPATH, "test_phylo-pers4genomes")
    # Create tree dir, with empty output file
    os.makedirs(treedir)
    out_grp = os.path.join(treedir, "TESTlaunch.grp.aln")
    open(out_grp, "w").close()
    dname = "TESTlaunch"
    quiet = False
    assert pal.launch_group_by_genome(all_genomes, alnfile, status, out_grp, dname, quiet) is True
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    assert tutil.compare_order_content(out_grp, exp_grp)
    assert "Grouping alignments per genome" in caplog.text


def test_launch_gbg_existsempty(caplog):
    """
    Giving an alignment file, and a list of genomes, and saying that we did not re-create this
    concatenated file. As the group by genome file already exists (even if empty...), it won't
    recreate it, just return None and keep it.
    """
    caplog.set_level(logging.DEBUG)
    all_genomes = ["GENO.1216.00002", "GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001"]
    alnfile = os.path.join(TESTPATH, "complete.cat.pers4genomes.aln")
    status = "OK"
    treedir = os.path.join(GENEPATH, "test_phylo-pers4genomes")
    # Create tree dir, with empty output file
    os.makedirs(treedir)
    out_grp = os.path.join(treedir, "TESTlaunch.grp.aln")
    open(out_grp, "w").close()
    dname = "TESTlaunch"
    quiet = False
    assert pal.launch_group_by_genome(all_genomes, alnfile, status, out_grp, dname, quiet) is True
    assert "Grouping alignments per genome" not in caplog.text
    with open(out_grp, "r") as outf:
        assert outf.readlines() == []


def test_launch_gbg_ok_notexist(caplog):
    """
    Giving an alignment file, and a list of genomes, and a tree directory, and saying that we
    did not re-create this concatenated file. As group by genome file does not exist,
    it will create it as if concatenation was just done, and return True
    """
    caplog.set_level(logging.DEBUG)
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    alnfile = os.path.join(TESTPATH, "complete.cat.pers4genomes.aln")
    status = "OK"
    treedir = os.path.join(GENEPATH, "test_phylo-pers4genomes")
    # Create tree dir
    os.makedirs(treedir)
    out_grp = os.path.join(treedir, "TESTlaunch.grp.aln")
    dname = "TESTlaunch"
    quiet = False
    assert pal.launch_group_by_genome(all_genomes, alnfile, status, out_grp, dname, quiet) is True
    exp_grp = os.path.join(EXPPATH, "exp_pers4genomes.grp.aln")
    assert tutil.compare_order_content(out_grp, exp_grp)
    assert "Grouping alignments per genome" in caplog.text


def test_launch_gbg_nogenome(caplog):
    """
    Giving an alignment file, and a list of genomes with one genome missing, check that it
    returns False, with an error message for the protein not having a corresponding genome,
    and does not create grp file
    """
    caplog.set_level(logging.DEBUG)
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001"]
    alnfile = os.path.join(TESTPATH, "complete.cat.pers4genomes.aln")
    status = "OK"
    treedir = os.path.join(GENEPATH, "test_phylo-pers4genomes")
    # Create tree dir
    os.makedirs(treedir)
    out_grp = os.path.join(treedir, "TESTlaunch.grp.aln")
    dname = "TESTlaunch"
    quiet = False
    assert pal.launch_group_by_genome(all_genomes, alnfile, status, treedir, dname, quiet) is False
    assert not os.path.isfile(out_grp)
    assert "Grouping alignments per genome" in caplog.text


def test_concat(caplog):
    """
    Given a list of families, and a directory where are alignment files, check that the files
    corresponding to the given families are concatenated as expected, and it returns "Done" and
    expected output filename
    """
    caplog.set_level(logging.DEBUG)
    # Prepare aldir with all needed alignment files
    aldir = os.path.join(GENEPATH, "test_concat_aldir")
    dname = "TEST-concat"
    prefix = os.path.join(aldir, dname)
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    orig_btr11 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.11.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    btr11 = os.path.join(aldir, dname + "-mafft-prt2nuc.11.aln")
    os.makedirs(aldir)
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    shutil.copyfile(orig_btr11, btr11)
    # Other parameters, and run concatenation
    fam_nums = [1, 8, 11]
    quiet = False
    output, mess = pal.concat_alignments(fam_nums, prefix, quiet)
    assert output == os.path.join(aldir, dname + "-complete.cat.aln")
    ref_concat = os.path.join(EXPPATH, "exp_concat_4genomes-fam1-8-11.aln")
    assert tutil.compare_order_content(output, ref_concat)
    assert mess == "Done"
    assert "Concatenating all alignment files" in caplog.text


def test_concat_quiet(caplog):
    """
    Given a list of families, and a directory where are alignment files, check that the files
    corresponding to the given families are concatenated as expected, and it returns "Done" and
    expected output filename
    """
    caplog.set_level(logging.DEBUG)
    # Prepare aldir with all needed alignment files
    aldir = os.path.join(GENEPATH, "test_concat_aldir")
    dname = "TESTconcat"
    prefix = os.path.join(aldir, dname)
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    orig_btr11 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.11.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    btr11 = os.path.join(aldir, dname + "-mafft-prt2nuc.11.aln")
    os.makedirs(aldir)
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    shutil.copyfile(orig_btr11, btr11)
    # Other parameters, and run concatenation
    fam_nums = [1, 8, 11]
    quiet = True
    output, mess = pal.concat_alignments(fam_nums, prefix, quiet)
    assert output == os.path.join(aldir, dname + "-complete.cat.aln")
    ref_concat = os.path.join(EXPPATH, "exp_concat_4genomes-fam1-8-11.aln")
    assert tutil.compare_order_content(output, ref_concat)
    assert mess == "Done"
    assert "Concatenating all alignment files" in caplog.text


def test_concat_noalignfile(caplog):
    """
    Given a list of families, and a directory where are alignment files except for 1 family,
    for which we do not have the alignment file
    """
    caplog.set_level(logging.DEBUG)
    # Prepare aldir with all needed alignment files
    aldir = os.path.join(GENEPATH, "test_concat_aldir")
    dname = "TESTconcat"
    prefix = os.path.join(aldir, dname)
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    orig_btr11 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.11.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    btr11 = os.path.join(aldir, dname + "-mafft-prt2nuc.11.aln")
    os.makedirs(aldir)
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    shutil.copyfile(orig_btr11, btr11)
    # Other parameters, and run concatenation
    fam_nums = [1, 8, 11, 15]
    quiet = False
    with pytest.raises(SystemExit):
        pal.concat_alignments(fam_nums, prefix, quiet)
    assert ("The alignment file "
            "test/data/align/generated_by_unit-tests/test_concat_aldir/"
            "TESTconcat-mafft-prt2nuc.15.aln "
            "does not exist. Please check the families you want, and their corresponding "
            "alignment files") in caplog.text


def test_concat_outexists(caplog):
    """
    Given a list of families, and a directory where are alignment files and the concatenated
    file, check that it keeps this outfile, and returns "OK" with a warning message
    """
    caplog.set_level(logging.DEBUG)
    # Prepare aldir with all needed alignment files
    aldir = os.path.join(GENEPATH, "test_concat_aldir")
    dname = "TESTconcat"
    prefix = os.path.join(aldir, dname)
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    orig_btr11 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.11.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    btr11 = os.path.join(aldir, dname + "-mafft-prt2nuc.11.aln")
    os.makedirs(aldir)
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    shutil.copyfile(orig_btr11, btr11)
    # Create empty concatenated file
    outempty = os.path.join(aldir, dname + "-complete.cat.aln")
    open(outempty, "w").close()
    # Other parameters, and run concatenation
    fam_nums = [1, 8, 11]
    quiet = False
    output, mess = pal.concat_alignments(fam_nums, prefix, quiet)
    assert output == outempty
    assert mess == "OK"
    with open(output, "r") as outf:
        assert outf.readlines() == []
    assert "Alignments already concatenated" in caplog.text
    assert ("Alignments already concatenated in "
            "test/data/align/generated_by_unit-tests/test_concat_aldir/"
            "TESTconcat-complete.cat.aln. "
            "Program will use it for next steps") in caplog.text


def test_postalign(caplog):
    """
    Test that when running post-alignment on a folder containing all expected alignment files,
    it creates concatenated alignments, and a folder Phylo with the alignments grouped by genome.
    """
    caplog.set_level(logging.DEBUG)
    # define parameters
    fam_nums = [1, 8, 11]
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    outdir = os.path.join(GENEPATH, "test_post-align")
    aldir = os.path.join(outdir, "aldir_post-align")
    os.makedirs(aldir)
    dname = "TESTpost"
    prefix = os.path.join(aldir, dname)
    quiet = False
    # Prepare aldir with all needed alignment files
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    orig_btr11 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.11.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    btr11 = os.path.join(aldir, dname + "-mafft-prt2nuc.11.aln")
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    shutil.copyfile(orig_btr11, btr11)
    # Run post-alignment
    pal.post_alignment(fam_nums, all_genomes, prefix, outdir, dname, quiet)
    # Check that concatenated file is created and with expected content
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    assert os.path.isfile(out_concat)
    ref_concat = os.path.join(EXPPATH, "exp_concat_4genomes-fam1-8-11.aln")
    assert tutil.compare_order_content(out_concat, ref_concat)
    # Check that grouped by genome file is created, with expected content
    treedir = os.path.join(outdir, "Phylo-" + dname)
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    assert os.path.isfile(out_grp)
    exp_grp = os.path.join(EXPPATH, "exp_grp_4genomes-fam1-8-11.aln")
    assert tutil.compare_order_content(out_grp, exp_grp)
    # check logs
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text


def test_postalign_missalign(caplog):
    """
    Test that when running post-alignment on a folder containing all expected alignment files
    except 1, it exits with an error message indicating the missing alignment file.
    """
    caplog.set_level(logging.DEBUG)
    # define parameters
    fam_nums = [1, 8, 11]
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    outdir = os.path.join(GENEPATH, "test_post-align_missalign")
    aldir = os.path.join(outdir, "aldir_post-align")
    os.makedirs(aldir)
    dname = "TESTpost"
    prefix = os.path.join(aldir, dname)
    quiet = False
    # Prepare aldir with all needed alignment files
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    # Run post-alignment
    with pytest.raises(SystemExit):
        pal.post_alignment(fam_nums, all_genomes, prefix, outdir, dname, quiet)
    assert ("The alignment file test/data/align/generated_by_unit-tests/"
            "test_post-align_missalign/aldir_post-align/TESTpost-mafft-prt2nuc.11.aln "
            "does not exist. Please check the families you want, and their corresponding "
            "alignment files") in caplog.text
    # Check that concatenated file is not created
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    assert not os.path.isfile(out_concat)


def test_postalign_missgenome(caplog):
    """
    Test that when running post-alignment on a folder containing all expected alignment files,
    but giving incomplete list of genomes, it exits with error message specifying protein which
    does not belong to any given genome name.
    """
    caplog.set_level(logging.DEBUG)
    # define parameters
    fam_nums = [1, 8, 11]
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001"]
    outdir = os.path.join(GENEPATH, "test_post-align")
    aldir = os.path.join(outdir, "aldir_post-align")
    os.makedirs(aldir)
    dname = "TESTpost"
    prefix = os.path.join(aldir, dname)
    quiet = False
    # Prepare aldir with all needed alignment files
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    orig_btr11 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.11.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    btr11 = os.path.join(aldir, dname + "-mafft-prt2nuc.11.aln")
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    shutil.copyfile(orig_btr11, btr11)
    # Run post-alignment
    pal.post_alignment(fam_nums, all_genomes, prefix, outdir, dname, quiet)
    # Check that concatenated file is created and with expected content
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    assert os.path.isfile(out_concat)
    ref_concat = os.path.join(EXPPATH, "exp_concat_4genomes-fam1-8-11.aln")
    assert tutil.compare_order_content(out_concat, ref_concat)
    # Check that grouped by genome file is not created
    treedir = os.path.join(outdir, "Phylo-" + dname)
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    assert not os.path.isfile(out_grp)
    # check logs
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text
    assert "An error occurred. We could not group sequences by genome." in caplog.text


def test_postalign_diffnbseq(caplog):
    """
    Test that when running post-alignment on a folder containing all expected alignment files,
    but giving incomplete list of genomes, it exits with error message specifying protein which
    does not belong to any given genome name.
    """
    caplog.set_level(logging.DEBUG)
    # define parameters
    fam_nums = [1, 8, 11]
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001"]
    outdir = os.path.join(GENEPATH, "test_post-align")
    aldir = os.path.join(outdir, "aldir_post-align")
    os.makedirs(aldir)
    dname = "TESTpost"
    prefix = os.path.join(aldir, dname)
    quiet = False
    # Prepare aldir with all needed alignment files
    orig_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    orig_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    orig_btr11 = os.path.join(TESTPATH, "error-mafft-prt2nuc.11.sup-seq.aln")
    btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    btr11 = os.path.join(aldir, dname + "-mafft-prt2nuc.11.aln")
    shutil.copyfile(orig_btr1, btr1)
    shutil.copyfile(orig_btr8, btr8)
    shutil.copyfile(orig_btr11, btr11)
    # Run post-alignment
    pal.post_alignment(fam_nums, all_genomes, prefix, outdir, dname, quiet)
    # Check that concatenated file is created and with expected content
    out_concat = os.path.join(aldir, dname + "-complete.cat.aln")
    assert os.path.isfile(out_concat)
    ref_concat = os.path.join(TESTPATH, "error-concat-sup-seq.aln")
    assert tutil.compare_order_content(out_concat, ref_concat)
    # Check that grouped by genome file is not created
    treedir = os.path.join(outdir, "Phylo-" + dname)
    out_grp = os.path.join(treedir, dname + ".grp.aln")
    assert not os.path.isfile(out_grp)
    # check logs
    assert "Concatenating all alignment files" in caplog.text
    assert "Grouping alignments per genome" in caplog.text
    assert "An error occurred. We could not group sequences by genome." in caplog.text
