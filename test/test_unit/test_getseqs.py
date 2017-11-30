#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the get_seqs submodule in align module
"""
import os
import shutil

import pytest

import genomeAPCAT.align_module.get_seqs as gseq
from genomeAPCAT import utils

# Define common variables
FASTA = os.path.join("test", "data", "pangenome", "test_files", "example_db", "Proteins",
                     "GEN2.1017.00001.prt")
ALDIR = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALDIR, "exp_files")
TESTPATH = os.path.join(ALDIR, "test_files")
DBPATH = os.path.join("test", "data", "pangenome", "test_files", "example_db")
LOGFILE_BASE = "test_getseqs-logs"


# Setup and teardown: create logger
def setup_module():
    """
    create logger at start of this test module
    """
    utils.init_logger(LOGFILE_BASE, 0, '', verbose=1)


def teardown_module():
    """
    Remove log files at the end of this test module
    """
    os.remove(LOGFILE_BASE + ".log")
    os.remove(LOGFILE_BASE + ".log.details")
    os.remove(LOGFILE_BASE + ".log.err")


# Start tests
def test_extract_noseq_out_given():
    """
    Test that when giving an open fasta file, an empty list of sequences to extract, and
    a open file to which extracted sequences must be written, the output file is empty.
    """
    to_extract = {}
    outfile = "test_noextract_out-given.prt"
    with open(FASTA, "r") as fasf, open(outfile, "w") as outf:
        gseq.extract_sequences(to_extract, fasf, outf=outf)
    with open(outfile, "r") as outf:
        assert outf.readlines() == []
    os.remove(outfile)


def test_extract_seq_out_given():
    """
    Test that when giving an open fasta file, a list of 3 sequences to extract, and
    a open file to which extracted sequences must be written, it writes the
    expected sequences to output file.
    """
    to_extract = ["GEN2.1017.00001.b0001_00001", "GEN2.1017.00001.b0002_00006",
                  "GEN2.1017.00001.i0004_00012"]
    outfile = "test_extract_out-given.prt"
    with open(FASTA, "r") as fasf, open(outfile, "w") as outf:
        gseq.extract_sequences(to_extract, fasf, outf=outf)
    exp_extracted = os.path.join(EXPPATH, "exp_extracted.prt")
    same_files(outfile, exp_extracted)
    os.remove(outfile)


def test_extract_seq_out_given_wrongname():
    """
    Test that when giving an open fasta file, a list of 4 sequences to extract, and
    a open file to which extracted sequences must be written, it writes the
    expected sequences to output file. If a sequence does not exist, it is just not written.
    """
    to_extract = {"GEN2.1017.00001.b0001_00001": "", "GEN2.1017.00001.b0002_00006": "",
                  "GEN2.1017.00001.i0004_00012": "", "toto": ""}
    outfile = "test_extract_out-given.prt"
    with open(FASTA, "r") as fasf, open(outfile, "w") as outf:
        gseq.extract_sequences(to_extract, fasf, outf=outf)
    exp_extracted = os.path.join(EXPPATH, "exp_extracted.prt")
    same_files(outfile, exp_extracted)
    os.remove(outfile)


def test_extract_seq_out_allsame():
    """
    Test that when giving an open fasta file, a list of 3 sequences to extract with a
    corresponding output file (same for all 3 proteins), it writes the
    expected sequences to output file.
    """
    out = "test_extract1.prt"
    to_extract = {"GEN2.1017.00001.b0001_00001": out,
                  "GEN2.1017.00001.b0002_00006": out,
                  "GEN2.1017.00001.i0004_00012": out}
    with open(FASTA, "r") as fasf:
        gseq.extract_sequences(to_extract, fasf, files_todo=[out])
    exp_extracted = os.path.join(EXPPATH, "exp_extracted.prt")
    same_files(out, exp_extracted)
    os.remove(out)


def test_extract_seq_out_different():
    """
    Test that when giving an open fasta file, a list of 3 sequences to extract with a
    corresponding output file for each, it writes the expected sequences to the expected output
    file.
    """
    out1 = "test_extract1.prt"
    out2 = "test_extract2.prt"
    to_extract = {"GEN2.1017.00001.b0001_00001": out1,
                  "GEN2.1017.00001.b0002_00006": out2,
                  "GEN2.1017.00001.i0004_00012": out1}
    with open(FASTA, "r") as fasf:
        gseq.extract_sequences(to_extract, fasf, files_todo=[out1, out2])
    exp_extracted1 = os.path.join(EXPPATH, "exp_extracted1.prt")
    exp_extracted2 = os.path.join(EXPPATH, "exp_extracted2.prt")
    same_files(out1, exp_extracted1)
    same_files(out2, exp_extracted2)
    os.remove(out1)
    os.remove(out2)


def test_extract_seq_out_different_notasked():
    """
    Test that when giving an open fasta file, a list of 3 sequences to extract with a
    corresponding output file for each, it writes the
    expected sequences to expected output file, only if this file is contained in 'files_todo'.
    If not, file is not created, and protein not extracted.
    """
    out1 = "test_extract1.prt"
    out2 = "test_extract2.prt"
    to_extract = {"GEN2.1017.00001.b0001_00001": out1,
                  "GEN2.1017.00001.b0002_00006": out2,
                  "GEN2.1017.00001.i0004_00012": out1}
    with open(FASTA, "r") as fasf:
        gseq.extract_sequences(to_extract, fasf, files_todo=[out1])
    exp_extracted1 = os.path.join(EXPPATH, "exp_extracted1.prt")
    same_files(out1, exp_extracted1)
    assert not os.path.isfile(out2)
    os.remove(out1)


def test_get_names_files():
    """
    Test that given an open tab file (containing 2 columns: name of sequence to extract,
    file where it must be extracted), it returns the expected dict of sequences to extract
    """
    tabfile = os.path.join(EXPPATH, "exp_getentry-gen-ESCO2.txt")
    with open(tabfile, "r") as tabf:
        toext = gseq.get_names_to_extract(tabf, outfile=None)
    exp_toext = {'ESCO2_00001': "Aldir/TEST6-current.1.gen",
                 'ESCO2_22': 'Aldir/TEST6-current.2.gen',
                 'ESCO2_456': 'Aldir/TEST6-current.4.gen',
                 'ESCO2_46': 'Aldir/TEST6-current.3.gen'}
    assert toext == exp_toext


def test_get_names_wrongformat(caplog):
    """
    Test that given an open tab file, containing only the sequence name (no file associated) fir
    1 sequence, it exits with the expected error message
    """
    tabfile = os.path.join(TESTPATH, "wrong_getentry.prt")
    with open(tabfile, "r") as tabf:
        with pytest.raises(SystemExit):
            gseq.get_names_to_extract(tabf, outfile=None)
    assert ("Your file test/data/align/test_files/wrong_getentry.prt does not contain an output "
            "filename for ESCO2_22. Please give an output filename for each sequence to "
            "extract, or give a general output filename where all sequences will be "
            "extracted.") in caplog.text


def test_get_names_out_and2columns():
    """
    Test that given an open tab file (containing 2 columns: name of sequence to extract,
    file where it must be extracted), and an output file name, it returns the expected dict of
    sequences to extract (all to extract to the same given output file, the 2nd column of tab
    file is ignored).
    """
    tabfile = os.path.join(EXPPATH, "exp_getentry-gen-ESCO2.txt")
    outfile = "test_getnames_2columns"
    with open(tabfile, "r") as tabf:
        toext = gseq.get_names_to_extract(tabf, outfile=outfile)
    exp_toext = {'ESCO2_00001': outfile,
                 'ESCO2_22': outfile,
                 'ESCO2_456': outfile,
                 'ESCO2_46': outfile}
    assert toext == exp_toext


def test_get_names_out_and1column():
    """
    Test that given an open tab file (containing only sequence names), and an output file name,
    it returns the expected dict of sequences to extract (all to extract to the same given
    output file).
    """
    tabfile = os.path.join(TESTPATH, "getentry_1column.txt")
    outfile = "test_getnames_1column"
    with open(tabfile, "r") as tabf:
        toext = gseq.get_names_to_extract(tabf, outfile=outfile)
    exp_toext = {'ESCO2_00001': outfile,
                 'ESCO2_22': outfile,
                 'ESCO2_456': outfile,
                 'ESCO2_46': outfile}
    assert toext == exp_toext


def test_get_genome_all_seqs():
    """
    Test that given a fasta file, and a tab file containing all sequences to extract, with the
    files to which it must be extracted, it extracts everything in the right file.
    """
    tabfile = os.path.join(TESTPATH, "getentry_all_2columns.txt")
    todo = ["file1.txt", "file2.txt"]
    gseq.get_genome_seqs(FASTA, tabfile, todo)
    for i in range(1, 3):
        outfile = "file{}.txt".format(i)
        exp_file = os.path.join(EXPPATH, "exp_extracted{}.prt".format(i))
        assert os.path.isfile(outfile)
        same_files(outfile, exp_file)
        os.remove(outfile)


def test_get_genome_seqs_outgiven_2cols():
    """
    Test that given a fasta file, and a tab file containing all sequences to extract, with the
    files to which it must be extracted, and an output file, it extracts all sequences to the
    same output file, ignoring the ones given in tab file
    """
    tabfile = os.path.join(TESTPATH, "getentry_all_2columns.txt")
    outfile = "fileout.txt"
    todo = []
    gseq.get_genome_seqs(FASTA, tabfile, todo, outfile)
    assert os.path.isfile(outfile)
    exp_file = os.path.join(EXPPATH, "exp_extracted.prt")
    same_files(outfile, exp_file)
    os.remove(outfile)


def test_get_genome_seqs_outgiven_1col():
    """
    Test that given a fasta file, and a tab file containing only all sequences to extract,
    (no filename), and an output file, it extracts all sequences to the same output file.
    """
    tabfile = os.path.join(TESTPATH, "getentry_all_1column.txt")
    outfile = "fileout.txt"
    todo = []
    gseq.get_genome_seqs(FASTA, tabfile, todo, outfile)
    assert os.path.isfile(outfile)
    exp_file = os.path.join(EXPPATH, "exp_extracted.prt")
    same_files(outfile, exp_file)
    os.remove(outfile)
    assert not os.path.isfile("file1.txt")
    assert not os.path.isfile("file2.txt")


def test_get_genome_seqs_1notasked():
    """
    Test that given a fasta file, and a tab file containing all sequences to extract, with the
    files to which it must be extracted, and only 1 of them in 'files_todo', it extracts only
    the proteins going to this file.
    """
    tabfile = os.path.join(TESTPATH, "getentry_all_2columns.txt")
    todo = ["file1.txt"]
    gseq.get_genome_seqs(FASTA, tabfile, todo)
    outfile = "file1.txt"
    assert os.path.isfile(outfile)
    exp_file = os.path.join(EXPPATH, "exp_extracted1.prt")
    same_files(outfile, exp_file)
    os.remove(outfile)
    assert not os.path.isfile("file2.txt")


def test_get_genome_seqs_exists(caplog):
    """
    Test that when the output file given aleady exists, it does not overwrite it, but returns a
    warning message saying that it will be used for next step.
    """
    tabfile = os.path.join(TESTPATH, "getentry_all_2columns.txt")
    outfile = "fileout.txt"
    open(outfile, "w").close()
    todo = []
    gseq.get_genome_seqs(FASTA, tabfile, todo, outfile)
    assert ("Sequences are already extracted in fileout.txt. This will be used for next step. If "
            "you want to re-extract all sequences, use option -F (or --force)") in caplog.text
    with open(outfile, "r") as outf:
        assert outf.readlines() == []
    os.remove(outfile)


def test_check_extract_empty():
    """
    Test that when no file is present, it returns all gen and prt files, for all families
    """
    all_fams = [1, 6, 10, 11, 13, 14]
    aldir = "Aldir"
    dname = "TESTcheckExtract"
    files_todo = gseq.check_existing_extract(all_fams, aldir, dname)
    assert len(files_todo) == 12
    gens = [os.path.join(aldir, "{}-current.{}.gen").format(dname, num) for num in all_fams]
    prts = [os.path.join(aldir, "{}-current.{}.prt").format(dname, num) for num in all_fams]
    assert len(gens + prts) == 12
    assert set(gens + prts) == set(files_todo)


def test_check_extract_prtfile():
    """
    Test that when there is only the prt files for 2 families in aldir (no gen file), and mafft
    file, prt2nuc file for a family + concatenated file, it removes all these files,
    and returns all prt and gen files for all families
    """
    all_fams = [1, 6, 10, 11, 13, 14]
    aldir = "Aldir-notempty"
    dname = "TESTcheckExtract"
    # Create prt files for families 6 and 13
    # Create mafft and prt2nuc files for family 6
    # Create concatenated file
    os.makedirs(aldir)
    prt1 = os.path.join(aldir, "{}-current.6.prt".format(dname))
    prt2 = os.path.join(aldir, "{}-current.13.prt".format(dname))
    mafft = os.path.join(aldir, "{}-mafft-align.6.aln".format(dname))
    btr = os.path.join(aldir, "{}-mafft-prt2nuc.6.aln".format(dname))
    concat = os.path.join(aldir, "{}-complete.cat.aln".format(dname))
    for outf in [prt1, prt2, mafft, btr, concat]:
        open(outf, "w").close()
    # Run check extract
    files_todo = gseq.check_existing_extract(all_fams, aldir, dname)
    assert len(files_todo) == 12
    gens = [os.path.join(aldir, "{}-current.{}.gen").format(dname, num) for num in all_fams]
    prts = [os.path.join(aldir, "{}-current.{}.prt").format(dname, num) for num in all_fams]
    assert len(gens + prts) == len(files_todo)
    assert set(gens + prts) == set(files_todo)
    for outf in [prt1, prt2, mafft, btr, concat]:
        assert not os.path.isfile(outf)
    shutil.rmtree(aldir)


def test_check_extract_files1genome():
    """
    Test that when there are already prt, gen, mafft and prt2nuc files for a family, it keeps them,
    and returns all prt and gen files for all other families (except the one already existing).
    + if concatenation file exists, it is removed (as we will realign all families except 1).
    """
    all_fams = [1, 6, 10, 11, 13, 14]
    aldir = "Aldir-notempty"
    dname = "TESTcheckExtract"
    # Create prt, gen, mafft and prt2nuc files for family 6 + concatenated file
    os.makedirs(aldir)
    prt = os.path.join(aldir, "{}-current.6.prt".format(dname))
    gen = os.path.join(aldir, "{}-current.6.gen".format(dname))
    mafft = os.path.join(aldir, "{}-mafft-align.6.aln".format(dname))
    btr = os.path.join(aldir, "{}-mafft-prt2nuc.6.aln".format(dname))
    concat = os.path.join(aldir, "{}-complete.cat.aln".format(dname))
    for outf in [prt, gen, mafft, btr, concat]:
        open(outf, "w").close()
    # Run check extract
    files_todo = gseq.check_existing_extract(all_fams, aldir, dname)
    assert len(files_todo) == 10
    fams_todo = [1, 10, 11, 13, 14]
    gens = [os.path.join(aldir, "{}-current.{}.gen").format(dname, num) for num in fams_todo]
    prts = [os.path.join(aldir, "{}-current.{}.prt").format(dname, num) for num in fams_todo]
    assert len(gens + prts) == len(files_todo)
    assert set(gens + prts) == set(files_todo)
    for outf in [prt, gen, mafft, btr]:
        assert os.path.isfile(outf)
    assert not os.path.isfile(concat)
    shutil.rmtree(aldir)


def test_check_extract_allexist():
    """
    Test that when there are already all prt and gen files for all families,
    it keeps them, and returns an empty list of files to extract.
    If concatenation file exists, it is also kept.
    """
    all_fams = [1, 6, 10, 11, 13, 14]
    aldir = "Aldir-notempty"
    dname = "TESTcheckExtract"
    # Create prt and gen files for all families
    os.makedirs(aldir)
    gens = [os.path.join(aldir, "{}-current.{}.gen").format(dname, num) for num in all_fams]
    prts = [os.path.join(aldir, "{}-current.{}.prt").format(dname, num) for num in all_fams]
    mafft = os.path.join(aldir, "{}-mafft-align.6.aln".format(dname))
    btr = os.path.join(aldir, "{}-mafft-prt2nuc.6.aln".format(dname))
    concat = os.path.join(aldir, "{}-complete.cat.aln".format(dname))
    for outf in gens + prts + [mafft, btr, concat]:
        open(outf, "w").close()
    # Run check extract
    files_todo = gseq.check_existing_extract(all_fams, aldir, dname)
    assert files_todo == []
    for outf in gens + prts + [mafft, btr, concat]:
        assert os.path.isfile(outf)
    shutil.rmtree(aldir)


def test_get_all_seqs(caplog):
    """
    Test that when giving a list of family numbers, and output directories are empty,
    it extracts all expected proteins and genes.
    => Default. empty output, give database and 2 families to extract and getentry files
    exist in Listdir
    """
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    dname = "TESTgetAllSeq"
    listdir = "Listdir"
    aldir = "Align"
    all_fams = [1, 6]
    quiet = False
    # Create listdir and aldir and put all getentry files in listdir
    os.makedirs(listdir)
    os.makedirs(aldir)
    ref_listdir = os.path.join(TESTPATH, "test_listdir")
    ref_aldir = os.path.join(EXPPATH, "exp_aldir")
    for gen in all_genomes:
        genome_gen = os.path.join(ref_listdir, "getentry-gen_{}".format(gen))
        genome_prt = os.path.join(ref_listdir, "getentry-prt_{}".format(gen))
        gen_out = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        prt_out = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        shutil.copyfile(genome_gen, gen_out)
        shutil.copyfile(genome_prt, prt_out)
    gseq.get_all_seqs(all_genomes, dname, DBPATH, listdir, aldir, all_fams, quiet)
    # For each family, check that prt and gen files exist, and their content
    for fam in all_fams:
        fam_prt = os.path.join(aldir, "{}-current.{}.prt".format(dname, fam))
        assert os.path.isfile(fam_prt)
        exp_fam_prt = os.path.join(ref_aldir, "current.{}.prt".format(fam))
        same_files(fam_prt, exp_fam_prt)
        fam_gen = os.path.join(aldir, "{}-current.{}.gen".format(dname, fam))
        assert os.path.isfile(fam_gen)
        exp_fam_gen = os.path.join(ref_aldir, "current.{}.gen".format(fam))
        same_files(fam_gen, exp_fam_gen)
    # Check logs
    assert "Extracting proteins and genes from all genomes" in caplog.text
    for gen in all_genomes:
        assert "Extracting proteins and genes from {}".format(gen) in caplog.text
    shutil.rmtree(listdir)
    shutil.rmtree(aldir)


def test_get_all_seqs_prt6(caplog):
    """
    Test that when giving a list of family numbers, and output directories contain only a prt
    file for 1 family, it removes this prt file and it extracts all expected proteins and genes.
    => Aldir with prt file for fam 6. Others as default
    """
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    dname = "TESTgetAllSeq"
    listdir = "Listdir"
    aldir = "Align"
    all_fams = [1, 6]
    quiet = False
    # Create listdir and aldir and put all getentry files in listdir
    os.makedirs(listdir)
    os.makedirs(aldir)
    ref_listdir = os.path.join(TESTPATH, "test_listdir")
    ref_aldir = os.path.join(EXPPATH, "exp_aldir")
    prt6 = os.path.join(aldir, "{}-current.6.prt".format(dname))
    # Create empty file for prt of family 6
    open(prt6, "w").close()
    for gen in all_genomes:
        genome_gen = os.path.join(ref_listdir, "getentry-gen_{}".format(gen))
        genome_prt = os.path.join(ref_listdir, "getentry-prt_{}".format(gen))
        gen_out = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        prt_out = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        shutil.copyfile(genome_gen, gen_out)
        shutil.copyfile(genome_prt, prt_out)
    assert os.path.isfile(prt6)
    gseq.get_all_seqs(all_genomes, dname, DBPATH, listdir, aldir, all_fams, quiet)
    # For each family, check that prt and gen files exist, and their content
    for fam in all_fams:
        fam_prt = os.path.join(aldir, "{}-current.{}.prt".format(dname, fam))
        assert os.path.isfile(fam_prt)
        exp_fam_prt = os.path.join(ref_aldir, "current.{}.prt".format(fam))
        same_files(fam_prt, exp_fam_prt)
        fam_gen = os.path.join(aldir, "{}-current.{}.gen".format(dname, fam))
        assert os.path.isfile(fam_gen)
        exp_fam_gen = os.path.join(ref_aldir, "current.{}.gen".format(fam))
        same_files(fam_gen, exp_fam_gen)
    # Check logs
    assert "Extracting proteins and genes from all genomes" in caplog.text
    for gen in all_genomes:
        assert "Extracting proteins and genes from {}".format(gen) in caplog.text
    shutil.rmtree(listdir)
    shutil.rmtree(aldir)


def test_get_all_seqs_prtgen6(caplog):
    """
    Test that when giving a list of family numbers, and output directories contain a prt and a gen
    file for 1 family, it extracts all expected proteins and genes for other families, but keeps
    the current file for family already having prt and gen
    + add mafft and prt2nuc files for this family, and check that they are not removed
    + add concatenate file, and check that it is removed
    => prt and gen files in Aldir for fam 6. Others as default
    """
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    dname = "TESTgetAllSeq"
    listdir = "Listdir"
    aldir = "Align"
    all_fams = [1, 6]
    quiet = False
    # Create listdir and aldir and put all getentry files in listdir
    os.makedirs(listdir)
    os.makedirs(aldir)
    ref_listdir = os.path.join(TESTPATH, "test_listdir")
    ref_aldir = os.path.join(EXPPATH, "exp_aldir")
    # Create empty files for prt, gen, mafft and prt2nuc files of family 6
    prt6 = os.path.join(aldir, "{}-current.6.prt".format(dname))
    gen6 = os.path.join(aldir, "{}-current.6.gen".format(dname))
    mafft6 = os.path.join(aldir, "{}-mafft-align.6.aln".format(dname))
    prt2nuc6 = os.path.join(aldir, "{}-mafft-prt2nuc.6.aln".format(dname))
    # Add concatenate file
    concat = os.path.join(aldir, "{}-complete.cat.aln".format(dname))
    for outf in [prt6, gen6, mafft6, prt2nuc6, concat]:
        open(outf, "w").close()
    for gen in all_genomes:
        genome_gen = os.path.join(ref_listdir, "getentry-gen_{}".format(gen))
        genome_prt = os.path.join(ref_listdir, "getentry-prt_{}".format(gen))
        gen_out = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        prt_out = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        shutil.copyfile(genome_gen, gen_out)
        shutil.copyfile(genome_prt, prt_out)
    gseq.get_all_seqs(all_genomes, dname, DBPATH, listdir, aldir, all_fams, quiet)
    # For family 1, check that prt and gen files exist and are as expected
    fam_prt = os.path.join(aldir, "{}-current.1.prt".format(dname))
    assert os.path.isfile(fam_prt)
    exp_fam_prt = os.path.join(ref_aldir, "current.1.prt")
    same_files(fam_prt, exp_fam_prt)
    fam_gen = os.path.join(aldir, "{}-current.1.gen".format(dname))
    assert os.path.isfile(fam_gen)
    exp_fam_gen = os.path.join(ref_aldir, "current.1.gen")
    same_files(fam_gen, exp_fam_gen)
    # For family 6 , check that all filesare present and empty
    for outf in [prt6, gen6, mafft6, prt2nuc6]:
        assert os.path.isfile(outf)
        with open(outf, "r") as out:
            assert out.readlines() == []
    # Check that concat file was removed
    assert not os.path.isfile(concat)
    # Check logs
    assert "Extracting proteins and genes from all genomes" in caplog.text
    for gen in all_genomes:
        assert "Extracting proteins and genes from {}".format(gen) in caplog.text
    shutil.rmtree(listdir)
    shutil.rmtree(aldir)


def test_get_all_seqs_allexist(caplog):
    """
    Test that when giving a list of family numbers, and all prt and gen extraction files for this
    family already exist in Aldir, it does not change them, and returns an warning message saying
    that all prt and gen are already extracted.
    + concat file present, and not removed
    """
    all_genomes = ["GEN2.1017.00001", "GEN4.1111.00001", "GENO.1017.00001", "GENO.1216.00002"]
    dname = "TESTgetAllSeq"
    listdir = "Listdir"
    aldir = "Align"
    all_fams = [1, 6]
    quiet = False
    # Create listdir and aldir and put all getentry files in listdir
    os.makedirs(listdir)
    os.makedirs(aldir)
    ref_listdir = os.path.join(TESTPATH, "test_listdir")
    ref_aldir = os.path.join(EXPPATH, "exp_aldir")
    # Create empty files for prt, gen, mafft and prt2nuc files of family 6
    created_files = []
    for fam in all_fams:
        prt = os.path.join(aldir, "{}-current.{}.prt".format(dname, fam))
        gen = os.path.join(aldir, "{}-current.{}.gen".format(dname, fam))
        mafft = os.path.join(aldir, "{}-mafft-align.{}.aln".format(dname, fam))
        prt2nuc = os.path.join(aldir, "{}-mafft-prt2nuc.{}.aln".format(dname, fam))
        created_files += [prt, gen, mafft, prt2nuc]
        for outf in [prt, gen, mafft, prt2nuc]:
            open(outf, "w").close()
    # Add concatenate file
    concat = os.path.join(aldir, "{}-complete.cat.aln".format(dname))
    open(concat, "w").close()
    for gen in all_genomes:
        genome_gen = os.path.join(ref_listdir, "getentry-gen_{}".format(gen))
        genome_prt = os.path.join(ref_listdir, "getentry-prt_{}".format(gen))
        gen_out = os.path.join(listdir, "{}-getEntry_gen_{}.txt".format(dname, gen))
        prt_out = os.path.join(listdir, "{}-getEntry_prt_{}.txt".format(dname, gen))
        shutil.copyfile(genome_gen, gen_out)
        shutil.copyfile(genome_prt, prt_out)
    gseq.get_all_seqs(all_genomes, dname, DBPATH, listdir, aldir, all_fams, quiet)
    # check that all files are present and empty
    for outf in created_files:
        assert os.path.isfile(outf)
        with open(outf, "r") as out:
            assert out.readlines() == []
    # Check that concat file was removed
    assert os.path.isfile(concat)
    # Check logs
    assert ("All extraction files already existing (see detailed log for " 
            "more information)") in caplog.text
    assert ("All prt and gene files for all families already exist. The program " 
            "will use them for the next step. If you want to re-extract a given " 
            "family, remove its prt and gen extraction files. If you want to " 
            "re-extract all families, use option -F (or --force).") in caplog.text
    shutil.rmtree(listdir)
    shutil.rmtree(aldir)


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
