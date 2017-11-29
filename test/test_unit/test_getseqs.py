#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the get_seqs submodule in align module
"""
import os

import pytest

import genomeAPCAT.align_module.get_seqs as gseq

# Define common variables
FASTA = os.path.join("test", "data", "pangenome", "test_files", "example_db", "Proteins",
                     "GEN2.1017.00001.prt")
ALDIR = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALDIR, "exp_files")
TESTPATH = os.path.join(ALDIR, "test_files")


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
