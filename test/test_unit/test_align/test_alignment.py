#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the alignment submodule in align module
"""
import logging
import os

import shutil

import genomeAPCAT.align_module.alignment as al


# Define common variables
ALDIR = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALDIR, "exp_files")
TESTPATH = os.path.join(ALDIR, "test_files")


def test_check_len():
    """
    Test that when giving an alignment file with all sequences having the same length,
    it returns this length, and the number of sequences aligned
    """
    aln_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    logger = logging.getLogger("test_alignment_check_len")
    fam = 1
    res = al.check_lens(aln_file, fam, logger)
    assert res[0] == 157
    assert res[1] == 4


def test_check_len_diff(caplog):
    """
    Test that when giving an alignment file with sequences not having all the same length,
    it returns False, and an error message
    """
    aln_file = os.path.join(TESTPATH, "test_alignment-mafft-diff-len.aln")
    logger = logging.getLogger("test_alignment_check_len")
    fam = 1
    assert not al.check_lens(aln_file, fam, logger)
    assert ("alignments for family 1 do not all have the same length. "
            "Lengths found are: {") in caplog.text
    assert "126" in caplog.text
    assert "154" in caplog.text
    assert "157" in caplog.text
    assert "}\n" in caplog.text


def test_check_btr():
    """
    Test that when giving a nucleotide alignment, and a number of sequences corresponding to the
    number of sequences in the alignment file, it returns True.
    """
    btr_file = os.path.join(TESTPATH, "test_alignment-btr-ok.aln")
    nbfal = 4
    logger = logging.getLogger("test_check_btr")
    message = "problem!"
    assert al.check_nb_seqs(btr_file, nbfal, logger, message) == 4


def test_check_btr_wrongnbfam(caplog):
    """
    Test that when giving a nucleotide alignment, and a number of sequences different from the
    number of sequences in the alignment file, it returns False, and an error message.
    """
    btr_file = os.path.join(TESTPATH, "test_alignment-btr-ok.aln")
    nbfal = 5
    logger = logging.getLogger("test_check_btr")
    message = "Problem with fam 1"
    assert not al.check_nb_seqs(btr_file, nbfal, logger, message)
    assert "Problem with fam 1 (4)" in caplog.text


def test_backtranslate(caplog):
    """
    Test that when giving an alignment file and a file with nucleic sequences, it generates a
    btr file, with the nucleotide alignments as expected, and returns True
    """
    num_fam = 1
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    btr_file = "test_btr_output.aln"
    nbfal = 4
    logger = logging.getLogger("test_backtranslate")
    # check that it returns true
    assert al.back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger) == 4
    # Check logs
    assert "Back-translating family 1" in caplog.text
    # Chek btr file content
    assert os.path.isfile(btr_file)
    exp_btr = os.path.join(EXPPATH, "exp_btr.1.aln")
    same_files(btr_file, exp_btr)
    os.remove(btr_file)


def test_backtranslate_wrongnbfam(caplog):
    """
    Test that when giving an alignment file and a file with nucleic sequences, it generates a
    btr file, with the nucleotide alignments as expected, but if we give a wrong number of
    expected sequences in btr file, it returns False with an error message
    """
    num_fam = 1
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    btr_file = "test_btr_output.aln"
    nbfal = 3
    logger = logging.getLogger("test_backtranslate")
    # check that it returns true
    assert not al.back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger)
    # Check logs
    assert "Back-translating family 1" in caplog.text
    assert ("fam 1: different number of proteins aligned in "
            "test/data/align/test_files/test_alignment-mafft-ok.aln (3) and genes "
            "back-translated in test_btr_output.aln (4)") in caplog.text
    # Chek btr file content
    assert os.path.isfile(btr_file)
    exp_btr = os.path.join(EXPPATH, "exp_btr.1.aln")
    same_files(btr_file, exp_btr)
    os.remove(btr_file)


def test_backtranslate_problem(caplog):
    """
    Test that when giving an alignment file and a wrong filename for the  nucleic sequences,
    it generates an empty btr file, and returns False, and an error message.
    """
    num_fam = 1
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1")
    btr_file = "test_btr_output.aln"
    nbfal = 4
    logger = logging.getLogger("test_backtranslate")
    # check that it returns true
    assert not al.back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger)
    # Check logs
    assert "Back-translating family 1" in caplog.text
    assert ("Problem while trying to backtranslate test/data/align/test_files/test_alignment-"
            "mafft-ok.aln to a nucleotide alignment") in caplog.text
    # Chek btr file content
    assert os.path.isfile(btr_file)
    with open(btr_file, "r") as btrf:
        assert btrf.readlines() == []
    os.remove(btr_file)


def test_check_mafft_align():
    """
    Test that when giving an alignment file, and a number of sequences equal to the number of
    sequences in the alignment file, it returns this number.
    """
    nbfal = 4
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    message = "problem!"
    logger = logging.getLogger("test_check_mafft_align")
    assert al.check_nb_seqs(mafft_file, nbfal, logger, message) == 4


def test_check_mafft_align_wrongnbfam(caplog):
    """
    Test that when giving an alignment file, and a number of sequences equal to the number of
    sequences in the alignment file, it returns this number.
    """
    nbfal = 3
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    message = "problem!"
    logger = logging.getLogger("test_check_mafft_align")
    assert not al.check_nb_seqs(mafft_file, nbfal, logger, message)
    assert 'problem! (4)' in caplog.text


def test_mafft_align(caplog):
    """
    Test that when giving a file containing extracted proteins, it aligns them as expected
    and returns True
    """
    num_fam = 1
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    mafft_file = "test_mafft_align.aln"
    nbfprt = 4
    logger = logging.getLogger("test_check_mafft_align")
    assert al.mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger) == 4
    assert "Aligning family 1" in caplog.text
    assert os.path.isfile(mafft_file)
    exp_mafft = os.path.join(EXPPATH, "exp_aldir", "mafft-align.1.aln")
    same_files(mafft_file, exp_mafft)
    os.remove(mafft_file)


def test_mafft_align_error(caplog):
    """
    Test that when giving a wrong file with the sequence to extract (non existing file),
    it returns false with the expected error message
    """
    num_fam = 1
    prt_file = "prtfile"
    mafft_file = "test_mafft_align.aln"
    nbfprt = 4
    logger = logging.getLogger("test_check_mafft_align")
    assert not al.mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger)
    assert "Aligning family 1" in caplog.text
    assert "Problem while trying to align fam 1" in caplog.text
    assert os.path.isfile(mafft_file)
    with open(mafft_file, "r") as btrf:
        assert btrf.readlines() == []
    os.remove(mafft_file)


def test_mafft_align_wrongnbfam(caplog):
    """
    Test that when giving a file containing extracted protein to align, it aligns them,
    but as we don't give the same number of expected aligned proteins, it returns False,
    with an error message
    """
    num_fam = 1
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    mafft_file = "test_mafft_align.aln"
    nbfprt = 3
    logger = logging.getLogger("test_check_mafft_align")
    assert not al.mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger)
    assert "Aligning family 1" in caplog.text
    assert ("fam 1: different number of proteins extracted in "
            "test/data/align/exp_files/exp_aldir/current.1.prt (3) and proteins aligned "
            "in test_mafft_align.aln (4)") in caplog.text
    assert os.path.isfile(mafft_file)
    exp_mafft = os.path.join(EXPPATH, "exp_aldir", "mafft-align.1.aln")
    same_files(mafft_file, exp_mafft)
    os.remove(mafft_file)


def test_check_extract(caplog):
    """
    Test that given the 3 files (prt, gen and miss) for a given family, having the expected number
    of sequences, it returns True
    """
    num_fam = 1
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    miss_file = "test_check_extract_miss-file.txt"
    ngenomes = 5
    logger = logging.getLogger("test_check_extract")
    with open(miss_file, "w") as missf:
        missf.write("Genome5")
    assert al.check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger) == 4
    assert "Checking extractions for family 1" in caplog.text
    os.remove(miss_file)


def test_check_extract_wrongnbmiss(caplog):
    """
    Test that given the 3 files: 4 proteins extracted in gen and prt, empty miss file,
    and saying that there are 5 genomes in the dataset, it returns False, with an error message
    indicating the wrong sum of prt_seqs + miss (4 instead of 5 expected)
    """
    num_fam = 1
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    miss_file = "test_check_extract_miss-file.txt"
    ngenomes = 5
    logger = logging.getLogger("test_check_extract")
    open(miss_file, "w").close()
    assert not al.check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger)
    assert "Checking extractions for family 1" in caplog.text
    assert ("fam 1: wrong sum of missing genomes (0) and prt extracted (4) for 5 genomes in "
            "the dataset.") in caplog.text
    os.remove(miss_file)


def test_check_extract_wrongnbgen(caplog):
    """
    Test that given the 3 files: 4 proteins extracted in gen, 5 in prt, empty miss file,
    and saying that there are 5 genomes in the dataset, it returns False, with an error message
    indicating the wrong sum of gen_seqs + miss (4 instead of 5 expected)
    """
    num_fam = 1
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    prt_file_ref = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    prt_file = "test_check_extract_prt-file.prt"
    shutil.copyfile(prt_file_ref, prt_file)
    miss_file = "test_check_extract_miss-file.txt"
    ngenomes = 5
    logger = logging.getLogger("test_check_extract")
    # Create empty miss file
    open(miss_file, "w").close()
    # Add a 5th sequence in prt file
    with open(prt_file, "a") as prtf:
        prtf.write(">fifth_sequence\nATTCGC\n")
    assert not al.check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger)
    assert "Checking extractions for family 1" in caplog.text
    assert ("fam 1: wrong sum of missing genomes (0) and gen extracted (4) for 5 genomes in "
            "the dataset.") in caplog.text
    os.remove(miss_file)
    os.remove(prt_file)


def test_family_align(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, it aligns all families, and backtranslates the alignment to
    nucleotides as expected, and returns the number of sequences aligned.
    """
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    mafft_file = "test_fam_align.8.aln"
    btr_file = "test_fam_align_btr.8.aln"
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == 3
    assert "Checking extractions for family 8" in caplog.text
    assert "Aligning family 8" in caplog.text
    assert "Back-translating family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    same_files(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    same_files(btr_file, exp_btr)
    os.remove(mafft_file)
    os.remove(btr_file)


def test_family_align_wrongextract(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 5, it returns False, error message for wrong sum of prt + miss,
    and removes existing mafft and btr files.
    """
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    mafft_file = "test_fam_align.8.aln"
    btr_file = "test_fam_align_btr.8.aln"
    open(mafft_file, "w").close()
    open(btr_file, "w").close()
    num_fam = 8
    ngenomes = 5
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) is False
    assert "Checking extractions for family 8" in caplog.text
    assert ("fam 8: wrong sum of missing genomes (1) and prt extracted (3) for 5 genomes in "
            "the dataset") in caplog.text
    # Check that mafft and btr files were removed
    assert not os.path.isfile(mafft_file)
    assert not os.path.isfile(btr_file)


def test_family_align_mafftok_nobtr(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, and the alignment file with the 3 sequences aligned,
    it back-translates it and returns 3, the number of sequences back-translated
    """
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create mafft file with all alignments
    mafft_file = "test_fam_align.8.aln"
    ref_mafft_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    shutil.copyfile(ref_mafft_file, mafft_file)
    btr_file = "test_fam_align_btr.8.aln"
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == 3
    assert "Checking extractions for family 8" in caplog.text
    assert "Back-translating family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    same_files(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    same_files(btr_file, exp_btr)
    os.remove(mafft_file)
    os.remove(btr_file)


def test_family_align_mafftempty_btrempty(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, and the alignment file with the 3 sequences aligned,
    it back-translates it and returns 3, the number of sequences back-translated
    """
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create empty mafft file
    mafft_file = "test_fam_align.8.aln"
    open(mafft_file, "w").close()
    # Create empty btr
    btr_file = "test_fam_align_btr.8.aln"
    open(btr_file, "w").close()
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == 3
    assert "Checking extractions for family 8" in caplog.text
    assert ("fam 8: Will redo alignment, because found a different number of proteins extracted "
            "in test/data/align/exp_files/exp_aldir-pers/current.8.prt (3) and proteins aligned "
            "in existing test_fam_align.8.aln (0)") in caplog.text
    assert "Aligning family 8" in caplog.text
    assert "Back-translating family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    same_files(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    same_files(btr_file, exp_btr)
    os.remove(mafft_file)
    os.remove(btr_file)


def test_family_align_mafftexists_emptybtr(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, and the alignment file with the 3 sequences aligned,
    and an empty btr file, it regenerates the back-translated file, and returns 3,
    the number of sequences back-translated.
    """
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create mafft file with all alignments
    mafft_file = "test_fam_align.8.aln"
    ref_mafft_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    shutil.copyfile(ref_mafft_file, mafft_file)
    # Create empty btr file
    btr_file = "test_fam_align_btr.8.aln"
    open(btr_file, "w").close()
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == 3
    assert "Checking extractions for family 8" in caplog.text
    assert ("fam 8: Will redo back-translation, because found a different number of "
            "proteins aligned in test_fam_align.8.aln (3) and genes back-translated in existing "
            "test_fam_align_btr.8.aln (0)") in caplog.text
    assert "Back-translating family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    same_files(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    same_files(btr_file, exp_btr)
    os.remove(mafft_file)
    os.remove(btr_file)


def test_family_align_mafftexists_btrok(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, and the alignment file with the 3 sequences aligned,
    and the already generated back-translated file, it returns "OK" as it did not generate any
    new file.
    """
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create mafft file with all alignments
    mafft_file = "test_fam_align.8.aln"
    ref_mafft_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    shutil.copyfile(ref_mafft_file, mafft_file)
    # Create btr file with alignments
    btr_file = "test_fam_align_btr.8.aln"
    ref_btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    shutil.copyfile(ref_btr_file, btr_file)
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == "OK"
    assert "Checking extractions for family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    same_files(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    same_files(btr_file, exp_btr)
    os.remove(mafft_file)
    os.remove(btr_file)

# TODO: no mafft file, btr file (empty), and problem while running mafft -> should return False,
# and btr file removed.


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
