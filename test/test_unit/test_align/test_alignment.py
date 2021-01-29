#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the alignment submodule in align module
"""
import logging
import os
import subprocess
import pytest
import shutil

import multiprocessing

import PanACoTA.align_module.alignment as al
import test.test_unit.utilities_for_tests as tutil


# Define common variables
ALDIR = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALDIR, "exp_files")
TESTPATH = os.path.join(ALDIR, "test_files")
GENEPATH = os.path.join(ALDIR, "generated_by_unit-tests")


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
    caplog.set_level(logging.DEBUG)
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
    number of sequences in the alignment file, it returns the number of sequences.
    """
    btr_file = os.path.join(TESTPATH, "test_alignment-btr-ok.aln")
    nbfal = 4
    logger = logging.getLogger("test_check_btr")
    message = "problem!"
    assert al.check_nb_seqs(btr_file, nbfal, logger, message) == 4


def test_check_btr_listsizes():
    """
    Test that when giving a nucleotide alignment, and a list of numbers, including 1
    corresponding to the number of sequences corresponding to the
    number of sequences in the alignment file, it returns the number of sequences in the
    alignment file.
    """
    btr_file = os.path.join(TESTPATH, "test_alignment-btr-ok.aln")
    nbfal = [3, 4, 5]
    logger = logging.getLogger("test_check_btr")
    message = "problem!"
    assert al.check_nb_seqs(btr_file, nbfal, logger, message) == 4


def test_check_btr_wrongnbfam(caplog):
    """
    Test that when giving a nucleotide alignment, and a number of sequences different from the
    number of sequences in the alignment file, it returns False, and an error message.
    """
    caplog.set_level(logging.DEBUG)
    btr_file = os.path.join(TESTPATH, "test_alignment-btr-ok.aln")
    nbfal = 5
    logger = logging.getLogger("test_check_btr")
    message = "Problem with fam 1"
    assert not al.check_nb_seqs(btr_file, nbfal, logger, message)
    assert "Problem with fam 1 (4)" in caplog.text


def test_check_btr_listwrongnbfam(caplog):
    """
    Test that when giving a nucleotide alignment, and a number of sequences different from the
    number of sequences in the alignment file, it returns False, and an error message.
    """
    caplog.set_level(logging.DEBUG)
    btr_file = os.path.join(TESTPATH, "test_alignment-btr-ok.aln")
    nbfal = [3, 5, 7]
    logger = logging.getLogger("test_check_btr")
    message = "Problem with fam 1"
    assert not al.check_nb_seqs(btr_file, nbfal, logger, message)
    assert "Problem with fam 1 (4)" in caplog.text


def test_backtranslate(caplog):
    """
    Test that when giving an alignment file and a file with nucleic sequences, it generates a
    btr file, with the nucleotide alignments as expected, and returns True
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    btr_file = os.path.join(GENEPATH, "test_btr_output.aln")
    nbfal = 4
    logger = logging.getLogger("test_backtranslate")
    # check that it returns true
    assert al.back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger) == 4
    # Check logs
    assert "Back-translating family 1" in caplog.text
    # Chek btr file content
    assert os.path.isfile(btr_file)
    exp_btr = os.path.join(EXPPATH, "exp_btr.1.aln")
    assert tutil.compare_file_content(btr_file, exp_btr)


def test_backtranslate_wrongnbfam(caplog):
    """
    Test that when giving an alignment file and a file with nucleic sequences, it generates a
    btr file, with the nucleotide alignments as expected, but if we give a wrong number of
    expected sequences in btr file, it returns False with an error message
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    btr_file = os.path.join(GENEPATH, "test_btr_output.aln")
    nbfal = 3
    logger = logging.getLogger("test_backtranslate")
    # check that it returns true
    assert not al.back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger)
    # Check logs
    assert "Back-translating family 1" in caplog.text
    assert ("fam 1: different number of proteins aligned in "
            "test/data/align/test_files/test_alignment-mafft-ok.aln (3) and genes back-translated in "
            "test/data/align/generated_by_unit-tests/test_btr_output.aln ") in caplog.text
    # Chek btr file content
    assert os.path.isfile(btr_file)
    exp_btr = os.path.join(EXPPATH, "exp_btr.1.aln")
    assert tutil.compare_file_content(btr_file, exp_btr)


def test_backtranslate_problem(caplog):
    """
    Test that when giving an alignment file and a wrong filename for the  nucleic sequences,
    it generates an empty btr file, and returns False, and an error message.
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    mafft_file = os.path.join(TESTPATH, "test_alignment-mafft-ok.aln")
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1")
    btr_file = os.path.join(GENEPATH, "test_btr_output.aln")
    nbfal = 4
    logger = logging.getLogger("test_backtranslate")
    # check that it returns true
    assert not al.back_translate(num_fam, mafft_file, gen_file, btr_file, nbfal, logger)
    # Check logs
    assert "Back-translating family 1" in caplog.text
    assert ("Problem while trying to backtranslate test/data/align/test_files/test_alignment-"
            "mafft-ok.aln to a nucleotide alignment") in caplog.text
    # Chek btr file content
    assert not os.path.isfile(btr_file)


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


def test_mafft_align(caplog):
    """
    Test that when giving a file containing extracted proteins, it aligns them as expected
    and returns True
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    mafft_file = os.path.join(GENEPATH, "test_mafft_align.aln")
    nbfprt = 4
    logger = logging.getLogger("test_check_mafft_align")
    assert al.mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger) == 4
    assert "Aligning family 1" in caplog.text
    assert os.path.isfile(mafft_file)
    exp_mafft = os.path.join(EXPPATH, "exp_aldir", "mafft-align.1.aln")
    assert tutil.compare_file_content(mafft_file, exp_mafft)


def test_mafft_align_error(caplog):
    """
    Test that when giving a wrong file with the sequence to extract (non existing file),
    it returns false with the expected error message
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    prt_file = "prtfile"
    mafft_file = os.path.join(GENEPATH, "test_mafft_align.aln")
    nbfprt = 4
    logger = logging.getLogger("test_check_mafft_align")
    assert not al.mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger)
    assert "Aligning family 1" in caplog.text
    assert "Problem while trying to align fam 1" in caplog.text
    assert not os.path.isfile(mafft_file)


def test_mafft_align_wrongnbfam(caplog):
    """
    Test that when giving a file containing extracted protein to align, it aligns them,
    but as we don't give the same number of expected aligned proteins, it returns False,
    with an error message
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    mafft_file = os.path.join(GENEPATH, "test_mafft_align.aln")
    nbfprt = 3
    logger = logging.getLogger("test_check_mafft_align")
    assert not al.mafft_align(num_fam, prt_file, mafft_file, nbfprt, logger)
    assert "Aligning family 1" in caplog.text
    assert ("fam 1: different number of proteins extracted in "
            "test/data/align/exp_files/exp_aldir/current.1.prt (3) and proteins aligned "
            "in test/data/align/generated_by_unit-tests/test_mafft_align.aln (4)") in caplog.text
    assert os.path.isfile(mafft_file)
    exp_mafft = os.path.join(EXPPATH, "exp_aldir", "mafft-align.1.aln")
    assert tutil.compare_file_content(mafft_file, exp_mafft)


def test_check_extract(caplog):
    """
    Test that given the 3 files (prt, gen and miss) for a given family, having the expected number
    of sequences, it returns True
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    miss_file = os.path.join(GENEPATH, "test_check_extract_miss-file.txt")
    ngenomes = 5
    logger = logging.getLogger("test_check_extract")
    with open(miss_file, "w") as missf:
        missf.write("Genome5")
    assert al.check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger) == 4
    assert "Checking extractions for family 1" in caplog.text


def test_check_extract_wrongnbmiss(caplog):
    """
    Test that given the 3 files: 4 proteins extracted in gen and prt, empty miss file,
    and saying that there are 5 genomes in the dataset, it returns False, with an error message
    indicating the wrong sum of prt_seqs + miss (4 instead of 5 expected)
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    prt_file = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    miss_file = os.path.join(GENEPATH, "test_check_extract_miss-file.txt")
    ngenomes = 5
    logger = logging.getLogger("test_check_extract")
    open(miss_file, "w").close()
    assert not al.check_extractions(num_fam, miss_file, prt_file, gen_file, ngenomes, logger)
    assert "Checking extractions for family 1" in caplog.text
    assert ("fam 1: wrong sum of missing genomes (0) and prt extracted (4) for 5 genomes in "
            "the dataset.") in caplog.text


def test_check_extract_wrongnbgen(caplog):
    """
    Test that given the 3 files: 4 proteins extracted in gen, 5 in prt, empty miss file,
    and saying that there are 5 genomes in the dataset, it returns False, with an error message
    indicating the wrong sum of gen_seqs + miss (4 instead of 5 expected)
    """
    caplog.set_level(logging.DEBUG)
    num_fam = 1
    gen_file = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    prt_file_ref = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    prt_file = os.path.join(GENEPATH, "test_check_extract_prt-file.prt")
    shutil.copyfile(prt_file_ref, prt_file)
    miss_file = os.path.join(GENEPATH, "test_check_extract_miss-file.txt")
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


def test_family_align(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, it aligns all families, and backtranslates the alignment to
    nucleotides as expected, and returns the number of sequences aligned.
    """
    caplog.set_level(logging.DEBUG)
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    mafft_file = os.path.join(GENEPATH, "test_fam_align.8.aln")
    btr_file = os.path.join(GENEPATH, "test_fam_align_btr.8.aln")
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
    assert tutil.compare_file_content(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    assert tutil.compare_file_content(btr_file, exp_btr)


def test_family_align_wrongextract(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 5, it returns False, error message for wrong sum of prt + miss,
    and removes existing mafft and btr files.
    """
    caplog.set_level(logging.DEBUG)
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    mafft_file = os.path.join(GENEPATH, "test_fam_align.8.aln")
    btr_file = os.path.join(GENEPATH, "test_fam_align_btr.8.aln")
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
    caplog.set_level(logging.DEBUG)
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create mafft file with all alignments
    btr_file = os.path.join(GENEPATH, "test_fam_align_btr.8.aln")
    mafft_file = os.path.join(GENEPATH, "test_fam_align.8.aln")
    ref_mafft_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    shutil.copyfile(ref_mafft_file, mafft_file)
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == 3
    assert "Checking extractions for family 8" in caplog.text
    assert "Back-translating family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    assert tutil.compare_file_content(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    assert tutil.compare_file_content(btr_file, exp_btr)


def test_family_align_mafftempty_btrempty(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, and the alignment file with the 3 sequences aligned,
    it back-translates it and returns 3, the number of sequences back-translated
    """
    caplog.set_level(logging.DEBUG)
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create empty mafft file
    mafft_file = os.path.join(GENEPATH, "test_fam_align.8.aln")
    open(mafft_file, "w").close()
    # Create empty btr
    btr_file = os.path.join(GENEPATH, "test_fam_align_btr.8.aln")
    open(btr_file, "w").close()
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == 3
    assert "Checking extractions for family 8" in caplog.text
    assert ("fam 8: Will redo alignment, because found a different number of proteins extracted "
            "in test/data/align/exp_files/exp_aldir-pers/current.8.prt (3) and proteins aligned "
            "in existing test/data/align/generated_by_unit-tests/test_fam_align.8.aln (0)") in caplog.text
    assert "Aligning family 8" in caplog.text
    assert "Back-translating family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    assert tutil.compare_file_content(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    assert tutil.compare_file_content(btr_file, exp_btr)


def test_family_align_mafftok_emptybtr(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, and the alignment file with the 3 sequences aligned,
    and an empty btr file, it re-does the backtranslation, and returns True
    """
    caplog.set_level(logging.DEBUG)
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create mafft file with all alignments
    mafft_file = os.path.join(GENEPATH, "test_fam_align.8.aln")
    ref_mafft_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    shutil.copyfile(ref_mafft_file, mafft_file)
    # Create empty btr file
    btr_file = os.path.join(GENEPATH, "test_fam_align_btr.8.aln")
    open(btr_file, "w").close()
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                               ngenomes, logger) == 3
    assert "Checking extractions for family 8" in caplog.text
    # Check content of mafft and btr files
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    assert tutil.compare_file_content(mafft_file, exp_mafft)
    assert tutil.compare_file_content(btr_file, exp_btr)


def test_family_align_mafftok_btrok(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, and the alignment file with the 3 sequences aligned,
    and the already generated back-translated file, it returns "OK" as it did not generate any
    new file.
    """
    caplog.set_level(logging.DEBUG)
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    # Create mafft file with all alignments
    mafft_file = os.path.join(GENEPATH, "test_fam_align.8.aln")
    ref_mafft_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    shutil.copyfile(ref_mafft_file, mafft_file)
    # Create btr file with alignments
    btr_file = os.path.join(GENEPATH, "test_fam_align_btr.8.aln")
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
    assert tutil.compare_file_content(mafft_file, exp_mafft)
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    assert tutil.compare_file_content(btr_file, exp_btr)


def test_family_align_nomafft_btrempty_errormafft(caplog):
    """
    Test that when giving prt file (3 extracted), gen file (3 extracted), miss file (1)
    and total nb genomes = 4, no alignment file, and an empty btr file, and make changes in
    path so that mafft crashes -> should return False, and remove btr and mafft files.
    """
    caplog.set_level(logging.DEBUG)
    orig_mafft = subprocess.check_output("which mafft".split()).decode().strip()
    temp_mafft = orig_mafft + "-orig"
    try:
        shutil.move(orig_mafft, temp_mafft)
    except PermissionError:
        pytest.skip("No permission to change mafft bin filename")
    prt_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    gen_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    mafft_file = os.path.join(GENEPATH, "test_fam_align.8.aln")
    # Create empty btr file
    btr_file = os.path.join(GENEPATH, "test_fam_align_btr.8.aln")
    open(btr_file, "w").close()
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_fam_align")
    assert not al.family_alignment(prt_file, gen_file, miss_file, mafft_file, btr_file, num_fam,
                                   ngenomes, logger)
    assert "Checking extractions for family 8" in caplog.text
    assert "Aligning family 8" in caplog.text
    assert ("command '>mafft --auto test/data/align/exp_files/exp_aldir-pers/current.8.prt' "
            "is not possible") in caplog.text
    # Check content of mafft and btr files
    assert not os.path.isfile(mafft_file)
    assert not os.path.isfile(btr_file)
    shutil.move(temp_mafft, orig_mafft)


def test_check_addmissing_ok():
    """
    Test that when giving a btr file with all sequences same length, and as many sequences as
    expected, it returns True
    """
    btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    num_fam = 8
    ngenomes = 3
    logger = logging.getLogger("test_check_add_missing")
    assert al.check_add_missing(btr_file, num_fam, ngenomes, logger) is True


def test_check_addmissing_difflen():
    """
    Test that when giving a btr file with sequences with different lengths, it returns False
    """
    btr_file = os.path.join(TESTPATH, "test_alignment-mafft-diff-len.aln")
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_check_add_missing")
    assert al.check_add_missing(btr_file, num_fam, ngenomes, logger) is False


def test_check_addmissing_missgenomes_prev():
    """
    Test that when giving a btr file with all sequences of same length, but not the same number
    of sequences than expected, and comes from previous run, it returns the length of alignment,
    without any error message.
    """
    btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_check_add_missing")
    assert al.check_add_missing(btr_file, num_fam, ngenomes, logger, prev=True) == 789


def test_check_addmissing_missgenomes_noprev(caplog):
    """
    Test that when giving a btr file with all sequences of same length, but not the same number
    of sequences than expected, and does not come from previous run, it returns the length of
    alignment, with an error message.
    """
    caplog.set_level(logging.DEBUG)
    btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    num_fam = 8
    ngenomes = 4
    logger = logging.getLogger("test_check_add_missing")
    assert al.check_add_missing(btr_file, num_fam, ngenomes, logger) == 789
    assert ("ERROR: family 8 contains 3 genomes in total instead of the 4 genomes in "
            "input.") in caplog.text


def test_add_missing_btrok():
    """
    Giving a btr_file with as many genomes as expected, and this btr file was just created (means
    that there is no missing genome in this family), it returns True
    """
    btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    # Empty miss file (no missing genome in this family)
    miss_file = "test_add_missing"
    num_fam = 8
    ngenomes = 3
    status1 = True
    logger = logging.getLogger("test_add_missing")
    assert al.add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, status1, logger)


def test_add_missing_btralreadyok(caplog):
    """
    Giving a btr_file with as many genomes as expected, and this btr file was not recreated (
    already existed from previous run), it returns "OK" with a warning message specifying that
    program will use already existing file.
    """
    caplog.set_level(logging.DEBUG)
    btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    # Empty miss file (no missing genome in this family)
    miss_file = "test_add_missing"
    num_fam = 8
    ngenomes = 3
    status1 = "OK"
    logger = logging.getLogger("test_add_missing")
    assert al.add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, status1, logger) == "OK"
    assert ("Alignment already done for family 8. The program will use it for next "
            "steps") in caplog.text


def test_add_missing_btr_difflen(caplog):
    """
    Giving a btr_file with as many genomes as expected, but different alignment lengths for the
    different sequences. It returns False, with an error message saying that there are different
    lengths in alignment.
    """
    caplog.set_level(logging.DEBUG)
    btr_file = os.path.join(TESTPATH, "test_alignment-mafft-diff-len.aln")
    # Empty miss file (no missing genome in this family)
    miss_file = "test_add_missing"
    num_fam = 8
    ngenomes = 4
    status1 = "OK"
    logger = logging.getLogger("test_add_missing")
    assert al.add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, status1, logger) is False
    assert ("alignments for family 8 do not all have the same length. Lengths found "
            "are: {") in caplog.text
    assert "126" in caplog.text
    assert "154" in caplog.text
    assert "157" in caplog.text
    assert "}\n" in caplog.text


def test_add_missing_btrmiss(caplog):
    """
    Giving a btr_file with all sequences with same lengths, but 1 sequence less that the
    expected number, + miss file with the name of the missing genome.
    Returns True, after adding missing genome to btr_file (overwrites it).
    """
    caplog.set_level(logging.DEBUG)
    ref_btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    # Copy ref_btr to another file, which will be overwritten.
    btr_file = os.path.join(GENEPATH, "test_add_missing_genomes.aln")
    shutil.copyfile(ref_btr_file, btr_file)
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    num_fam = 8
    ngenomes = 4
    status1 = "OK"
    logger = logging.getLogger("test_add_missing")
    assert al.add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, status1, logger) is True
    assert "Adding missing genomes for family 8" in caplog.text
    final_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    assert tutil.compare_file_content(btr_file, final_btr)


def test_add_missing_btrmiss2(caplog):
    """
    Giving a btr_file with all sequences with same lengths, but 2 sequences less that the
    expected number, + miss file with the name of 1 missing genome.
    Returns False, after adding 1 missing genome to btr_file. + error message saying that a
    genome is missing in btr_file.
    """
    caplog.set_level(logging.DEBUG)
    ref_btr_file = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    # Copy ref_btr to another file, which will be overwritten.
    btr_file = os.path.join(GENEPATH, "test_add_missing_genomes.aln")
    shutil.copyfile(ref_btr_file, btr_file)
    miss_file = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    num_fam = 8
    ngenomes = 5
    status1 = True
    logger = logging.getLogger("test_add_missing")
    assert al.add_missing_genomes(btr_file, miss_file, num_fam, ngenomes, status1, logger) is False
    assert "Adding missing genomes for family 8" in caplog.text
    assert ("ERROR: family 8 contains 4 genomes in total instead of the 5 genomes in "
            "input") in caplog.text
    final_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    assert tutil.compare_file_content(btr_file, final_btr)


def test_handle_family_true():
    """
    Giving an aldir with correct prt, gen and miss files, it creates mafft and btr files as
    expected, and returns True
    """
    aldir = os.path.join(GENEPATH, "aldir")
    prefix = os.path.join(aldir, "TESThandlefam")
    num_fam = 8
    ngenomes = 4
    q = multiprocessing.Manager().Queue()
    # Create aldir, and put prt, gen and miss files in it
    ref_prt = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_miss = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    cur_prt = os.path.join(aldir, "TESThandlefam-current.8.prt")
    cur_gen = os.path.join(aldir, "TESThandlefam-current.8.gen")
    cur_miss = os.path.join(aldir, "TESThandlefam-current.8.miss.lst")
    # Create aldir inside GENEPATH, empty directory created for this test
    os.makedirs(aldir)
    shutil.copyfile(ref_prt, cur_prt)
    shutil.copyfile(ref_gen, cur_gen)
    shutil. copyfile(ref_miss, cur_miss)
    args = (prefix, num_fam, ngenomes, q)
    assert al.handle_family(args) is True
    cur_mafft = os.path.join(aldir, "TESThandlefam-mafft-align.8.aln")
    cur_btr = os.path.join(aldir, "TESThandlefam-mafft-prt2nuc.8.aln")
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    assert tutil.compare_order_content(cur_mafft, exp_mafft)
    assert tutil.compare_order_content(cur_btr, exp_btr)
    q.put(None)
    assert "Checking extractions for family 8" in q.get().message
    assert "Aligning family 8" in q.get().message
    assert ("mafft --auto test/data/align/generated_by_unit-tests/aldir/"
            "TESThandlefam-current.8.prt") in q.get().message
    assert "Back-translating family 8" in q.get().message
    assert "Adding missing genomes for family 8" in q.get().message
    assert not q.get()


def test_handle_family_true_nomiss():
    """
    Giving an aldir with correct prt, gen and miss files (no missing file in this family),
    it creates mafft and btr files as expected, and returns True. No need to add missing genomes
    as there are not in this family.
    """
    aldir = os.path.join(GENEPATH, "aldir")
    prefix = os.path.join(aldir, "TESThandlefam")
    num_fam = 8
    ngenomes = 3
    q = multiprocessing.Manager().Queue()
    # Create aldir, and put prt, gen and miss files in it
    ref_prt = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    cur_prt = os.path.join(aldir, "TESThandlefam-current.8.prt")
    cur_gen = os.path.join(aldir, "TESThandlefam-current.8.gen")
    cur_miss = os.path.join(aldir, "TESThandlefam-current.8.miss.lst")
    os.makedirs(aldir)
    shutil.copyfile(ref_prt, cur_prt)
    shutil.copyfile(ref_gen, cur_gen)
    open(cur_miss, "w").close()
    args = (prefix, num_fam, ngenomes, q)
    assert al.handle_family(args)
    cur_mafft = os.path.join(aldir, "TESThandlefam-mafft-align.8.aln")
    cur_btr = os.path.join(aldir, "TESThandlefam-mafft-prt2nuc.8.aln")
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-btr.8.aln")
    assert tutil.compare_order_content(cur_mafft, exp_mafft)
    assert tutil.compare_order_content(cur_btr, exp_btr)
    q.put(None)
    assert "Checking extractions for family 8" in q.get().message
    assert "Aligning family 8" in q.get().message
    assert ("mafft --auto test/data/align/generated_by_unit-tests/aldir/"
            "TESThandlefam-current.8.prt") in q.get().message
    assert "Back-translating family 8" in q.get().message
    assert not q.get()


def test_handle_family_emptyaln_true():
    """
    Giving an aldir with correct prt, gen, miss files, and an empty mafft file, it creates mafft
    and btr files as expected, and returns True
    """
    aldir = os.path.join(GENEPATH, "aldir")
    prefix = os.path.join(aldir, "TESThandlefam")
    num_fam = 8
    ngenomes = 4
    q = multiprocessing.Manager().Queue()
    # Create aldir, and put prt, gen and miss files in it
    ref_prt = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_miss = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    cur_prt = os.path.join(aldir, "TESThandlefam-current.8.prt")
    cur_gen = os.path.join(aldir, "TESThandlefam-current.8.gen")
    cur_miss = os.path.join(aldir, "TESThandlefam-current.8.miss.lst")
    cur_mafft = os.path.join(aldir, "TESThandlefam-mafft-align.8.aln")
    # Put real gen, prt and miss files to aldir, create empty mafft file
    os.makedirs(aldir)
    shutil.copyfile(ref_prt, cur_prt)
    shutil.copyfile(ref_gen, cur_gen)
    shutil. copyfile(ref_miss, cur_miss)
    open(cur_mafft, "w").close()
    args = (prefix, num_fam, ngenomes, q)
    assert al.handle_family(args)
    cur_btr = os.path.join(aldir, "TESThandlefam-mafft-prt2nuc.8.aln")
    exp_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    exp_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    assert tutil.compare_order_content(cur_mafft, exp_mafft)
    assert tutil.compare_order_content(cur_btr, exp_btr)
    q.put(None)
    assert "Checking extractions for family 8" in q.get().message
    assert ("fam 8: Will redo alignment, because found a different number of proteins "
            "extracted in "
            "test/data/align/generated_by_unit-tests/aldir/TESThandlefam-current.8.prt (3) "
            "and proteins aligned in existing "
            "test/data/align/generated_by_unit-tests/aldir/TESThandlefam-mafft-align.8.aln "
            "(0)") in q.get().message
    assert "Aligning family 8" in q.get().message
    assert ("mafft --auto test/data/align/generated_by_unit-tests/aldir/"
            "TESThandlefam-current.8.prt") in q.get().message
    assert "Back-translating family 8" in q.get().message
    assert "Adding missing genomes for family 8" in q.get().message
    assert not q.get()


def test_handle_family_emptybtr_true():
    """
    Giving an aldir with correct prt, gen, miss and mafft files, and an empty btr files,
    it redoes back-translation and returns True
    """
    aldir = os.path.join(GENEPATH, "aldir")
    prefix = os.path.join(aldir, "TESThandlefam")
    num_fam = 8
    ngenomes = 4
    q = multiprocessing.Manager().Queue()
    # Create aldir, and put prt, gen and miss files in it
    ref_prt = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_miss = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    ref_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    ref_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    cur_prt = os.path.join(aldir, "TESThandlefam-current.8.prt")
    cur_gen = os.path.join(aldir, "TESThandlefam-current.8.gen")
    cur_miss = os.path.join(aldir, "TESThandlefam-current.8.miss.lst")
    cur_mafft = os.path.join(aldir, "TESThandlefam-mafft-align.8.aln")
    cur_btr = os.path.join(aldir, "TESThandlefam-mafft-prt2nuc.8.aln")
    # Put real gen, prt and miss files to aldir, create empty mafft file
    os.makedirs(aldir)
    shutil.copyfile(ref_prt, cur_prt)
    shutil.copyfile(ref_gen, cur_gen)
    shutil.copyfile(ref_miss, cur_miss)
    shutil.copyfile(ref_mafft, cur_mafft)
    open(cur_btr, "w").close()
    args = (prefix, num_fam, ngenomes, q)
    assert al.handle_family(args)
    assert tutil.compare_order_content(cur_mafft, ref_mafft)
    assert tutil.compare_order_content(cur_btr, ref_btr)
    q.put(None)
    assert "Checking extractions for family 8" in q.get().message
    assert ("fam 8: Will redo back-translation, because found a different number of proteins "
            "aligned in "
            "test/data/align/generated_by_unit-tests/aldir/TESThandlefam-mafft-align.8.aln (3) "
            "and genes back-translated in "
            "existing test/data/align/generated_by_unit-tests/"
            "aldir/TESThandlefam-mafft-prt2nuc.8.aln (0)") in q.get().message
    assert "Back-translating family 8" in q.get().message
    assert "Adding missing genomes for family 8" in q.get().message
    assert not q.get()


def test_handle_family_already_ok():
    """
    Giving an aldir with correct prt, gen, miss, mafft and btr files, it returns "OK",
    not creating any new file, but with a warning message saying that already existing files are
    used for next steps.
    """
    aldir = os.path.join(GENEPATH, "aldir")
    prefix = os.path.join(aldir, "TESThandlefam")
    num_fam = 8
    ngenomes = 4
    q = multiprocessing.Manager().Queue()
    # Create aldir, and put prt, gen, miss, mafft and btr files in it
    ref_prt = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_miss = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    ref_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    ref_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    cur_prt = os.path.join(aldir, "TESThandlefam-current.8.prt")
    cur_gen = os.path.join(aldir, "TESThandlefam-current.8.gen")
    cur_miss = os.path.join(aldir, "TESThandlefam-current.8.miss.lst")
    cur_mafft = os.path.join(aldir, "TESThandlefam-mafft-align.8.aln")
    cur_btr = os.path.join(aldir, "TESThandlefam-mafft-prt2nuc.8.aln")
    os.makedirs(aldir)
    shutil.copyfile(ref_prt, cur_prt)
    shutil.copyfile(ref_gen, cur_gen)
    shutil.copyfile(ref_miss, cur_miss)
    shutil.copyfile(ref_mafft, cur_mafft)
    shutil.copyfile(ref_btr, cur_btr)
    args = (prefix, num_fam, ngenomes, q)
    assert al.handle_family(args) == "OK"
    q.put(None)
    assert "Checking extractions for family 8" in q.get().message
    assert ("Alignment already done for family 8. The program will use it for next "
            "steps") in q.get().message
    assert not q.get()


def test_handle_family_wrongextract():
    """
    Giving an aldir with correct prt, gen, miss, mafft and btr files, but says that there are 5
    genomes in dataset (instead of 4), it returns False, with an error message saying that
    extractions are wrong.
    """
    aldir = os.path.join(GENEPATH, "aldir")
    prefix = os.path.join(aldir, "TESThandlefam")
    num_fam = 8
    ngenomes = 5
    q = multiprocessing.Manager().Queue()
    # Create aldir, and put prt, gen, miss, mafft and btr files in it
    ref_prt = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_miss = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    ref_mafft = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    ref_btr = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    cur_prt = os.path.join(aldir, "TESThandlefam-current.8.prt")
    cur_gen = os.path.join(aldir, "TESThandlefam-current.8.gen")
    cur_miss = os.path.join(aldir, "TESThandlefam-current.8.miss.lst")
    cur_mafft = os.path.join(aldir, "TESThandlefam-mafft-align.8.aln")
    cur_btr = os.path.join(aldir, "TESThandlefam-mafft-prt2nuc.8.aln")
    os.makedirs(aldir)
    shutil.copyfile(ref_prt, cur_prt)
    shutil.copyfile(ref_gen, cur_gen)
    shutil.copyfile(ref_miss, cur_miss)
    shutil.copyfile(ref_mafft, cur_mafft)
    shutil.copyfile(ref_btr, cur_btr)
    args = (prefix, num_fam, ngenomes, q)
    assert not al.handle_family(args)
    assert not os.path.isfile(cur_mafft)
    assert not os.path.isfile(cur_btr)
    q.put(None)
    assert "Checking extractions for family 8" in q.get().message
    assert ("fam 8: wrong sum of missing genomes (1) and prt extracted (3) for 5 genomes in the "
            "dataset.") in q.get().message
    assert not q.get()


def test_align_all_true(caplog):
    """
    Giving aldir with prt, gen and miss files for families 1 and 8, as well as concat file (
    empty). It should create mafft and prt2nuc files for families 1 and 8, remove concat
    file, and return True.
    """
    caplog.set_level(logging.DEBUG)
    aldir = os.path.join(GENEPATH, "aldir")
    dname = "TESTalign-all"
    prefix = os.path.join(aldir, dname)
    all_fams = [1, 8]
    ngenomes = 4
    quiet = False
    threads = 1
    # Create aldir, and put all prt, gen and miss files for families 1 and 8
    os.makedirs(aldir)
    ref_prt8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_miss8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    ref_prt1 = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    ref_gen1 = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    cur_prt8 = os.path.join(aldir, dname + "-current.8.prt")
    cur_gen8 = os.path.join(aldir, dname + "-current.8.gen")
    cur_miss8 = os.path.join(aldir, dname + "-current.8.miss.lst")
    cur_prt1 = os.path.join(aldir, dname + "-current.1.prt")
    cur_gen1 = os.path.join(aldir, dname + "-current.1.gen")
    cur_miss1 = os.path.join(aldir, dname + "-current.1.miss.lst")
    refs = [ref_gen8, ref_gen1, ref_prt8, ref_prt1, ref_miss8]
    todos = [cur_gen8, cur_gen1, cur_prt8, cur_prt1, cur_miss8]
    for f1, f2 in zip(refs, todos):
        shutil.copyfile(f1, f2)
    # empty miss file for family 1
    open(cur_miss1, "w").close()
    assert os.path.isfile(cur_miss1)
    # create empty concat-file, to check that it is removed
    concat = os.path.join(aldir, dname + "-complete.cat.aln")
    open(concat, "w").close()
    assert al.align_all_families(prefix, all_fams, ngenomes, dname, quiet, threads)
    # Check output files
    out_mafft1 = os.path.join(aldir, dname + "-mafft-align.1.aln")
    out_btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    out_mafft8 = os.path.join(aldir, dname + "-mafft-align.8.aln")
    out_btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    exp_mafft1 = os.path.join(EXPPATH, "exp_aldir", "mafft-align.1.aln")
    exp_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    exp_mafft8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    exp_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    assert tutil.compare_order_content(out_mafft1, exp_mafft1)
    assert tutil.compare_order_content(out_btr1, exp_btr1)
    assert tutil.compare_order_content(out_mafft8, exp_mafft8)
    assert tutil.compare_order_content(out_btr8, exp_btr8)
    assert not os.path.isfile(concat)
    # Check logs
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    assert "Checking extractions for family 1" in caplog.text
    assert "Aligning family 1" in caplog.text
    assert "Back-translating family 1" in caplog.text
    assert "Checking extractions for family 8" in caplog.text
    assert "Aligning family 8" in caplog.text
    assert "Back-translating family 8" in caplog.text
    assert "Adding missing genomes for family 8" in caplog.text


def test_align_all_exists_true(caplog):
    """
    Giving aldir with prt, gen, miss, mafft and prt2nuc files for families 1 and 8, as well as
    concat file (empty). It should keep concat file, and return True.
    """
    caplog.set_level(logging.DEBUG)
    aldir = os.path.join(GENEPATH, "aldir")
    dname = "TESTalign-all"
    prefix = os.path.join(aldir, dname)
    all_fams = [1, 8]
    ngenomes = 4
    quiet = False
    threads = 2
    # Create aldir, and put all prt, gen and miss files for families 1 and 8
    os.makedirs(aldir)
    ref_prt8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_miss8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.miss.lst")
    ref_prt1 = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    ref_gen1 = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    exp_mafft1 = os.path.join(EXPPATH, "exp_aldir", "mafft-align.1.aln")
    exp_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    exp_mafft8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-align.8.aln")
    exp_btr8 = os.path.join(EXPPATH, "exp_aldir-pers", "mafft-prt2nuc.8.aln")
    cur_prt8 = os.path.join(aldir, dname + "-current.8.prt")
    cur_gen8 = os.path.join(aldir, dname + "-current.8.gen")
    cur_miss8 = os.path.join(aldir, dname + "-current.8.miss.lst")
    cur_prt1 = os.path.join(aldir, dname + "-current.1.prt")
    cur_gen1 = os.path.join(aldir, dname + "-current.1.gen")
    cur_miss1 = os.path.join(aldir, dname + "-current.1.miss.lst")
    out_mafft1 = os.path.join(aldir, dname + "-mafft-align.1.aln")
    out_btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    out_mafft8 = os.path.join(aldir, dname + "-mafft-align.8.aln")
    out_btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    refs = [ref_gen8, ref_gen1, ref_prt8, ref_prt1, ref_miss8, exp_mafft1, exp_mafft8, exp_btr1,
            exp_btr8]
    todos = [cur_gen8, cur_gen1, cur_prt8, cur_prt1, cur_miss8, out_mafft1, out_mafft8, out_btr1,
             out_btr8]
    for f1, f2 in zip(refs, todos):
        shutil.copyfile(f1, f2)
    # empty miss file for family 1
    open(cur_miss1, "w").close()
    # create empty concat-file, to check that it is removed
    concat = os.path.join(aldir, dname + "-complete.cat.aln")
    open(concat, "w").close()
    assert al.align_all_families(prefix, all_fams, ngenomes, dname, quiet, threads)
    # Check output files
    assert tutil.compare_order_content(out_mafft1, exp_mafft1)
    assert tutil.compare_order_content(out_btr1, exp_btr1)
    assert tutil.compare_order_content(out_mafft8, exp_mafft8)
    assert tutil.compare_order_content(out_btr8, exp_btr8)
    assert os.path.isfile(concat)
    with open(concat, "r") as conf:
        assert conf.readlines() == []
    # Check logs
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    assert "Checking extractions for family 1" in caplog.text
    assert ("Alignment already done for family 1. The program will use it for next "
            "steps") in caplog.text
    assert "Checking extractions for family 8" in caplog.text
    assert ("Alignment already done for family 8. The program will use it for next "
            "steps") in caplog.text


def test_align_all_false(caplog):
    """
    Giving aldir with prt, gen, miss files for families 1 and 8, as well as
    concat file (empty). But for family 8, miss file is empty: total number of genomes does not
    correspond. It should return error message for fam 8, info logs for fam 1, remove concat file,
    and return False
    """
    caplog.set_level(logging.DEBUG)
    aldir = os.path.join(GENEPATH, "aldir")
    dname = "TESTalign-all"
    prefix = os.path.join(aldir, dname)
    all_fams = [1, 8]
    ngenomes = 4
    quiet = False
    threads = 1
    # Create aldir, and put all prt, gen and miss files for families 1 and 8
    os.makedirs(aldir)
    ref_prt8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.prt")
    ref_gen8 = os.path.join(EXPPATH, "exp_aldir-pers", "current.8.gen")
    ref_prt1 = os.path.join(EXPPATH, "exp_aldir", "current.1.prt")
    ref_gen1 = os.path.join(EXPPATH, "exp_aldir", "current.1.gen")
    cur_prt8 = os.path.join(aldir, dname + "-current.8.prt")
    cur_gen8 = os.path.join(aldir, dname + "-current.8.gen")
    cur_miss8 = os.path.join(aldir, dname + "-current.8.miss.lst")
    cur_prt1 = os.path.join(aldir, dname + "-current.1.prt")
    cur_gen1 = os.path.join(aldir, dname + "-current.1.gen")
    cur_miss1 = os.path.join(aldir, dname + "-current.1.miss.lst")
    refs = [ref_gen8, ref_gen1, ref_prt8, ref_prt1]
    todos = [cur_gen8, cur_gen1, cur_prt8, cur_prt1]
    for f1, f2 in zip(refs, todos):
        shutil.copyfile(f1, f2)
    # empty miss file for family 1 and 8
    open(cur_miss1, "w").close()
    open(cur_miss8, "w").close()
    # create empty concat-file, to check that it is removed
    concat = os.path.join(aldir, dname + "-complete.cat.aln")
    open(concat, "w").close()
    assert not al.align_all_families(prefix, all_fams, ngenomes, dname, quiet, threads)
    # Check output files
    out_mafft1 = os.path.join(aldir, dname + "-mafft-align.1.aln")
    out_btr1 = os.path.join(aldir, dname + "-mafft-prt2nuc.1.aln")
    out_mafft8 = os.path.join(aldir, dname + "-mafft-align.8.aln")
    out_btr8 = os.path.join(aldir, dname + "-mafft-prt2nuc.8.aln")
    exp_mafft1 = os.path.join(EXPPATH, "exp_aldir", "mafft-align.1.aln")
    exp_btr1 = os.path.join(EXPPATH, "exp_aldir", "mafft-prt2nuc.1.aln")
    assert tutil.compare_order_content(out_mafft1, exp_mafft1)
    assert tutil.compare_order_content(out_btr1, exp_btr1)
    assert not os.path.isfile(out_mafft8)
    assert not os.path.isfile(out_btr8)
    assert not os.path.isfile(concat)
    # Check logs
    assert ("Starting alignment of all families: protein alignment, back-translation to "
            "nucleotides, and add missing genomes in the family") in caplog.text
    assert "Checking extractions for family 1" in caplog.text
    assert "Aligning family 1" in caplog.text
    assert "Back-translating family 1" in caplog.text
    assert "Checking extractions for family 8" in caplog.text
    assert ("fam 8: wrong sum of missing genomes (0) and prt extracted (3) for 4 genomes in the "
            "dataset.") in caplog.text
