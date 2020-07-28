#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the pan_to_pergenome submodule in align module
"""
import os
import shutil
import logging

import pytest

import PanACoTA.align_module.pan_to_pergenome as p2p
import test.test_unit.utilities_for_tests as tutil

ALPATH = os.path.join("test", "data", "align")
EXPPATH = os.path.join(ALPATH, "exp_files")
TESTPATH = os.path.join(ALPATH, "test_files")
GENEPATH = os.path.join(ALPATH, "generated_by_unit-tests")


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


ALL_PROTS = {"ESCO1": {"ESCO1_00001": '1',
                       "ESCO1_00002": '4'},
             "ESCO2": {"ESCO2_00001": '1',
                       "ESCO2_22": '2',
                       "ESCO2_456": '4',
                       "ESCO2_46": '3'},
             "ESCO3": {"ESCO3_1": '2',
                       "ESCO3_12": '1',
                       "ESCO3_4564": '3',
                       "ESCO3_00123": '4',
                       "ESCO3_8": '2'},
             "ESCO4": {"ESCO4_00001": '1',
                       "ESCO4_00002": '4',
                       "ESCO4_00003": '3',
                       "ESCO4_00004": '2',
                       "ESCO4_00006": '4'},
             "ESCO5": {"ESCO5_1": '1',
                       "ESCO5_2": '3',
                       "ESCO5_3": '2',
                       "ESCO5_4": '4',
                       "ESCO5_5": '2'},
             "ESCO6": {"ESCO6_1": '4',
                       "ESCO6_2": '3',
                       "ESCO6_3": '1'}}
FAM_GENOMES = {'1': ["ESCO1", "ESCO2", "ESCO3", "ESCO4", "ESCO5", "ESCO6"],
               '2': ["ESCO2", "ESCO3", "ESCO4", "ESCO5"],
               '3': ["ESCO2", "ESCO3", "ESCO4", "ESCO5", "ESCO6"],
               '4': ["ESCO1", "ESCO2", "ESCO3", "ESCO4", "ESCO5", "ESCO6"]}
SEVERAL = {'1': [],
           '2': ["ESCO3", "ESCO5"],
           '3': [],
           '4': ["ESCO4"]}
ALL_GENOMES = ["ESCO1", "ESCO2", "ESCO3", "ESCO4", "ESCO5", "ESCO6"]


def test_get_per_genome(caplog):
    """
    Test that when giving a persistent genome file and a list of genomes,
    it creates all expected files in output/Listdir
    """
    caplog.set_level(logging.DEBUG)
    pers = os.path.join("test", "data", "persgenome", "exp_files", "exp_pers-floor-mixed.txt")
    list_gen = os.path.join(TESTPATH, "listfile.txt")
    dname = "TEST-all-gembase"
    outdir = os.path.join(GENEPATH, "test_get_per_genome")
    all_genomes, aldir, listdir, fams = p2p.get_per_genome(pers, list_gen, dname, outdir)
    assert ("Reading PersGenome and constructing lists of missing genomes "
            "in each family") in caplog.text
    exp_al = os.path.join(outdir, "Align-TEST-all-gembase")
    exp_list = os.path.join(outdir, "List-TEST-all-gembase")
    assert exp_al == aldir
    assert os.path.isdir(aldir)
    assert exp_list == listdir
    assert os.path.isdir(listdir)
    exp_genomes = ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002", "GENO.1216.00003"]
    assert all_genomes == exp_genomes
    exp_fams = ['1', '3', '5', '8', '10', '11', '12']
    assert set(fams) == set(exp_fams)


def test_prot_per_strain():
    """
    Test parser of persistent genome file
    """
    pers = os.path.join("test", "data", "persgenome", "exp_files", "exp_pers-floor-mixed.txt")
    all_prots, fams_genomes, several = p2p.proteins_per_strain(pers)
    exp_several = {'1': ["GENO.1216.00002"],
                   '3': [],
                   '5': [],
                   '8': [],
                   '10': [],
                   '11': [],
                   '12': ["GENO.1216.00003"]}
    assert several == exp_several
    exp_fams = {'1': ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002", "GENO.1216.00003"],
                '3': ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002", "GENO.1216.00003"],
                '5': ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002", "GENO.1216.00003"],
                '8': ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002"],
                '10': ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002"],
                '11': ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002"],
                '12': ["GEN4.1111.00001", "GENO.0817.00001", "GENO.1216.00002", "GENO.1216.00003"]}
    assert fams_genomes == exp_fams
    exp_prots = {"GEN4.1111.00001": {"GEN4.1111.00001.b0001_00001": '1',
                                     "GEN4.1111.00001.b0001_00009": '3',
                                     "GEN4.1111.00001.i0001_00002": '5',
                                     "GEN4.1111.00001.i0001_00007": '8',
                                     "GEN4.1111.00001.i0001_00004": '10',
                                     "GEN4.1111.00001.i0001_00005": '11',
                                     "GEN4.1111.00001.i0001_00008": '12'
                                     },
                 "GENO.0817.00001": {"GENO.0817.00001.b0001_00002": '1',
                                     "GENO.0817.00001.b0002_00011": '3',
                                     "GENO.0817.00001.b0002_00003": '5',
                                     "GENO.0817.00001.i0002_00009": '8',
                                     "GENO.0817.00001.i0002_00004": '10',
                                     "GENO.0817.00001.i0002_00005": '11',
                                     "GENO.0817.00001.i0002_00010": '12'
                                     },
                 "GENO.1216.00002": {"GENO.1216.00002.b0001_00001": '1',
                                     "GENO.1216.00002.i0001_00002": '1',
                                     "GENO.1216.00002.b0002_00010": '3',
                                     "GENO.1216.00002.i0001_00003": '5',
                                     "GENO.1216.00002.b0001_00008": '8',
                                     "GENO.1216.00002.i0001_00005": '10',
                                     "GENO.1216.00002.i0001_00006": '11',
                                     "GENO.1216.00002.b0002_00009": '12'
                                     },
                 "GENO.1216.00003": {"GENO.1216.00003.i0001_00003": '1',
                                     "GENO.1216.00003.i0001_01010": '3',
                                     "GENO.1216.00003.i0080_00010": '5',
                                     "GENO.1216.00003.i0001_00004": '12',
                                     "GENO.1216.00003.i0001_01000": '12'
                                     }
                 }
    assert all_prots == exp_prots


def test_prot_per_strain_member_bis(caplog):
    """
    Test parser of persistent genome file when a same member is in 2 different families
    """
    caplog.set_level(logging.DEBUG)
    pers = os.path.join(TESTPATH, "pers_genome_member-bis.txt")
    all_prots, fams_genomes, several = p2p.proteins_per_strain(pers)
    assert "problem: ESCO2_2 already exists, in family 5. Conflict with family 32" in caplog.text
    exp_several = {'1': [], '5': [], '12': [], '32': []}
    assert several == exp_several
    exp_fams = {'1': ["ESCO_1", "ESCO2", "ESCO3", "ESCO4"], '5': ["ESCO_1", "ESCO2", "ESCO4"],
                '12': ["ESCO_1", "ESCO2", "ESCO3"], '32': ["ESCO_1", "ESCO2", "ESCO3", "ESCO4"]}
    assert fams_genomes == exp_fams
    exp_prots = {"ESCO_1": {"ESCO_1_1": '1', "ESCO_1_2": '5', "ESCO_1_3": '12', "ESCO_1_4": '32'},
                 "ESCO2": {"ESCO2_1": '1', "ESCO2_2": '32', "ESCO2_3": '12'},
                 "ESCO3": {"ESCO3_1": '1', "ESCO3_3": '12', "ESCO3_4": '32'},
                 "ESCO4": {"ESCO4_1": '1', "ESCO4_3": '5', "ESCO4_4": '32'}}
    assert all_prots == exp_prots


def test_get_genomes():
    """
    Test parser of list of genomes
    """
    lstfile = os.path.join("test", "data", "annotate", "exp_files", "results_test_func-default",
                           "LSTINFO-list_genomes-func-test-default.lst")
    all_genomes = p2p.get_all_genomes(lstfile)
    assert all_genomes == ["ESCO.1015.00001", "ESCO.1116.00002", "GENO.1015.00001"]


def test_write_getentry():
    """
    Test that when giving a list of genomes with their persistent gene names,
    it creates all expected files.
    """
    listdir = os.path.join(GENEPATH, "Listdir")
    aldir = os.path.join(GENEPATH, "Aldir")
    # Create align folder
    os.makedirs(listdir)
    dname = "TEST6"
    p2p.write_getentry_files(ALL_PROTS, SEVERAL, listdir, aldir, dname, ALL_GENOMES)
    # Check creation and content of all files
    genfiles = [os.path.join(listdir, "{}-getEntry_gen_ESCO{}.txt".format(dname, num)) for num in
                range(1, 7)]
    expgens = [os.path.join(EXPPATH, "exp_getentry-gen-ESCO{}.txt".format(num)) for num in
               range(1, 7)]
    for fexp, fout in zip(expgens, genfiles):
        print(fexp, fout)
        assert tutil.compare_file_content(fexp, fout)
    prtfiles = [os.path.join(listdir, "{}-getEntry_prt_ESCO{}.txt".format(dname, num)) for num in
                range(1, 7)]
    expprts = [os.path.join(EXPPATH, "exp_getentry-prt-ESCO{}.txt".format(num)) for num in
               range(1, 7)]
    for fexp, fout in zip(expprts, prtfiles):
        assert tutil.compare_file_content(fexp, fout)


def test_write_getentry_error(caplog):
    """
    Test that when giving a list of genomes with their persistent gene names,
    but for 2 genomes, there is no persistent gene, it exists, with an error message
    """
    caplog.set_level(logging.DEBUG)
    all_prots = {"ESCO1": {"ESCO1_00001": '1',
                           "ESCO1_00002": '4'},
                 "ESCO2": {"ESCO2_00001": '1',
                           "ESCO2_22": '2',
                           "ESCO2_456": '4',
                           "ESCO2_46": '3'},
                 "ESCO3": {"ESCO3_1": '2',
                           "ESCO3_12": '1',
                           "ESCO3_4564": '3',
                           "ESCO3_00123": '4',
                           "ESCO3_8": '2'},
                 "ESCO6": {"ESCO6_1": '4',
                           "ESCO6_2": '3',
                           "ESCO6_3": '1'}}
    several = {'1': [],
               '2': ["ESCO3"],
               '3': [],
               '4': []}
    listdir = os.path.join(GENEPATH, "Listdir")
    aldir = os.path.join(GENEPATH, "Aldir")
    # Create align folder
    os.makedirs(listdir)
    dname = "TEST6"
    with pytest.raises(SystemExit):
        p2p.write_getentry_files(all_prots, several, listdir, aldir, dname, ALL_GENOMES)
    assert ("There is not any protein for genome ESCO4 in any family! The program will close, "
            "please fix this problem to be able to run the alignments") in caplog.text
    assert ("There is not any protein for genome ESCO5 in any family! The program will close, "
            "please fix this problem to be able to run the alignments") in caplog.text
    # Check creation and content of all files
    genfiles = [os.path.join(listdir, "{}-getEntry_gen_ESCO{}.txt".format(dname, num)) for num in
                list(range(1, 4)) + [6]]
    expgens = [os.path.join(EXPPATH, "exp_getentry-gen-ESCO{}.txt".format(num)) for num in
               list(range(1, 4)) + [6]]
    for fexp, fout in zip(expgens, genfiles):
        assert tutil.compare_file_content(fexp, fout)
    prtfiles = [os.path.join(listdir, "{}-getEntry_prt_ESCO{}.txt".format(dname, num)) for num in
                list(range(1, 4)) + [6]]
    expprts = [os.path.join(EXPPATH, "exp_getentry-prt-ESCO{}.txt".format(num)) for num in
               list(range(1, 4)) + [6]]
    for fexp, fout in zip(expprts, prtfiles):
        assert tutil.compare_file_content(fexp, fout)


def test_write_genome_prt_exists():
    """
    Test that when only prt file exists, it overwrites it and generates
    expected prt and gen files
    """
    listdir = os.path.join(GENEPATH, "Listdir")
    aldir = os.path.join(GENEPATH, "Aldir")
    # Create align folder
    os.makedirs(listdir)
    dname = "test_write_genome"
    strain = "ESCO4"
    members = ALL_PROTS[strain]

    # Create prt file
    fileprt = os.path.join(listdir, f"{dname}-getEntry_prt_ESCO4.txt")
    with open(fileprt, "w") as prtf:
        prtf.write("Wrong prt file\n")
    p2p.write_genome_file(listdir, aldir, dname, strain, members, SEVERAL)

    # Check creation of files and content
    expprt = os.path.join(EXPPATH, "exp_getentry-prt-ESCO4_write-prt.txt")
    assert tutil.compare_file_content(fileprt, expprt)
    filegen = os.path.join(listdir, f"{dname}-getEntry_gen_ESCO4.txt")
    expgen = os.path.join(EXPPATH, "exp_getentry-gen-ESCO4_write-prt.txt")
    assert tutil.compare_file_content(expgen, filegen)


def test_write_genome_gen_exists():
    """
    Test that when only gen file exists, it overwrites it and generates
    expected prt and gen files
    """
    listdir = os.path.join(GENEPATH, "Listdir")
    aldir = os.path.join(GENEPATH, "Aldir")
    # Create align folder
    os.makedirs(listdir)
    dname = "test_write_genome"
    strain = "ESCO4"
    member4 = ALL_PROTS[strain]
    # Create prt file
    filegen = os.path.join(listdir, f"{dname}-getEntry_gen_ESCO4.txt")
    with open(filegen, "w") as genf:
        genf.write("Wrong gen file\n")
    p2p.write_genome_file(listdir, aldir, dname, strain, member4, SEVERAL)

    # Check creation of files and content
    fileprt = os.path.join(listdir, f"{dname}-getEntry_prt_ESCO4.txt")
    expprt = os.path.join(EXPPATH, "exp_getentry-prt-ESCO4_write-prt.txt")
    assert tutil.compare_file_content(fileprt, expprt)
    expgen = os.path.join(EXPPATH, "exp_getentry-gen-ESCO4_write-prt.txt")
    assert tutil.compare_file_content(filegen, expgen)


def test_write_genome_gen_prt_exist(caplog):
    """
    Test that when gen and prt files already exist, it does not do anything.
    Those files will be used for next steps.
    """
    caplog.set_level(logging.DEBUG)
    listdir = os.path.join(GENEPATH, "Listdir")
    aldir = os.path.join(GENEPATH, "Aldir")
    # Create align folder
    os.makedirs(listdir)
    dname = "TEST6"
    strain = "ESCO4"
    member4 = ALL_PROTS[strain]
    # Create gen and prt files
    filegen = os.path.join(listdir, "{}-getEntry_gen_ESCO4.txt".format(dname))
    with open(filegen, "w") as genf:
        genf.write("Wrong gen file\n")
    fileprt = os.path.join(listdir, "{}-getEntry_prt_ESCO4.txt".format(dname))
    with open(fileprt, "w") as prtf:
        prtf.write("Wrong prt file\n")
    p2p.write_genome_file(listdir, aldir, dname, strain, member4, SEVERAL)

    # Check log
    assert ("For genome ESCO4, "
            "test/data/align/generated_by_unit-tests/Listdir/TEST6-getEntry_prt_ESCO4.txt and "
            "test/data/align/generated_by_unit-tests/Listdir/TEST6-getEntry_gen_ESCO4.txt "
            "already exist. The program will use them to extract "
            "proteins and genes. If you prefer to rewrite them, use option "
            "-F (or --force).".format(fileprt, filegen)) in caplog.text

    # Check content of prt and gen has not changed
    with open(fileprt, "r") as prtf:
        lines = prtf.readlines()
        assert lines == ["Wrong prt file\n"]
    with open(filegen, "r") as prtf:
        lines = prtf.readlines()
        assert lines == ["Wrong gen file\n"]


def test_write_genome():
    """
    Test that given a genome, it writes the list of its proteins
    and genes in expected files.
    """
    listdir = os.path.join(GENEPATH, "Listdir")
    aldir = os.path.join(GENEPATH, "Aldir")
    # Create align folder
    os.makedirs(listdir)
    dname = "test_write_genome"
    strain = "ESCO4"
    members = ALL_PROTS[strain]
    p2p.write_genome_file(listdir, aldir, dname, strain, members, SEVERAL)

    # Check creation of files and content
    fileprt = os.path.join(listdir, f"{dname}-getEntry_prt_ESCO4.txt")
    expprt = os.path.join(EXPPATH, "exp_getentry-prt-ESCO4_write-prt.txt")
    assert tutil.compare_file_content(fileprt, expprt)
    filegen = os.path.join(listdir, f"{dname}-getEntry_gen_ESCO4.txt")
    expgen = os.path.join(EXPPATH, "exp_getentry-gen-ESCO4_write-prt.txt")
    assert tutil.compare_file_content(filegen, expgen)


def test_write_missing():
    """
    Test that given families with genomes present, genomes with several numbers and
    list of all genomes, it returns, for each family, the genomes which will not
    be considered.
    """
    dname = "test_write_missing"
    p2p.write_missing_genomes(FAM_GENOMES, SEVERAL, ALL_GENOMES, GENEPATH, dname)

    exp_res = [None, [], ["ESCO1", "ESCO3", "ESCO5", "ESCO6"], ["ESCO1"], ["ESCO4"]]
    for num in range(1,5):
        miss_file = os.path.join(GENEPATH, f"{dname}-current.{num}.miss.lst")
        assert os.path.isfile(miss_file)
        assert tutil.compare_file_to_list(miss_file, exp_res[num])
