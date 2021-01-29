#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import PanACoTA.utils as utils
import test.test_unit.utilities_for_tests as utilities
import pytest
import os
import logging
import shutil
import matplotlib

matplotlib.use('AGG')

# Define variables used by several tests
DATA_DIR = os.path.join("test", "data", "annotate")
TEST_DIR = os.path.join(DATA_DIR, "test_files")
BASELINE_DIR = os.path.abspath(os.path.join(DATA_DIR, "exp_files", "baseline"))
GENEPATH = os.path.join(DATA_DIR, "generated_by_unit-tests")
LOGFILE_BASE = "test_prokka"
LOGFILES = [LOGFILE_BASE + ext for ext in [".log", ".log.debug", ".log.details", ".log.err"]]

@pytest.fixture(autouse=True)
def setup_teardown_module():
    """
    Remove log files at the end of this test module
    """
    # Init logger to level detail (15)
    # utils.init_logger(LOGFILE_BASE, logging.DEBUG, 'test_utils', verbose=1)
    os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH)
    print("teardown")


# Start tests
def test_check_install():
    """
    Try to run prokka, which is installed, and check that there is no problem
    """
    assert utils.check_installed("prokka")


def test_check_install_error():
    """
    Try to run a command which does not exist, and check that it closes the program
    with exit code 1
    """
    assert not utils.check_installed("plop false command...")


@pytest.mark.mpl_image_compare(baseline_dir=BASELINE_DIR, tolerance=5, backend="agg")
def test_plot_dist():
    """
    Plot a given distribution, and check that output is as expected
    """
    logger = logging.getLogger("test_utils")
    values = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 10]
    limit = 3
    res_dir = os.path.join("test", "data", "annotate")
    # reffile = os.path.join("test", "data", "annotate", "exp_files", "res_plot_distr.png")
    title = "Distribution test"
    text = "Max L90 ="
    myfig = utils.plot_distr(values, limit, title, text, logger)
    return myfig


def test_skipped_prokka_prodigal(caplog):
    """
    Test that when the list of skipped genomes (because of prokka run) is not empty,
    it writes the right message.
    """
    caplog.set_level(logging.DEBUG)
    skipped = ["toto", "genome", "genome2"]
    utils.write_warning_skipped(skipped)
    assert ("prokka had problems while annotating some genomes, or did not "
            "find any gene. Hence, they are not "
            "formatted, and absent from your output database. Please look at the "
            "current error log (<output_directory>/PanACoTA-annotate_list_genomes[-date].log.err) "
            "to get more information on the problems. "
            "Here are those genomes:") in caplog.text
    assert "\n\t- toto\n\t- genome\n\t- genome2" in caplog.text
    utils.write_warning_skipped(skipped, prodigal_only=True)
    assert ("prodigal had problems while annotating some genomes, or did not "
            "find any gene. Hence, they are not "
            "formatted, and absent from your output database. Please look at the "
            "current error log (<output_directory>/PanACoTA-annotate_list_genomes[-date].log.err) "
            "to get more information on the problems. "
            "Here are those genomes:") in caplog.text
    assert "\n\t- toto\n\t- genome\n\t- genome2" in caplog.text


def test_skipped_format(caplog):
    """
    Test that when the list of skipped genomes (format step could not run) is not empty,
    it writes the right message.
    """
    caplog.set_level(logging.DEBUG)
    skipped_format = ["toto", "genome", "genome2"]
    utils.write_warning_skipped(skipped_format, do_format=True)
    assert ("Some genomes were annotated by prokka, but could not be formatted, "
            "and are hence absent from your output database. Please look at "
            "'<output_directory>/PanACoTA-annotate_list_genomes[-date].log.err' and .details "
            "files to get more information about why they could not be formatted") in caplog.text
    assert ("\n\t- toto\n\t- genome\n\t- genome2\n") in caplog.text


def test_write_discarded(caplog):
    """
    Test that the list of discarded genomes is written as expected.
    """
    caplog.set_level(logging.DEBUG)
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join(DATA_DIR, "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], gpaths[4], 9876546, 6, 2]
               }
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3"]
    list_file = os.path.join("titi", "toto", "list_genomes.txt")
    utils.write_genomes_info(genomes, kept_genomes, list_file, GENEPATH, qc=False)
    outfile = os.path.join(GENEPATH, "discarded-list_genomes.lst")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_test_write_discard.lst")
    assert "2 genomes were discarded" in caplog.text
    # There is no order in the discarded file. So, just check that the lines
    # written are as expected.
    assert utilities.compare_file_content(outfile, exp_file)


def test_write_discarded_1genome(caplog):
    """
    Test that the list of discarded genomes is written as expected.
    """
    caplog.set_level(logging.DEBUG)
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join(DATA_DIR, "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], gpaths[4], 9876546, 6, 2]
               }
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3", "genome2"]
    list_file = os.path.join("titi", "toto", "list_genomes.txt")
    utils.write_genomes_info(genomes, kept_genomes, list_file, GENEPATH, qc=False)
    outfile = os.path.join(GENEPATH, "discarded-list_genomes.lst")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_test_write_discard_1genome.lst")
    assert "1 genome was discarded" in caplog.text
    # There is no order in the discarded file. So, just check that the lines
    # written are as expected.
    assert utilities.compare_file_content(outfile, exp_file)


def test_write_discarded_qc(caplog):
    """
    Test that the list with information on all genomes when we run with QC only is
    written as expected
    """
    caplog.set_level(logging.DEBUG)
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join(DATA_DIR, "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], gpaths[4], 9876546, 6, 2]
               }
    kept_genomes = []
    list_file = os.path.join(DATA_DIR, "list_genomes.txt")
    utils.write_genomes_info(genomes, kept_genomes, list_file, GENEPATH, qc=True)
    outfile = os.path.join(GENEPATH, "ALL-GENOMES-info-list_genomes.lst")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_test_write_info_qc.lst")
    assert ("Writing information on genomes in "
            "test/data/annotate/generated_by_unit-tests/ALL-GENOMES-info-list_genomes.lst") in caplog.text
    # There is no order in the discarded file. So, just check that the lines
    # written are as expected.
    assert utilities.compare_file_content(outfile, exp_file)


def test_write_discarded_empty():
    """
    Test that when the list of genomes is empty, but the list of kept-genomes is not
    (should never happen...), it writes only the header of discarded lst file.
    """
    genomes = {}
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3"]
    list_file = os.path.join(DATA_DIR, "list_genomes.txt")
    utils.write_genomes_info(genomes, kept_genomes, list_file, GENEPATH)
    outfile = os.path.join(GENEPATH, "discarded-list_genomes.lst")
    with open(outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "orig_name\tto_annotate\tgsize\tnb_conts\tL90\n"


def test_write_discarded_all_kept(caplog):
    """
    Test that when all genomes are kept, the discarded lst file only contains the
    header line.
    """
    caplog.set_level(logging.DEBUG)
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join(DATA_DIR, "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], gpaths[4], 9876546, 6, 2]
               }
    kept_genomes = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome3", "genome2", "genome1"]
    list_file = os.path.join(DATA_DIR, "list_genomes.txt")
    utils.write_genomes_info(genomes, kept_genomes, list_file, GENEPATH)
    outfile = os.path.join(GENEPATH, "discarded-list_genomes.lst")
    with open(outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "orig_name\tto_annotate\tgsize\tnb_conts\tL90\n"
    assert "0 genome was discarded" in caplog.text


def test_write_lstinfo():
    """
    Test that lstinfo file is written as expected.
    """
    gnames = ["H299_H561.fasta", "B2_A3_5.fasta-problems", "genome1", "genome2", "genome3"]
    gpaths = [os.path.join(DATA_DIR, "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["toto.0417", gpaths[0], gpaths[0], 12656, 3, 1],
               gnames[1]: ["toto.0417", gpaths[1], gpaths[1], 456464645, 5, 1],
               gnames[2]: ["toto.0417", gpaths[2], gpaths[2], 4564855, 156, 40],
               gnames[3]: ["toto.0417", gpaths[3], gpaths[3], 6549, 16, 8],
               gnames[4]: ["toto.0417", gpaths[4], gpaths[4], 9876546, 6, 2]
               }
    list_file = os.path.join("toto", "list_genomes.txt")
    utils.write_lstinfo(list_file, genomes, GENEPATH)
    outfile = os.path.join(GENEPATH, "LSTINFO-list_genomes.lst")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_test_write_lstinfo.lst")
    assert utilities.compare_order_content(outfile, exp_file)


def test_write_lstinfo_nogenome():
    """
    Test that when there is no genome fully annotated, lstinfo contains
    only header.
    """
    genomes = {}
    list_file = os.path.join("toto", "list_genomes.txt")
    utils.write_lstinfo(list_file, genomes, GENEPATH)
    outfile = os.path.join(GENEPATH, "LSTINFO-list_genomes.lst")
    with open(outfile, "r") as outf:
        all_lines = outf.readlines()
        assert len(all_lines) == 1
        assert all_lines[0] == "gembase_name\torig_name\tto_annotate\tgsize\tnb_conts\tL90\n"


def test_sort_gene():
    """
    Test that genomes are sorted by species first, and then by strain number.
    """
    # genomes = {genome_orig, [gembase, path, gsize, nbcont, L90]}
    genomes = ["genome.0417.00010", "toto.0417.00010", "genome1.0417.00002",
               "genome.0417.00015", "totn.0417.00010", "genome.0417.00009",
               "genome.0517.00001", "toto.0417.00011"]
    sorted_genomes = sorted(genomes, key=utils.sort_genomes_by_name)
    exp = ["genome.0517.00001", "genome.0417.00009", "genome.0417.00010",
           "genome.0417.00015", "genome1.0417.00002", "totn.0417.00010",
           "toto.0417.00010", "toto.0417.00011"]
    assert sorted_genomes == exp


def test_sort_gene_tuple():
    """
    Test that genomes are sorted by species first, and then by strain number.
    """
    # genomes = {genome_orig, [gembase, path, gsize, nbcont, L90]}
    genomes = {"name1": ["genome.0417.00010", "path/to/genome", 123456, 50, 3],
               "name2": ["toto.0417.00010", "path/to/genome", 123456, 50, 3],
               "name3": ["genome1.0417.00002", "path/to/genome", 123456, 50, 3],
               "name4": ["genome.0417.00015", "path/to/genome", 123456, 50, 3],
               "name5": ["totn.0417.00010", "path/to/genome", 123456, 50, 3],
               "name6": ["genome.0417.00009", "path/to/genome", 123456, 50, 3],
               "name7": ["genome.0517.00001", "path/to/genome", 123456, 50, 3],
               "name8": ["toto.0417.00011", "path/to/genome", 123456, 50, 3], }
    sorted_genomes = sorted(genomes.items(), key=utils.sort_genomes_by_name)
    exp = [("name7", genomes["name7"]), ("name6", genomes["name6"]),
           ("name1", genomes["name1"]), ("name4", genomes["name4"]),
           ("name3", genomes["name3"]), ("name5", genomes["name5"]),
           ("name2", genomes["name2"]), ("name8", genomes["name8"])]
    assert sorted_genomes == exp


def test_sort_gene_noformat():
    """
    Test that when genomes are not in the gembase format, it returns
    genomes in alphabetical order
    """
    # genomes = {genome_orig, [gembase, path, gsize, nbcont, L90]}
    genomes = {"name1": ["genome", "path/to/genome", 123456, 50, 3],
               "name2": ["toto", "path/to/genome", 123456, 50, 3],
               "name3": ["genome1   ", "path/to/genome", 123456, 50, 3],
               "name4": ["genome1bis", "path/to/genome", 123456, 50, 3],
               "name5": ["mygenome", "path/to/genome", 123456, 50, 3],
               "name6": ["agenome_nogembase", "path/to/genome", 123456, 50, 3],
               "name7": ["agenome_nogembase2", "path/to/genome", 123456, 50, 3],
               "name8": ["toto.0417.00011", "path/to/genome", 123456, 50, 3], }
    sorted_genomes = sorted(genomes.items(), key=utils.sort_genomes_by_name)
    exp = [("name6", genomes["name6"]), ("name7", genomes["name7"]),
           ("name1", genomes["name1"]), ("name3", genomes["name3"]),
           ("name4", genomes["name4"]), ("name5", genomes["name5"]),
           ("name2", genomes["name2"]), ("name8", genomes["name8"])]
    assert sorted_genomes == exp


def test_sort_genome_l90():
    """
    Test that when genomes are not in the gembase format, it returns
    genomes in alphabetical order
    """
    # genomes = {genome_orig, [gembase, path, gsize, nbcont, L90]}
    genomes = {"name1": ["genome", "path/to/genome", 123456, 1, 2],
               "name2": ["toto", "path/to/genome", 123456, 100, 4],
               "name3": ["genome1   ", "path/to/genome", 123456, 10, 2],
               "name4": ["genome1bis", "path/to/genome", 123456, 11, 2],
               "name5": ["mygenome", "path/to/genome", 123456, 8, 3],
               "name6": ["agenome_nogembase", "path/to/genome", 123456, 3, 1],
               "name7": ["agenome_nogembase2", "path/to/genome", 123456, 4, 1],
               "name8": ["toto.0417.00011", "path/to/genome", 123456, 50, 5], }
    sorted_genomes = sorted(genomes.items(), key=utils.sort_genomes_l90_nbcont)
    exp = [("name6", genomes["name6"]), ("name7", genomes["name7"]),
           ("name1", genomes["name1"]), ("name3", genomes["name3"]),
           ("name4", genomes["name4"]), ("name5", genomes["name5"]),
           ("name2", genomes["name2"]), ("name8", genomes["name8"])]
    assert sorted_genomes == exp


def test_sort_proteins():
    """
    Test that when a list of proteins is given, it returns the list sorted by :
    - species
    - strain number
    - protein number
    """
    proteins = ["ESCO.0417.00010.i0001_12354", "ESCA.0617.00001.i0001_005",
                "ESCO.0517.00001.i0001_12354", "ESCO.0417.00010.i0001_1",
                "TOTO.0417.00001.i0001_12", "TOTO.0717.00001.i0008_9"]
    sorted_prot = sorted(proteins, key=utils.sort_proteins)
    exp = ["ESCA.0617.00001.i0001_005", "ESCO.0517.00001.i0001_12354",
           "ESCO.0417.00010.i0001_1", "ESCO.0417.00010.i0001_12354",
           "TOTO.0717.00001.i0008_9", "TOTO.0417.00001.i0001_12"]
    assert sorted_prot == exp


def test_sort_proteins_other_format():
    """
    Test that when a list of proteins is given, it returns the list sorted by :
    - species
    - strain number
    - protein number
    for all proteins where this information is available,
    and then, the other proteins, sorted by :
    - species (string before "_")
    - protein number (num after "_")
    """
    proteins = ["ESCO.0417.00010.i0001_12354", "ESCO.0617.00001.i0001_005", "my_protein_0002",
                "ESCO.0517.00001.i0001_12354", "my_protein_0001", "ESCO.0417.00010.i0001_1",
                "a_my_prot_15", "a_my_prot_2",
                "TOTO.0417.00001.i0001_12", "TOTO.0717.00001.i0008_9"]
    sorted_prot = sorted(proteins, key=utils.sort_proteins)
    exp = ["ESCO.0617.00001.i0001_005", "ESCO.0517.00001.i0001_12354",
           "ESCO.0417.00010.i0001_1", "ESCO.0417.00010.i0001_12354",
           "TOTO.0717.00001.i0008_9", "TOTO.0417.00001.i0001_12",
           "a_my_prot_2", "a_my_prot_15", "my_protein_0001", "my_protein_0002"]
    assert sorted_prot == exp


def test_sort_proteins_error_format1(caplog):
    """
    Test that when a protein name does not follow the format <alpha_num>_<num>,
    it gives an error.
    """
    caplog.set_level(logging.DEBUG)
    proteins = ["ESCO.0417.00010.i0001_12354", "ESCO.0617.00001.i0001_005",
                "ESCO.0517.00001.i0001_12354", "error-protein"]
    with pytest.raises(SystemExit):
        sorted(proteins, key=utils.sort_proteins)
    assert ("ERROR: Protein error-protein does not have the required format. It must contain, "
            "at least <alpha-num>_<num_only>, and at best "
            "<name>.<date>.<strain_num>.<contig_info>_<prot_num>. "
            "Please change its name.") in caplog.text


def test_read_genomes_nofile(caplog):
    """
    Test that when the genome list file provided does not exist, it
    ends the program with an error message
    """
    logger = logging.getLogger("default")
    with pytest.raises(SystemExit):
        utils.read_genomes("toto.txt", "TOTO", "0417", "db/path", "tmppath", logger)
    assert "ERROR: Your list file " in caplog.text
    assert "toto.txt" in caplog.text
    assert "does not exist. Please provide a list file." in caplog.text
    assert "Ending program." in caplog.text


def test_read_genomes_wrongname():
    """
    Test that when the list file contains only genome names which do not exist,
    it returns an empty list of genomes to annotate/format.
    """
    logger = logging.getLogger("default")
    name = "ESCO"
    date = "0417"
    dbpath = os.path.join(DATA_DIR, "genomes")
    tmppath = "tmppath"
    list_file = os.path.join(DATA_DIR, "test_files", "list_genomes-wrongNames.txt")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath, logger)
    assert genomes == {}


def test_read_genomes_ok(caplog):
    """
    Test that when the list file contains genomes existing, it returns the expected list
    of genomes
    """
    logger = logging.getLogger("default")
    name = "ESCO"
    date = "0417"
    dbpath = os.path.join(DATA_DIR, "genomes")
    tmppath = "tmppath"
    list_file = os.path.join(DATA_DIR, "test_files", "list_genomes.lst")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath, logger)
    exp = {"A_H738.fasta": ["ESCO.0417"],
           "B2_A3_5.fasta-split5N.fna-short-contig.fna": ["ESCO.0417"],
           "H299_H561.fasta": ["ABCD.0417"], "genome2.fasta": ["TOTO.0417"],
           "genome3.fasta": ["ESCO.0512"], "genome4.fasta": ["TOTO.0417"],
           "genome5.fasta": ["TOTO.0114"]}
    assert exp == genomes
    assert "genome.fst genome file does not exist. It will be ignored." in caplog.text


def test_read_genomes_errors(caplog):
    """
    Test that when the list file contains errors in name and date provided,
    it returns the expected errors, and the expected genome list.
    """
    logger = logging.getLogger("default")
    name = "ESCO"
    date = "0417"
    dbpath = os.path.join(DATA_DIR, "genomes")
    tmppath = "tmppath"
    list_file = os.path.join(DATA_DIR, "test_files", "list_genomes-errors.txt")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath, logger)
    exp = {"A_H738.fasta": ["ESCO.0417"],
           "B2_A3_5.fasta-split5N.fna-short-contig.fna": ["ESCO.0417"],
           "H299_H561.fasta": ["ABCD.0417"], "genome2.fasta": ["TOTO.0417"],
           "genome3.fasta": ["ESCO.0512"], "genome4.fasta": ["ESCO.0417"],
           "genome5.fasta": ["ESCO.0417"]}
    assert genomes == exp
    assert ("Invalid name/date given for genome A_H738.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default name (ESCO) and date (0417) will "
            "be used.") in caplog.text
    assert ("Invalid name abc given for genome B2_A3_5.fasta-split5N.fna-short-contig.fna. Only put "
           "4 alphanumeric characters in your date and name. For "
           "this genome, the default name (ESCO) will "
           "be used.") in caplog.text
    assert ("Invalid date 152 given for genome H299_H561.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default date (0417) will "
            "be used.") in caplog.text
    assert ("Invalid date 1-03 given for genome genome2.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default date (0417) will "
            "be used.") in caplog.text
    assert ("genome.fst genome file does not exist. "
            "It will be ignored.") in caplog.text
    assert ("Invalid name a/b2 given for genome genome3.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default name (ESCO) will "
            "be used.") in caplog.text
    assert ("Invalid name #esc given for genome genome5.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default name (ESCO) will "
            "be used.") in caplog.text
    assert ("Invalid date 1_14 given for genome genome5.fasta. Only put "
            "4 alphanumeric characters in your date and name. For "
            "this genome, the default date (0417) will "
            "be used.") in caplog.text


def test_read_genomes_multi_files(caplog):
    """
    Test that when the list file contains several filenames for 1 same genome,
    it returns the expected genome list, the expected errors (when some genome
    files do not exist) and the expected concatenated files.
    """
    logger = logging.getLogger("default")
    name = "ESCO"
    date = "0417"
    tmppath = os.path.join(GENEPATH, "tmppath")
    os.mkdir(tmppath)
    dbpath = os.path.join(DATA_DIR, "genomes")
    list_file = os.path.join(DATA_DIR, "test_files", "list_genomes-multi-files.txt")
    genomes = utils.read_genomes(list_file, name, date, dbpath, tmppath, logger)
    exp = {"A_H738.fasta-all.fna": ["ESCO.0417"],
           "H299_H561.fasta-all.fna": ["ABCD.0417"], "genome2.fasta": ["TOTO.0417"],
           "genome3.fasta": ["ESCO.0512"], "genome4.fasta": ["TOTO.0417"],
           "genome5.fasta": ["TOTO.0114"]}
    assert exp == genomes
    # Check error messages
    assert ("genome.fna genome file does not exist. Its file will be ignored "
            "when concatenating ['H299_H561.fasta', 'genome6.fasta', 'genome.fna']") in caplog.text
    assert "genome.fst genome file does not exist. It will be ignored." in caplog.text
    assert ("toto.fst genome file does not exist. Its file will be ignored "
            "when concatenating ['toto.fst', 'toto.fasta', 'genome.fst']") in caplog.text
    assert ("toto.fasta genome file does not exist. Its file will be ignored "
            "when concatenating ['toto.fst', 'toto.fasta', 'genome.fst']") in caplog.text
    assert ("genome.fst genome file does not exist. Its file will be ignored "
            "when concatenating ['toto.fst', 'toto.fasta', 'genome.fst']") in caplog.text
    assert ("None of the genome files in ['toto.fst', 'toto.fasta', 'genome.fst'] exist. "
            "This genome will be ignored.") in caplog.text
    # Check that files were concatenated as expected
    concat1 = os.path.join(tmppath, "A_H738.fasta-all.fna")
    exp_concat1 = os.path.join(dbpath, "A_H738-and-B2_A3_5.fna")
    concat2 = os.path.join(tmppath, "H299_H561.fasta-all.fna")
    exp_concat2 = os.path.join(dbpath, "H299_H561-and-genome6.fna")
    assert utilities.compare_order_content(concat1, exp_concat1)
    assert utilities.compare_order_content(concat2, exp_concat2)


def test_read_genomes_info_nofile(caplog):
    """
    Read lstinfo file and get all genomes information
    Check that when the file does not exist, it exits with appropriate error message
    """
    caplog.set_level(logging.DEBUG)
    name = "ESCO"
    with pytest.raises(SystemExit):
        utils.read_genomes_info("toto.txt", name)
    assert 'Reading given information on your genomes in toto.txt' in caplog.text
    assert ("ERROR: The info file toto.txt that you gave does not exist. "
            "Please provide the right path/name for this file.") in caplog.text
    assert "Ending program" in caplog.text


def test_read_genomes_info_wrongheader(caplog):
    """
    Read lstinfo file and get all genomes information
    When wrong header (not all required columns), exits and appropriate error message
    """
    caplog.set_level(logging.DEBUG)
    name = "ESCO"
    lstinfo_file = os.path.join(DATA_DIR, "test_files", "lstinfo-wrong-header.lst")
    with pytest.raises(SystemExit):
        utils.read_genomes_info(lstinfo_file, name)
    assert ("Reading given information on your genomes in "
            "test/data/annotate/test_files/lstinfo-wrong-header.lst") in caplog.text
    assert ("ERROR: It seems that your info file test/data/annotate/test_files/lstinfo-wrong-header.lst "
            "does not have a header, "
            "or this header does not have, at least, the required columns tab separated: ") in caplog.text
    assert("to_annotate, gsize nb_conts and L90 (in any order).\nEnding program.") in caplog.text


def test_read_genomes_info_noheader(caplog):
    """
    Read lstinfo file and get all genomes information
    When no header, exits and appropriate error message
    """
    caplog.set_level(logging.DEBUG)
    name = "ESCO"
    lstinfo_file = os.path.join(DATA_DIR, "test_files", "lstinfo-no-header.lst")
    with pytest.raises(SystemExit):
        utils.read_genomes_info(lstinfo_file, name)
    assert ("Reading given information on your genomes in "
            "test/data/annotate/test_files/lstinfo-no-header.lst") in caplog.text
    assert ("ERROR: It seems that your info file test/data/annotate/test_files/lstinfo-no-header.lst "
            "does not have a header, "
            "or this header does not have, at least, the required columns tab separated: ") in caplog.text
    assert("to_annotate, gsize nb_conts and L90 (in any order).\nEnding program.") in caplog.text


def test_read_genomes_info_not_int(caplog):
    """
    Read lstinfo file and get all genomes information
    When a nbcont column is not an int for 1 genome, writes error message and ignores the genome
    """
    caplog.set_level(logging.DEBUG)
    name = "ESCO"
    lstinfo_file = os.path.join(DATA_DIR, "test_files", "lstinfo-not-int-nbcont.lst")
    genomes = utils.read_genomes_info(lstinfo_file, name)
    assert ("Reading given information on your genomes in "
            "test/data/annotate/test_files/lstinfo-not-int-nbcont.lst") in caplog.text
    assert 'Found 2 genomes in total' in caplog.text
    assert ("For genome A_H738-and-B2_A3_5, at least one of your columns 'gsize', "
            "'nb_conts' or 'L90' contains a non numeric value. This genome will be ignored") in caplog.text
    exp = {"genome1.fasta":
            ["genome1", "test/data/annotate/genomes/genome1.fasta",
             "test/data/annotate/genomes/genome1.fasta", 800, 6, 5],
           "genome7.fasta":
            ["genome7", "test/data/annotate/genomes/genome7.fasta",
             "test/data/annotate/genomes/genome7.fasta", 79705, 80, 65]}
    assert genomes == exp


def test_read_genomes_info_not_all_filled(caplog):
    """
    Read lstinfo file and get all genomes information
    1 column not filled for at least 1 genome: exits and write appropriate error message
    """
    caplog.set_level(logging.DEBUG)
    name = "ESCO"
    lstinfo_file = os.path.join(DATA_DIR, "test_files", "lstinfo-not-all-filled.lst")
    with pytest.raises(SystemExit):
        utils.read_genomes_info(lstinfo_file, name)
    assert ("Reading given information on your genomes in "
            "test/data/annotate/test_files/lstinfo-not-all-filled.lst") in caplog.text
    assert ("ERROR: Check that all fields of test/data/annotate/test_files/lstinfo-not-all-filled.lst "
            "are filled in each line (can be 'NA')") in caplog.text


def test_read_genomes_info_no_path(caplog):
    """
    Read lstinfo file and get all genomes information
    When a genome in lstinfo does not have its corresponding file, write appropriate error message
    and ignores genome
    """
    caplog.set_level(logging.DEBUG)
    name = "ESCO"
    lstinfo_file = os.path.join(DATA_DIR, "test_files", "lstinfo-miss-1genome.lst")
    genomes = utils.read_genomes_info(lstinfo_file, name)
    assert ("data/annotate/genomes/A_H738-and-B2_A3_5.fna genome file does not exist. "
            "This genome will be ignored") in caplog.text
    assert 'Found 2 genomes in total' in caplog.text
    exp = {"genome1.fasta":
            ["genome1", "test/data/annotate/genomes/genome1.fasta",
             "test/data/annotate/genomes/genome1.fasta", 800, 6, 5],
           "genome7.fasta":
            ["genome7", "test/data/annotate/genomes/genome7.fasta",
             "test/data/annotate/genomes/genome7.fasta", 79705, 80, 65]}
    assert genomes == exp


def test_read_genomes_info_no_genomes(caplog):
    """
    Read lstinfo file and get all genomes information
    When no genome in lstinfo correspond to existing paths, exits and appropriate error message
    """
    caplog.set_level(logging.DEBUG)
    name = "ESCO"
    lstinfo_file = os.path.join(DATA_DIR, "test_files", "lstinfo-no-genome.lst")
    with pytest.raises(SystemExit):
        utils.read_genomes_info(lstinfo_file, name)
    assert ("Reading given information on your genomes in "
            "test/data/annotate/test_files/lstinfo-no-genome.lst") in caplog.text
    assert ("No genome listed in test/data/annotate/test_files/lstinfo-no-genome.lst "
            "was found.") in caplog.text


def test_read_genomes_info_ok(caplog):
    """
    Read lstinfo file and get all genomes information and returns it as expected
    """
    caplog.set_level(logging.DEBUG)
    list_file = os.path.join(DATA_DIR, "test_files", "lstinfo.lst")
    name = "ESCO"
    genomes = utils.read_genomes_info(list_file, name)
    assert ("Reading given information on your genomes in "
            "test/data/annotate/test_files/lstinfo.lst") in caplog.text
    assert 'Found 3 genomes in total' in caplog.text
    exp = {"genome1.fasta":
            ["genome1", "test/data/annotate/genomes/genome1.fasta",
             "test/data/annotate/genomes/genome1.fasta", 800, 6, 6],
           "A_H738-and-B2_A3_5.fna":
            ["A_H738-and-B2_A3_5", "test/data/annotate/genomes/A_H738-and-B2_A3_5.fna",
             "test/data/annotate/genomes/A_H738-and-B2_A3_5.fna", 7000, 78, 5],
           "genome7.fasta":
            ["genome7", "test/data/annotate/genomes/genome7.fasta",
             "test/data/annotate/genomes/genome7.fasta", 79705, 80, 65]}
    assert genomes == exp


def test_read_genomes_info_date_ok(caplog):
    """
    Read lstinfo file and get all genomes information
    """
    caplog.set_level(logging.DEBUG)
    list_file = os.path.join(DATA_DIR, "test_files", "lstinfo.lst")
    name = "ESCO"
    date = "0720"
    genomes = utils.read_genomes_info(list_file, name, date=date)
    assert ("Reading given information on your genomes in "
            "test/data/annotate/test_files/lstinfo.lst") in caplog.text
    assert 'Found 3 genomes in total' in caplog.text
    exp = {"test/data/annotate/genomes/genome1.fasta":
            ["ESCO.0720", "test/data/annotate/genomes/genome1.fasta",
             "test/data/annotate/genomes/genome1.fasta", 800, 6, 6],
           "test/data/annotate/genomes/A_H738-and-B2_A3_5.fna":
            ["ESCO.0720", "test/data/annotate/genomes/A_H738-and-B2_A3_5.fna",
             "test/data/annotate/genomes/A_H738-and-B2_A3_5.fna", 7000, 78, 5],
           "test/data/annotate/genomes/genome7.fasta":
            ["ESCO.0720", "test/data/annotate/genomes/genome7.fasta",
             "test/data/annotate/genomes/genome7.fasta", 79705, 80, 65]}
    assert genomes == exp


def test_check_format():
    """
    test that format is ok or not
    """
    assert utils.check_format("ESC1")
    assert utils.check_format("1234")
    assert not utils.check_format("12")
    assert not utils.check_format("ESC*")


def test_check_resdirlst(caplog):
    """
    Test that when the result directory already contains .lst files in LSTINFO,
    program ends with an error message.
    """
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(GENEPATH, "LSTINFO"))
    open(os.path.join(GENEPATH, "LSTINFO", "toto.lst"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(GENEPATH)
    assert ("ERROR: Your output directory already has .lst files in the "
            "LSTINFO folder. Provide another result directory, or remove the "
            "files in this one.") in caplog.text


def test_check_resdirprt(caplog):
    """
    Test that when the result directory already contains .prt files in Proteins,
    program ends with an error message.
    """
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(GENEPATH, "Proteins"))
    open(os.path.join(GENEPATH, "Proteins", "toto.prt"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(GENEPATH)
    assert ("ERROR: Your output directory already has .prt files in the "
            "Proteins folder. Provide another result directory, or remove the "
            "files in this one.") in caplog.text


def test_check_resdirgen(caplog):
    """
    Test that when the result directory already contains .gen files in Genes,
    program ends with an error message.
    """
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(GENEPATH, "Genes"))
    open(os.path.join(GENEPATH, "Genes", "toto.gen"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(GENEPATH)
    assert ("ERROR: Your output directory already has .gen files in the "
            "Genes folder. Provide another result directory, or remove the "
            "files in this one.") in caplog.text


def test_check_resdirrep(caplog):
    """
    Test that when the result directory already contains .fna files in Replicons,
    program ends with an error message.
    """
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(GENEPATH, "Replicons"))
    open(os.path.join(GENEPATH, "Replicons", "toto.fna"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(GENEPATH)
    assert ("ERROR: Your output directory already has .fna files in the "
            "Replicons folder. Provide another result directory, or remove the "
            "files in this one.") in caplog.text


def test_check_resdirgff(caplog):
    """
    Test that when the result directory already contains .fna files in Replicons,
    program ends with an error message.
    """
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(GENEPATH, "gff3"))
    open(os.path.join(GENEPATH, "gff3", "toto.gff"), "w").close()
    with pytest.raises(SystemExit):
        utils.check_out_dirs(GENEPATH)
    assert ("ERROR: Your output directory already has .gff files in the "
            "gff3 folder. Provide another result directory, or remove the "
            "files in this one.") in caplog.text


def test_check_resdirotherext():
    """
    Test that when the result directory contains txt files in Replicons dir,
    there is no problem.
    """
    # Create output directory with a lst file in LSTINFO
    os.makedirs(os.path.join(GENEPATH, "Replicons"))
    open(os.path.join(GENEPATH, "Replicons", "toto.txt"), "w").close()
    utils.check_out_dirs(GENEPATH)


def test_check_resdirnodir():
    """
    Test that when the result directory does not already exist, there is no problem.
    """
    utils.check_out_dirs("totoresdir")


def test_run_cmd_error_noquit(caplog):
    """
    Test that when we try to run a command which does not exist, it returns an error message,
    but does not exit the program (eof=False).
    """
    cmd = "toto"
    error = "error trying to run toto"
    assert utils.run_cmd(cmd, error) == 1
    assert "error: command '>toto' is not possible" in caplog.text


def test_run_cmd_error_noquit_logger(caplog):
    """
    Test that when we try to run a command which does not exist, it returns an error message,
    but does not exit the program (eof=False). With a given logger where error is written.
    """
    cmd = "toto"
    logger = logging.getLogger("default")
    error = "error trying to run toto"
    assert utils.run_cmd(cmd, error, logger=logger) == 1
    assert "error: command '>toto' is not possible" in caplog.text


def test_run_cmd_error_quit(caplog):
    """
    Test that when we try to run a command which does not exist, it returns an error message,
    and exits the program (eof=True)
    """
    cmd = "toto"
    error = "error trying to run toto"
    with pytest.raises(SystemExit):
        utils.run_cmd(cmd, error, eof=True)
    assert "command '>toto' is not possible" in caplog.text


def test_run_cmd_retcode_non0(caplog):
    """
    Test that when the command fails, it returns a non-zero int, writes an error message,
    but does not quit the program (eof=False by default).
    """
    cmd = "prokka"
    error = "error trying to run prodigal"
    assert utils.run_cmd(cmd, error).returncode != 0
    assert error in caplog.text


def test_run_cmd_retcode_non0_quit(caplog):
    """
    Test that when the command fails, it returns a non-zero int, writes an error message,
    and exits the program (eof=True).
    """
    cmd = "prokka"
    error = "error trying to run prodigal"
    with pytest.raises(SystemExit):
        utils.run_cmd(cmd, error, eof=True)
    assert error in caplog.text


def test_run_cmd_error_stderrfile(caplog):
    """
    Test that when we try to run a command which does not exist, and direct its output to
    a file instead of stderr, we have the expected error written in the given error file.
    """
    cmd = "toto"
    error = "error trying to run toto"
    outfile = open(os.path.join(GENEPATH, "stderr_run_cmd.txt"), "w")
    assert utils.run_cmd(cmd, error, stderr=outfile) == 1
    assert "error: command '>toto' is not possible" in caplog.text
    outfile.close()


def test_run_cmd_error_stdoutfile(caplog):
    """
    Test that when we try to run a command which does not exist, and direct its output to
    a file instead of stderr, we have the expected error written in the given error file.
    """
    cmd = "toto"
    error = "error trying to run toto"
    outfile = open(os.path.join(GENEPATH, "stdout_run_cmd.txt"), "w")
    assert utils.run_cmd(cmd, error, stdout=outfile) == 1
    assert "error: command '>toto' is not possible" in caplog.text
    outfile.close()


def test_rename_contigs():
    """
    From a given sequence, rename all its contigs with the given gembase name + a number,
    and save the output sequence to the given res_path.
    Check that the output file is as expected.
    """
    logger = logging.getLogger("default")
    gpath = os.path.join(DATA_DIR, "genomes", "H299_H561.fasta")
    gembase_name = "ESCO.0216.00005"
    outfile = os.path.join(GENEPATH, "H299_H561.fasta-short-contig.fna")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_H299_H561-ESCO00005.fna")
    contigs, sizes = utils.get_genome_contigs_and_rename(gembase_name, gpath, outfile, logger)
    assert contigs == {"H561_S27":"ESCO.0216.00005.0001",
                       "H561_S28":"ESCO.0216.00005.0002",
                       "H561_S29":"ESCO.0216.00005.0003"}
    assert sizes == {"ESCO.0216.00005.0001":3480,
                     "ESCO.0216.00005.0002":7080,
                     "ESCO.0216.00005.0003":2583}
    assert utilities.compare_order_content(outfile, exp_file)


def test_rename_contigs_empty_fasta(caplog):
    """
    From an empty fasta file, ask to rename contigs. Should return empty tuple of dict, and error
    message.
    """
    logger = logging.getLogger("default")
    gpath = os.path.join(GENEPATH, "empty.fasta")
    # Create empty fasta file
    open(gpath, "w").close()
    gembase_name = "ESCO.0216.00005"
    outfile = os.path.join(GENEPATH, "H299_H561.fasta-short-contig.fna")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_H299_H561-ESCO00005.fna")
    contigs, sizes = utils.get_genome_contigs_and_rename(gembase_name, gpath, outfile, logger)
    assert contigs == {}
    assert sizes == {}
    # Check that outfile is empty
    with open(outfile, "r") as of:
        assert len(of.readlines()) == 0
    # Check error message
    assert ("Your genome test/data/annotate/generated_by_unit-tests/empty.fasta does "
            "not contain any sequence, or is not in fasta format") in caplog.text


def test_rename_contigs_non_fasta(caplog):
    """
    From an file which is not empty, but not in fasta format (no headers starting with '>'),
    check that contigs, sizes and outfile are empty, and error message is ok
    """
    logger = logging.getLogger("default")
    gpath = os.path.join(TEST_DIR, "non-fasta.seq")
    # Create empty fasta file
    open(gpath, "w").close()
    gembase_name = "ESCO.0216.00005"
    outfile = os.path.join(GENEPATH, "H299_H561.fasta-short-contig.fna")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_H299_H561-ESCO00005.fna")
    contigs, sizes = utils.get_genome_contigs_and_rename(gembase_name, gpath, outfile, logger)
    assert contigs == {}
    assert sizes == {}
    # Check that outfile is empty
    with open(outfile, "r") as of:
        assert len(of.readlines()) == 0
    # Check error message
    assert ("Your genome test/data/annotate/test_files/non-fasta.seq does "
            "not contain any sequence, or is not in fasta format") in caplog.text


def test_rename_contigs_duplicate(caplog):
    """
    From a given sequence, there are 2 contigs named "contig2". Stops and returns false
    """
    logger = logging.getLogger("default")
    gpath = os.path.join(DATA_DIR, "genomes", "genome-duplicated-header.fasta")
    gembase_name = "ESCO.0216.00005"
    outfile = os.path.join(GENEPATH, "genome_dup_error.fna")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_H299_H561-ESCO00005.fna")
    contigs, sizes = utils.get_genome_contigs_and_rename(gembase_name, gpath, outfile, logger)
    assert not contigs
    assert not sizes
    with open(outfile, "r") as of:
        assert of.readline().startswith(">ESCO.0216.00005.0001")
        of.readline() # skip sequence
        assert of.readline().startswith(">ESCO.0216.00005.0002")
        of.readline() # skip sequence
        assert of.readline().startswith(">ESCO.0216.00005.0003")
    assert ("several contigs have the same name contig2 in test/data/annotate/genomes/"
            "genome-duplicated-header.fasta.") in caplog.text


def test_rename_contigs_duplicate_last(caplog):
    """
    The last contig of the sequence has the same name as a previous contig. Stops and returns false
    """
    logger = logging.getLogger("default")
    gpath = os.path.join(DATA_DIR, "genomes", "genome-duplicated-header-last.fasta")
    gembase_name = "ESCO.0216.00005"
    outfile = os.path.join(GENEPATH, "genome_dup_error.fna")
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_H299_H561-ESCO00005.fna")
    contigs, sizes = utils.get_genome_contigs_and_rename(gembase_name, gpath, outfile, logger)
    assert not contigs
    assert not sizes
    with open(outfile, "r") as of:
        assert of.readline().startswith(">ESCO.0216.00005.0001")
        of.readline() # skip sequence
        assert of.readline().startswith(">ESCO.0216.00005.0002")
        of.readline() # skip sequence
        assert of.readline().startswith(">ESCO.0216.00005.0003")
    assert ("several contigs have the same name contig2 in test/data/annotate/genomes/"
            "genome-duplicated-header-last.fasta.") in caplog.text


def test_cat_nobar():
    """
    Check that when cat is called on a list of several files, the output file
    contains what is expected (concatenation of content of all input files)
    """
    genomes = ["genome1.fasta", "genome_long_header.fst", "genome6.fasta"]
    list_files = [os.path.join(DATA_DIR, "genomes", gen) for gen in genomes]
    outfile = os.path.join(GENEPATH, "test_catfile_nobar.txt")
    utils.cat(list_files, outfile)
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_test_cat-genomes.fst")
    assert utilities.compare_file_content(outfile, exp_file)


def test_cat_bar():
    """
    Check that when cat is called on a list of several files, the output file
    contains what is expected (concatenation of content of all input files)
    """
    genomes = ["genome1.fasta", "genome_long_header.fst", "genome6.fasta"]
    list_files = [os.path.join(DATA_DIR, "genomes", gen) for gen in genomes]
    outfile = os.path.join(GENEPATH, "test_catfile-bar.txt")
    title = "test cat progressbar"
    utils.cat(list_files, outfile, title=title)
    exp_file = os.path.join(DATA_DIR, "exp_files", "res_test_cat-genomes.fst")
    assert utilities.compare_file_content(outfile, exp_file)


def test_detail():
    """
    Check value returned for detail level
    """
    assert utils.detail_lvl() == 15


def test_grep():
    """
    Check that it returns all lines containing the given pattern
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = os.path.join(GENEPATH, "filein.txt")
    with open(filein, "w") as ff:
        ff.write(lines)
    pattern = "[0-9]+toto[.]txt [a-z]{6}"
    lines_grep = utils.grep(filein, pattern)
    exp = ["1toto.txt otherletters",
           "123toto.txt otherletters"]
    assert exp == lines_grep


def test_grep_count():
    """
    Check that it returns the number of lines containing the given pattern
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = os.path.join(GENEPATH, "filein.txt")
    with open(filein, "w") as ff:
        ff.write(lines)
    pattern = "[0-9]+toto[.]txt [a-z]{6}"
    lines_grep = utils.grep(filein, pattern, counts=True)
    assert lines_grep == 2


def test_count_lines():
    """
    Test that it returns the number of lines in given file
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = os.path.join(GENEPATH, "filein.txt")
    with open(filein, "w") as ff:
        ff.write(lines)
    nbline = utils.count(filein)
    assert nbline == 5


def test_count_words():
    """
    Test that it returns the number of words in given file
    """
    lines = ("1toto.txt otherletters\n"
             "123toto.txtotherletters\n"
             "123toto.txt other letters\n"
             "123toto:txt otherletters\n"
             "123toto.txt otherletters\n")
    filein = os.path.join(GENEPATH, "filein.txt")
    with open(filein, "w") as ff:
        ff.write(lines)
    nbword = utils.count(filein, get="words")
    assert nbword == 10


def test_count_error(caplog):
    """
    test that when we want to count something else than 'lines' or 'words', it returns an error
    """
    with pytest.raises(SystemExit):
        utils.count("filein", get="letters")
    assert "Choose what you want to count among ['lines', 'words']" in caplog.text


def test_save_bin():
    """
    Test that python objects are correctly saved into binary file
    """
    obj1 = {1: "toto", "key": 455, "key2": "val"}
    obj2 = [1, 2, 5, "toto", "plop"]
    obj3 = "a string"
    objects = [obj1, obj2, obj3]
    fileout = os.path.join(GENEPATH, "fileout.txt")
    utils.save_bin(objects, fileout)
    import _pickle as pickle
    with open(fileout, "rb") as binf:
        obj = pickle.load(binf)
    assert objects == obj


def test_load_bin():
    """
    Test that python objects are correctly recovered from binary file
    """
    fileout = os.path.join("test", "data", "utils", "objects_saved.bin")
    found = utils.load_bin(fileout)
    obj1 = {1: "toto", "key": 455, "key2": "val"}
    obj2 = [1, 2, 5, "toto", "plop"]
    obj3 = "a string"
    objects = [obj1, obj2, obj3]
    assert found == objects


def test_list_to_str():
    """
    Test that the given list is returned as expected string
    """
    inlist = [1, 2, "toto", {1: 2, "1": 5}, 1e-6]
    tostr = utils.list_to_str(inlist, sep=" ")
    exp1 = "1 2 toto {1: 2, '1': 5} 1e-06\n"
    exp2 = "1 2 toto {'1': 5, 1: 2} 1e-06\n"
    assert tostr == exp1 or tostr == exp2


def test_list_to_str_default():
    """
    Test that the given list is returned as expected string
    """
    inlist = [1, 2, "toto", {1: 2, "1": 5}, 1e-6]
    tostr = utils.list_to_str(inlist)
    exp1 = "1\t2\ttoto\t{1: 2, '1': 5}\t1e-06\n"
    exp2 = "1\t2\ttoto\t{'1': 5, 1: 2}\t1e-06\n"
    assert tostr == exp1 or tostr == exp2


def test_write_list():
    """
    test that given a list, it writes all its elements, 1 by line, in given outfile
    """
    outfile = os.path.join(GENEPATH, "toto.txt")
    inlist = [1, 2, "toto", {1: 2, "1": 5}, 1e-6]
    utils.write_list(inlist, outfile)
    with open(outfile, "r") as of:
        lines = of.readlines()
    assert len(lines) == 5
    assert "1\n" in lines
    assert "2\n" in lines
    assert "toto\n" in lines
    assert "{1: 2, '1': 5}\n" in lines
    assert "1e-06\n" in lines


def test_remove_exits():
    """
    Test that a given file is removed if it exists
    """
    infile = os.path.join(GENEPATH, "toto.txt")
    open(infile, "w").close()
    assert os.path.isfile(infile)
    utils.remove(infile)
    assert not os.path.isfile(infile)


def test_remove_not_exist():
    """
    Test that removing a file which does not exist brings no error
    """
    infile = os.path.join(GENEPATH, "toto.txt")
    assert not os.path.isfile(infile)
    utils.remove(infile)
    assert not os.path.isfile(infile)
