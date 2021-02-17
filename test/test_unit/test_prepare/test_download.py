#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the download_genomes_func submodule in prepare module
"""
import os
import logging
import glob
import shutil
import pytest

import PanACoTA.prepare_module.download_genomes_func as downg


DATA_TEST_DIR = os.path.join("test", "data", "prepare")
GENEPATH = os.path.join(DATA_TEST_DIR, "generated_by_unit-tests")


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
    if os.path.isdir(GENEPATH):
        content = os.listdir(GENEPATH)
        for f in content:
            assert f.startswith(".fuse")
    else:
        os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")


def test_to_database():
    """
    Test that all fna.gz files are uncompressed and moved to a created Database_init folder
    """
    out_dir = os.path.join(DATA_TEST_DIR, "genomes")
    nb_gen, db_init_dir = downg.to_database(out_dir, "refseq")
    db_dir = os.path.join(DATA_TEST_DIR, "genomes", "Database_init")
    assert os.path.isdir(db_dir)
    files_all = glob.glob(os.path.join(db_dir, "*"))
    files_fna = glob.glob(os.path.join(db_dir, "*.fna"))
    # Check that there are only 3 files in result database
    assert len(files_all) == len(files_fna)
    # And that those files are .fna files
    assert len(files_fna) == 3
    # Check that we have as many genomes as expected, and that the output database has the
    # expected name
    assert nb_gen == 3
    assert db_init_dir == db_dir
    assert os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))

    shutil.rmtree(db_dir)


def test_to_database_nofolder_refseq(caplog):
    """
    Test behavior when the folder that should contain refseq downloaded genomes does not exist
    -> should exit with error message
    """
    caplog.set_level(logging.DEBUG)
    with pytest.raises(SystemExit):
        downg.to_database(GENEPATH, "genbank")

    assert "ERROR" in caplog.text
    assert ("The folder containing genomes downloaded from NCBI genbank "
            "(test/data/prepare/generated_by_unit-tests/genbank/bacteria) "
            "does not exist.") in caplog.text
    assert ("Check that you really downloaded sequences (fna.gz) and that they are "
            "in this folder") in caplog.text


def test_to_database_nofolder_per_genome(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, but there are no folders inside
    -> should exit with error message
    """
    outdir = os.path.join(GENEPATH, "out-refseq")
    empty_dir = os.path.join(outdir, "refseq", "bacteria")
    os.makedirs(empty_dir)
    caplog.set_level(logging.DEBUG)
    with pytest.raises(SystemExit):
        downg.to_database(outdir, "refseq")
    # Check error message is as expected
    assert "ERROR" in caplog.text
    assert ("The folder supposed to contain genomes downloaded from NCBI refseq "
            "(test/data/prepare/generated_by_unit-tests/out-refseq/refseq/bacteria) "
            "exists but is empty") in caplog.text
    assert ("Check that you really downloaded sequences (fna.gz)") in caplog.text

    # Same with genbank
    outdir = os.path.join(GENEPATH, "out-genbank")
    empty_dir = os.path.join(outdir, "genbank", "bacteria")
    os.makedirs(empty_dir)
    caplog.set_level(logging.DEBUG)
    with pytest.raises(SystemExit):
        downg.to_database(outdir, "genbank")
    # Check error message is as expected
    assert "ERROR" in caplog.text
    assert ("The folder supposed to contain genomes downloaded from NCBI genbank "
            "(test/data/prepare/generated_by_unit-tests/out-genbank/genbank/bacteria) "
            "exists but is empty") in caplog.text
    assert ("Check that you really downloaded sequences (fna.gz)") in caplog.text


def test_to_database_1empty_genome_folder(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, there are subfolders inside,
    but 1 of them is empty: warning message informing that this genome will be ignored
    """
    caplog.set_level(logging.DEBUG)
    out_dir = os.path.join(GENEPATH, "1empty_genome_folder")
    refseq_dir = os.path.join(DATA_TEST_DIR, "genomes")
    # Copy content of refseq in genomes test data to output folder that will be used
    shutil.copytree(refseq_dir, out_dir)

    # Empty 1 directory: move its file to 'out_dir'
    to_remove = os.path.join(out_dir, "refseq", "bacteria", "ACOR003", "ACOR003.0519.fna.gz")
    os.remove(to_remove)
    # Run to_database
    nb_gen, db_dir = downg.to_database(out_dir, "refseq")
    assert nb_gen == 2
    assert db_dir == os.path.join(out_dir, "Database_init")

    # Check that a warning message was raised, indicating that genome is ignored
    assert "WARNING" in caplog.text
    assert ("Problem with genome in ACOR003: no compressed fasta file downloaded. "
            "This genome will be ignored.") in caplog.text
    assert not os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR001.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))


def test_to_database_several_genomes(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, there are subfolders inside,
    but 1 of them contains more than 1 genome: warning message informing that this
    genome will be ignored
    """
    out_dir = os.path.join(GENEPATH, "genomes")
    refseq_dir = os.path.join(DATA_TEST_DIR, "genomes")
    # Copy content of refseq in genomes test data to output folder that will be used
    shutil.copytree(refseq_dir, out_dir)

    # Create a new gz file in one of the genome directories
    to_create_filename = "ACOR002.0519.bis.fna.gz" # Name of file that must be created
    to_fill_dir = "ACOR002" # Directory containing file to create
    to_create_path = os.path.join(out_dir, "refseq", "bacteria", to_fill_dir, to_create_filename)
    # Create empty gz file
    open(to_create_path, "w").close()

    # Run to_database, and check that only 2 genomes were considered
    nb_gen, db_dir = downg.to_database(out_dir, "refseq")
    assert nb_gen == 2
    assert db_dir == os.path.join(out_dir, "Database_init")

    # Check that a warning message was raised, indicating that genome is ignored
    caplog.set_level(logging.DEBUG)
    assert "WARNING" in caplog.text
    assert ("Problem with genome in ACOR002: several compressed fasta files found. "
            "This genome will be ignored.") in caplog.text
    assert not os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR001.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))


def test_to_database_1genome_wrong_format(caplog):
    """
    Test behavior when the folder refseq/bacteria exists, there is 1 genome per subfolder,
    but 1 genome cannot be unzipped
    """
    # out_dir = os.path.join(DATA_TEST_DIR, "genomes")
    # gz_genomes_folder = os.path.join(out_dir, "refseq", "bacteria")

    out_dir = os.path.join(GENEPATH, "genomes")
    refseq_dir = os.path.join(DATA_TEST_DIR, "genomes")
    # Copy content of refseq in genomes test data to output folder that will be used
    shutil.copytree(refseq_dir, out_dir)

    # Name of directory directly containing the original gz file
    to_corrupt_dir = "ACOR001"
    to_corrupt_filename = "ACOR001.0519.fna.gz"
    to_corrupt_path = os.path.join(out_dir, "refseq", "bacteria",  to_corrupt_dir,
                                   to_corrupt_filename)
    # Create fake gz file (txt file)
    false_gz = open(to_corrupt_path, "w")
    false_gz.write("This is not a gz file")
    false_gz.close()

    # Run to_database
    nb_gen, db_dir = downg.to_database(out_dir, "refseq")
    assert nb_gen == 2
    assert db_dir == os.path.join(out_dir, "Database_init")

    # Check that a error message was raised, indicating that genome is ignored
    caplog.set_level(logging.DEBUG)
    assert "ERROR" in caplog.text
    assert ("Error while trying to uncompress "
            "test/data/prepare/generated_by_unit-tests/genomes/Database_init/ACOR001.0519.fna.gz. "
            "This genome will be ignored") in caplog.text
    # Check that there are only 2 files in the database, and that they correspond
    # to uncompressed gz files
    list_db = os.listdir(db_dir)
    assert len(list_db) == 2
    assert not os.path.isfile(os.path.join(db_dir, to_corrupt_filename))
    assert os.path.isfile(os.path.join(db_dir, "ACOR002.0519.fna"))
    assert os.path.isfile(os.path.join(db_dir, "ACOR003.0519.fna"))


def test_download_specify_level(caplog):
    """
    Test that, given a taxid, and a species name,
    it downloads genomes in .gz, and uncompress them in the
    db folder (which is named as expected)

    We cannot compare log, as it is already catched by NCBI_genome_download
    """
    caplog.set_level(logging.INFO)

    species_linked = "Acetobacter_orleanensis"
    section = "refseq"
    NCBI_species = "Acetobacter orleanensis"
    NCBI_species_taxid = "104099"
    NCBI_taxid = ""
    outdir = os.path.join(GENEPATH, "test_download_refseq")
    threads = 1
    levels = ""

    db_dir, nb_gen = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                                outdir, threads)
    # Check path to uncompressed files is as expected
    assert db_dir == os.path.join(outdir, "Database_init")
    # Check number of genomes downloaded. We cannot know the exact value, as it is updated
    # everyday. But in nov. 2019, there are 4 genomes. So, there must be at least those 4 genomes
    assert nb_gen >= 4
    # And that db_dir exists and contains nb_gen files
    assert os.path.isdir(db_dir)
    assert len(os.listdir(db_dir)) == nb_gen

    # Check that assembly summary file wwas donwloaded as expected
    sum_file = os.path.join(outdir, "assembly_summary-Acetobacter_orleanensis.txt" )
    assert os.path.isfile(sum_file)
    # Check number of genomes in summary file, and how many with scaffold or complete
    # assembly level -> will check that when asking only for those levels, we get the same number
    other = 0
    scaf = 0
    comp = 0
    with open(sum_file, "r") as sf:
        sf.readline()  # skip header
        for line in sf:
            if "complete" in line.split("\t")[13].lower():
                comp += 1
            elif "scaffold" in line.split("\t")[13].lower():
                scaf += 1
            else:
                other += 1
    assert other + scaf + comp == nb_gen

    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) >= 4
    # Check log giving species name + species taxid
    assert 'Downloading all genomes for NCBI species = Acetobacter orleanensis (NCBI_species_taxid = 104099)' in caplog.text

    # Re-run, but only asking for complete and scaffold
    outdir2 = os.path.join(GENEPATH, "test_download_refseq_only-scaf")
    levels2 = "scaffold,complete"
    db_dir2, nb_gen2 = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid,
                                                  levels2, outdir2, threads)
    assert scaf + comp == nb_gen2
    assert db_dir2 == os.path.join(outdir2, "Database_init")
    # Check log giving species name + species taxid + levels given
    assert ("Downloading all genomes for NCBI species = Acetobacter orleanensis "
            "(NCBI_species_taxid = 104099). (Only those assembly levels: scaffold,complete)") in caplog.text


def test_download_only_spetaxid(caplog):
    """
    Test that, given a species taxid, it downloads all genomes of the species in .gz, and uncompress them in the
    db folder (which is named as expected)

    We cannot compare log, as it is already catched by NCBI_genome_download
    """
    caplog.set_level(logging.INFO)
    species_linked = "toto"
    section = "refseq"
    NCBI_species = None
    NCBI_species_taxid = "104099"
    NCBI_taxid = ""
    outdir = os.path.join(GENEPATH, "test_download_refseq_noSpe")
    threads = 1
    levels = ""

    db_dir, nb_gen = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                                outdir, threads)

    # Check path to uncompressed files is as expected
    assert db_dir == os.path.join(outdir, "Database_init")
    # Check number of genomes downloaded. We cannot know the exact value, as it is updated everyday. But in nov. 2019, there are 4 genomes. So, there must be at least those 4 genomes
    assert nb_gen >= 4
    # And that db_dir exists and contains nb_gen files
    assert os.path.isdir(db_dir)
    assert len(os.listdir(db_dir)) == nb_gen
    # Check log giving only species taxid
    assert "Downloading all genomes for  NCBI_species_taxid = 104099" in caplog.text

    # Check that assembly summary file was donwloaded as expected
    sum_file = os.path.join(outdir, "assembly_summary-toto.txt" )
    assert os.path.isfile(sum_file)

    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) >= 4


def test_download_taxid_and_spetaxid(caplog):
    """
    Test that, given a species taxid and a taxid, it downloads only the genome(s) corresponding to taxid (intersection)

    We cannot compare log, as it is already catched by NCBI_genome_download
    """
    caplog.set_level(logging.INFO)
    species_linked = "toto-spe"
    section = "refseq"
    NCBI_species = None
    NCBI_species_taxid = "104099"
    NCBI_taxid = "1231342"
    levels = ""
    threads = 1
    outdir2 = os.path.join(GENEPATH, "test_download_refseq_noSpeandSpecific")
    db_dir2, nb_gen2 = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                                  outdir2, threads)

    # Check path to uncompressed files is as expected
    assert db_dir2 == os.path.join(outdir2, "Database_init")
    # Check number of genomes downloaded. We cannot know the exact value, as it is updated everyday. But in nov. 2019, there are 4 genomes. So, there must be at least those 4 genomes
    assert nb_gen2 == 1
    # And that db_dir exists and contains nb_gen files
    assert os.path.isdir(db_dir2)
    assert len(os.listdir(db_dir2)) == 1
    # Check log giving only species taxid
    assert "Downloading all genomes for  NCBI_species_taxid = 104099 (and NCBI_taxid = 1231342)" in caplog.text

    # Check that assembly summary file was donwloaded as expected
    sum_file = os.path.join(outdir2, "assembly_summary-toto-spe.txt" )
    assert os.path.isfile(sum_file)

    # Check that the NCBI_genome_download output directory exists
    ngd_outdir2 = os.path.join(outdir2, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir2)
    assert len(os.listdir(ngd_outdir2)) == 1


def test_download_taxid_and_spename(caplog):
    """
    Test that, given a taxid and a species name, it downloads only the genome(s) corresponding to taxid (intersection)

    We cannot compare log, as it is already catched by NCBI_genome_download
    """
    caplog.set_level(logging.INFO)
    species_linked = "aceor"
    section = "refseq"
    NCBI_species = "Acetobacter orleanensis"
    NCBI_species_taxid = ""
    NCBI_taxid = "1231342"
    levels = ""
    threads = 1
    outdir = os.path.join(GENEPATH, "test_download_refseq_noSpeandSpecific")
    db_dir, nb_gen = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                                outdir, threads)

    # Check path to uncompressed files is as expected
    assert db_dir == os.path.join(outdir, "Database_init")
    # Check number of genomes downloaded. We cannot know the exact value, as it is updated everyday. But in nov. 2019, there are 4 genomes. So, there must be at least those 4 genomes
    assert nb_gen == 1
    # And that db_dir exists and contains nb_gen files
    assert os.path.isdir(db_dir)
    assert len(os.listdir(db_dir)) == 1
    # Check log giving only species taxid
    assert "Downloading all genomes for NCBI species = Acetobacter orleanensis (and NCBI_taxid = 1231342)" in caplog.text

    # Check that assembly summary file was donwloaded as expected
    sum_file = os.path.join(outdir, "assembly_summary-aceor.txt" )
    assert os.path.isfile(sum_file)

    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) == 1


def test_download_specific_strain(caplog):
    """
    Test that, given a taxid of a specific strain, it only downloads this one

    We cannot compare log, as it is already catched by NCBI_genome_download
    """
    caplog.set_level(logging.INFO)
    species_linked = "toto"
    NCBI_species = None
    section = "refseq"
    NCBI_species_taxid = ""
    NCBI_taxid = "1123862"
    outdir = os.path.join(GENEPATH, "test_download_refseq_specific")
    threads = 1
    levels = ""

    db_dir, nb_gen = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                                outdir, threads)

    # Check path to uncompressed files is as expected
    assert db_dir == os.path.join(outdir, "Database_init")
    # Check number of genomes downloaded. We cannot know the exact value, as it is updated everyday. But in nov. 2019, there are 4 genomes. So, there must be at least those 4 genomes
    assert nb_gen == 1
    # And that db_dir exists and contains nb_gen files
    assert os.path.isdir(db_dir)
    assert len(os.listdir(db_dir)) == nb_gen

    # Check that assembly summary file was donwloaded as expected
    sum_file = os.path.join(outdir, "assembly_summary-toto.txt" )
    assert os.path.isfile(sum_file)

    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    # And that it contains folders
    assert os.path.isdir(ngd_outdir)
    assert len(os.listdir(ngd_outdir)) == 1

    # Check log giving only specific taxid
    assert "Downloading all genomes for  NCBI_taxid = 1123862" in caplog.text


def test_download_2taxid(caplog):
    """
    Give a taxid of a subspecies and a taxid of a specific strain. Should download all genomes
    of the subspecies + specific strain.
    If only the subspecies taxid, the specific strain is not downloaded.
    """
    caplog.set_level(logging.INFO)
    species_linked = "salmo"
    section = "refseq"
    NCBI_species = None
    NCBI_species_taxid = ""
    # 913079 is the subspecies Salmonella enterica subsp. enterica serovar Mississippi
    # 1212561  is the strain Salmonella enterica subsp. enterica serovar Mississippi strain 2010K-1406
    NCBI_taxid = "913079,1212561"
    outdir = os.path.join(GENEPATH, "test_download_refseq_2taxid")
    threads = 1
    levels = ""
    db_dir, nb_gen = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                                outdir, threads)

    # Check path to uncompressed files is as expected
    assert db_dir == os.path.join(outdir, "Database_init")
    # Check number of genomes downloaded. We cannot know the exact value, as it is updated everyday. But in nov. 2019, there are 4 genomes. So, there must be at least those 4 genomes
    assert nb_gen >= 13
    # And that db_dir exists and contains nb_gen files
    assert os.path.isdir(db_dir)
    assert len(os.listdir(db_dir)) == nb_gen
    # Check log giving only species taxid
    assert "Downloading all genomes for  NCBI_taxid = 913079,1212561" in caplog.text

    # Check that assembly summary file was donwloaded as expected
    sum_file = os.path.join(outdir, "assembly_summary-salmo.txt" )
    assert os.path.isfile(sum_file)

    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    assert os.path.isdir(ngd_outdir)


    # Redo, without the specific strain taxid. Should download the same -1 (not the specific strain)
    species_linked = "salmo"
    section = "refseq"
    NCBI_species = None
    NCBI_species_taxid = ""
    # 913079 is the subspecies Salmonella enterica subsp. enterica serovar Mississippi
    NCBI_taxid_1 = "913079"
    outdir_1 = os.path.join(GENEPATH, "test_download_refseq_2taxid_1")
    threads = 1
    levels = ""
    db_dir_1, nb_gen_1 = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid_1, levels,
                                                    outdir_1, threads)
    assert nb_gen == nb_gen_1 + 1
    assert "Downloading all genomes for  NCBI_taxid = 913079" in caplog.text


def test_download_refseq_vs_genbank(caplog):
    """
    Give a taxid of a subspecies download strains from refseq, and then from genbank.
    Currently, no strains in refseq, and 2 in genbank.
    39831 = Klebsiella pneumoniae subsp. rhinoscleromatis
    Later, there can be some in refseq, but always at least 2 more in genbank
    """
    caplog.set_level(logging.INFO)
    species_linked = "refseq-genbank"
    section = "refseq"
    NCBI_species = None
    NCBI_species_taxid = ""
    NCBI_taxid = "39831"
    outdir = os.path.join(GENEPATH, "test_download_refseq_genbank")
    levels = ""
    threads = 1

    # With refseq, no genome found
    with pytest.raises(SystemExit):
        downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                   outdir, threads)

    # Check path to uncompressed files does not exist
    assert not os.path.isdir(os.path.join(outdir, "Database_init"))

    # Check that the NCBI_genome_download output directory was not created
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    assert not os.path.isdir(ngd_outdir)

    # Check logs
    assert "ERROR" in caplog.text
    assert ("Could not download genomes. Check that you gave valid NCBI taxid and/or "
            "NCBI species name. If you gave both, check that given taxID and name really "
            "correspond to the same species.") in caplog.text

    # REDO with genbank instead of refseq
    section = "genbank"
    outdir2 = os.path.join(GENEPATH, "test_download_genbank")
    db_dir, nb_gen = downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                                outdir2, threads)

    # Check path to uncompressed files is as expected
    assert db_dir == os.path.join(outdir2, "Database_init")
    # Check number of genomes downloaded. We cannot know the exact value, as it is updated everyday. But in nov. 2019, there are 4 genomes. So, there must be at least those 4 genomes
    assert nb_gen >= 2
    # And that db_dir exists and contains nb_gen files
    assert os.path.isdir(db_dir)
    assert len(os.listdir(db_dir)) == nb_gen
    # Check log giving only species taxid
    assert "Downloading all genomes for  NCBI_taxid = 39831" in caplog.text
    # Check that assembly summary file was donwloaded as expected
    sum_file = os.path.join(outdir2, "assembly_summary-refseq-genbank.txt" )
    assert os.path.isfile(sum_file)
    # Check that the NCBI_genome_download output directory exists
    ngd_outdir = os.path.join(outdir2, "genbank", "bacteria")
    assert os.path.isdir(ngd_outdir)


def test_download_wrongTaxID(caplog):
    """
    Test that, when a non existing taxid is given, it exits (with error message)

    We cannot compare log, as it is already catched by NCBI_genome_download
    """
    species_linked = "Acetobacter_orleanensis"
    NCBI_species = None
    section = "refseq"
    NCBI_species_taxid = "10409"
    NCBI_taxid = ""
    outdir = os.path.join(GENEPATH, "test_download_refseq_wrongTaxID")
    threads = 1
    levels = ""
    with pytest.raises(SystemExit):
        downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                   outdir, threads)

    # Check path to uncompressed files does not exist
    assert not os.path.isdir(os.path.join(outdir, "Database_init"))

    # Check that the NCBI_genome_download output directory was not created
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    assert not os.path.isdir(ngd_outdir)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert "ERROR" in caplog.text
    assert ("Could not download genomes. Check that you gave valid NCBI taxid and/or "
            "NCBI species name. If you gave both, check that given taxID and name really "
            "correspond to the same species.") in caplog.text

    # Check that output directory was not created
    assert not os.path.isdir(outdir)


def test_download_diffSpeTaxID(caplog):
    """
    Test that, when a taxID and a species name are given, but those 2 elements do not
    match with the same genomes, it exists with error message

    We cannot compare log, as it is already catched by NCBI_genome_download
    """
    species_linked = "Acetobacter_orleanensis"
    section = "refseq"
    NCBI_species = "Acetobacter fabarum"
    NCBI_species_taxid = "104099"
    NCBI_taxid = ""
    outdir = os.path.join(GENEPATH, "test_download_refseq_wrongTaxID")
    threads = 1
    levels = ""
    with pytest.raises(SystemExit):
        downg.download_from_ncbi(species_linked, section, NCBI_species, NCBI_species_taxid, NCBI_taxid, levels,
                                   outdir, threads)

    # Check path to uncompressed files does not exist
    assert not os.path.isdir(os.path.join(outdir, "Database_init"))

    # Check that the NCBI_genome_download output directory was not created
    ngd_outdir = os.path.join(outdir, "refseq", "bacteria")
    assert not os.path.isdir(ngd_outdir)

    # Check logs
    caplog.set_level(logging.DEBUG)
    assert "ERROR" in caplog.text
    assert ("Could not download genomes. Check that you gave valid NCBI taxid and/or "
            "NCBI species name. If you gave both, check that given taxID and name really "
            "correspond to the same species.") in caplog.text

    # Check that output directory was not created
    assert not os.path.isdir(outdir)
