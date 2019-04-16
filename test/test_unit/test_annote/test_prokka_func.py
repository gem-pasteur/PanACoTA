#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for utils.py
"""

import os
import shutil
import logging
import genomeAPCAT.annote_module.prokka_functions as pfunc
import genomeAPCAT.utils as utils

logfile_base = "panacota"

# Define methods and variables shared by several tests
def my_logger():
    """
    logger given to function called by a subprocess
    """
    import multiprocessing
    m = multiprocessing.Manager()
    q = m.Queue()
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    return q, logging.getLogger('process')


def teardown_function(function):
    """
    Remove logfiles when test is done
    """
    if os.path.isfile(logfile_base + ".log"):
        os.remove(logfile_base + ".log")
    if os.path.isfile(logfile_base + ".log.err"):
        os.remove(logfile_base + ".log.err")
    if os.path.isfile(logfile_base + ".log.details"):
        os.remove(logfile_base + ".log.details")


# Start tests
def test_count_tbl():
    """
    Count the different features found in the tbl file, and return
    nbcont, nbCDS, nbGene, nbCRISPR
    """
    tblfile = os.path.join("test", "data", "annotate", "test_files", "original_name.fna-prokkaRes",
                           "prokka_out_for_test.tbl")
    ncont, ncds, ngene, ncris = pfunc.count_tbl(tblfile)
    assert ncont == 7
    assert ncds == 13
    assert ngene == 15
    assert ncris == 2


def test_count_headers():
    """
    Count how many sequences there are in the given multi-fasta file
    """
    seqfile = os.path.join("test", "data", "annotate", "genomes", "genome4.fasta")
    nb = pfunc.count_headers(seqfile)
    assert nb == 5


def test_check_prokka_no_outdir():
    """
    Test that prokka returns the right error message when output directory does not exist
    """
    logger = my_logger()
    outdir = "toto"
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont, logger[1])
    q = logger[0]
    assert q.qsize() == 1
    msg = "Prokka could not run properly. Look at prokka.log for more information."
    assert q.get().message == msg


def test_check_prokka_notbl():
    """
    Check that check_prokka returns false when a tbl file is missing, and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_notbl")
    name = "prokka_out_for_test-misstbl"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = "prokka.log"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-misstbl original_name-error.fna: no .tbl file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_sevtbl():
    """
    Check that check_prokka returns false when there is more than 1 tbl file,
    and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_notbl")
    name = "prokka_out_for_test-misstbl"
    gpath = "path/to/nogenome/original_name-error.fna"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + "2.tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = "prokka.log"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-misstbl original_name-error.fna: several .tbl files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_nofaa():
    """
    Check that check_prokka returns false when a faa file is missing, and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_nofaa")
    name = "prokka_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missfaa original_name.fna: no .faa file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_sevfaa():
    """
    Check that check_prokka returns false when there is more than 1 faa file,
    and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_nofaa")
    name = "prokka_out_for_test-missfaa"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + "2.faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missfaa original_name.fna: several .faa files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_noffn():
    """
    Check that check_prokka returns false when a ffn file is missing, and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_noffn")
    name = "prokka_out_for_test-missffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missffn original_name.fna: no .ffn file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_sevffn():
    """
    Check that check_prokka returns false when there is more than 1 ffn file,
    and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_noffn")
    name = "prokka_out_for_test-missffn"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + "2.ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missffn original_name.fna: several .ffn files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_nogff():
    """
    Check that check_prokka returns false when a ffn file is missing, and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_noffn")
    name = "prokka_out_for_test-missgff"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-missgff original_name.fna: no .gff file"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_sevgff():
    """
    Check that check_prokka returns false when there is more than 1 ffn file,
    and an error message
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "out_test_noffn")
    name = "prokka_out_for_test-sevgff"
    os.makedirs(out_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + "2.gff"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = "prokka_out_for_test-sevgff original_name.fna: several .gff files"
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_wrong_cont():
    """
    Check that check_prokka returns an error message when the number of contigs in tbl
    file is not as expected
    """
    logger = my_logger()
    outdir = os.path.join("test", "data", "annotate", "test_files",
                          "original_name.fna-prokkaRes")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 10
    assert not pfunc.check_prokka(outdir, logf, name, gpath, nbcont, logger[1])
    msg = ("prokka_out_for_test original_name.fna: no matching number of contigs; "
           "nbcontig=10; in tbl =7")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg


def test_check_prokka_wrong_tbl_cds():
    """
    Check that check_prokka returns an error message when the number of CDS in tbl
    file is different from the number of headers in faa file
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "res_checkProkkaWrongTbl")
    os.makedirs(out_dir)
    name = "prokka_out_for_test-wrongCDS"
    tblfile = os.path.join("test", "data", "annotate", "test_files", name + ".tbl")
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))

    shutil.copyfile(tblfile, os.path.join(out_dir, name + ".tbl"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg1 = ("prokka_out_for_test-wrongCDS original_name.fna: "
            "no matching number of proteins between tbl and faa; "
            "faa=13; in tbl =12")
    msg2 = ("prokka_out_for_test-wrongCDS original_name.fna: "
            "no matching number of genes between tbl and ffn; "
            "ffn=17; in tbl =14genes 2CRISPR")
    q = logger[0]
    assert q.qsize() == 2
    assert q.get().message == msg1
    assert q.get().message == msg2
    shutil.rmtree(out_dir)


def test_check_prokka_wrong_tbl_crispr():
    """
    Check that check_prokka returns an error message when the number of headers in ffn
    file is different from the number of CDS + CRISPR in tbl file (1CRISPR in tbl, 2 in ffn)
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "res_checkProkkaWrongCRISPR")
    os.makedirs(out_dir)
    name = "prokka_out_for_test-wrongtblCRISP"
    tblfile = os.path.join("test", "data", "annotate", "test_files", name + ".tbl")
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(out_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    shutil.copyfile(tblfile, os.path.join(out_dir, name + ".tbl"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert not pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    msg = ("prokka_out_for_test-wrongtblCRISP original_name.fna: "
           "no matching number of genes between tbl and ffn; "
           "ffn=17; in tbl =15genes 1CRISPR")
    q = logger[0]
    assert q.qsize() == 1
    assert q.get().message == msg
    shutil.rmtree(out_dir)


def test_check_prokka_tbl_crispr_newversion():
    """
    Check that check_prokka does not return an error message when the number of headers in ffn
    file is equal to the number of CDS in tbl file (1CRISPR in tbl, 0 in ffn), but
    does not contain the CRISPRs found in tbl
    As the new version of prokka (1.12) does not put crisprs in .ffn
    """
    logger = my_logger()
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    out_dir = os.path.join("test", "data", "annotate", "res_checkProkkaWrongCRISPRnewversion")
    os.makedirs(out_dir)
    name = "prokka_out_for_test-wrongtblCRISPnewversion"
    ffnfile = os.path.join("test", "data", "annotate", "test_files", name + ".ffn")
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".tbl"),
                    os.path.join(out_dir, name + ".tbl"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(out_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(out_dir, name + ".gff"))
    shutil.copyfile(ffnfile, os.path.join(out_dir, name + ".ffn"))
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert pfunc.check_prokka(out_dir, logf, name, gpath, nbcont, logger[1])
    shutil.rmtree(out_dir)


def test_check_prokka_ok():
    """
    Check that everything is ok with prokka results (tbl, faa and ffn files exist,
    and number of CDS, CRISPR and genes correspond between them)
    """
    logger = my_logger()
    outdir = os.path.join("test", "data", "annotate", "test_files", "original_name.fna-prokkaRes")
    name = "prokka_out_for_test"
    logf = "prokka.log"
    gpath = "path/to/nogenome/original_name.fna"
    nbcont = 7
    assert pfunc.check_prokka(outdir, logf, name, gpath, nbcont, logger[1])


def test_run_prokka_out_exists_ok():
    """
    Test that when the output directory already exists, and files inside are OK,
    run_prokka returns True, with a warning message indicating that prokka did not rerun.
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    gpath = "path/to/nogenome/original_name.fna"
    outdir = os.path.join("test", "data", "annotate", "test_files")
    cores_prokka = 1
    name = "prokka_out_for_test"
    force = False
    nbcont = 7
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont, logger[0])
    assert pfunc.run_prokka(arguments)

    q = logger[0]
    assert q.qsize() == 4
    # start annotating :
    assert q.get().message.startswith("Start annotating")
    # warning prokka results folder exists:
    assert q.get().message.startswith("Prokka results folder already exists.")
    # Results in result folder are ok
    assert q.get().message.startswith("Prokka did not run again, formatting step used already "
                                      "generated results of Prokka in ")
    # End annotation:
    assert q.get().message.startswith("End annotating")


def test_run_prokka_out_exists_error():
    """
    Test that when the output directory already exists, and 1 file is missing,
    run_prokka returns False, and writes the warning message saying that prokka did not
    rerun, + the warning message for the missing file(s).
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    ori_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name.fna-prokkaRes")
    ori_name = "prokka_out_for_test"
    new_prok_dir = os.path.join("test", "data", "annotate", "test_files",
                                "original_name-error-prokkaRes")
    name = "prokka_out_for_test-wrongCDS"
    os.makedirs(new_prok_dir)
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".ffn"),
                    os.path.join(new_prok_dir, name + ".ffn"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".faa"),
                    os.path.join(new_prok_dir, name + ".faa"))
    shutil.copyfile(os.path.join(ori_prok_dir, ori_name + ".gff"),
                    os.path.join(new_prok_dir, name + ".gff"))
    gpath = "path/to/nogenome/original_name-error"
    outdir = os.path.join("test", "data", "annotate", "test_files")
    cores_prokka = 1
    force = False
    nbcont = 7
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont, logger[0])
    assert not pfunc.run_prokka(arguments)
    q = logger[0]
    assert q.qsize() == 4
    # start annotating :
    assert q.get().message.startswith("Start annotating")
    # warning prokka results folder exists:
    assert q.get().message == "Prokka results folder already exists."
    # error, no tbl file
    msg = "prokka_out_for_test-wrongCDS original_name-error: no .tbl file"
    assert q.get().message == msg
    # warning, files in outdir are not as expected
    assert q.get().message.startswith("Problems in the files contained in your already existing "
                                      "output dir ")
    shutil.rmtree(new_prok_dir)


def test_run_prokka_out_exists_force():
    """
    Test that when the output directory already exists with wrong files, but force is on,
    prokka is rerun and outputs the right files
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    gpath = os.path.join("test", "data", "annotate", "genomes", "H299_H561.fasta")
    outdir = os.path.join("test", "data", "annotate")
    out_prokdir = os.path.join(outdir, "H299_H561.fasta-prokkaRes")
    name = "test_runprokka_H299"
    # Put empty tbl, faa, ffn files in prokka output dir, to check that they are overridden
    os.makedirs(out_prokdir)
    open(os.path.join(out_prokdir, name + ".tbl"), "w").close()
    open(os.path.join(out_prokdir, name + ".faa"), "w").close()
    open(os.path.join(out_prokdir, name + ".ffn"), "w").close()
    cores_prokka = 5
    force = True
    nbcont = 3
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont, logger[0])
    assert pfunc.run_prokka(arguments)
    # Check content of tbl, ffn and faa files
    exp_dir = os.path.join("test", "data", "annotate", "exp_files",
                           "H299_H561.fasta-short-contig.fna-prokkaRes",
                           "test_runprokka_H299")
    out_tbl = os.path.join(out_prokdir, name + ".tbl")
    out_faa = os.path.join(out_prokdir, name + ".faa")
    out_ffn = os.path.join(out_prokdir, name + ".ffn")
    assert os.path.isfile(out_tbl)
    with open(exp_dir + ".tbl", "r") as expf, open(out_tbl, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    assert os.path.isfile(out_faa)
    with open(exp_dir + ".faa", "r") as expf, open(out_faa, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    assert os.path.isfile(out_ffn)
    with open(exp_dir + ".ffn", "r") as expf, open(out_ffn, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    shutil.rmtree(out_prokdir)
    q = logger[0]
    assert q.qsize() == 3
    assert q.get() .message.startswith("Start annotating")
    assert q.get() .message == ("Prokka results folder already exists, but removed because "
                                "--force option used")
    assert q.get() .message.startswith("End annotating")
    os.remove(os.path.join(outdir, "H299_H561.fasta-prokka.log"))


def test_run_prokka_out_doesnt_exist():
    """
    Test that when the output directory does not exist, it creates it, and runs prokka
    with all expected outfiles
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    gpath = os.path.join("test", "data", "annotate", "genomes", "H299_H561.fasta")
    prok_dir = os.path.join("test", "data", "annotate")
    out_dir = os.path.join(prok_dir, "H299_H561.fasta-prokkaRes")
    cores_prokka = 5
    name = "test_runprokka_H299"
    force = False
    nbcont = 3
    arguments = (gpath, prok_dir, cores_prokka, name, force, nbcont, logger[0])
    assert pfunc.run_prokka(arguments)
    # Check content of tbl, ffn and faa files
    exp_dir = os.path.join("test", "data", "annotate", "exp_files",
                           "H299_H561.fasta-short-contig.fna-prokkaRes",
                           "test_runprokka_H299")
    out_tbl = os.path.join(out_dir, name + ".tbl")
    out_faa = os.path.join(out_dir, name + ".faa")
    out_ffn = os.path.join(out_dir, name + ".ffn")
    out_gff = os.path.join(out_dir, name + ".gff")
    assert os.path.isfile(out_tbl)
    with open(exp_dir + ".tbl", "r") as expf, open(out_tbl, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    assert os.path.isfile(out_faa)
    with open(exp_dir + ".faa", "r") as expf, open(out_faa, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    assert os.path.isfile(out_ffn)
    with open(exp_dir + ".ffn", "r") as expf, open(out_ffn, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    assert os.path.isfile(out_gff)
    with open(exp_dir + ".gff", "r") as expf, open(out_gff, "r") as outf:
        for line_exp, line_out in zip(expf, outf):
            assert line_exp == line_out
    q = logger[0]
    assert q.qsize() == 2
    assert q.get().message.startswith("Start annotating")
    assert q.get().message.startswith("End annotating")
    shutil.rmtree(out_dir)
    os.remove(os.path.join(prok_dir, "H299_H561.fasta-prokka.log"))


def test_run_prokka_out_problem_running():
    """
    Check that when a problem occurs while trying to run prokka, run_prokka returns False,
    and the error message indicating to read in the log why it couldn't run
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    gpath = os.path.join("test", "data", "annotate", "genomes", "H299 H561.fasta")
    outdir = os.path.join("test", "data", "annotate")
    cores_prokka = 5
    name = "test_runprokka_H299-error"
    force = False
    nbcont = 3
    logf = os.path.join(outdir, "H299 H561.fasta-prokka.log")
    arguments = (gpath, outdir, cores_prokka, name, force, nbcont, logger[0])
    assert not pfunc.run_prokka(arguments)
    q = logger[0]
    assert q.qsize() == 2
    assert q.get().message.startswith("Start annotating")
    assert q.get().message == ("Prokka could not run properly. Look at {} for more "
                               "information.").format(logf)
    os.remove(logf)



def test_run_all_1by1():
    """
    Check that when running with 3 threads (not parallel), prokka runs as expected,
    and returns True for each genome
    Start and end must be ordered: (start1, end1, start2, end2) or (start2, end2, start1, end1)
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join("test", "data", "annotate", "genomes", genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join("test", "data", "annotate", "genomes", genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, 456464645, 1, 465]}
    threads = 3
    force = False
    prok_folder = os.path.join("test", "data", "annotate")
    res_dir = os.path.join("test", "data", "annotate")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
    #  arguments = [(genomes[g][1], prok_folder, cores_prokka, genomes[g][0],
                  # force, genomes[g][3], q)
    assert final[genome1]
    assert final[genome2]
    shutil.rmtree(os.path.join(prok_folder, genome1 + "-prokkaRes"))
    shutil.rmtree(os.path.join(prok_folder, genome2 + "-prokkaRes"))
    os.remove(os.path.join(res_dir, genome1 + "-prokka.log"))
    os.remove(os.path.join(res_dir, genome2 + "-prokka.log"))
    q = logger[0]
    assert q.qsize() == 5
    assert q.get().message == 'Annotating all genomes with prokka'
    # Messages for start and end annotation of the different genomes
    message_start_annot1 = ("Start annotating test_runall_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_start_annot2 = ("Start annotating test_runall_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    message_end_annot1 = ("End annotating test_runall_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_end_annot2 = ("End annotating test_runall_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    qget = q.get().message
    assert qget == message_start_annot1 or message_start_annot2
    if qget == message_start_annot1:
        # Ending annotation of first genome (same genome as started because running 1by1)
        assert q.get().message == message_end_annot1
    else:
        assert q.get().message == message_end_annot2
    qget2 = q.get().message
    assert qget2 == message_start_annot1 or message_start_annot2
    if qget2 == message_start_annot2:
        # Ending annotation of first genome (same genome as started because running 1by1)
        assert q.get().message == message_end_annot2
    else:
        assert q.get().message == message_end_annot1


def test_run_all_parallel_more_threads():
    """
    Check that there is no problem when running with more threads than genomes (each genome
    uses nb_threads/nb_genome threads)
    Start and end are not necessarily in the same order (ex: start1, start2, end2, end1)
    """
    logger = my_logger()
    utils.init_logger(logfile_base, 0, '')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    genome1 = "H299_H561.fasta"
    gpath1 = os.path.join("test", "data", "annotate", "genomes", genome1)
    genome2 = "A_H738.fasta"
    gpath2 = os.path.join("test", "data", "annotate", "genomes", genome2)
    genomes = {genome1: ["test_runall_1by1_1", gpath1, 12656, 3, 0],
               genome2: ["test_runall_1by1_2", gpath2, 456464645, 1, 465]}
    threads = 8
    force = False
    prok_folder = os.path.join("test", "data", "annotate")
    res_dir = os.path.join("test", "data", "annotate")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
    assert final[genome1]
    assert final[genome2]
    shutil.rmtree(os.path.join(prok_folder, genome1 + "-prokkaRes"))
    shutil.rmtree(os.path.join(prok_folder, genome2 + "-prokkaRes"))
    os.remove(os.path.join(res_dir, genome1 + "-prokka.log"))
    os.remove(os.path.join(res_dir, genome2 + "-prokka.log"))
    q = logger[0]
    assert q.qsize() == 5
    assert q.get().message == 'Annotating all genomes with prokka'
    message_start_annot1 = ("Start annotating test_run_all_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_start_annot2 = ("Start annotating test_run_all_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    # Starting annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_start_annot1 or message_start_annot2
    message_end_annot1 = ("End annotating test_run_all_1by1_1 test/data/annotate/genomes/"
                            "H299_H561.fasta")
    message_end_annot2 = ("End annotating test_run_all_1by1_2 test/data/annotate/genomes/"
                            "A_H738.fasta")
    # Ending annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_end_annot1 or message_end_annot2
    # Starting annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_start_annot1 or message_start_annot2
    # Ending annotation of first genome (can be genome1 or genome2, because dicts are unordered)
    assert q.get().message == message_end_annot1 or message_end_annot2


def test_run_all_parallel_less_threads():
    """
    Check that there is no problem when running with less threads than genomes (each genomes
    uses 2 threads)
    Genomes H299 and A_H738 should run well, but genomes genome* have problems (no CDS found),
    so check_prokka should return false.
    """
    utils.init_logger(logfile_base, 0, '')
    # genomes = {genome: [name, gpath, size, nbcont, l90]}
    gnames = ["H299_H561.fasta", "A_H738.fasta", "genome1.fasta", "genome2.fasta", "genome3.fasta"]
    gpaths = [os.path.join("test", "data", "annotate", "genomes", name) for name in gnames]
    genomes = {gnames[0]: ["test_runall_1by1_1", gpaths[0], 12656, 3, 1],
               gnames[1]: ["test_runall_1by1_2", gpaths[1], 456464645, 1, 1],
               gnames[2]: ["test_runall_1by1_2", gpaths[2], 456464645, 4, 1],
               gnames[3]: ["test_runall_1by1_2", gpaths[3], 456464645, 3, 1],
               gnames[4]: ["test_runall_1by1_2", gpaths[4], 456464645, 1, 1]
               }
    threads = 4
    force = False
    prok_folder = os.path.join("test", "data", "annotate")
    res_dir = os.path.join("test", "data", "annotate")
    final = pfunc.run_prokka_all(genomes, threads, force, prok_folder)
    assert final[gnames[0]]
    assert final[gnames[1]]
    assert not final[gnames[2]]
    assert not final[gnames[3]]
    assert not final[gnames[4]]
    for name in gnames:
        shutil.rmtree(os.path.join(prok_folder, name + "-prokkaRes"))
        os.remove(os.path.join(res_dir, name + "-prokka.log"))

