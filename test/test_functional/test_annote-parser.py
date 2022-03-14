#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of 'annotate' subcommand
"""
import pytest
import argparse
import time
import os
import shutil

from PanACoTA.subcommands import annotate as annot
DBDIR = os.path.join("test", "data", "annotate")
GENEPATH = os.path.join(DBDIR, "generated_by_unit-tests")


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


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert "[-d DB_PATH] -r RES_PATH [-l LIST_FILE] [-n NAME] [-Q]" in err
    assert "[--info FROM_INFO] [--prodigal] [--small] [--l90 L90]" in err
    assert "[--nbcont NBCONT] [--cutn CUTN] [--date DATE] [--tmp TMPDIR]" in err
    assert "[--annot_dir ANNOTDIR] [-F] [--threads THREADS] [-v]" in err
    assert "[-q] [-h]" in err
    assert "the following arguments are required: -r" in err


def test_parser_noname(capsys):
    """
    Test that when the script is called without any name for the genomes not -Q option,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-r respath".split())
    _, err = capsys.readouterr()
    assert ("You must specify your genomes dataset name in 4 characters with "
            "'-n name' option (type -h for more information). Or, if you do not want "
            "to annotate and format your genomes but just to run quality control, use "
            "option '-Q") in err


# def test_parser_wrongname(capsys):
#     """
#     Test that when the script is called with a genome name with more than 4 characters,
#     it returns an error message
#     """
#     parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
#     annot.build_parser(parser)
#     with pytest.raises(SystemExit):
#         annot.parse(parser, "list_file -d dbpath -r respath -n genome".split())
#     _, err = capsys.readouterr()
#     assert ("The genome name must contain 4 characters. For example, this name can "
#             "correspond to the 2 first letters of genus, and 2 first letters of "
#             "species, e.g. ESCO for Escherichia Coli.") in err


def test_parser_negative_cont(capsys):
    """
    Test that when the script is called with a limit of contig number <0,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-d dbpath -r respath -n g123 --nbcont -5".split())
    _, err = capsys.readouterr()
    assert "The maximum number of contigs allowed must be a positive number." in err


def test_parser_high_cont(capsys):
    """
    Test that when the script is called with a negative limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-Q -r respath -n g123 --nbcont 10005".split())
    _, err = capsys.readouterr()
    assert "We do not support genomes with more than 9999 contigs." in err


def test_parser_wrong_cont(capsys):
    """
    Test that when the script is called with a non integer limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-Q -r respath -n g123 --nbcont 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --nbcont: invalid int value: 10.5" in err


def test_parser_wrongl90(capsys):
    """
    Test that when the user does not give an int for the l90 limit, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-n TOTO -r respath -n g123 --l90 l90".split())
    _, err = capsys.readouterr()
    assert "argument --l90: invalid int value: 'l90'" in err


def test_parser_wrong_cut(capsys):
    """
    Test that when the user does not give an int for the cutN value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-Q -d dbpath -r respath -n g123 --cutn 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --cutn: invalid int value: '10.5'" in err


def test_parser_wrong_thread(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-Q -d dbpath -r respath -n g123 --threads 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --threads threads: invalid int value: 10.5" in err


def test_parser_wrong_date(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-Q -d dbpath -r respath -n g123 --date 417".split())
    _, err = capsys.readouterr()
    assert ("The date must contain 4 characters. Usually, it contains 4 digits, "
            "corresponding to the month (2 digits) and year (2 digits).") in err


def test_parser_q_and_v(capsys):
    """
    Test that when the user wants both quiet and verbose option, it gives an error message
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-d dbpath -r respath -n g123 -q -v".split())
    _, err = capsys.readouterr()
    assert ("Choose between a verbose output (-v) or a quiet output (-q). "
            "You cannot have both.") in err


def test_parser_default():
    """
    Test that when run with the minimum required arguments, all default values are
    as expected.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, "-r respath -n g123 -l list_genomes -d dbpath".split())
    assert options.list_file == "list_genomes"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "g123"
    assert options.l90 == 100
    assert options.nbcont == 999
    assert options.cutn == 5
    assert options.threads == 1
    assert options.date == time.strftime("%m%y")
    assert not options.force
    assert not options.qc_only
    assert not options.from_info
    assert not options.prodigal_only


def test_parser_values():
    """
    Test that values for L90, nbcontig, cutn, threads, date are taken into account
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, ("-l list_file -d dbpath -r respath -n g123 --l90 2 "
                                   "--nbcont 10 --cutn 0 --threads 2 --date toto "
                                   "--prodigal -F").split())
    assert options.list_file == "list_file"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "g123"
    assert options.l90 == 2
    assert options.nbcont == 10
    assert options.cutn == 0
    assert options.threads == 2
    assert options.date == "toto"
    assert options.force
    assert not options.qc_only
    assert not options.from_info
    assert options.prodigal_only


def test_parser_wrongforce(capsys):
    """
    Test that when run with '-F' option + a value, it returns an error message.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-d dbpath -r respath -n g123 -F 10".split())
    _, err = capsys.readouterr()
    assert "unrecognized arguments: 10" in err


def test_parser_qc():
    """
    Test that when run with '-Q' option (for QC only) and no name given for the genome, it
    is set to "NONE"
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, "-l list_file -d dbpath -r respath -Q".split())
    assert options.list_file == "list_file"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "NONE"
    assert options.l90 == 100
    assert options.nbcont == 999
    assert options.cutn == 5
    assert options.threads == 1
    assert options.date == time.strftime("%m%y")
    assert not options.force
    assert options.qc_only


def test_parser_info_cutn(capsys):
    """
    Test that when run with --info and --cutn x : error message
    If we run from info file, will not touch the sequences.
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-l list_file -d dbpath -r respath -n name "
                            "--info infofile --cutn 10".split())
    _, err = capsys.readouterr()
    assert ("If you provide a list of genomes with their calculated L90 and number of contigs, "
            "PanACoTA will use the given sequences as is. It will not cut them. So, you cannot "
            "use both --cutn and --info.") in err

def test_info_and_lstfile(capsys):
    """
    Test that there is an error message if user gives both -l infofile and --info LSTINFO
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-d dbpath -r respath -n name --nbcont 20 -l toto --info info".split())
    _, err = capsys.readouterr()
    assert ("Either you want to annotate raw sequences (name of files in '-l infofile') "
            "which will first go through the QC process, "
            "OR you already did QC on your sequences and just want to annotate them "
            "(information on those sequences in '--info LSTINFO-file'). "
            "Please choose one of these 2 possibilities.") in err


def test_parser_noinfo_nolist(capsys):
    """
    Test that when run without --info nor -l : error message
    Must provide one of them
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-d dbpath -r respath -n name --nbcont 20".split())
    _, err = capsys.readouterr()
    assert ("You must provide a list of genomes to annotate. Either raw genomes "
            "(see -l option), or genomes with quality information (see --info option).") in err


def test_parser_noinfo_nodbpath(capsys):
    """
    Test that when run without --info nor -d : error message
    Needs a path to genomes to annotate!
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-r respath -n name --nbcont 20 -l listgenomes".split())
    _, err = capsys.readouterr()
    assert ("You must provide a path to your database genome sequences (-d <db_path>). "
            "If you already have a LSTINFO file, it contains this db_path. Use it "
            "with --info <lstinfo file> option.") in err


def test_parser_info_dbpath(capsys):
    """
    Test that when run with both --info and -d : error message
    Must know which genomes to annotate between the 2
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-r respath -n name --info file -d dbpath".split())
    _, err = capsys.readouterr()
    assert ("If you run from your LSTINFO file, this one already contains the path "
            "of genomes to annotate. Remove -d <db_path> option.") in err


def test_parser_small_noprodigal(capsys):
    """
    Test that when run with both --small but do not ask to use prodigal, it returns error
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    with pytest.raises(SystemExit):
        annot.parse(parser, "-r respath -n name --small".split())
    _, err = capsys.readouterr()
    assert("You cannot use --small option with prokka. "
           "Either use prodigal, or remove this option") in err


def test_parser_filter(capsys):
    """
    Test that warnings are written (when will split l90 and/or nbcont)
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, "-l list_file -d dbpath -r respath -Q --l90 10".split())
    assert options.list_file == "list_file"
    assert options.db_path == "dbpath"
    assert options.res_path == "respath"
    assert options.name == "NONE"
    assert options.l90 == 10
    assert options.nbcont == 999
    assert options.cutn == 5
    assert options.threads == 1
    assert options.date == time.strftime("%m%y")
    assert not options.force
    assert options.qc_only
    stdout, _ = capsys.readouterr()
    assert(" !! Your genomes will be filtered, and only the ones with 'L90' <= 10 and 'number of contigs' < 999 "
           "will be kept. If you want to change those thresholds, use '--l90' and '--nbcont' options.") in stdout
    assert("! Your genomes will be split when sequence contains at least 5'N' in a row. "
           "If you want to change this threshold, see --cutn option.") in stdout



def test_parser_nosplit(capsys):
    """
    Test that warnings are written (when will split l90 and/or nbcont)
    """
    parser = argparse.ArgumentParser(description="Annotate all genomes", add_help=False)
    annot.build_parser(parser)
    options = annot.parse(parser, "--info infofile -r respath -Q".split())
    assert not options.list_file
    assert not options.db_path
    assert options.from_info == "infofile"
    assert options.res_path == "respath"
    assert options.name == "NONE"
    assert options.l90 == 100
    assert options.nbcont == 999
    assert options.cutn == 5
    assert options.threads == 1
    assert options.date == time.strftime("%m%y")
    assert not options.force
    assert options.qc_only
    stdout, _ = capsys.readouterr()
    assert (" !! Your sequences will be used as is by PanACoTA. Be sure you already split your sequences at each row of X 'N' if needed.") in stdout
    assert ("PanACoTA will use the values (L90, nbcont) given in your info file. It will ignore the genomes for which those values are incorrect. It will also ignore genomes with more than 999 contigs.")
