#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of 'prepare' subcommand
"""
import argparse
import pytest

from PanACoTA.subcommands import prepare


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Prepare genomes", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "".split())
    _, err = capsys.readouterr()
    print(err)
    assert "error: " in err
    assert ("As you did not put the '--norefseq' nor the '-M' option, it means that you want "
            "to download refseq (or genbank) genomes. But you did not provide any information, so PanACoTA "
            "cannot guess which species you want to download. Specify NCBI_taxid (-t)") in err
    assert ("NCBI species taxid (-T) and/or NCBI_species (-g) to download, "
            "or add one of the 2 options (--norefseq or -M) "
            "if you want to skip the 'download step'.") in err


def test_cutn_noint(capsys):
    """
    Test that when user is giving a number of 'N' from which to cut which is:
     - not a number
     - not an int
     - <0
    it gives error message
    """
    parser = argparse.ArgumentParser(description="Prepare genomes", add_help=False)
    prepare.build_parser(parser)
    # Not a number
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--cutn ten".split())
    _, err = capsys.readouterr()
    assert "error: argument --cutn: invalid int value: 'ten'" in err
    # Not an int
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--cutn 1.5".split())
    _, err = capsys.readouterr()
    assert "error: argument --cutn: invalid int value: '1.5'" in err
    # Negative number
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--cutn -5".split())
    _, err = capsys.readouterr()
    assert "error: argument --cutn must be a positive integer: invalid int value: '-5'" in err


def test_l90_noint(capsys):
    """
    Test that when user is giving a number for max L90 which is not valid:
     - not a number
     - not an int
    it gives error message
    """
    parser = argparse.ArgumentParser(description="Prepare genomes", add_help=False)
    prepare.build_parser(parser)
    # Not a number
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--l90 ten".split())
    _, err = capsys.readouterr()
    assert "error: argument --l90: invalid int value: 'ten'" in err
    # Not an int
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--l90 1.5".split())
    _, err = capsys.readouterr()
    assert "error: argument --l90: invalid int value: '1.5'" in err


def test_parser_negative_cont(capsys):
    """
    Test that when the script is called with a limit of contig number <0,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--nbcont -5".split())
    _, err = capsys.readouterr()
    assert "The maximum number of contigs allowed must be a positive number." in err


def test_parser_high_cont(capsys):
    """
    Test that when the script is called with a negative limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--nbcont 10005".split())
    _, err = capsys.readouterr()
    assert "We do not support genomes with more than 9999 contigs." in err


def test_parser_wrong_cont(capsys):
    """
    Test that when the script is called with a non integer limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--nbcont 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --nbcont: invalid int value: 10.5" in err


def test_parser_wrong_level(capsys):
    """
    Test that when the script is called with a non integer limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-t 1234 -o toto -l toto".split())
    _, err = capsys.readouterr()
    assert ("Please choose between available assembly levels: 'all', 'complete', "
            "'chromosome', 'scaffold', 'contig'. If several levels, provide a "
            "comma-separated list. Invalid value: 'toto'") in err


def test_parser_wrong_level_notcomma(capsys):
    """
    Test that when the script is called with a non integer limit of contig number,
    it returns an error message
    """
    parser = argparse.ArgumentParser(description="prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-t 1234 -o outdir -l complete.scaffold".split())
    _, err = capsys.readouterr()
    assert ("Please choose between available assembly levels: 'all', 'complete', "
            "'chromosome', 'scaffold', 'contig'. If several levels, provide a "
            "comma-separated list. Invalid value: 'complete.scaffold'") in err


def test_max_mash_dist(capsys):
    """
    Test that when user is giving a number for max_dist which is not valid:
     - not a number
     - > 1
     - <0
    it gives error message
    """
    parser = argparse.ArgumentParser(description="Prepare genomes", add_help=False)
    prepare.build_parser(parser)
    # Not a number
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--max_dist ten".split())
    _, err = capsys.readouterr()
    assert "error: mash distance: invalid float value: 'ten'" in err
    # > 1
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--max_dist 1.5".split())
    _, err = capsys.readouterr()
    assert "error: mash distance must be between 0 and 1: invalid value: '1.5'" in err
    # < 0
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--max_dist -0.5".split())
    _, err = capsys.readouterr()
    assert "error: mash distance must be between 0 and 1: invalid value: '-0.5'" in err


def test_min_mash_dist(capsys):
    """
    Test that when user is giving a number for max_dist which is not valid:
     - not a number
     - > 1
     - <0
    it gives error message
    """
    parser = argparse.ArgumentParser(description="Prepare genomes", add_help=False)
    prepare.build_parser(parser)
    # Not a number
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--min_dist ten".split())
    _, err = capsys.readouterr()
    assert "error: mash distance: invalid float value: 'ten'" in err
    # > 1
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--min_dist 1.5".split())
    _, err = capsys.readouterr()
    assert "error: mash distance must be between 0 and 1: invalid value: '1.5'" in err
    # < 0
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--min_dist -0.5".split())
    _, err = capsys.readouterr()
    assert "error: mash distance must be between 0 and 1: invalid value: '-0.5'" in err


def test_min_sup_max(capsys):
    '''
    Test that we get an error message if min_dist > max_dist
    '''
    parser = argparse.ArgumentParser(description="Prepare genomes", add_help=False)
    prepare.build_parser(parser)
    # Not a number
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--min_dist 0.9 --max_dist=0.8 --norefseq -o toto".split())
    _, err = capsys.readouterr()
    assert "min_dist (0.9) cannot be higher than max_dist (0.8)" in err


def test_parser_wrong_thread(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-p 10.5".split())
    _, err = capsys.readouterr()
    assert "argument --threads threads: invalid int value: 10.5" in err
    # Negative number of threads
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-p -1".split())
    _, err = capsys.readouterr()
    assert ("Please provide a positive number of threads (or 0 for all threads): "
            "Invalid value: -1") in err


def test_parser_more_threads(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    import multiprocessing
    nb_cpu = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-p 30000".split())
    _, err = capsys.readouterr()
    assert (f"You have {nb_cpu} threads on your computer, you cannot ask for more: "
            "invalid value: 30000") in err


def test_parser_all_threads(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    import multiprocessing
    nb_cpu = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    options = prepare.parse(parser, "-p 0 --norefseq -o toto".split())
    assert options.parallel == nb_cpu
    assert options.norefseq == True
    assert options.only_mash == False


def test_parse_missing_arg(capsys):
    """
    running prepare without NCBI info nor mash_only nor norefseq -> error asking one of those
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-p 1".split())
    _, err = capsys.readouterr()
    assert ("As you did not put the '--norefseq' nor the '-M' option, it means that you want "
            "to download refseq (or genbank) genomes. But you did not provide any information, so PanACoTA "
            "cannot guess which species you want to download. Specify NCBI_taxid (-t)") in err
    assert ("NCBI species taxid (-T) and/or NCBI_species (-g) to download, "
            "or add one of the 2 options (--norefseq or -M) "
            "if you want to skip the 'download step'.") in err


def test_norefseq_nooutdir(capsys):
    """
    Try running without refseq, but not giving an output directory
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "--norefseq".split())
    _, err = capsys.readouterr()
    assert ("You must provide an output directory, where your results will be saved.") in err


def test_onlymash_noinfo(capsys):
    """
    Try running without refseq, but not giving an output directory
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-M".split())
    _, err = capsys.readouterr()
    assert ("If you want to run only Mash filtering steps, please give the info file with "
            "the required information (see '--info' option") in err


def test_onlymash_nooutdir(capsys):
    """
    Try running without refseq, but not giving an output directory
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-M --info toto ".split())
    _, err = capsys.readouterr()
    assert ("If you want to run only Mash filtering steps, please give the output "
            "directory where you want to save your results (see '-o' option)") in err


def test_verbose_quiet(capsys):
    """
    Try running without refseq, but not giving an output directory
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-q -vv -M --info toto -o outdir".split())
    _, err = capsys.readouterr()
    assert ("Choose between a verbose output (-v) or a quiet output (-q). "
            "You cannot have both.") in err


def test_parser_nospecies(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    options = prepare.parse(parser, "-T 1234".split())
    assert not options.norefseq
    assert not options.only_mash
    assert options.ncbi_species_taxid == "1234"
    assert options.ncbi_taxid == ""
    assert options.ncbi_species_name == ""
    out, err = capsys.readouterr()
    assert ("WARNING: you did not provide a species name ('-g species' option) "
            "nor an output directory ('-o outdir'). "
            "All files will be downloaded in a folder called with the NCBI species "
            "taxid 1234 instead of the species name.") in out


def test_parser_nospecies_nospeid(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    options = prepare.parse(parser, "-t 1234".split())
    assert not options.norefseq
    assert not options.only_mash
    assert options.ncbi_species_taxid == ""
    assert options.ncbi_taxid == "1234"
    assert options.ncbi_species_name == ""
    out, err = capsys.readouterr()
    assert ("WARNING: you did not provide a species name ('-g species' option) "
            "nor a species taxid ('-T spetaxid') nor an output directory ('-o outdir'). "
            "All files will be downloaded in a folder called with the NCBI "
            "taxid 1234.") in out


def test_parser_default_cutn(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    options = prepare.parse(parser, "-t 1234 -o outdir -g species".split())
    assert not options.norefseq
    assert not options.only_mash
    assert options.ncbi_taxid == "1234"
    assert options.ncbi_species_name == "species"
    out, err = capsys.readouterr()
    assert ("!! Your genomes will be split when sequence contains at "
            "least 5'N' in a row. If you want to change this threshold, use "
            "'--cutn n' option (n=0 if you do not want to cut)") in out


def test_parser_error_section_name(capsys):
    """
    Test when we give another value than refseq or genbank to -s section option
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    with pytest.raises(SystemExit):
        prepare.parse(parser, "-t 1234 -o outdir -s species --cutn 1".split())
    _, err = capsys.readouterr()
    assert ("argument -s: invalid choice: 'species' (choose from 'refseq', 'genbank')") in err


def test_parser_default_l90_nb_cont(capsys):
    """
    Test that when the user does not give an int for the threads value, it returns an
    error message.
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    options = prepare.parse(parser, "-t 1234 -o outdir -g species --cutn 1".split())
    assert not options.norefseq
    assert not options.only_mash
    assert options.ncbi_taxid == "1234"
    assert options.ncbi_species_name == "species"
    out, err = capsys.readouterr()
    assert ("!! Your genomes will be filtered, and only the ones with 'L90' <= 100 "
            "and 'number of contigs' < 999 will be kept. If you want to change those "
            "thresholds, use '--l90' and '--nbcont' options.") in out


def test_parser_info_notonlymash(capsys):
    """
    Giving an info file, but not asking for only_mash -> useless info file
    """
    parser = argparse.ArgumentParser(description="Prepare", add_help=False)
    prepare.build_parser(parser)
    options = prepare.parse(parser, "-T 1234 -o outdir -g species --cutn 1 --info toto".split())
    assert not options.norefseq
    assert not options.only_mash
    assert options.ncbi_species_taxid == "1234"
    assert options.ncbi_species_name == "species"
    out, err = capsys.readouterr()
    assert ("!! You gave an info file (--info option), but did not ask to run only Mash "
            "step (-M option). Your info file will be ignored (and renamed with '.back' "
            "at the end), and another one will be created with the new calculated values.") in out
