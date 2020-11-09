#!/usr/bin/env python3
# coding: utf-8

"""
Functional tests for the parser of tree subcommand
"""
import argparse
import pytest

from PanACoTA.subcommands import tree


def test_parser_noarg(capsys):
    """
    Test that when the script is called without any argument, an error message appears,
    indicating the required arguments.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "".split())
    _, err = capsys.readouterr()
    assert "usage: " in err
    assert "-a ALIGNMENT -o OUTDIR" in err
    assert "[-s {fasttree,fastme,quicktree,iqtree,iqtree2}] [-b BOOT]" in err
    assert "[--threads THREADS] [-m MODEL]" in err
    assert "[-B] [--mem MEMORY" in err
    assert "[-v]" in err
    assert "[-q] [-h]" in err
    assert "the following arguments are required: -a, -o" in err


def test_parser_threadnotint(capsys):
    """
    Test that when the number of threads given is not an int, it returns an error message
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align -o outdir --threads 1.2".split())
    _, err = capsys.readouterr()
    assert "argument --threads threads: invalid int value: 1.2" in err


def test_parser_thread_toomany(capsys):
    """
    Test that when the number of threads given is higher than the total number of threads,
    it returns the expected error message
    """
    import multiprocessing
    nb = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align --threads {} -o outdir".format(nb + 3).split())
    _, err = capsys.readouterr()
    assert ("You have {} threads on your computer, you cannot ask for more: "
            "invalid value: {}".format(nb, nb+3)) in err


def test_parser_thread_neg(capsys):
    """
    Test that when the number of threads given is a negative number, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align --threads -5".split())
    _, err = capsys.readouterr()
    assert ("Please provide a positive number of threads (or 0 for all threads): "
            "Invalid value: -5") in err


def test_parser_boot_notenough(capsys):
    """
    Test that when the number of threads given is a negative number, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align -b 10 -o outdir".split())
    _, err = capsys.readouterr()
    assert ("With IQtree, number of replicates for bootstraps must be >= 1000") in err


def test_parser_quicktree_parallel(capsys):
    """
    Test that when soft is quicktree, and e ask for more than 1 thread, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align -s quicktree --threads 2 -o outdir".split())
    _, err = capsys.readouterr()
    assert ("You cannot run quicktree with multiple threads. Choose another software, or remove "
            "the --threads option.") in err


def test_parser_quicktree_model(capsys):
    """
    Test that when soft is quicktree, and we ask for a specific model, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-o outdir -a align -s quicktree -m F84".split())
    _, err = capsys.readouterr()
    assert ("Quicktree only runs the NJ algorithm. You cannot choose a DNA substitution "
            "model.") in err


def test_parser_quicktree_writeboot(capsys):
    """
    Test that when soft is quicktree, and we ask for a specific model, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-o outdir -a align -s quicktree -B".split())
    _, err = capsys.readouterr()
    assert "'-B' option is only available with FastME and IQtree." in err


def test_parser_fastme_wrongmodel(capsys):
    """
    Test that when soft is fastme, and we ask for GTR model, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align -o outdir -s fastme -m GTR".split())
    _, err = capsys.readouterr()
    assert ("GTR is not an available model for fastme. Please choose an available DNA model (see "
            "-h for more details)") in err


def test_parser_fasttree_wrongmodel(capsys):
    """
    Test that when soft is fasttree, and we ask for RY model, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align -o outdir -s fasttree -m RY".split())
    _, err = capsys.readouterr()
    assert ("RY is not an available model for fasttree. Please choose an available DNA model (see "
            "-h for more details)") in err


def test_parser_iqtree_wrongmodel(capsys):
    """
    Test that when soft is fasttree, and we ask for RY model, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-a align -o outdir -s iqtree -m toto".split())
    _, err = capsys.readouterr()
    assert ("toto is not an available model for iqtree. Please choose an available DNA model (see "
            "-h for more details)") in err


def test_parser_fasttree_writeboot(capsys):
    """
    Test that when soft is fasttree, and we ask to write bootstraps, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-o outdir -a align -s fasttree -B".split())
    _, err = capsys.readouterr()
    assert "'-B' option is only available with FastME and IQtree" in err


def test_parser_fasttree_fast(capsys):
    """
    Test that when soft is fasttree, and we ask to write bootstraps, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-o outdir -a align -s fasttree -fast".split())
    _, err = capsys.readouterr()
    assert "-fast option is available only for IQtree, and not compatible with '-B' and '-b' options (bootstraps)" in err


def test_parser_fastme_memory(capsys):
    """
    Test that when soft is fastme, and we ask for a memory amount, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-o outdir -a align -s fastme --mem 10GB".split())
    _, err = capsys.readouterr()
    assert "'--mem' option is only available for IQtree" in err


def test_parser_iqtree_boot_fast(capsys):
    """
    Test that when soft is fasttree, and we ask to write bootstraps, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-o outdir -a align -b 1000 -fast".split())
    _, err = capsys.readouterr()
    assert "-fast option is available only for IQtree, and not compatible with '-B' and '-b' options (bootstraps)" in err


def test_parser_iqtree_writeboot_fast(capsys):
    """
    Test that when soft is fasttree, and we ask to write bootstraps, it returns the expected
    error message.
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    with pytest.raises(SystemExit):
        tree.parse(parser, "-o outdir -a align -B -fast".split())
    _, err = capsys.readouterr()
    assert "-fast option is available only for IQtree, and not compatible with '-B' and '-b' options (bootstraps)" in err


def test_parser_default():
    """
    Test that when giving only align file, it returns the expected values for all other arguments
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    args = tree.parse(parser, "-a align -o outdir".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "iqtree2"
    assert args.model == "GTR"
    assert args.write_boot is False
    assert args.threads == 1
    assert args.verbose == 0
    assert args.quiet == False


def test_parser_default_fastme():
    """
    Test that when giving only align file, and to use fastme, it returns the expected values for
    all other arguments (with default DNA model)
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    args = tree.parse(parser, "-a align -o outdir -s fastme".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "fastme"
    assert args.model == "T"
    assert args.write_boot is False
    assert args.threads == 1
    assert args.verbose == 0
    assert args.quiet == False


def test_parser_default_fasttree():
    """
    Test that when giving only align file, and to use fastme, it returns the expected values for
    all other arguments (with default DNA model)
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    args = tree.parse(parser, "-a align -o outdir -s fasttree".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "fasttree"
    assert args.model == "-gtr"
    assert args.write_boot is False
    assert args.threads == 1
    assert args.verbose == 0
    assert args.quiet == False


def test_parser_all_threads():
    """
    Test that when giving alignment file, and 0 to use all threads, it returns expected values
    for all arguments
    """
    import multiprocessing
    nb = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    args = tree.parse(parser, "-a align -o outdir --threads 0".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "iqtree2"
    assert args.model == "GTR"
    assert args.write_boot is False
    assert args.threads == nb
    assert args.verbose == 0
    assert args.quiet == False


def test_parser_threads_ok():
    """
    Test that when giving align file and a correct number of threads (positive int), it returns the
    expected values for all arguments
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    import multiprocessing
    nb = multiprocessing.cpu_count()
    args = tree.parse(parser, f"-a align -o outdir --threads {nb}".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "iqtree2"
    assert args.model == "GTR"
    assert args.write_boot is False
    assert args.threads == nb
    assert args.verbose == 0
    assert args.quiet == False


def test_parser_fastme_modelkey():
    """
    Test that when giving align file, using fastme, with p-distance model, it returns expected
    values for all arguments
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    args = tree.parse(parser, "-a align -o outdir -s fastme -m p-distance".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "fastme"
    assert args.model == "p"
    assert args.write_boot is False
    assert args.threads == 1
    assert args.verbose == 0
    assert args.quiet == False


def test_parser_fastme_modelval():
    """
    Test that when giving align file, using fastme, with 'Y' model, it returns expected
    values for all arguments
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    args = tree.parse(parser, "-a align -o outdir -s fastme -m Y".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "fastme"
    assert args.model == "Y"
    assert args.write_boot is False
    assert args.threads == 1
    assert args.verbose == 0
    assert args.quiet == False


def test_parser_fasttree_jc():
    """
    Test that when giving align file, using fasttree, with JC model, it returns expected
    values for all arguments
    """
    parser = argparse.ArgumentParser(description="Tree", add_help=False)
    tree.build_parser(parser)
    args = tree.parse(parser, "-a align -o outdir -s fasttree -m JC".split())
    assert args.alignment == "align"
    assert args.boot is None
    assert args.outdir == "outdir"
    assert args.soft == "fasttree"
    assert args.model == ""
    assert args.write_boot is False
    assert args.threads == 1
    assert args.verbose == 0
    assert args.quiet == False
