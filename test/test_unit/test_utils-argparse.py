#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the functions of utils-argparse.py dealing with checking arguments given to argparse
"""

import logging
import os
import pytest
import shutil
import argparse
import multiprocessing

import PanACoTA.utils_argparse as autils

CONFFILE = os.path.join("test", "data", "utils", "configfile.ini")

def test_gen_name():
    """
    Test that, when giving a gene name, if it does not have 4 characters (letters or num),
    it returns an error message
    """
    assert autils.gen_name("TOTO") == "TOTO"
    assert autils.gen_name("1234") == "1234"
    assert autils.gen_name("T1O2") == "T1O2"
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.gen_name("a long name")
    assert ("The genome name must contain 4 characters. For example, this name "
            "can correspond to the 2 first letters of genus, and 2 first letters of "
            "species, e.g. ESCO for Escherichia Coli") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.gen_name("-gdd")
    assert ("The genome name must contain 4 characters. For example, this name can correspond "
            "to the 2 first letters of genus, and 2 first letters of species, e.g. "
            "ESCO for Escherichia Coli") in str(err.value)


def test_date_name(capsys):
    """
    Test that when a date is given, it returns expected message
    """
    assert autils.date_name("0920") == "0920"
    assert autils.date_name("se20") == "se20"
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.date_name("september 2020")
    assert ("The date must contain 4 characters. Usually, it contains 4 digits, "
            "corresponding to the month (2 digits) and year (2 digits)") in str(err.value)


def test_get_date():
    """
    test that it returns current date
    """
    import time
    t = time.strftime("%m%y")
    assert autils.get_date() == t


def test_cont_num():
    """
    Test that given value of contig number is valid
    """
    assert autils.cont_num(10) == 10
    with pytest.raises(argparse.ArgumentTypeError) as raised_err:
        autils.cont_num("-2")
    assert ("The maximum number of contigs allowed must be a "
            "positive number") in str(raised_err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.cont_num("10000")
    assert ("We do not support genomes with more than 9999 contigs") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.cont_num("a")
    assert ("argument --nbcont: invalid int value: a") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        autils.cont_num("1.1")
    assert ("argument --nbcont: invalid int value: 1.") in str(err.value)


def test_thread_num():
    """
    Test that given number of threads is as expected
    """
    assert autils.thread_num("1") == 1
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num("1.1")
    assert ("argument --threads threads: invalid int value: 1.1") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num("a")
    assert ("argument --threads threads: invalid int value: a") in str(err.value)
    try:
        nb_cpu = len(os.sched_getaffinity(0))
    except AttributeError:
        nb_cpu = multiprocessing.cpu_count()
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num(str(nb_cpu*2))
    assert (f"You have {nb_cpu} threads on your computer, you cannot ask for more: invalid value: "
            f"{nb_cpu*2}") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.thread_num("-1")
    assert ("Please provide a positive number of threads (or 0 for all threads): "
            "Invalid value: -1") in str(err.value)
    assert autils.thread_num(0) == nb_cpu


def test_positive_int():
    """
    Test checking that given argument is a positive integer
    """
    assert autils.positive_int("1") == 1
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.positive_int("1.1")
    assert ("argument --cutn: invalid int value: '1.1'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.positive_int("-1")
    assert ("error: argument --cutn must be a positive integer: "
            "invalid int value: '-1'") in str(err.value)


def test_mash_dist():
    """
    Test checking that given value is ok for a mash distance
    """
    assert autils.mash_dist("0.05") == 0.05
    assert autils.mash_dist("1e-4") == 0.0001
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("1.1.1")
    assert ("error: mash distance: invalid float value: '1.1.1'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("one")
    assert ("error: mash distance: invalid float value: 'one'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("1.000001")
    assert ("error: mash distance must be between 0 and 1: "
            "invalid value: '1.000001'") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.mash_dist("-1e-4")
    assert ("error: mash distance must be between 0 and 1: "
            "invalid value: '-0.0001'") in str(err.value)


def test_percentage():
    """
    Test value given for parameter -t tol:
    - string corresponding to a float or an int -> returns float value
    - a non digital param for a percentage raises appropriate error
    - a negative number -> appropriate error
    - number>1 -> appropriate error
    """
    assert autils.percentage("0.5") == 0.5
    assert type(autils.percentage("0.5")) == float
    assert autils.percentage("1e-1") == 0.1
    assert autils.percentage("1") == 1.0
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.percentage("one")
    assert ("argument -t tol: invalid float value: one") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.percentage(-0.5)
    assert ("The minimum %% of genomes required in a family to be persistent must "
            "be in [0, 1]. Invalid value: -0.5") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.percentage("1.1")
    assert ("The minimum %% of genomes required in a family to be persistent must "
            "be in [0, 1]. Invalid value: 1.1") in str(err.value)


def test_perc_id():
    """
    Same as test_percentage, but for parameter -i percentage_id
    """
    assert autils.perc_id("0.5") == 0.5
    assert type(autils.perc_id("0.5")) == float
    assert autils.perc_id("1e-1") == 0.1
    assert autils.perc_id("1") == 1.0
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.perc_id("one")
    assert ("argument -i percentage_id: invalid float value: one") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.perc_id(-0.5)
    assert ("The minimum %% of identity must be in [0, 1]. Invalid value: -0.5") in str(err.value)
    with pytest.raises(argparse.ArgumentTypeError) as err:
        a = autils.perc_id("1.1")
    assert ("The minimum %% of identity must be in [0, 1]. Invalid value: 1.1") in str(err.value)


def test_conf_parser_init_empty(capsys):
    """
    test class Conf_all_parser init when no config file or empty config file
    """
    # No conf file, no section
    c = autils.Conf_all_parser("")
    assert c.sec_dicts == {}
    # No conf file, empty section list
    c = autils.Conf_all_parser("", [])
    assert c.sec_dicts == {}
    # No conf file, sections
    c = autils.Conf_all_parser("", ["sec1", "sec2"])
    assert c.sec_dicts == {"sec1": {}, "sec2": {}}

    confdir = os.path.join("test", "data", "generated_by_utils")
    os.makedirs(confdir)
    conffile = os.path.join(confdir, "conf.ini")
    open(conffile, "w").close()
    # Empty conffile, no section
    c = autils.Conf_all_parser(conffile)
    assert c.sec_dicts == {}
    # Empty conffile, sections
    c = autils.Conf_all_parser(conffile, ["sec1", "sec2"])
    assert c.sec_dicts == {"sec1": {}, "sec2": {}}
    shutil.rmtree(confdir)


def test_conf_parser_init_notexisting(capsys):
    """
    test class Conf_all_parser init when given configfile does not exist
    """
    with pytest.raises(SystemExit):
        autils.Conf_all_parser("not_existing_conffile.ini")
    out, err = capsys.readouterr()
    assert "Error: config file not_existing_conffile.ini not found." in out


def test_conf_parser_init_duplicate(capsys):
    """
    test class Conf_all_parser configfile contains a duplicated parameter
    -> raises error
    """
    conffile = os.path.join("test", "data", "utils", "configfile-dup.ini")
    with pytest.raises(SystemExit):
        autils.Conf_all_parser(conffile)
    out, err = capsys.readouterr()
    assert "option 'bool_y' in section 'sec_bool' already exists" in out


def test_conf_parser_init():
    """
    Test config parser with a given config file. Check value of defaults etc.
    """
    # configfile but no section
    c = autils.Conf_all_parser(CONFFILE)
    assert c["sec1"]["toto"] == "parameter"
    assert c["sec1"]["param1"] == "10"
    assert c["sec2"]["param1"] == "3"
    assert c["sec3"]["param1"] == "3"
    assert c["sec1"]["param2"] == "10"
    assert c["sec2"]["param2"] == ''
    assert c["sec3"]["param2"] == "10"
    assert c.sec_dicts == {}

    # configfile and 2 sections among 3 existing in configfile
    c = autils.Conf_all_parser(CONFFILE, ["sec1", "sec3"])
    assert c["sec1"]["toto"] == "parameter"
    assert c["sec1"]["param1"] == "10"
    assert c["sec2"]["param1"] == "3"
    assert c["sec2"]["myval"] == "parameter"
    assert c["sec3"]["param1"] == "3"
    assert c["sec1"]["param2"] == "10"
    assert c["sec2"]["param2"] == ''
    assert c["sec3"]["param2"] == "10"
    assert c.sec_dicts == {"sec1": {"param1": "10", "param2": "10", 'sec1p': "",
                                    "toto": "parameter", "titi": "my value"},
                           "sec3": {"param1": "3", "param2": "10"}}

    # configfile 2 sections, 1 not in configfile
    c = autils.Conf_all_parser(CONFFILE, ["sec1", "sec4"])
    assert c["sec1"]["toto"] == "parameter"
    assert c["sec1"]["param1"] == "10"
    assert c["sec2"]["param1"] == "3"
    assert c["sec4"]["param1"] == "3"  # created sec4 with default parameters
    assert c["sec4"]["param2"] == "10"  # created sec4 with default parameters
    assert not "toto" in c["sec4"]  # created sec4 with default parameters
    assert c.sec_dicts == {"sec1": {"param1": "10", "param2": "10", 'sec1p': "",
                                    "toto": "parameter", "titi": "my value"},
                           "sec4": {}}  # But sec4 dict is empty (no param given in configfile)


def test_conf_parser_init_noCleanStr():
    """
    Test config parser with a given config file. Do not clean str (keep surounding "")
    """
    # configfile but no section
    cfg_file = os.path.join("test", "data", "utils", "configfile-str.ini")
    c = autils.Conf_all_parser(cfg_file, clean_str=False)
    assert c["sec2"]["param1"] == "3"
    assert c["sec2"]["param2"] == ''
    assert c["sec2"]["myval"] == '"parameter"'
    assert c["sec2"]["param3"] == 'myparameter'
    assert c["sec3"]["param1"] == "3"
    assert c["sec3"]["param2"] == "10"
    assert c["sec_bool"]["bool num_false"] == "0"
    assert c["sec_bool"]["bool_t"] == "ON"
    assert c["sec_bool"]["bool_f"] == '"off"'
    assert c.sec_dicts == {}


def test_conf_parser_clean_str(capsys):
    """
    Check that surounding "" are removed
    """
    cfg_file = os.path.join("test", "data", "utils", "configfile-str.ini")
    c = autils.Conf_all_parser(cfg_file, ["sec2", "sec_bool"])
    assert c.get_section_dict("sec2") == {"param2": "", "myval": "parameter", 
                                          'param3': "myparameter",
                                          "param1": "3", "paramint": "2"}
    assert c.get_section_dict("sec_bool") == {"bool num_false": "0", "bool_num_true": "1", 
                                              "bool_f": "off", "bool_t": "ON", "bool_n": "no",
                                              "bool_y": "YES", 
                                              "bool_false": "FalSe", "bool_true": "tRUE",
                                              "param1": "3", "param2": "10"}


def test_conf_parser_get_section(capsys):
    """
    Test get dict of values for a given section
    """
    c = autils.Conf_all_parser(CONFFILE, ["sec1", "sec4"])
    # Get sec1 dict
    assert c.get_section_dict("sec1") == {"param1": "10", "param2": "10", 'sec1p': "",
                                          "toto": "parameter", "titi": "my value"}
    assert c.get_section_dict("sec4") == {}
    # Try to get sec2 dict, but does not exist
    with pytest.raises(SystemExit):
        c.get_section_dict("sec2")
    out, err = capsys.readouterr()
    assert ('No section sec2 in test/data/utils/configfile.ini') in out


def test_conf_parser_add_default(capsys):
    """
    Test add default parameters to config parser.
    If parameter given already exists, do nothing. If does not exist, add it with given value
    """
    c = autils.Conf_all_parser(CONFFILE, ["sec1", "sec4"])
    defaults = {"param1": "55", "defparam": "default value"}
    # Add default parameters to sec1
    c.add_default(defaults, "sec1")
    assert c.get_section_dict("sec1") == {"param1": "10", "param2": "10", 'sec1p': "",
                                          "toto": "parameter", "titi": "my value",
                                          "defparam": "default value"}
    assert c.get_section_dict("sec4") == {}
    c["sec1"]["param1"] == "10"
    c["sec1"]["defparam"] == "default value"
    c["sec4"]["param1"] == "3"   # sec4 has default parameter found in configfile

    # Add default parameters to sec4 (sec1 already have them)
    c.add_default(defaults, "sec4")
    assert c.get_section_dict("sec1") == {"param1": "10", "param2": "10", 'sec1p': "",
                                          "toto": "parameter", "titi": "my value",
                                          "defparam": "default value"}
    assert c.get_section_dict("sec4") == {"param1": "55", "defparam": "default value"}
    c["sec1"]["param1"] == "10"
    c["sec1"]["defparam"] == "default value"
    c["sec4"]["param1"] == "55"


def test_conf_parser_update(capsys):
    """
    Test update section
    Update current parameters with values given
    """
    c1 = autils.Conf_all_parser(CONFFILE, ["sec1", "sec4"])
    update = {"param1": "55", "defparam": "default value"}
    # Add default parameters to sec1
    c1.update(update, "sec1")
    assert c1.get_section_dict("sec1") == {"param1": "55", "param2": "10", 'sec1p': "",
                                          "toto": "parameter", "titi": "my value",
                                          "defparam": "default value"}
    assert c1.get_section_dict("sec4") == {}
    c1["sec1"]["param1"] == "10"
    c1["sec1"]["defparam"] == "default value"
    c1["sec4"]["param1"] == "3"   # sec4 has default parameter found in configfile

    # Add default parameters to sec1
    c2 = autils.Conf_all_parser(CONFFILE, ["sec1", "sec4"])
    c2.update(update, "sec4")
    assert c2.get_section_dict("sec1") == {"param1": "10", "param2": "10", 'sec1p': "",
                                          "toto": "parameter", "titi": "my value"}
    assert c2.get_section_dict("sec4") == {"param1": "55", "defparam": "default value"}
    c2["sec1"]["param1"] == "10"
    c2["sec4"]["param1"] == "55"   # sec4 has default parameter found in configfile


def test_conf_parser_setbool(capsys):
    """
    try to convert a given parameter to a boolean
    """
    c1 = autils.Conf_all_parser(CONFFILE, ["sec_bool"])
    # 0/1 to False/True
    c1.set_boolean("sec_bool", "bool num_false")
    assert c1["sec_bool"]["bool num_false"] == "0"
    assert c1.sec_dicts["sec_bool"]["bool num_false"] == False
    c1.set_boolean("sec_bool", "bool_num_true")
    assert c1["sec_bool"]["bool_num_true"] == "1"
    assert c1.sec_dicts["sec_bool"]["bool_num_true"] == True

    # off/on to False/True
    c1.set_boolean("sec_bool", "bool_f")
    assert c1["sec_bool"]["bool_f"] == "off"
    assert c1.sec_dicts["sec_bool"]["bool_f"] == False
    c1.set_boolean("sec_bool", "bool_t")
    assert c1["sec_bool"]["bool_t"] == "ON"
    assert c1.sec_dicts["sec_bool"]["bool_t"] == True

    # no/yes to False/True
    c1.set_boolean("sec_bool", "bool_n")
    assert c1["sec_bool"]["bool_n"] == "no"
    assert c1.sec_dicts["sec_bool"]["bool_n"] == False
    c1.set_boolean("sec_bool", "bool_y")
    assert c1["sec_bool"]["bool_y"] == "YES"
    assert c1.sec_dicts["sec_bool"]["bool_y"] == True

    # false/true to False/True
    c1.set_boolean("sec_bool", "bool_false")
    assert c1["sec_bool"]["bool_false"] == "FalSe"
    assert c1.sec_dicts["sec_bool"]["bool_false"] == False
    c1.set_boolean("sec_bool", "bool_true")
    assert c1["sec_bool"]["bool_true"] == "tRUE"
    assert c1.sec_dicts["sec_bool"]["bool_true"] == True

    # With a boolean called with "", not cleaned
    cfg_file = os.path.join("test", "data", "utils", "configfile-str.ini")
    c_str = autils.Conf_all_parser(cfg_file, ["sec_bool"], clean_str=False)
    with pytest.raises(SystemExit):
        c_str.set_boolean("sec_bool", "bool_f")
    out, err = capsys.readouterr()
    assert ('ERROR: bool_f must be a boolean. Wrong value: "off".') in out

    # With a boolean called with "" but cleaned
    cfg_file = os.path.join("test", "data", "utils", "configfile-str.ini")
    c_str = autils.Conf_all_parser(cfg_file, ["sec_bool"])
    c_str.set_boolean("sec_bool", "bool_f")
    assert c1["sec_bool"]["bool_f"] == 'off'
    assert c1.sec_dicts["sec_bool"]["bool_f"] == False

    # error
    with pytest.raises(SystemExit):
        c1.set_boolean('sec1', "param1")
    out, err = capsys.readouterr()
    assert ('ERROR: param1 must be a boolean. Wrong value: 10.') in out


def test_conf_parser_setint(capsys):
    """
    try to convert a given parameter to an int
    """
    c1 = autils.Conf_all_parser(CONFFILE, ["sec1"])
    c1.set_int("sec1", "param1")
    assert c1["sec1"]["param1"] == "10"
    assert c1.sec_dicts["sec1"]["param1"] == 10

    # With an int called with "" but cleaned
    cfg_file = os.path.join("test", "data", "utils", "configfile-str.ini")
    c_str = autils.Conf_all_parser(cfg_file, ["sec2"])
    assert c_str.sec_dicts["sec2"]["paramint"] == "2"
    c_str.set_int("sec2", "paramint")
    assert c_str["sec2"]["paramint"] == "2"
    assert c_str.sec_dicts["sec2"]["paramint"] == 2

    with pytest.raises(SystemExit):
        c1.set_int("sec1", "toto")
    out, err = capsys.readouterr()
    assert ('ERROR: toto must be an int. Wrong value: parameter.') in out


def test_conf_parser_setfloat(capsys):
    """
    try to convert a given parameter to a float
    """
    c1 = autils.Conf_all_parser(CONFFILE, ["sec_float"])
    # float with decimal
    c1.set_float("sec_float", "float_param")
    assert c1["sec_float"]["float_param"] == "0.12"
    assert c1.sec_dicts["sec_float"]["float_param"] == 0.12

    # exp format
    c1.set_float("sec_float", "float_param2")
    assert c1["sec_float"]["float_param2"] == "1e5"
    assert c1.sec_dicts["sec_float"]["float_param2"] == 100000

    with pytest.raises(SystemExit):
        c1.set_float("sec1", "toto")
    out, err = capsys.readouterr()
    assert ('ERROR: toto must be a float. Wrong value: parameter.') in out

    # With a fload called with "" but cleaned
    cfg_file = os.path.join("test", "data", "utils", "configfile-str.ini")
    c_str = autils.Conf_all_parser(cfg_file, ["sec_float"])
    assert c_str.sec_dicts["sec_float"]["float_param"] == "0.12"
    c_str.set_float("sec_float", "float_param")
    assert c_str["sec_float"]["float_param"] == "0.12"
    assert c_str.sec_dicts["sec_float"]["float_param"] == 0.12

    # With a fload called with "" not cleaned
    cfg_file = os.path.join("test", "data", "utils", "configfile-str.ini")
    c_str = autils.Conf_all_parser(cfg_file, ["sec_float"], clean_str=False)
    assert c_str.sec_dicts["sec_float"]["float_param"] == '"0.12"'
    with pytest.raises(SystemExit):
        c_str.set_float("sec_float", "float_param")
    out, err = capsys.readouterr()
    assert ('ERROR: float_param must be a float. Wrong value: "0.12".') in out
