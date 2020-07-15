#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the persistent_functions submodule in corepers module
"""
import os
import logging

import PanACoTA.corepers_module.persistent_functions as persf
import test.test_unit.utilities_for_tests as tutils


PERS_PATH = os.path.join("test", "data", "persgenome")
EXP_PATH = os.path.join(PERS_PATH, "exp_files")
FAMS_BY_STRAIN = \
    {'1': {'GEN4.1111.00001': ['GEN4.1111.00001.b0001_00001'],
           'GENO.0817.00001': ['GENO.0817.00001.b0001_00002'],
           'GENO.1216.00002': ['GENO.1216.00002.b0001_00001', 'GENO.1216.00002.i0001_00002'],
           'GENO.1216.00003': ['GENO.1216.00003.i0001_00003']
           },
     '10': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00004'],
            'GENO.0817.00001': ['GENO.0817.00001.i0002_00004'],
            'GENO.1216.00002': ['GENO.1216.00002.i0001_00005']
            },
     '11': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00005'],
            'GENO.0817.00001': ['GENO.0817.00001.i0002_00005'],
            'GENO.1216.00002': ['GENO.1216.00002.i0001_00006']
            },
     '12': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00008'],
            'GENO.0817.00001': ['GENO.0817.00001.i0002_00010'],
            'GENO.1216.00002': ['GENO.1216.00002.b0002_00009'],
            'GENO.1216.00003': ['GENO.1216.00003.i0001_00004', 'GENO.1216.00003.i0001_01000']
            },
     '13': {'GENO.1216.00002': ['GENO.1216.00002.b0003_00011'],
            'GENO.0817.00001': ['GENO.0817.00001.i0002_00010']},
     '14': {'GENO.1216.00002': ['GENO.1216.00002.b0003_00012']},
     '2': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00003']},
     '3': {'GEN4.1111.00001': ['GEN4.1111.00001.b0001_00009'],
           'GENO.0817.00001': ['GENO.0817.00001.b0002_00011'],
           'GENO.1216.00002': ['GENO.1216.00002.b0002_00010'],
           'GENO.1216.00003': ['GENO.1216.00003.i0001_01010']
           },
     '4': {'GENO.0817.00001': ['GENO.0817.00001.b0001_00001']},
     '5': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00002'],
           'GENO.0817.00001': ['GENO.0817.00001.b0002_00003'],
           'GENO.1216.00002': ['GENO.1216.00002.i0001_00003'],
           'GENO.1216.00003': ['GENO.1216.00003.i0080_00010']
           },
     '6': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00006'],
           'GENO.0817.00001': ['GENO.0817.00001.i0002_00006', 'GENO.0817.00001.i0002_00007'],
           'GENO.1216.00002': ['GENO.1216.00002.i0001_00007']
           },
     '7': {'GENO.0817.00001': ['GENO.0817.00001.i0002_00008']},
     '8': {'GEN4.1111.00001': ['GEN4.1111.00001.i0001_00007'],
           'GENO.0817.00001': ['GENO.0817.00001.i0002_00009'],
           'GENO.1216.00002': ['GENO.1216.00002.b0001_00008']
           },
     '9': {'GENO.1216.00002': ['GENO.1216.00002.i0001_00004']}
     }

FAMILIES = {'1': ['GEN4.1111.00001.b0001_00001', 'GENO.0817.00001.b0001_00002',
                  'GENO.1216.00002.b0001_00001', 'GENO.1216.00002.i0001_00002',
                  'GENO.1216.00003.i0001_00003'],
            '10': ['GEN4.1111.00001.i0001_00004', 'GENO.0817.00001.i0002_00004',
                   'GENO.1216.00002.i0001_00005'],
            '11': ['GEN4.1111.00001.i0001_00005', 'GENO.0817.00001.i0002_00005',
                   'GENO.1216.00002.i0001_00006'],
            '12': ['GEN4.1111.00001.i0001_00008', 'GENO.0817.00001.i0002_00010',
                   'GENO.1216.00002.b0002_00009', 'GENO.1216.00003.i0001_00004',
                   'GENO.1216.00003.i0001_01000'],
            '13': ['GENO.1216.00002.b0003_00011', 'GENO.0817.00001.i0002_00010'],
            '14': ['GENO.1216.00002.b0003_00012'],
            '2': ['GEN4.1111.00001.i0001_00003'],
            '3': ['GEN4.1111.00001.b0001_00009', 'GENO.0817.00001.b0002_00011',
                  'GENO.1216.00002.b0002_00010', 'GENO.1216.00003.i0001_01010'],
            '4': ['GENO.0817.00001.b0001_00001'],
            '5': ['GEN4.1111.00001.i0001_00002', 'GENO.0817.00001.b0002_00003',
                  'GENO.1216.00002.i0001_00003', 'GENO.1216.00003.i0080_00010'],
            '6': ['GEN4.1111.00001.i0001_00006', 'GENO.0817.00001.i0002_00006',
                  'GENO.0817.00001.i0002_00007', 'GENO.1216.00002.i0001_00007'],
            '7': ['GENO.0817.00001.i0002_00008'],
            '8': ['GEN4.1111.00001.i0001_00007', 'GENO.0817.00001.i0002_00009',
                  'GENO.1216.00002.b0001_00008'],
            '9': ['GENO.1216.00002.i0001_00004']
            }


def test_uniq_mems():
    """
    Test that it returns True when there is only 1 member of each genome in the given family
    """
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain3": ["member3"]}
    assert persf.uniq_members(family)


def test_nonuniq_mems():
    """
    Test that it returns False when there is more than 1 member in a genome
    """
    family = {"strain1": ["member", "member-bis"],
              "strain2": ["member2"],
              "strain3": ["member3"]}
    assert not persf.uniq_members(family)


def test_zero_mem(caplog):
    """
    Test that when there are no member in 1 family, and 1 in all others, it
    returns True but with a warning message.
    """
    caplog.set_level(logging.DEBUG)
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain3": ["member3"],
              "strain4": []}
    assert persf.uniq_members(family)
    assert "Problem, no members for strain4!" in caplog.text


def test_zero_mem_notuniq(caplog):
    """
    Test that when there are no member in 1 family, and several in 1 other family, it
    returns False + a warning message.
    """
    caplog.set_level(logging.DEBUG)
    family = {"strain1": [],
              "strain2": ["member2"],
              "strain3": ["member3", "hello"],
              "strain4": ["my_member"]}
    assert not persf.uniq_members(family)
    assert "Problem, no members for strain1!" in caplog.text

def test_mixed():
    """
    Test that it returns true when there is exactly 1 member in 3 genomes/4 (and minimum
    asked is 3), and several members in the other strain
    """
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain3": ["member3"],
              "strain4": ["member4", "member4bis"]}
    assert persf.mixed_family(family, 3)


def test_mixed_empty(caplog):
    """
    Test that when there is exactly 1 member in 3 genomes / 4 (and min 3 asked),
    and 0 in the other genomes, it returns True.
    """
    caplog.set_level(logging.DEBUG)
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain2b": [],
              "strain3b": [],
              "strain4": [],
              "strain3": ["member3"],
              "strain5": []}
    assert persf.mixed_family(family, 3)
    assert "Problem, no members for strain2b" in caplog.text
    assert "Problem, no members for strain3b" in caplog.text
    assert "Problem, no members for strain4" in caplog.text


def test_not_mixed():
    """
    Test that it returns false when there are less than the required
    number of genomes with exactly 1 member
    """
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain3": ["member3"],
              "strain4": ["member4", "member4bis"]}
    assert not persf.mixed_family(family, 4)


def test_mixed_more():
    """
    Test that it returns true when there are more than the required
    number of genomes with exactly 1 member
    """
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain3": ["member3"],
              "strain4": ["member4"],
              "strain5": ["member5"]}
    assert persf.mixed_family(family, 3)


def test_write_pers():
    """
    Test that output file is written as expected
    """
    fams = {9: ["member_3", "member_12", "other_member_2"],
            3: ["member_10", "member_100", "member_1"],
            10: ["member_1", "member_2", "member_3"],
            1: ["my_protein_3", "my_protein_12", "my_protein_2"],
            5: ["ESCO.1216.00003.i001_01001", "SAEN.0215.00003.i009_00001",
                "ESCO.1017.00003.b001_00001", "ESCO.0812.00002.i002_02000",
                "ESCO.0812.00003.i002_02000"]}
    outfile = "test-persistent_families.txt"
    persf.write_persistent(fams, outfile)
    expfile = os.path.join(EXP_PATH, "exp_persgenome1.txt")
    assert tutils.compare_order_content(outfile, expfile)
    os.remove(outfile)


def test_get_core_strict(caplog):
    """
    Getting a core genome (4 genomes, all having exactly 1 member)
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['3', '5']}
    assert exp_fams == fams
    assert ("The core genome contains 2 families, each one having "
            "exactly 4 members, from the 4 different genomes.") in caplog.text


def test_get_core_multi(caplog):
    """
    Getting a multi core genome (4 genomes, having at least 1 member)
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4, multi=True)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['3', '5', '1', '12']}
    assert exp_fams == fams
    assert ("The persistent genome contains 4 families with members present in "
            "at least 4 different genomes (100% of the total number of genomes)") in caplog.text


def test_get_99pers_floor_strict(caplog):
    """
    Getting a strict persistent at floor(99%) -> at least 3 genomes with 1 member, others
    absent.
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4, tol=0.99, floor=True)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['3', '5', '8', '10', '11']}
    assert exp_fams == fams
    assert ("The persistent genome contains 5 families, each one "
            "having exactly 1 member from at least 99.0% of the 4 different genomes "
            "(that is 3 genomes). The other genomes are absent from the family.") in caplog.text


def test_get_99pers_floor_mixed(caplog):
    """
    Getting a mixed persistent at floor(99%) -> at least 3 genomes with 1 member, others
    anything
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4, tol=0.99, floor=True, mixed=True)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['1', '3', '5', '8', '10',
                                                                       '11', '12']}
    assert exp_fams == fams
    assert ("The persistent genome contains 7 families, each one having exactly 1 member from at least "
            "99.0% of the genomes (3 genomes). In the remaining "
            "1.0% genomes, there can be 0, 1 or several members.") in caplog.text


def test_get_99pers_floor_multi(caplog):
    """
    Getting a multi persistent genome at floor(99%) -> at least 3 genomes (any number of members)
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4, tol=0.99, floor=True, multi=True)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['1', '3', '5', '6', '8',
                                                                       '10', '11', '12']}
    assert exp_fams == fams
    assert ("The persistent genome contains 8 families with members present in "
            "at least 3 different genomes (99.0% of the total number of genomes).") in caplog.text


def test_get_99pers_strict(caplog):
    """
    Getting a persistent genome at 99% (ceil) -> 4 genomes with exactly 1member
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4, tol=0.99)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['3', '5']}
    assert exp_fams == fams
    assert ("The persistent genome contains 2 families, each one "
            "having exactly 1 member from at least 99.0% of the 4 different genomes "
            "(that is 4 genomes). The other genomes are absent from the family.") in caplog.text


def test_get_99pers_mixed(caplog):
    """
    Getting a mixed persistent genome at 99% (ceil) -> 4 genomes with exactly 1member
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4, tol=0.99, mixed=True)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['3', '5']}
    assert exp_fams == fams
    assert ("The persistent genome contains 2 families, each one having exactly 1 member from at least "
            "99.0% of the genomes (4 genomes). In the remaining "
            "1.0% genomes, there can be 0, 1 or several members.") in caplog.text


def test_get_99pers_multi(caplog):
    """
    Getting a multi persistent genome at 99% (ceil) -> 3 genomes with exactly 1member,
    other with anything
    """
    caplog.set_level(logging.DEBUG)
    fams = persf.get_pers(FAMS_BY_STRAIN, FAMILIES, 4, tol=0.99, multi=True)
    exp_fams = {num: mems for num, mems in FAMILIES.items() if num in ['1', '3', '5', '12']}
    assert exp_fams == fams
    assert ("The persistent genome contains 4 families with members present in "
            "at least 4 different genomes (99.0% of the total number of genomes).") in caplog.text
