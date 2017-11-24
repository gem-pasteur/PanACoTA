#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the persistent_functions submodule in corepers module
"""

import genomeAPCAT.corepers_module.persistent_functions as persf

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
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain3": ["member3"],
              "strain4": []}
    assert persf.uniq_members(family)
    assert "problem, no members for strain4!" in caplog.text


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
    and 0 in the other genome, it returns True, and a warning message
    """
    family = {"strain1": ["member"],
              "strain2": ["member2"],
              "strain2b": [],
              "strain3b": [],
              "strain4": [],
              "strain3": ["member3"],
              "strain4": []}
    assert persf.mixed_family(family, 3)
    assert "problem, no members for strain" in caplog.text
