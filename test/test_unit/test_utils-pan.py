#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the utils_pangenome submodule of genomeAPCAT
"""


from genomeAPCAT import utils_pangenome as upan


def test_read_gene():
    """
    Check that when reading a given gene name, it extracts expected information
    """
    gene = "ESCO.1016.00012.i012_00015"
    num = "1"
    fams_by_strain = {"1": {}}
    all_strains = []
    upan.read_gene(gene, num, fams_by_strain, all_strains)
    assert all_strains == ["ESCO.1016.00012"]
    assert fams_by_strain == {"1": {"ESCO.1016.00012": [gene]}}


def test_read_gene_other_format():
    """
    Check that when reading a gene name which does not have the gembase
    format, it still returns the right information
    """
    gene = "my_gene-name_other_0001"
    num = "1"
    fams_by_strain = {"1": {}}
    all_strains = []
    upan.read_gene(gene, num, fams_by_strain, all_strains)
    assert all_strains == ["my_gene-name_other"]
    assert fams_by_strain == {"1": {"my_gene-name_other": [gene]}}


def test_read_gene_strain_known():
    """
    Check that when reading a gene name, and its corresponding
    strain already exists in fams_by_strain and all_strains, it just adds
    this new gene to fams_by_strain, and does nothing to all_strains
    """
    gene = "ESCO.1016.00012.i012_00015"
    num = "1"
    fams_by_strain = {"1": {"ESCO.1016.00012": ["ESCO.1016.00012.i001_01"]}}
    all_strains = ["ESCO.1016.00012"]
    upan.read_gene(gene, num, fams_by_strain, all_strains)
    assert all_strains == ["ESCO.1016.00012"]
    assert fams_by_strain == {"1": {"ESCO.1016.00012": ["ESCO.1016.00012.i001_01",
                                                        gene]}}

