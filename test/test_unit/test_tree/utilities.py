#!/usr/bin/env python3
# coding: utf-8

"""
Functions used by tree unit tests
"""
import numbers

from Bio import Phylo


def is_tree_lengths(treefile):
    """
    Check that given tree is in newick format with branch lengths only (no bootstrap)

    Parameters
    ----------
    treefile: str
        Path to file containing tree to check

    Returns
    -------
    bool
    """
    mytree = Phylo.read(treefile, "newick")
    root = False
    for elem in mytree.find_clades():
        # Check that each clade has a branch length and no bootstrap
        if not isinstance(elem.branch_length, float):
            if not root:
                root = True
            else:
                print("No branch length")
                return False
        if elem.confidence is not None:
            print("Bootstrap value")
            return False
    return True


def is_tree_bootstrap(treefile):
        """
        Check that given tree is in newick format with branch lengths and bootstrap values

        Parameters
        ----------
        treefile: str
            Path to file containing tree to check

        Returns
        -------
        bool
        """
        mytree = Phylo.read(treefile, "newick")
        root = False
        for elem in mytree.find_clades():
            # Check that each clade has a branch length
            if not isinstance(elem.branch_length, numbers.Real):
                if not root:
                    root = True
                    continue
                else:
                    print("No branch length")
                    return False
            # If it is a leaf, check that no bootstrap value. Otherwise, check that bootstrap value
            if elem.name is None:
                if not isinstance(elem.confidence, numbers.Real):
                    print("No bootstrap", elem.confidence)
                    return False
            else:
                if elem.confidence is not None:
                    print("Bootstrap for leaf")
                    return False
        return True
