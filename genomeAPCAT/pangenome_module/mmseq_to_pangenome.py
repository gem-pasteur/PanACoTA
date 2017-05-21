#!/usr/bin/env python
# coding: utf-8

"""
After clustering with mmseq, and converting the result to a tsv file, we have a
file with 2 columns, corresponding to all couples which are in the same cluster. The first column
is the representative sequence of the cluster, and the second column is another sequence in the cluster (can be the representative itself also).

We then want to convert that into a "pangenome-like" file:
1 line per family: numfamily and then all IDs of proteins in the family

@author: gem
November 2016
"""
import os


def main(tsvfile, outfile=None):
    """
    Converts the given tsv file into a pangenome file
    """
    if not outfile:
        outpath = os.path.dirname(tsvfile)
        base = os.path.basename(tsvfile)
        outfile = os.path.join(outpath, "PanGenome-" + base + ".lst")
    clusters = mmseq_tsv_to_clusters(tsvfile)
    clusters_to_file(clusters, outfile)


def mmseq_tsv_to_clusters(mmseq):
    """
    Reads the output of mmseq as a tsv file, and converts it to a python dist
    """
    clusters = {}  # {representative: [all members]}
    with open(mmseq, "r") as mmsf:
        for line in mmsf:
            repres, other = line.strip().split()
            if repres in clusters:
                clusters[repres].append(other)
            else:
                clusters[repres] = [repres]
    return clusters


def clusters_to_file(clust, fileout):
    """
    Write all clusters to a file
    """
    with open(fileout, "w") as fout:
        num = 1
        for _, fam in clust.items():
            fout.write(str(num))
            for mem in sorted(fam, key=lambda x: (int(x.split(".")[2]), int(x.split("_")[1]))):
                fout.write(" " + mem)
            fout.write("\n")
            num += 1

