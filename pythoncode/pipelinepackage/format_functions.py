#!/usr/bin/env python3
# coding: utf-8

"""
Functions to:
- convert prokka tbl file to our tab file
- convert prokka ffn and faa headers to our format

@author gem
April 2017
"""

def tbl2lst(tblfile, lstfile):
    """
    Read prokka tbl file, and convert it to the lst file.

    * prokka tbl file format:
    ```
    >Feature contig_name
    start end type
            [EC_number text]
            [gene text]
            inference ab initio prediction:Prodigal:2.6
            [inference text]
            locus_tag test
            product text
    ```
    where `type` can be CDS, tRNA, rRNA, etc ...
    lines between [] are not always present

    * lst file format:
    ```
    start end strand type locus gene_name | product | EC_number | inference 2
    ```
    with the same types as prokka file, and strain is C (complement) or D (direct)
    locus is: `<genome_name>.<i or b><contig_num>_<protein_num>`

    """
    crispr_num = 1
    start, end = -1, -1
    strand, gtype, genome, cont_num, locus_num = [""] * 5
    gene_name, product, ecnum, inf2 = ["NA"] * 4
    cont_loc = "b"
    with open(tblfile, "r") as tblf, open(lstfile, "w") as lstf:
        for line in tblf:
            elems = line.strip().split("\t")
            # New feature (= new contig)
            if line.startswith(">Feature"):
                # if there was something in the previous contig, print last gene of
                # previous contig (with contigLoc=b), and reinitiate for new contig/genes
                if start != -1 and end != -1:
                    # if last gene was a crispr
                    if gtype == "CRISPR":
                        locus_num = "CRISPR" + str(crispr_num)
                        gene_name = "crispr"
                        product = "crispr-array"
                        crispr_num += 1
                    cont_loc = "b"
                    locus_name = "{}.{}{}_{}".format(genome, cont_loc,
                                                     str(cont_num).zfill(4),
                                                     str(locus_num).zfill(5))
                    more_info = "| {} | {} | {}".format(product, ecnum, inf2)
                    lst_line = "\t".join([str(start), str(end), strain, gtype,
                                          locus_name, gene_name, more_info])
                    lstf.write(lst_line + "\n")

                    # init for next feature
                    start, end = -1, -1
                    strand, gtype, genome, cont_num, locus_num = [""] * 5
                    gene_name, product, ecnum, inf2 = ["NA"] * 4
                    cont_loc = "b"
                contig = line.strip().split()[-1]
                genome = ".".join(contig.split(".")[:-1])
                cont_num = contig.split(".")[-1]
            # Line indicating position of gene
            if len(elems) == 3:
                # if not first gene of the contig, print it
                if start != -1 and end != -1:
                    if gtype == "repeat_region":
                        gtype = "CRISPR"
                        locus_num = "CRISPR" + str(crispr_num)
                        gene_name = "crispr"
                        product = "crispr-array"
                        crispr_num += 1
                    locus_name = "{}.{}{}_{}".format(genome, cont_loc,
                                                     cont_num.zfill(4),
                                                     locus_num.zfill(5))
                    more_info = "| {} | {} | {}".format(product, ecnum, inf2)
                    lst_line = "\t".join([str(start), str(end), strain, gtype,
                                          locus_name, gene_name, more_info])
                    lstf.write(lst_line + "\n")
                    gtype, strand, locus_num = [""] * 3
                    start, end = -1, -1
                    gene_name, product, ecnum, inf2 = ["NA"] * 4
                    cont_loc = "i"
                start, end, gtype = int(elems[0]), int(elems[1]), elems[2]
                # Get strain of gene
                if end < start:
                    start, end = end, start
                    strain = "C"
                else:
                    strain = "D"
            # line indicating locus_tag
            if "locus_tag" in elems[0] and len(elems) == 2:
                locus_num = elems[-1].split("_")[-1]
            if "gene" in elems[0] and len(elems) == 2:
                gene_name = elems[1]
            if "product" in elems[0] and len(elems) == 2:
                product = elems[1]
            if "EC_number" in elems[0] and len(elems) == 2:
                ecnum = elems[1]
            if ("inference" in elems[0] and "ab initio prediction:Prodigal" not in line
                and len(elems) == 2):
                inf2 = elems[1]
        # Write last gene:
        if start != -1 and end != -1:
            # if last gene was a crispr
            if gtype == "repeat_region":
                gtype = "CRISPR"
                locus_num = "CRISPR" + str(crispr_num)
                gene_name = "crispr"
                product = "crispr-array"
                crispr_num += 1
            cont_loc = "b"
            locus_name = "{}.{}{}_{}".format(genome, cont_loc,
                                             str(cont_num).zfill(4),
                                             str(locus_num).zfill(5))
            more_info = "| {} | {} | {}".format(product, ecnum, inf2)
            lst_line = "\t".join([str(start), str(end), strain, gtype,
                                  locus_name, gene_name, more_info])
            lstf.write(lst_line + "\n")

if __name__ == '__main__':
    import sys
    tbl2lst(sys.argv[1], sys.argv[2])
