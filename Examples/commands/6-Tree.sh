# run example after running all other steps (annotate, pan, core, align)
PanACoTA tree -a Examples/5-align/Phylo-GENO3_1/GENO3_1.nucl.grp.aln -o Examples/6-tree

# run only tree step
PanACoTA tree -a Examples/input_files/tree-input/GENO3_1.nucl.grp.aln -o Examples/6-tree-alone