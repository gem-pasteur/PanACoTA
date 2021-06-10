# core genome
PanACoTA corepers -p Examples/3-pangenome/PanGenome-GENO3.All.prt-clust-0.8-mode1.lst -o Examples/4-corepers

# strict persistent at 95%
PanACoTA corepers -p Examples/3-pangenome/PanGenome-GENO3.All.prt-clust-0.8-mode1.lst -o Examples/4-corepers0.95 -t 0.95

# run only coregenome step from example
PanACoTA corepers -p Examples/input_files/core-input/PanGenome-example.lst -o Examples/4-corepers-alone