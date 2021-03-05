# run example (after running annotate, pan and core steps)
PanACoTA align -c Examples/4-corepers/PersGenome_<pangenome-filename>_1.lst -l Examples/2-res-prokka/LSTINFO-list_genomes.lst -n GENO3_1 -d Examples/2-res-prokka -o Examples/5-align

# only run align step of example
PanACoTA align -c Examples/input_files/align-input/coregenome-example.lst -l Examples/input_files/pan-input/LSTINFO-list_genomes.lst -n GENO3_1 -d Examples/input_files/pan-input -o Examples/5-align