# after running the example command-line for annotate
PanACoTA pangenome -l Examples/2-res-prokka/LSTINFO-list_genomes.lst -n GENO3 -d Examples/2-res-prokka/Proteins -i 0.8 -o Examples/3-pangenome

# to run only pangenome step of the example
PanACoTA pangenome -l Examples/input_files/pan-input/LSTINFO-list_genomes.lst -n GENO3 -d Examples/input_files/pan-input/Proteins -i 0.8 -o Examples/3-pangenome-alone
