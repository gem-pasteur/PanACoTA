# Only QC
PanACoTA annotate -d Examples/genomes_init -l Examples/input_files/list_genomes.lst -r Examples/2-res-QC -Q

# With prokka
PanACoTA annotate -d Examples/genomes_init -r Examples/2-res-prokka -l Examples/input_files/list_genomes.lst -n GENO --l90 3 --nbcont 10


# With prodigal
PanACoTA annotate -d Examples/genomes_init -r Examples/2-res-prodigal -l Examples/input_files/list_genomes.lst -n GENO  --l90 3 --prodigal --small 
