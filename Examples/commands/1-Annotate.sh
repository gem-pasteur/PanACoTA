# With prokka
PanACoTA annotate -d Examples/genomes -r Examples/1-res-Annotate-prokka Examples/input_files/list_genomes.lst -n EXAM

# With prodigal
PanACoTA annotate -d Examples/genomes_init -l Examples/input_files/list_genomes.lst -n EXAM -r Examples/1-res-Annotate-prodigal --prodigal --small
