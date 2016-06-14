#! /bin/bash

if [ $# != 4 ]
then
	echo "usage: $0 <oriname> <gemname> <respath> <scriptdir>

    - genomes_LSTINFO: LSTINFO file containing list of genomes with original name, gembase name, etc.
    - db_path: path to the database containing all genome sequences and prokka result folders
    - res_path: path to put the gembase formatted database (will contain the folders Replicons, Genes, Proteins etc.
    - scriptdir: folder where all scripts are."


	exit
fi

oriname=$1
gemname=$2
respath=$3   # will contain folders LSTINFO, Genes, Proteins and Replicons
scriptdir=$4

mkdir -p $respath/LSTINFO
mkdir -p $respath/Replicons
mkdir -p $respath/Genes
mkdir -p $respath/Proteins

# Generate lstinfo file
awk -f $scriptdir/tblToLst.awk $oriname-prokka11Res/$gemname.tbl > $respath/LSTINFO/$gemname.lst
# copy input multifasta file to <gembaseName>.fna in Replicons
cp $oriname"-gembase.fna" $respath/Replicons/$gemname.fna
# generate .gen (in Genes) and .prt (in Proteins)
$(which python) $scriptdir/create_prt_gen.py -p $oriname"-prokka11Res/"$gemname -l $respath
