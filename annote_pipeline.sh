#!/bin/bash

show_help() {
cat << EOF
Usage: $0 -l LSTINFO_file -d dbpath -r respath [-m email]

	-l LSTINFO_file: File with 1 genome per line, 2 columns. First column with gembase name of the genome, 2nd column with the current name of the genome (with file extension). The gembase name is: GGSS.mmyy.nnnnn with GGSS the 2 first letters of gender and 2 first letters of species, mmyy the mounth and year, nnnnn the strain number (with trainling 0s if needed).
	-d dbpath: path to the folder containing all genome sequences to annotate.
	-r respath: path to the folder where all results must be generated (folders LSTINFO, Genes, Replicons, Proteins)
	[-m email]: give here your email address to be notified when all is finished.


Output:
	- In your given respath, you will find 4 folders: LSTINFO (information on each genome, with gene annotations), Genes (nuc. gene sequences), Proteins (aa proteins sequences), Replicons (input sequences but with formatted headers).
	- In the database, folders with prokka results will be created for each input genome.
	- In your given respath, a file called log-<LSTINFO_file>-<current_date>.out will be generated. You can find there all logs: problems during annotation (hence no formatting step ran), and problems during formatting step. All steps start with a '*', and problems start with a ' - '.

EOF
}

# parse arguments
while getopts ":l:d:r:mh" opt; do
  case $opt in
    l)
		lstinfo=$OPTARG
      	;;
    d)
		dbpath=$OPTARG
		;;
    r)
		respath=$OPTARG
		;;
	e)
		email=$OPTARG
		;;
  	h)
  		show_help
  		exit 1
  		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
     	exit 1
      	;;
    \?)
     	echo "Invalid option: -$OPTARG" >&2
     	exit 1
      	;;
  	esac
done

# Check that required arguments are given
if [ -z $lstinfo ]; then
	echo "'-l LSTINFO_file' required. Use -h for more information on usage."
	exit 1
fi
if [ -z $dbpath ]; then
	echo "'-d dbpath' required. Use -h for more information on usage."
	exit 1
fi
if [ -z $respath ]; then
	echo "'-r respath' required. Use -h for more information on usage."
	exit 1
fi

lstfile=$(basename $lstinfo)
nbtask=`wc -l < $lstinfo`
date=`date +"%d-%m-%y"`

echo "-> LSTINFO file: " $lstfile
echo "-> Database in: " $dbpath
echo "-> number of genomes: " $nbtask

if [[ "${dbpath:0:1}" != / && "${dbpath:0:2}" != ~[/a-z] ]]; then
	dbpath=$(pwd)/$dbpath
fi

if [[ "${lstinfo:0:1}" != / && "${lstinfo:0:2}" != ~/ ]]; then
	lstinfo=$(pwd)/$lstinfo
fi

scriptdir=`pwd`

source /local/gensoft/adm/etc/profile.d/modules.sh
module load Python/2.7.11

# rename contigs with gembase name for all genomes and get genome info (nb contigs and size)
qsub -q gem -N prep-$lstfile -cwd -S $(which python) prepare_sequences.py -d $dbpath -l $lstinfo

# For each genome:
#
#	- annotate with prokka
#	- check prokka run
#	- if prokka run ok, translate output to gembase format
#	- check gembase format generated
qsub -t 1-$nbtask -q gem -hold_jid prep-$lstfile -N prokka-gembase_$lstfile-$date -wd $dbpath -pe thread 2 prokka_array.sh $lstinfo $scriptdir $respath

if [ -z $email ]; then
	qsub -q gem -N post-$lstfile -cwd -hold_jid prokka-gembase_$lstfile-$date post_pipeline.sh $dbpath $lstfile
else
	qsub -q gem -N post-$lstfile -M $email -m e -cwd -hold_jid prokka-gembase_$lstfile-$date post_pipeline.sh $dbpath $lstfile
fi