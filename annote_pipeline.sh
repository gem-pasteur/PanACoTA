#!/bin/bash



show_help() {
cat << EOF
Usage: $0 -l LSTINFO_file -d dbpath -r respath [-m email -f -c cpus]

	-l LSTINFO_file: File with 1 genome per line, 2 columns, no header. First column with gembase name of the genome, 2nd column with the current name of the genome (with file extension). The gembase name is: GGSS.mmyy.nnnnn with GGSS the 2 first letters of gender and 2 first letters of species, mmyy the month and year, nnnnn the strain number (with trainling 0s if needed).
	-d dbpath: path to the folder containing all genome sequences to annotate.
	-r respath: path to the folder where all results must be generated (folders LSTINFO, Genes, Replicons, Proteins)
	[-m email]: give here your email address to be notified when all is finished.
	[-f]: add this option if you want to run prokka even if the result folder already exists. Otherwise, if the prokka folder exists, the pipeline will run the formatting step with the already generated results. Note that this will be applied to all genomes having a result folder. If you want to rerun prokka only on a specific genome, remove its result folder before running ths script without the '-f' option.
	[-c cpus]: number of cpus to use to run prokka. By default, 2 cpus are used.


Output:
	- In your given respath, you will find 4 folders: LSTINFO (information on each genome, with gene annotations), Genes (nuc. gene sequences), Proteins (aa proteins sequences), Replicons (input sequences but with formatted headers).
	- In the database, folders with prokka results will be created for each input genome. If errors are generated during prokka step, you can look at the log file to see what was wrong.
	- In your given respath, a file called log-<LSTINFO_file>-<current_date>.out will be generated. You can find there all logs: problems during annotation (hence no formatting step ran), and problems during formatting step. All steps start with a '*', and problems start with a ' - '.
	- In your given respath, a file called err-<LSTINFO_file>-<current_date>.err will be generated, containing information on errors occured. If this file is empty, then annotation and formatting steps finished without any problem for all genomes.


EOF
}

# parse arguments
while getopts "l:d:r:m:fc:h" opt; do
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
	m)
		email=$OPTARG
		;;
	f)
		force="--force"
		;;
	c)
		cpus=$OPTARG
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
if [ -z $cpus ]; then
	cpus=2
fi

# Get given parameters
lstfile=$(basename $lstinfo)
nbtask=`wc -l < $lstinfo`
dateinit=`date +"%d-%m-%y.%H-%M"`

echo "-> LSTINFO file: " $lstfile
echo "-> Database in: " $dbpath
echo "-> number of genomes: " $nbtask

# Put absolute path to database
if [[ "${dbpath:0:1}" != / && "${dbpath:0:2}" != ~[/a-z] ]]; then
	dbpath=$(pwd)/$dbpath
fi

# Put absolute path to LSTINFO file
if [[ "${lstinfo:0:1}" != / && "${lstinfo:0:2}" != ~/ ]]; then
	lstinfo=$(pwd)/$lstinfo
fi

# Put absolute path to respath
if [[ "${respath:0:1}" != / && "${respath:0:2}" != ~[/a-z] ]]; then
        respath=$(pwd)/$respath
fi

# Add execution rights to all scripts
scriptdir=`pwd`/scripts
chmod u+x $scriptdir/*

# Load python
source /local/gensoft/adm/etc/profile.d/modules.sh
module load Python/2.7.11

# rename contigs with gembase name for all genomes and get genome info (nb contigs and size)
qsub -q gem -N prep-$lstfile -cwd -S $(which python) $scriptdir/prepare_sequences.py -d $dbpath -l $lstinfo

# For each genome:
#
#	- annotate with prokka
#	- check prokka run
#	- if prokka run ok, translate output to gembase format
#	- check gembase format generated
qsub -t 1-$nbtask -q gem -hold_jid prep-$lstfile -N prokka-gembase_$lstfile-$dateinit -wd $dbpath -pe thread $cpus $scriptdir/prokka_array.sh $lstinfo-complete.lst $scriptdir $respath $dateinit $cpus $force

post_options="-q gem -N post-$lstfile-$dateinit -cwd -hold_jid prokka-gembase_$lstfile-$dateinit -o $respath/post-treatment.out -e $respath/post-treatment.err"
if [ -z $email ]; then
	qsub $post_options $scriptdir/post_pipeline.sh $dbpath $lstfile $dateinit $respath
else
	qsub $post_options -M $email -m e $scriptdir/post_pipeline.sh $dbpath $lstfile $dateinit $respath
fi
