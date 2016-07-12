#!/bin/bash

lstinfo=$1
scriptdir=$2
respath=$3
dateinit=$4

module load hmmer/3.1b1 aragorn/1.2.36 barrnap minced/0.1.6 blast+ prodigal infernal/1.1 ncbi_toolbox/20151127 signalp
module load prokka/1.11

vals=$(cat $lstinfo | head -n $(( ${SGE_TASK_ID} + 1 )) | tail -n 1)
vars=($vals)
gemname=${vars[0]}
oriname=${vars[1]}
nbcont=${vars[2]}
logfile=log-$gemname-$dateinit.out

if [ ! $gemname == "gembase_name" ]; then
	echo "--------------------------" > $logfile
	echo "*" `date +"%d/%m/%y - %T"` "prokka: task" ${SGE_TASK_ID}: $gemname $oriname $nbcont "contigs" >> $logfile
	prokka --outdir $oriname-prokka11Res --cpus 2 --prefix $gemname $oriname"-gembase.fna"
	# --centre to generate clean contig names (length <= 20)
	echo "*" `date +"%d/%m/%y - %T"` "checking prokka results" >> $logfile
	$scriptdir/check_prokka_run.sh $oriname $gemname $nbcont $logfile
	logsize=`wc -l < $logfile`
	# check that prokka11Res folder exists

	if (( $logsize == 3 )); then
		echo "*" `date +"%d/%m/%y - %T"` "Generate gembase files" >> $logfile
		$scriptdir/generate_gembase.sh $oriname $gemname $respath $scriptdir
		echo "*" `date +"%d/%m/%y - %T"` "Check gembase results" >> $logfile
		$scriptdir/check_gembase.sh $oriname $gemname $respath $logfile
	fi
else
	echo "prokka: task" ${SGE_TASK_ID}: $gemname $oriname $nbcont
fi
