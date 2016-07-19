#!/bin/bash

lstinfo=$1
scriptdir=$2
respath=$3
dateinit=$4
cpus=$5
force=$6


module load hmmer/3.1b1 aragorn/1.2.36 barrnap minced/0.1.6 blast+ prodigal infernal/1.1 ncbi_toolbox/20151127 signalp
module load prokka/1.11

vals=$(cat $lstinfo | head -n $(( ${SGE_TASK_ID} + 1 )) | tail -n 1)
vars=($vals)
gemname=${vars[0]}
oriname=${vars[1]}
nbcont=${vars[2]}
logfile=log-$gemname-$dateinit.out
errfile=err-$gemname-$dateinit.err

# Check that the oriname file exists. Otherwise, put a message in errfile and exit.
if [ ! -f $oriname-gembase.fna ]; then
	echo -e $oriname-gembase.fna "(from original file '$oriname') not found. Annotation and formatting pipeline cannot run.\n---------------" | tee -a $errfile -a $logfile
	exit 1
fi

if [ ! $gemname == "gembase_name" ]; then
	echo "*" `date +"%d/%m/%y - %T"` "prokka: task" ${SGE_TASK_ID}: $gemname $oriname $nbcont "contigs" >> $logfile
	# If prokka res folder already exists, write WARNING because it will use already generated results
#	if [[ -d $oriname-prokka11Res ]] && [[ -z $force ]]; then
#		warning="/!\\ `date +"%d/%m/%y - %T"` $gemname $oriname  WARNING: prokka results folder already exists. Prokka did not run again, formatting step used already generated results of Prokka in $oriname-prokka11Res. If you want to re-run prokka, first remove this result folder, or use '-f' option if you want to rerun prokka for all genomes."
#	fi
	prokka --outdir $oriname-prokka11Res --cpus $cpus $force --prefix $gemname $oriname"-gembase.fna"
	# --centre to generate clean contig names (length <= 20)
	echo "*" `date +"%d/%m/%y - %T"` "checking prokka results" >> $logfile
	$scriptdir/check_prokka_run.sh $oriname $gemname $nbcont $logfile $errfile

	# If annotation step ran well, continue with formatting
	logsize=`wc -l < $logfile`
	if (( $logsize == 2 )); then
		echo "*" `date +"%d/%m/%y - %T"` "Generate gembase files" >> $logfile
		$scriptdir/generate_gembase.sh $oriname $gemname $respath $scriptdir
		echo "*" `date +"%d/%m/%y - %T"` "Check gembase results" >> $logfile
		$scriptdir/check_gembase.sh $oriname $gemname $respath $logfile $errfile
	fi
	if [[ ! -z $warning ]]; then
		echo $warning | tee -a $logfile -a $errfile
	fi
	echo "--------------------------" >> $logfile
else
	echo "prokka: task" ${SGE_TASK_ID}: $gemname $oriname $nbcont
fi

errors=`wc -l < $errfile`
if (( $errors > 0 )); then
    echo "-----------------------" >> $errfile
fi
