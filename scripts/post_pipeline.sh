#!/bin/bash

if [ $# != 4 ]
then
	echo "usage: $0 <dbpath> <lstfile> <date_init> <respath>"
	exit
fi

# Get parameters
dbpath=$1
lstfile=$2
dateinit=$3
respath=$4

# Get current date and hour
date=`date +"%d-%m-%y.%H-%M"`

# Concatenate all out files into one
cat $dbpath/log-*-$dateinit.out > $respath/log-$lstfile-$date.out
# Concatenate all error files into one
cat $dbpath/err-*-$dateinit.err > $respath/err-$lstfile-$date.err

# Remove temporary files
rm prep-$lstfile*  # cluster logs of pre-annotation step
rm $dbpath/prokka-gembase_$lstfile*  # cluster logs of annotation steps
rm $dbpath/log-*-$dateinit.out  # individual output files
rm $dbpath/err-*-$dateinit.err  # individual error files

