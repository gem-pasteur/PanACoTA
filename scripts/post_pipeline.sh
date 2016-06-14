#!/bin/bash

#$ -o $4/post-treatment.out
#$ -e $4/post-treatment.err

if [ $# != 4 ]
then
	echo "usage: $0 <dbpath> <lstfile> <date_init> <respath>"
	exit
fi

dbpath=$1
lstfile=$2
dateinit=$3
respath=$4

date=`date +"%d-%m-%y.%H-%M"`

cat $dbpath/log-*-$dateinit.out > $respath/log-$lstfile-$date.out
rm prep-$lstfile*
# rm $dbpath/prokka-gembase_$lstfile*
rm $dbpath/log-*-$dateinit.out
