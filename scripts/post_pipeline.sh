#!/bin/bash

if [ $# != 3 ]
then
	echo "usage: $0 <dbpath> <lstfile> <date_init>"
	exit
fi

dbpath=$1
lstfile=$2
dateinit=$3

date=`date +"%d-%m-%y.%H-%M"`

cat $dbpath/log-*-$dateinit.out > log-$lstfile-$date.out
rm prep-$lstfile*
# rm $dbpath/prokka-gembase_$lstfile*
rm $dbpath/log-*-$dateinit.out
