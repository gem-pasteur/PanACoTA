#!/bin/bash

if [ $# != 2 ]
then
	echo "usage: $0 <dbpath> <lstfile>"
	exit
fi

dbpath=$1
lstfile=$2

date=`date +"%d-%m-%y"`

cat $dbpath/log-* > log-$lstfile-$date.out
rm prep-$lstfile*
# rm $dbpath/prokka-gembase_$lstfile*
rm $dbpath/log-*.out
