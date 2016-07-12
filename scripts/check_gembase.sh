#!/bin/bash

if (( $# != 5 )); then
	echo "usage: $0 <oriname> <gemname> <respath> <logfile> <errfile>
    "

    exit 1
fi

oriname=$1
gemname=$2
respath=$3
logfile=$4
errfile=$5

echo $gemname $oriname
if [ ! -f $respath"/LSTINFO/"$gemname.lst ]; then
    echo " - "$gemname": no .lst file in LSTINFO" >> $logfile
    echo " - "$gemname": check gembase: no .lst file in LSTINFO" >> $errfile
fi
if [ ! -f $respath/Replicons/$gemname.fna ]; then
    echo " - "$gemname": no .fna file in Replicons" >> $logfile
    echo " - "$gemname": check gembase: no .fna file in Replicons" >> $errfile
fi
if [ ! -f $respath"/Genes/"$gemname.gen ]; then
    echo " - "$gemname": no .gen file in Genes" >> $logfile
    echo " - "$gemname": check gembase: no .gen file in Genes" >> $errfile
fi
if [ ! -f $respath"/Proteins/"$gemname.prt ]; then
    echo " - "$gemname": no .prt file in Proteins" >> $logfile
    echo " - "$gemname": check gembase: no .prt file in Proteins" >> $errfile
fi
nbcrispr=`grep -c "repeat_region" $oriname-prokka11Res/$gemname.tbl`
genetab=`grep -c "locus_tag" $oriname-prokka11Res/$gemname.tbl`
genelst=`cat $respath"/LSTINFO/"$gemname.lst | wc -l`
# check number of lines in .lst = number of genes/crisprs in .tbl
if (( $nbcrispr + $genetab != $genelst )); then
	echo " - "$gemname ": error tbl to lst in total number" >> $logfile
    echo " - "$gemname ": check gembase: error tbl to lst in total number" >> $errfile
fi
# check number of crispr in .lst = number of crispr in .tbl
crisprlst=`grep -c "crispr-array" $respath"/LSTINFO/"$gemname.lst`
if (( $nbcrispr != crisprlst )); then
	echo " - "$gemname ": error in crispr number" >> $logfile
    echo " - "$gemname ": check gembase: error in crispr number" >> $errfile
fi
# check Replicons/.fna == _multifasta.fst
diffrep=`diff $oriname"-gembase.fna" $respath/Replicons/$gemname.fna`
if [ ! -z $diffrep ]; then
    echo " - "$gemname": error in replicons." >> $logfile
    echo " - "$gemname": check gembase: error in replicons." >> $errfile
fi
# check Proteins/.prt number of proteins = number of CDS in .tbl
cdstbl=`grep -c "CDS" $oriname-prokka11Res/$gemname.tbl`
cdsprt=`grep -c "^>" $respath"/Proteins/"$gemname.prt`
if (( $cdsprt != $cdstbl )); then
    echo " - "$gemname": error in number of proteins in .prt prt:"$cdsprt" tbl:"$cdstbl >> $logfile
    echo " - "$gemname": check gembase: error in number of proteins in .prt prt:"$cdsprt" tbl:"$cdstbl >> $errfile
fi
# check Genes/.gen number of features = number of features in tbl (genes+crispr)
gengen=`grep -c "^>" $respath"/Genes/"$gemname.gen`
if (( $gengen != $genetab + $nbcrispr )); then
    echo " - "$gemname": error in total number of genes in .gen gen:"$gengen" tbl:"$nbcrispr+$genetab >> $logfile
    echo " - "$gemname": check gembase: error in total number of genes in .gen gen:"$gengen" tbl:"$nbcrispr+$genetab >> $errfile
fi

errors=`wc -l < $errfile`
if (( $errors > 0 )); then
    echo "-----------------------" >> $errfile
fi
