#!/bin/bash
if (( $# != 4 )); then
	echo "usage: $0 <oriname> <gemname> <nbcont> <logfile>"
        exit 1
fi

oriname=$1
gemname=$2
nbcont=$3
logfile=$4

echo "running check prokka": $oriname $gemname

# check that prokka11Res folder exists
if [ ! -d $oriname-"prokka11Res" ]; then
    echo " - "$gemname $oriname": no directory prokka11Res" >> $logfile
# if exists, check inside folder
else
	# TODO: check if .tbl .faa and .ffn exist
	if [ ! -f $oriname-"prokka11Res"/$gemname.tbl ]; then
    	echo " - "$gemname $oriname": no .tbl file" >> $logfile
    fi
    if [ ! -f $oriname-"prokka11Res"/$gemname.faa ]; then
    	echo " - "$gemname $oriname": no .faa file" >> $logfile
    fi
    if [ ! -f $oriname-"prokka11Res"/$gemname.ffn ]; then
    	echo " - "$gemname $oriname": no .ffn file" >> $logfile
    fi
    # check that the nb of contigs explored by prokka (in tbl) corresponds to the nb of contigs (LSTINFO)
    nbfound=`grep "^>" $oriname-"prokka11Res/"$gemname.tbl | wc -l`
    if (( $nbfound != $nbcont )); then
        echo " - "$gemname $oriname": no matching number of contigs; nbcontig="$nbcont"; in tbl ="$nbfound >> $logfile
    fi
    # check that the number of proteins (CDS) found in tbl is the same as in the faa file
    protfaa=`grep -c "^>" $oriname"-prokka11Res/"$gemname.faa`
    prottab=`grep -c "CDS" $oriname"-prokka11Res/"$gemname.tbl`
    if (( $prottab != $protfaa )); then
        echo " - "$gemname $oriname": no matching number of proteins between tbl and faa; tbl="$prottab"; faa="$protfaa >> $logfile
    fi
    # check that the total number of genes in tbl is the same as in ffn file
    genetab=`grep -c "locus_tag" $oriname"-prokka11Res/"$gemname.tbl`
    crisprtab=`grep -c "repeat_region" $oriname"-prokka11Res/"$gemname.tbl`
    geneffn=`grep -c "^>" $oriname"-prokka11Res/"$gemname.ffn`
    if (( $genetab + $crisprtab != $geneffn )); then
        echo " - "$gemname": no matching number of genes between tbl and ffn; tbl="$genetab $crisprtab"; faa="$geneffn >> $logfile
    fi
fi