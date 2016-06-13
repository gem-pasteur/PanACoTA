#!/usr/bin/awk

BEGIN{
	FS="\t"  # tab separated input file
	OFS=FS   # tab separated output file
	start=0  # start position on contig
	end=0  # end position on contig
	strand=""  # strand on which is the gene: D or C
	type=""  # tRNA, CDS, rRNA etc.
	prefix=""  # SAEN.1015.<strainID>
	contigID=0  # on 3 digits
	contigLoc="b"  # i (inner: gene inside contig) or b (border: first or last gene of contig). First will be border
	genID=0  # on 5 digits, given by prokka
	genName="NA"  # gene name, given by prokka, NA if not given
	product="NA"  # product description given by Prokka, NA if not given
	ECnum="NA"  # EC_number given by prokka, NA if not given
	inf2="NA"  # 2nd line of "inference" given by Prokka, NA if not given
	crisprNm=1  # crispr ID to give if there are crispr found
	lastFeatureNR=0  # last line number with ">Feature". If a new ">Feature" is the line just after the last one, the last contig did not contain any gene.
}

# New contig
/>Feature/{
	# if not first contig, print last gene of previous contig (with contigLoc=b), and reinitiate for new contig/genes
	# and if last feature contained at least one gene (current line != lastFeatureNR+1)
	if(prefix != ""){
		# print last
		if (type == "CRISPR"){
			genID="CRISPR"crisprNm
			crisprNm=crisprNm+1
			genName="crispr"
			product="crispr-array"
		}
		if ( NR != lastFeatureNR + 1){
			contigLoc = "b"
			print start, end, strand, type, prefix contigLoc contigID "_" genID, genName, "| " product " | " ECnum " | " inf2
		}
		# init for next
		start=0
		end=0
		strand=""
		type=""
		prefix=""
		contigID=0
		genID=0
		genName="NA"
		product="NA"
		ECnum="NA"
		inf2="NA"
	}
	lastFeatureNR=NR
	split($0, a, " ")  # separate ">Feature" and <contig_name>
	split(a[2], tab, ".")  # split each parts of <contig_name> : SAEN, 1015, <strainID>, c<contigID>
	# prefix = SAEN.1015.<strainID>.
	for(i=1;i<=3;i++){
        prefix = prefix tab[i] "."
    }
    #  get contig ID from c<contigID>
    contigID=tab[4]
}

# New line indicating new gene (start, end, type, direction)
NF == 3{
	# if not first gene of contig, print last one, and reinitiate for next genes
	if(start != 0 && end != 0){
		if (type == "CRISPR"){
			genID="CRISPR"crisprNm
			crisprNm=crisprNm+1
			genName="crispr"
			product="crispr-array"
		}
		print start, end, strand, type, prefix contigLoc contigID "_" genID, genName, "| " product " | " ECnum " | " inf2
		# initialize for new gene
		start=0
		end=0
		strand=""
		type=""
		genID=0
		genName="NA"
		product="NA"
		ECnum="NA"
		inf2="NA"
		contigLoc="i"
	}

	# collect info for new gene
	if ($1 < $2){
		start=$1
		end=$2
		strand="D"
	}else{
		start=$2
		end=$1
		strand="C"
	}
	if($3 == "repeat_region"){
		type="CRISPR"
	}else{
		type=$3
	}
	# print start, end, strand, type, prefix contigLoc contigID, "end!!" 
}

# EC_number
NF == 5 && $4=="EC_number"{
	ECnum=$5
	gsub("[|]", " ", ECnum)
}

# gene name
NF == 5 && $4=="gene"{
	genName=$5
}

# inference info line 2
/inference/ && ! /ab initio prediction:Prodigal:2.60/{
	inf2=$5
	gsub("[|]", " ", inf2)
}

# gene ID
NF == 5 && $4=="locus_tag"{
	split($5, a, "_")
	genID=a[2]
}

# product
NF == 5 && $4=="product"{
	product=$5
	gsub("[|]", " ", product)
}

# print last gene
END{
	contigLoc = "b"
	if (type == "CRISPR"){
		genID="CRISPR"crisprNm
		crisprNm=crisprNm+1
		genName="crispr"
		product="crispr-array"
	}
	# if last one corresponds to a gene, and not a contig name (with no gene found)
    if (NR != lastFeatureNR){
        print start, end, strand, type, prefix contigLoc contigID "_" genID, genName, "| " product " | " ECnum " | " inf2
    }
}


