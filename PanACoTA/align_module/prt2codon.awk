# 1st file is the aln PRT file in fasta format, 2nd arg is the genes in 
# fasta format

# current: name of current protein in alignment PRT file or genes file
# elt = sequence number in alignment PRT file
# protseq: tab[current] with full PRT sequence of 'current' protein ('current' is the protein ID)
# genseq: tab[current] with full nuc sequence of 'current' protein ('current' is the gene ID)
# id: match between sequence number (elt) and sequence id (current)

BEGIN{
	elt =0
}

# new sequence
/^>/{
	# if sequence is protein
	if (FILENAME == ARGV[1]){
		current = substr($1, 2)  # get id of prot
		id[elt] = current  # id: tab with prot ids
		protseq[current] = ""
		elt++
	} # if sequence is gene
	else {
		current = substr($1, 2)
		genseq[current] = ""  
	}
	next
}

# Reading Prot sequence line
(FILENAME == ARGV[1]){
	protseq[current]= protseq[current] $1
	next
}

# Reading Gene sequence line
{
genseq[current] = genseq[current] $1
}

END{
	for (i=0; i<elt; i++){
		seq =""
		len = length(protseq[id[i]])
		prtpointr=1
		print ">" id[i] # "  " genseq[id[i]]

		for (j=1; j<=len; j++){
			if (substr(protseq[id[i]], j, 1)!="-"){
				seq = seq substr(genseq[id[i]], prtpointr*3-2, 3)
#				print substr(genseq[id[i]], prtpointr*3-2, 3)
				prtpointr++
			}
			else {
				seq = seq "---"
#				print "---"
			}
		}
		len = length(seq)

		for (j=1; j<=len; j+=60)
			print substr(seq, j, 60)
	}
}

