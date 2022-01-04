# BLAST-to-BED conversion with awk
# Peter Innes
# 2.26.21

if [ $# -lt 3 ] 
	then
		echo "
		This program converts BLAST results (in outfmt 6) to a six-column BED format for use in downstream analyses (e.g. extracting nucleotide sequences from a reference fasta).
		Required options: 
		
		[-i] input file with BLAST results
		[-o] output file
		[-n] gene/sequence name	
		"
	else
		while getopts ":i:o:n:" option; do 
		case $option in
		i) infile=$OPTARG;;
		o) outfile=$OPTARG;;
		n) name=$OPTARG;;
		esac
		done
		
        # have to subract 1 from the BED start position because bed is 0-based while BLAST is 1-based.
		awk 'BEGIN {OFS="\t"}; {
		if ($9 > $10)
			print $2, $10-1, $9, "'"$name"'", $3, "-";
		else
			print $2, $9-1, $10, "'"$name"'", $3, "+";
		}' $infile > $outfile
fi
