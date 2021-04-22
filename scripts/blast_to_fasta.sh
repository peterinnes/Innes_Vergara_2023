# take coordinates from BLAST hits and extract nucleotide sequence from the referance fasta
# Peter Innes
# last updated 4.22.21

if [ $# -lt 4 ] 
	then
		echo "
		This program extracts nucleotide sequences from a reference fasta based on coordinates from BLAST results,
		using bedtools and the BED file format. 
		
		Required options: 
		
		[-i] BLAST results in outfmt 7
		[-r] reference fasta file
		[-n] reference assembly name
		[-g] gene/sequence name	
		"
	else
		while getopts ":i:r:n:g:" option; do 
		case $option in
		i) blast_input=$OPTARG;;
		r) ref_fasta=$OPTARG;;
		n) ref_name=$OPTARG;;
		g) gene_name=$OPTARG;;
		esac
		done
		
		# First step is to convert BLAST format to BED format, using awk. Important to make note if gene is on  plus or minus strand
		awk 'BEGIN {OFS="\t"}; {
		if ($9 > $10)
			print $2, $10, $9, "'"$gene_name"'", $3, "-";
		else
			print $2, $9, $10, "'"$gene_name"'", $3, "+";
		}' $blast_input > ${gene_name}_$ref_name.exons.bed

		# Next, use bedtools getfasta to extract nucleotide sequence. Use the -s flag to force strandedness (reverse complement if on minus strand)
		bedtools getfasta \
		-fi $ref_fasta \
		-bed ${gene_name}_$ref_name.exons.bed \
		-s \
		-fo ${gene_name}_$ref_name.exons.fasta
fi
