# take coordinates from BLAST hits and extract nucleotide sequence from the referance fasta
# Peter Innes
# last updated:

if [ $# -lt 4 ] 
	then
		echo "
		This program extracts nucleotide sequences from a reference fasta based on coordinates from BLAST results,
		using bedtools and the BED file format. 
		
		Required options: 
		
		[-i] BLAST results in outfmt 6
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
		  
		# First step is to convert BLAST format to BED format, using awk. Important to make note if gene is on  plus or minus strand.
        # Also, subtract 1 from BED start position because BLAST is 0-based while BED is 1-based.
		awk 'BEGIN {OFS="\t"}; {
		if ($9 > $10)
			print $2, $10-1, $9, "'"$gene_name"'", $3, "-";
		else
			print $2, $9-1, $10, "'"$gene_name"'", $3, "+";
		}' $blast_input > exons_${gene_name}_$ref_name.bed

		# Next, use bedtools getfasta to extract nucleotide sequence. Use the -s flag to force strandedness (reverse complement if on minus strand)
		/data0/CannabisGenomics2021/Everybody/petah_pahka/software/bedtools2/bin/fastaFromBed \
		-fi $ref_fasta \
		-bed exons_${gene_name}_$ref_name.bed \
		-s \
		-fo exons_${gene_name}_$ref_name.fasta
        
        # Finally, merge the exons
        union -sequence exons_${gene_name}_$ref_name.fasta -outseq mRNA_${gene_name}_$ref_name.fasta
fi
