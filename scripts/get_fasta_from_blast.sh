# take coordinates from BLAST hits and extract nucleotide sequence from the referance fasta
# Peter Innes

if [ $# -lt 3 ] 
	then
		echo "
		This program extracts nucleotide sequences from a reference fasta based on coordinates from BLAST results,
		using bedtools and the BED file format. 
		
		Required options: 
		
		[-i] table with all blast results files and corresponding gene name, assembly name, and  assembly paths (tab separated)
		[-o] output directory
        [-b] directory with blast results
        "
	else
		while getopts ":i:o:b:" option; do 
		case $option in
		i) input=$OPTARG;;
        o) outdir=$OPTARG;;
        b) blastdir=$OPTARG;;
		esac
		done
		
        cat $input | while read -r blast gene_name ref_name ref_fasta; do
            echo "gene is $gene_name, ref is $ref_fasta..."
            
		        # First step is to convert BLAST format to BED format, using awk. Important to make note if gene is on  plus or minus strand.
            # Also, subtract 1 from BED start position because BLAST is 0-based while BED is 1-based.
		    
            awk 'BEGIN {OFS="\t"}; {
		    if ($9 > $10)
		    	print $2, $10-1, $9, "'"$gene_name"'", $3, "-";
		    else
		    	print $2, $9-1, $10, "'"$gene_name"'", $3, "+";
		    }' $blastdir$blast > ${outdir}exons_${gene_name}_$ref_name.bed

		    # Next, use bedtools getfasta to extract nucleotide sequence. Use the -s flag to force strandedness (reverse complement if on minus strand)
		    /data0/CannabisGenomics2021/Everybody/petah_pahka/software/bedtools2/bin/fastaFromBed \
		    -fi $ref_fasta \
		    -bed ${outdir}exons_${gene_name}_$ref_name.bed \
		    -s \
		    -fo ${outdir}exons_${gene_name}_$ref_name.fasta
            
        # Finally, merge the exons. Use sed to get rid of funky pipe character in chromosome names
        # because it causes union -sequence to fail. 
        cat ${outdir}exons_${gene_name}_$ref_name.fasta | sed -e 's/|//g' \
        | union -sequence stdin -outseq ${outdir}CDS_${gene_name}_$ref_name.fasta
        
        echo ">${gene_name}_${ref_name}" > temp
        tail -n +2 ${outdir}CDS_${gene_name}_$ref_name.fasta >> temp
        mv temp ${outdir}CDS_${gene_name}_$ref_name.fasta
        done
fi
