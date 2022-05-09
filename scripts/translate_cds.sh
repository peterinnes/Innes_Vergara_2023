# script to translate nucleotide sequence into amino acid sequence using EMBOSS transeq
out_dir=/data0/CannabisGenomics2021/Everybody/petah_pahka/getfasta_results

cd $out_dir
for fasta in `ls CDS*`; do
    gene=$(echo $fasta | sed 's/\.fasta//')
    transeq -trim -sequence $fasta -outseq aa_${gene}.fasta
done
    
