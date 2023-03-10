# script to translate nucleotide sequence into amino acid sequence using EMBOSS transeq
run_dir=~/CsPTx_genomics/data/getfasta_results/blast_exons/

cd $run_dir
for fasta in `ls CDS*`; do
    gene=$(echo $fasta | sed 's/\.fasta//')
    transeq -trim -sequence $fasta -outseq aa_${gene}.fasta
done
    
