##This directory contains code for mining Cannabis assemblies for the geranylpyrophosphate:olivetolate geranyltransferase enzyme (GOT) aka CBGA Synthase.
#Steps:

1. blast CsPT1, CsPT4, CsPT7 exons to all 11 assemblies (9 cannabis and 2 hops) with the script blastn_exons_CsPTx_all_assemblies.sh.

2. Then we manually parse exon blast results to pull out separate copies of the gene (if present), based on positions of exons. 

4. Finally, to get a fasta sequence for each gene from its corresponding reference genome, we must convert blast output to bed format, then pull fasta sequence from assembly using bedtools. This is done with 'get_fasta_from_blast.sh' using 'CsPTx_exon_get_fasta_input_table.tsv' as input.

4.1. For some reason, the chromosome names for the C24 and USO31 assemblies are not recognized by bedtools getfasta and so we have to use python get the the fasta sequences from these two genomes. This is done with 'bed_to_fasta.py.' using 'CsPTx_bed_to_fasta_python_input_table.csv' as input (i.e. commandline argument)


