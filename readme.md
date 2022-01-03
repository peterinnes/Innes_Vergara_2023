##This directory contains code for mining Cannabis assemblies for the geranylpyrophosphate:olivetolate geranyltransferase enzyme (GOT) aka CBGA Synthase. 
#Steps:
#1. blast CsPT1-7 genes to all 11 assemblies (9 cannabis and 2 hops) with the script blastn_all_assemblies.sh.
#2. For each blast results file, manually "clean" the hits. This is necessary because all these genes have 10 to 12 exons. There appear to be some duplicates of some of the exons, so we get rid of these and just keep the set of exons that align to the same region and comprise a complete single copy of the GOT gene. 
#3. Use the script trim_blast_exon_positions.R to trim coordinates of the blast hits (exons) so that exons don't overlap and cause redundant sequence. 
#4. Finally, to get a fasta sequence for each gene from its corresponding reference genome, we must convert the trimmed single copy blast output to bed format, then pull fasta sequence from assembly using bedtools. This is done with 'blast_to_fasta.sh'


