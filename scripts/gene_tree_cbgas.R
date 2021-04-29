library(msa)
library(Biostrings)
library(ape)
library(phangorn)


c24 <- readAAStringSet("C24/aa_CBGAs_C24.mRNA.fasta")
cannatonic <- readAAStringSet("Cannatonic/aa_CBGAs_Cannatonic.mRNA.fasta")
wild <- readAAStringSet("C_sativa_isolate_JL/aa_CBGAs_vs_Gao.mRNA.fasta")
cs10 <- readAAStringSet("cs10/aa_CBGAs_cs10.mRNA.fasta")
finola <- readAAStringSet("Finola/aa_CBGAs_Finola.mRNA.fasta")
jamaican_lion <- readAAStringSet("Jamaican_Lion/aa_CBGAs_Jamaican_Lion.mRNA.fasta")
pbbk <- readAAStringSet("PBBK_steepAsm/aa_CBGAs_PBBK_steepAsm.mRNA.fasta")
purp_kush <- readAAStringSet("Purple_Kush/aa_CBGAs_Purple_Kush.mRNA.fasta")
uso31 <- readAAStringSet("USO31/aa_CBGAs_USO31.mRNA.fasta")
ref <- readAAStringSet("aa_geranyltransferase.fasta")

seqs <- AAStringSet(c(c24, cannatonic, wild, cs10, finola, jamaican_lion, pbbk, purp_kush, uso31, ref))
head(seqs)
names(seqs) <- c("c24", "cannatonic", "wild (Gao et al)", "cs10", "finola", "jamaican_lion", "pbbk", "purp_kush", "uso31", "reference (patent)")
?msa
alignment <- msa(seqs,method='Muscle')
print(alignment)

AAStr <- as(alignment, "AAStringSet")
writeXStringSet(AAStr, file="aa_alignment.fasta")

# Transform the alignment file that you just generated into the 'phyDat' format that the phangorn package requires
cbgas_alignment <- read.phyDat("aa_alignment.fasta", format="fasta",type="AA")

dm <- dist.ml(cbgas_alignment)
head(dm)
# Infer a UPGMA tree based on the pairwise differences between sequences. UPGMA 
# stands for Unweighted Pair Group Method with Arithmetic Mean. This is a hierarchical 
# clustering method, meaning that pairs of taxa are clustered into a higher level cluster, 
# which is clustered to the next most similar cluster, and so on.
cbgas_UPGMA <- upgma(dm)
cbgas_NJ <- NJ(dm) #make a neighbor joining tree for fun. This is an alternative distance matrix-based tree building method.

cbgas_upgma_tree <- plot.phylo(cbgas_UPGMA, main="UPGMA gene tree of \n geranylpyrophosphate:olivetolate geranyltransferase \n aka CBGA Synthase")

