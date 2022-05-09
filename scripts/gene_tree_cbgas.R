library(msa)
library(Biostrings)
library(ape)
library(phangorn)

c24 <- readAAStringSet("C24/CBGAs_C24.mRNA.fasta")
cannatonic <- readAAStringSet("Cannatonic/CBGAs_Cannatonic.mRNA.fasta")
wild <- readAAStringSet("C_sativa_isolate_JL/CBGAs_vs_Gao.mRNA.fasta")
cs10 <- readAAStringSet("cs10/CBGAs_cs10.mRNA.fasta")
finola <- readAAStringSet("Finola/CBGAs_Finola.mRNA.fasta")
jamaican_lion <- readAAStringSet("Jamaican_Lion/CBGAs_Jamaican_Lion.mRNA.fasta")
pbbk <- readAAStringSet("PBBK_steepAsm/CBGAs_PBBK_steepAsm.mRNA.fasta")
purp_kush <- readAAStringSet("Purple_Kush/CBGAs_Purple_Kush.mRNA.fasta")
uso31 <- readAAStringSet("USO31/CBGAs_USO31.mRNA.fasta")
ref <- readAAStringSet("geranyltransferase.fasta")

seqs <- AAStringSet(c(c24, cannatonic, wild, cs10, finola, jamaican_lion, pbbk, purp_kush, uso31, ref))
head(seqs)
names(seqs) <- c("c24", "cannatonic", "wild (Gao et al)", "cs10", "finola", "jamaican_lion", "pbbk", "purp_kush", "uso31", "reference (patent)")
?msa
alignment <- msa(seqs,method='Muscle')
print(alignment)

AAStr <- as(alignment, "AAStringSet")
writeXStringSet(AAStr, file="alignment.fasta")

# Transform the alignment file that you just generated into the 'phyDat' format that the phangorn package requires
cbgas_alignment <- read.phyDat("alignment.fasta", format="fasta",type="AA")

dm <- dist.ml(cbgas_alignment)
head(dm)
# Infer a UPGMA tree based on the pairwise differences between sequences. UPGMA 
# stands for Unweighted Pair Group Method with Arithmetic Mean. This is a hierarchical 
# clustering method, meaning that pairs of taxa are clustered into a higher level cluster, 
# which is clustered to the next most similar cluster, and so on.
cbgas_UPGMA <- upgma(dm)
cbgas_NJ <- NJ(dm) #make a neighbor joining tree for fun. This is an alternative distance matrix-based tree building method.

cbgas_upgma_tree <- plot.phylo(cbgas_UPGMA, main="UPGMA gene tree of \n geranylpyrophosphate:olivetolate geranyltransferase \n aka CBGA Synthase")

#### Dn/Ds with ape, started 10.25.21 ####
# I think intraspecific Dn/Ds may not be valid? See Kryazhimskiy and Plotkin 2008
c24_nucl <- readDNAStringSet("C24/CBGAs_C24.mRNA.fasta")
cannatonic_nucl <- readDNAStringSet("Cannatonic/CBGAs_Cannatonic.mRNA.fasta")
wild_nucl <- readDNAStringSet("C_sativa_isolate_JL/CBGAs_C_sativa_Isolate_JL.mRNA.fasta")
cs10_nucl <- readDNAStringSet("cs10/CBGAs_cs10.mRNA.fasta")
finola_nucl <- readDNAStringSet("Finola/CBGAs_Finola.mRNA.fasta")
jamaican_lion_nucl <- readDNAStringSet("Jamaican_Lion/CBGAs_Jamaican_Lion.mRNA.fasta")
purp_kush_nucl <- readDNAStringSet("Purple_Kush/CBGAs_Purple_Kush.mRNA.fasta")
uso31_nucl <- readDNAStringSet("USO31/CBGAs_USO31.mRNA.fasta")
ref_nucl <- readDNAStringSet("geranyltransferase.fasta")
pbbk_nucl <- readDNAStringSet("PBBK_steepAsm/CBGAs_PBBK_steepAsm.mRNA.fasta")

dnds_seqs <- DNAStringSet(c(c24_nucl, cannatonic_nucl, wild_nucl, cs10_nucl,
                            finola_nucl, jamaican_lion_nucl, purp_kush_nucl, uso31_nucl, ref_nucl, pbbk_nucl))

names(dnds_seqs) <- c("c24", "cannatonic", "putative wild", "cs10", "finola", "jamaican_lion", "pbbk", "purp_kush", "uso31", "reference (patent)")
?msa
dnds_alignment <- as.DNAbin(msa(dnds_seqs,method='Muscle'))
dnds_results <- dnds(dnds_alignment)

dnds_results[which(dnds_results=="Inf")] <- NA
mean(dnds_results, na.rm = T) # 0.419

#### alignment and tree with all CsPTx genes (DNA sequence) ####
CsPTx_seqs <- readDNAStringSet("getfasta_results/merged_CDS_CsPTx.fasta")
CsPTx_alignment <- msa(CsPTx_seqs, method = "ClustalOmega")

CsPTx_Str <- as(CsPTx_alignment, "DNAStringSet")
writeXStringSet(CsPTx_Str, file="getfasta_results/merged_CDS_CsPTx.fasta.ClustalOmega")
CsPTx_phyDat <- read.phyDat("getfasta_results/merged_CDS_CsPTx.fasta.ClustalOmega", format = "fasta", type = "DNA")
dm <- dist.ml(CsPTx_phyDat)
CsPTx_NJ <- NJ(dm)
CsPTx_tree <- plot.phylo(CsPTx_NJ, main="NJ gene tree of \n Cannabis geranylpyrophosphate:olivetolate \n geranyltransferase genes")

# ML phylogeny
fit <- pml(CsPTx_NJ, CsPTx_phyDat)
print(fit)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
?plotBS
