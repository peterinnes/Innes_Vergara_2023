# This code is used for trimming/correcting exon positions
# in BLAST results for genes CsPT1-7 in each of the 9 Cannabis assemblies.
# These genes all have 10+ exons, with introns that are occasionally very large (15kb+).
# Often, it happens that the end of one exon and the beginning of the next are similar enough
# that there is a small bp overlap when comparing the query end and query start positions of exon i and exon i+1.
# In these cases, we need to edit the correponding subject positions such that we do not extract 
# additional nucleotides in downstream steps.

# Input: file with filepaths of all the blast results to be processed.
# Each individual blast results file should comprise a single gene, and a single copy of that gene. 
# i.e. if a gene has multiple copies in a particular assembly (e.g. CsPT4 appears three times in cs10 assembly) 
# Then the copies need to be split into their own files. 
# Further more, each file should contain only the exons/blast hit that comprise the gene. 
# Therefore, a manual processing step that involves deleting erroneous blast hits/exons, 
# as well as splitting gene copies into separate file, is required before running this code
# NOTE that gaps resulting from insertions in the subject (i.e. the particular assembly that we are searching through) are not accounted for here and will require manual editing. An example of this is found at the boundary of exons 5 and 6 for the CsPT1 gene, in most of the assemblies. 
# 
# Outputs of this script are used in the sript get_fasta_from_blast.sh, which first coverts the blast outfmt 6 files used here to the bed format, and then uses bedtools getfasta to retrieve the nucleotide sequence from assemblies. 

# IMPORTANT!!!! 'cleaned' blast results files must be sorted by the query start position (column 7)
#

system("rm ~/student_folders/Everybody/petah_pahka/list_of_blast_results_files.txt")
system("for file in blast_results/cleaned/*.cleaned.blast; do realpath $file >> ~/student_folders/Everybody/petah_pahka/list_of_blast_results_files.txt; done")
blast_files <- readLines("list_of_blast_results_files.txt")

for ( file in blast_files ){
  
  file_name <- strsplit(file, "/")[[1]][8]
  blast <- read.table(file)
  
  cat("gene is:",file_name, "\n")
  for( i in 1:(nrow(blast)-1) ){
    # skip exon (do not edit its positions) if it is a duplicate of the next exon.
    skip <- FALSE
    if( blast[i,7]==blast[i+1,7] ){
      skip <- TRUE
      print(skip)
    } 
    
    if ( !skip ){
      # check if plus or minus strand
      if( blast[i,9] < blast[i,10] ){
        strand <- "plus"
      } else {
        strand <- "minus"
      }
      
      # check for overlap, in query position, between end of preceding (current) exon
      # and the start of the trailing (next) exon
      overlap <- blast[i,8] - blast[i+1,7] + 1 #add one because inclusive
      PID_preceding <- blast[i,3]
      PID_trailing <- blast[i+1,3]
      cat("overlap between exon", i, "and exon", i+1, "is:", overlap, "\n")
      if( overlap>0 ){ #if an overlap exists
        # check which exon has higher PID. We will edit the position
        # of the one with lower PID. If PID is the same, edit end of preceding exon
        if( PID_preceding <= PID_trailing ) {
          # edit the query end position of preceeding (current) exon
          blast[i,8] <- blast[i,8] - overlap
          # edit the corresponding subject end position
          if ( strand=="plus" ){
            blast[i,10] <- blast[i,10] - overlap
          } else {
            blast[i,10] <- blast[i,10] + overlap
          }
        } else if( PID_preceding > PID_trailing ) {
          # edit the query start position of trailing (next) exon
          blast[i+1,7] <- blast[i+1,7] + overlap
          # edit the corresponding subject start position
          if ( strand=="plus" ){
            blast[i+1,9] <- blast[i+1,9] + overlap
          } else {
            blast[i+1,9] <- blast[i+1,9] - overlap
          }
        } 
      } else if( overlap<0 ) {
        cat("no overlap between exon", i, "and exon", i+1, "\n")
      }
    }  
  }
  write.table(blast, file=paste("blast_results/cleaned_and_trimmed/",file_name,".trimmed", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
}

#### Manual edits to blast coordinates to get rid of spurious gaps at exon boundaries ####
# Note we are only updating subject positions in the blast results; leaving query pos as they are
# Fixing exon 5 / 6 boundary in CsPT1

cs10_CsPT1 <- read.table("blast_results/cleaned_and_trimmed/cs10_CsPT1.cleaned.blast.trimmed") 
cs10_CsPT1[6,9] <- cs10_CsPT1[6,9] - 1 #minus strand
write.table(cs10_CsPT1, file = "blast_results/cleaned_and_trimmed/cs10_CsPT1.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Cannatonic_CsPT1 <- read.table("blast_results/cleaned_and_trimmed/Cannatonic_CsPT1.cleaned.blast.trimmed") 
Cannatonic_CsPT1[6,9] <- Cannatonic_CsPT1[6,9] - 1 #minus strand
write.table(Cannatonic_CsPT1, file = "blast_results/cleaned_and_trimmed/Cannatonic_CsPT1.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Jamaican_Lion_CsPT1 <- read.table("blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT1.cleaned.blast.trimmed") 
Jamaican_Lion_CsPT1[6,9] <- Jamaican_Lion_CsPT1[6,9] + 1 #plus strand
write.table(Jamaican_Lion_CsPT1, file = "blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT1.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

PBBK_CsPT1 <- read.table("blast_results/cleaned_and_trimmed/PBBK_CsPT1.cleaned.blast.trimmed") 
PBBK_CsPT1[6,9] <- PBBK_CsPT1[6,9] + 1 #plus strand
write.table(PBBK_CsPT1, file = "blast_results/cleaned_and_trimmed/PBBK_CsPT1.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Purple_Kush_CsPT1 <- read.table("blast_results/cleaned_and_trimmed/Purple_Kush_CsPT1.cleaned.blast.trimmed") 
Purple_Kush_CsPT1[6,9] <- Purple_Kush_CsPT1[6,9] + 1 #plus strand
write.table(Purple_Kush_CsPT1, file = "blast_results/cleaned_and_trimmed/Purple_Kush_CsPT1.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

USO31_CsPT1 <- read.table("blast_results/cleaned_and_trimmed/USO31_CsPT1_2021-11-22.cleaned.blast.trimmed") 
USO31_CsPT1[6,9] <- USO31_CsPT1[6,9] - 1 #minus strand
write.table(USO31_CsPT1, file = "blast_results/cleaned_and_trimmed/USO31_CsPT1_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

# Fixes to CsPT4? All are missing the first 7 bp so add these; Finola CsPT4 and cs10 CsPT4.2 have additional issues.
Finola_CsPT4 <- read.table("blast_results/cleaned_and_trimmed/Finola_CsPT4_2021-11-22.cleaned.blast.trimmed")
Finola_CsPT4[3,9] <- Finola_CsPT4[3,9] + 1 #plus strand
Finola_CsPT4[1,9] <- Finola_CsPT4[1,9] + 7 #minus strand
write.table(Finola_CsPT4, file = "blast_results/cleaned_and_trimmed/Finola_CsPT4_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Cannatonic_CsPT4.1 <- read.table("blast_results/cleaned_and_trimmed/Cannatonic_CsPT4.1_2021-11-22.cleaned.blast.trimmed")
Cannatonic_CsPT4.1[1,9] <- Cannatonic_CsPT4.1[1,9] - 7 #plus strand
write.table(Cannatonic_CsPT4.1, file="blast_results/cleaned_and_trimmed/Cannatonic_CsPT4.1_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Cannatonic_CsPT4.2 <- read.table("blast_results/cleaned_and_trimmed/Cannatonic_CsPT4.2_2021-11-22.cleaned.blast.trimmed")
Cannatonic_CsPT4.2[1,9] <- Cannatonic_CsPT4.2[1,9] + 7 #minus strand
write.table(Cannatonic_CsPT4.2, file="blast_results/cleaned_and_trimmed/Cannatonic_CsPT4.2_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

cs10_CsPT4.1 <- read.table("blast_results/cleaned_and_trimmed/cs10_CsPT4.1_2021-11-22.cleaned.blast.trimmed")
cs10_CsPT4.1[1,9] <- cs10_CsPT4.1[1,9] - 7 #plus strand
write.table(cs10_CsPT4.1, file="blast_results/cleaned_and_trimmed/cs10_CsPT4.1_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

cs10_CsPT4.2 <- read.table("blast_results/cleaned_and_trimmed/cs10_CsPT4.2_2021-11-22.cleaned.blast.trimmed")
cs10_CsPT4.2[1,9] <- cs10_CsPT4.2[1,9] - 7 #plus strand
cs10_CsPT4.2[5,10] <- cs10_CsPT4.2[5,10] + 1 #plus strand, lengthen exon 5 by 1 bp
cs10_CsPT4.2[7,9] <- cs10_CsPT4.2[7,9] - 5 #plus strand, lengthen start of exon 7 by 5bp
cs10_CsPT4.2[7,10] <- cs10_CsPT4.2[7,10] + 1 #plus strand, lengthen end of exon 7 by 1 bp
write.table(cs10_CsPT4.2, file="blast_results/cleaned_and_trimmed/cs10_CsPT4.2_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")
cs

cs10_CsPT4.3 <- read.table("blast_results/cleaned_and_trimmed/cs10_CsPT4.3_2021-11-22.cleaned.blast.trimmed")
cs10_CsPT4.3[1,9] <- cs10_CsPT4.3[1,9] - 7 #plus strand
write.table(cs10_CsPT4.3, file="blast_results/cleaned_and_trimmed/cs10_CsPT4.3_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Jamaican_Lion_CsPT4.1 <- read.table("blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT4.1_2021-11-22.cleaned.blast.trimmed")
Jamaican_Lion_CsPT4.1[1,9] <- Jamaican_Lion_CsPT4.1[1,9] + 7 #minus strand
write.table(Jamaican_Lion_CsPT4.1, file="blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT4.1_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Jamaican_Lion_CsPT4.2 <- read.table("blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT4.2_2021-11-22.cleaned.blast.trimmed")
Jamaican_Lion_CsPT4.2[1,9] <- Jamaican_Lion_CsPT4.2[1,9] + 7 #minus strand
write.table(Jamaican_Lion_CsPT4.2, file="blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT4.2_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Jamaican_Lion_CsPT4.3 <- read.table("blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT4.3_2021-11-22.cleaned.blast.trimmed")
Jamaican_Lion_CsPT4.3[1,9] <- Jamaican_Lion_CsPT4.3[1,9] + 7 #minus strand
write.table(Jamaican_Lion_CsPT4.3, file="blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT4.3_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Carmagnola_CsPT4 <- read.table("blast_results/cleaned_and_trimmed/C24_CsPT4_2021-11-22.cleaned.blast.trimmed")
Carmagnola_CsPT4[1,9] <- Carmagnola_CsPT4[1,9] - 7 #plus strand
write.table(Carmagnola_CsPT4, file = "blast_results/cleaned_and_trimmed/C24_CsPT4_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

isolate_JL_CsPT4 <- read.table("blast_results/cleaned_and_trimmed/isolate_JL_CsPT4_2021-11-22.cleaned.blast.trimmed")
isolate_JL_CsPT4[1,9] <- isolate_JL_CsPT4[1,9] - 7 #plus strand
write.table(isolate_JL_CsPT4, file = "blast_results/cleaned_and_trimmed/isolate_JL_CsPT4_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Purple_Kush_CsPT4 <- read.table("blast_results/cleaned_and_trimmed/Purple_Kush_CsPT4_2021-11-22.cleaned.blast.trimmed")
Purple_Kush_CsPT4[1,9] <- Purple_Kush_CsPT4[1,9] + 7 #minus strand
write.table(Purple_Kush_CsPT4, file = "blast_results/cleaned_and_trimmed/Purple_Kush_CsPT4_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

USO31_CsPT4 <- read.table("blast_results/cleaned_and_trimmed/USO31_CsPT4_2021-11-22.cleaned.blast.trimmed")
USO31_CsPT4[1,9] <- USO31_CsPT4[1,9] - 7 #plus strand
write.table(USO31_CsPT4, file = "blast_results/cleaned_and_trimmed/USO31_CsPT4_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")


# Fixes to CsPT7
cs10_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/cs10_CsPT7_2021-11-22.cleaned.blast.trimmed")
cs10_CsPT7[6,9] <- cs10_CsPT7[6,9] - 1 #minus strand
write.table(cs10_CsPT7, file="blast_results/cleaned_and_trimmed/cs10_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Carmagnola_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/C24_CsPT7_2021-11-22.cleaned.blast.trimmed")
Carmagnola_CsPT7[6,9] <- Carmagnola_CsPT7[6,9] - 1 #minus strand
write.table(Carmagnola_CsPT7, file="blast_results/cleaned_and_trimmed/Carmagnola_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Cannatonic_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/Cannatonic_CsPT7_2021-11-22.cleaned.blast.trimmed")
Cannatonic_CsPT7[6,9] <- Cannatonic_CsPT7[6,9] - 1 #minus strand
write.table(Cannatonic_CsPT7, file="blast_results/cleaned_and_trimmed/Cannatonic_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Finola_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/Finola_CsPT7_2021-11-22.cleaned.blast.trimmed")
Finola_CsPT7[6,9] <- Finola_CsPT7[6,9] - 1 #minus strand
write.table(Finola_CsPT7, file="blast_results/cleaned_and_trimmed/Finola_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Jamaican_Lion_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT7_2021-11-22.cleaned.blast.trimmed")
Jamaican_Lion_CsPT7[6,9] <- Jamaican_Lion_CsPT7[6,9] + 1 #plus strand
write.table(Jamaican_Lion_CsPT7, file="blast_results/cleaned_and_trimmed/Jamaican_Lion_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

PBBK_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/PBBK_CsPT7_2021-11-22.cleaned.blast.trimmed")
PBBK_CsPT7[6,9] <- PBBK_CsPT7[6,9] + 1 #plus strand
write.table(PBBK_CsPT7, file="blast_results/cleaned_and_trimmed/PBBK_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

Purple_Kush_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/Purple_Kush_CsPT7_2021-11-22.cleaned.blast.trimmed")
Purple_Kush_CsPT7[6,9] <- Purple_Kush_CsPT7[6,9] + 1 #plus strand
write.table(Purple_Kush_CsPT7, file="blast_results/cleaned_and_trimmed/Purple_Kush_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")

USO31_CsPT7 <- read.table("blast_results/cleaned_and_trimmed/USO31_CsPT7_2021-11-22.cleaned.blast.trimmed")
USO31_CsPT7[6,9] <- USO31_CsPT7[6,9] - 1 #minus strand
write.table(USO31_CsPT7, file="blast_results/cleaned_and_trimmed/USO31_CsPT7_2021-11-22.cleaned.blast.trimmed", row.names = F, col.names = F, quote = F, sep = "\t")




