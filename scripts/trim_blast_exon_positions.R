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

#IMPORTANT!!!! 'cleaned' blast results files must be sorted by the query start position (column 7)

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
