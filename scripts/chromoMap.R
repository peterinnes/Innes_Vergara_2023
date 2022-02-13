#### Plot positions of CsPT1,4, and 7 in the cs10 genome ####
install.packages("chromoMap")
library(chromoMap)
library(dplyr)
library(tibble)

chrom_file <- read.table("chromoMap/chromosome_file.txt")
# second column of the chrom file is the start position
chrom_file[2] <- rep(1,10)

# lengths (end position) of each chromosome, 1 - 9, plus X 
chrom_file[3] <- c(101209240, 96346938, 94670641, 91913879, 88181582, 79335105, 71238074, 64622176, 61561104, 104987320)

x_chrom <- chrom_file[10,]
names(x_chrom) <- NULL
x_chrom[2] <- 62000000
x_chrom[3] <- 62500000
# read-on bed files containing positions of CsPT1,4, and 7

CsPT1_bed <- read.table("getfasta_results/exons_CsPT1_cs10.bed")
CsPT4_bed <- read.table("getfasta_results/exons_CsPT4.1_cs10.bed")
CsPT7_bed <- read.table("getfasta_results/exons_CsPT7_cs10.bed")

# join all three
CsPTx_annots <- rbind(CsPT1_bed,CsPT4_bed,CsPT7_bed) %>%
  select(1:4) %>%
  rownames_to_column()

names(CsPTx_annots) <- NULL

chromoMap(list(x_chrom), list(CsPTx_annots), segment_annotation = T)


