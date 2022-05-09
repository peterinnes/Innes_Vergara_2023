#### Plot positions of CsPT1,4, and 7 in the cs10 genome ####
install.packages("chromoMap")
#BiocManager::install("biomaRt")
#BiocManager::install("Sushi")
devtools::install_github("greymonroe/genemodel")


library(chromoMap) #this is the main package we are using, to draw gene location on chromosome
library(genemodel) #package to draw exon structure of genes
#library(Sushi) #a different package to makes plots of exon structure
library(dplyr)
library(tidyr)
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


CsPT4.1_bed <- read.table("getfasta_results/exons_CsPT4.1_cs10.bed")
CsPT4.2_bed <- read.table("getfasta_results/exons_CsPT4.2_cs10.bed")
CsPT4.3_bed <- read.table("getfasta_results/exons_CsPT4.3_cs10.bed")
CsPT1_bed <- read.table("getfasta_results/exons_CsPT1_cs10.bed")
CsPT7_bed <- read.table("getfasta_results/exons_CsPT7_cs10.bed")

# join all 5
CsPTx_annots <- rbind(CsPT4.1_bed,CsPT4.2_bed,CsPT4.3_bed,CsPT1_bed,CsPT7_bed) %>%
  select(1:4) %>%
  rownames_to_column()
names(CsPTx_annots) <- NULL

chromoMap(list(x_chrom), list(CsPTx_annots), segment_annotation = T,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("red", "orange", "yellow", "green", "blue")))

#### make gene models for each gene using genemodel package ####
# NOTE: I haven't gotten this to work yet but I think it's just because R Studio needs to be updated on Keepers
# need our gene/exon coordiates in a different format
CsPT4.1_gm <- data.frame(CsPT4.1_bed) %>%
  unite("coordinates", 2:3, sep = "-") %>%
  mutate(type="coding_region") %>%
  select(type, coordinates)

genemodel.plot(CsPT4.1_gm, 62197791, 62216905, "forward")

plot(CsPT4.1_bed$start, CsPT4.1_bed$stop)
