#### Plot positions of CsPT1,4, and 7 in the cs10 genome ####
library(tidyr)
library(chromoMap)
library(GenomicRanges)
library(Gviz)

chrom_file <- read.table("chromoMap/chromosome_file.txt")
temp <- chrom_file[c(1,3)] 
write.table(x = temp, file = "data/annotations/cs10_genome_file.txt", sep = "\t", row.names = F, quote = F, col.names = F)
# second column of the chrom file is the start position
chrom_file[2] <- rep(1,10)

# lengths (end position) of each chromosome, 1 - 9, plus X 
chrom_file[3] <- c(101209240, 96346938, 94670641, 91913879, 88181582, 79335105, 71238074, 64622176, 61561104, 104987320)

x_chrom <- data.frame(chrom="ChrX", start=62150000, end=62450000)

# read-on bed files containing positions of CsPT1,4, and 7
CsPT4.1_bed <- read.table("data/getfasta_results/blast_exons/exons_CsPT4.1_cs10.edit.bed")
CsPT4.2_bed <- read.table("data/getfasta_results/blast_exons/exons_CsPT4.2_cs10.edit.bed")
CsPT4.3_bed <- read.table("data/getfasta_results/blast_exons/exons_CsPT4.3_cs10.bed")
CsPT1_bed <- read.table("data/getfasta_results/blast_exons/exons_CsPT1_cs10.bed")
CsPT7_bed <- read.table("data/getfasta_results/blast_exons/exons_CsPT7_cs10.bed")

# join start and end positions of all 5 genes
CsPTx_annots <- data.frame(name=c("CsPT1","CsPT7","CsPT4.1","CsPT4.2","CsPT4.3"),
                           chrom=rep("ChrX", 5), start=c(62381540,62390802,
                                                         62197784,62265622,
                                                         62334343), 
                           end=c(62388633,62395691,62216905,62277840,62348916))

#CsPTx_annots <- rbind(CsPT4.1_bed,CsPT4.2_bed,CsPT4.3_bed,CsPT1_bed,CsPT7_bed) %>%
#  dplyr::select(c(4,1,2,3)) #%>%
#  #tibble::rownames_to_column()

#CsPTx_annots[2] <- "ChrX"
CsPTx_annots[5] <- CsPTx_annots[1]
names(CsPTx_annots) <- NULL

#purple: #9260a7
#yellow: #f0ab31
#blue: #55c0e4
#green: #60b84e
#grey: #a7a4a6
  
chromoMap(list(x_chrom), list(CsPTx_annots), segment_annotation = T,
          labels = T,
          label_angle = -65,
          chr_length = 5,
          chr_width = 30,
          chr_color = "#605e60",
          text_font_size = 14,
          label_font = 12,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("#60b84e", "#f0ab31", "#55c0e4", "#55c0e4", "#55c0e4")),
          export.options = T)

#### Daniela's code for plotting OAC ####
oac_chrom <- chrom_file[9,]
names(oac_chrom) <- NULL
oac_chrom[2] <- 4200000 #sets the chunk of the chromosome we want to plot. I made it a larger region so it's not such a hugely different scale than the other chromoPlots. 
oac_chrom[3] <- 4400000

oac_bed <- read.table("oac_cs10_23052022R.bed")
oac_annots <- (oac_bed) %>%
  select(c(7,1,2,3,4)) #%>%
  #rownames_to_column()
names(oac_annots) <- NULL
# For some reason, chromoMap won't plot just a single gene/single label. So we need to differentiate the two exons. You can change in illustrator if you want.
oac_annots[1] <- c("OAC exon_1", "OAC exon_2")
oac_annots[5] <- c("OAC exon_1", "OAC exon_2")

chromoMap(list(oac_chrom), list(oac_annots), segment_annotation = T,
          labels = T,
          label_angle = -65,
          chr_length = 5,
          chr_width = 30,
          chr_color = "#605e60",
          text_font_size = 14,
          label_font = 12,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("#9260a7", "#9260a7")), #use same color for both exons
          export.options = T) #for easy export to PNG

#### plot gene models with Gviz package ####
options(ucscChromosomeNames = FALSE)
CsPT4.1 <- read.table("data/getfasta_results/blast_exons/exons_CsPT4.1_cs10.edit.bed", col.names = c("chromosome", "start", "end", "symbol", "score", "strand")) %>%
  mutate(width=end-start+1, feature="protein_coding", gene="CsPT4.1", exon=seq(1, 10, 1), transcript="CsPT4.1", symbol="CsPT4.1", color="blue") %>%
  dplyr::select(chromosome, start, end, width, strand, feature, gene, exon, transcript, symbol, color)

CsPT4.2 <- read.table("data/getfasta_results/blast_exons/exons_CsPT4.2_cs10.edit.bed", col.names = c("chromosome", "start", "end", "symbol", "score", "strand")) %>%
  mutate(width=end-start+1, feature="protein_coding", gene="CsPT4.2", exon=seq(1, 10, 1), transcript="CsPT4.2", symbol="CsPT4.2", color="blue") %>%
  dplyr::select(chromosome, start, end, width, strand, feature, gene, exon, transcript, symbol, color)

CsPT4.3 <- read.table("data/getfasta_results/blast_exons/exons_CsPT4.3_cs10.bed", col.names = c("chromosome", "start", "end", "symbol", "score", "strand")) %>%
  mutate(width=end-start+1, feature="protein_coding", gene="CsPT4.3", exon=seq(1, 10, 1), transcript="CsPT4.3", symbol="CsPT4.3", color="blue") %>%
  dplyr::select(chromosome, start, end, width, strand, feature, gene, exon, transcript, symbol, color)

CsPT1 <- read.table("data/getfasta_results/blast_exons/exons_CsPT1_cs10.bed", col.names = c("chromosome", "start", "end", "symbol", "score", "strand")) %>%
  mutate(width=end-start+1, feature="protein_coding", gene="CsPT1", exon=seq(1, 10, 1), transcript="CsPT1", symbol="CsPT1", color="green") %>%
  dplyr::select(chromosome, start, end, width, strand, feature, gene, exon, transcript, symbol, color)

CsPT7 <- read.table("data/getfasta_results/blast_exons/exons_CsPT7_cs10.bed", col.names = c("chromosome", "start", "end", "symbol", "score", "strand")) %>%
  mutate(width=end-start+7, feature="protein_coding", gene="CsPT7", exon=seq(7, 70, 7), transcript="CsPT7", symbol="CsPT7", color="yellow") %>%
  dplyr::select(chromosome, start, end, width, strand, feature, gene, exon, transcript, symbol, color)

#granges <- makeGRangesFromDataFrame(CsPT4.1_bed, start.field="start", end.field="end", strand.field="strand", seqnames.field=c("seq_id"), keep.extra.columns = T)

grtrack_CsPT4.1 <- GeneRegionTrack(CsPT4.1, name="CsPT4.1",
                                   feature = as.vector(CsPT4.1$color),
                                   blue="#55c0e4")
grtrack_CsPT4.2 <- GeneRegionTrack(CsPT4.2, name="CsPT4.2",
                                   feature = as.vector(CsPT4.2$color),
                                   blue="#55c0e4")
grtrack_CsPT4.3 <- GeneRegionTrack(CsPT4.3, name="CsPT4.3",
                                   feature = as.vector(CsPT4.3$color),
                                   blue="#55c0e4")
grtrack_CsPT1 <- GeneRegionTrack(CsPT1, name="CsPT1",
                                 feature = as.vector(CsPT1$color),
                                 green="#60b84e")
grtrack_CsPT7 <- GeneRegionTrack(CsPT7, name="CsPT7",
                                 feature = as.vector(CsPT7$color),
                                 yellow="#f0ab31")
#### plot the gene models ####
gtrack <- GenomeAxisTrack()

pdf("figures/gene_model_CsPT1.pdf", width = 5.175, height = 1.725)
plotTracks(list(gtrack,grtrack_CsPT1), from=62375087, to=62395087)
dev.off()

pdf("figures/gene_model_CsPT7.pdf", width = 5.175, height = 1.725)
plotTracks(list(gtrack,grtrack_CsPT7), from = 62383247, to = 62403247)
dev.off()

pdf("figures/gene_model_CsPT4.1.pdf", width = 5.175, height = 1.725)
plotTracks(list(gtrack, grtrack_CsPT4.1), from = 62197000 , to = 62217000)
dev.off()

pdf("figures/gene_model_CsPT4.2.pdf", width = 5.175, height = 1.725)
plotTracks(list(gtrack, grtrack_CsPT4.2), from = 62261622, to = 62281622)
dev.off()

pdf("figures/gene_model_CsPT4.3.pdf", width = 5.175, height = 1.725)
plotTracks(list(gtrack,grtrack_CsPT4.3), from =62331629, to = 62351629 )
dev.off()





