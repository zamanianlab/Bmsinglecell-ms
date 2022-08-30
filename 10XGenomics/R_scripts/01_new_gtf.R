# determining the ideal extension length for each gene for maximal mapping rates across the transcriptome. 

# Need to increase vector memory to 100GB prior to beginning or you will get nothing done. 
options(future.globals.maxSize = 8000 * 1024^2)

# install packages
install.packages("hdf5r")
install.packages("readbitmap")

# read in libraries
library(tidyverse)
library(DropletUtils)
library(ggplot2)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(rhdf5)
library(data.table)

# set working directory
setwd("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY")

# set directories for each raw ext data set
data.dir.0 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_0bp_ext/filtered_feature_bc_matrix/")
data.dir.50 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_50bp_ext/filtered_feature_bc_matrix/")
data.dir.100 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_100bp_ext/filtered_feature_bc_matrix/")
data.dir.150 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_150bp_ext/filtered_feature_bc_matrix/")
data.dir.200 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_200bp_ext/filtered_feature_bc_matrix/")
data.dir.250 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_250bp_ext/filtered_feature_bc_matrix/")
data.dir.300 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_300bp_ext/filtered_feature_bc_matrix/")
data.dir.400 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_400bp_ext/filtered_feature_bc_matrix/")
data.dir.500 <- ("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_500bp_ext/filtered_feature_bc_matrix/")

# read in filtered data
filt_0 <- Read10X(data.dir.0)
filt_50 <- Read10X(data.dir.50)
filt_100 <- Read10X(data.dir.100)
filt_150 <- Read10X(data.dir.150)
filt_200 <- Read10X(data.dir.200)
filt_250 <- Read10X(data.dir.250)
filt_300 <- Read10X(data.dir.300)
filt_400 <- Read10X(data.dir.400)
filt_500 <- Read10X(data.dir.500)

# create Seurat objects
bma_0 <- Seurat::CreateSeuratObject(counts = filt_0)
bma_50 <- Seurat::CreateSeuratObject(counts = filt_50)
bma_100 <- Seurat::CreateSeuratObject(counts = filt_100)
bma_150 <- Seurat::CreateSeuratObject(counts = filt_150)
bma_200 <- Seurat::CreateSeuratObject(counts = filt_200)
bma_250 <- Seurat::CreateSeuratObject(counts = filt_250)
bma_300 <- Seurat::CreateSeuratObject(counts = filt_300)
bma_400 <- Seurat::CreateSeuratObject(counts = filt_400)
bma_500 <- Seurat::CreateSeuratObject(counts = filt_500)


# write a csv file containing the count matrix from the Seurat object
bma_0_cts <- write.table(as.matrix(GetAssayData(object = bma_0, slot = "counts")), 
            '~/Desktop/bma_0_counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

bma_50_cts <- write.table(as.matrix(GetAssayData(object = bma_50, slot = "counts")), 
                         '~/Desktop/bma_50_counts.csv', 
                         sep = ',', row.names = T, col.names = T, quote = F)

bma_100_cts <- write.table(as.matrix(GetAssayData(object = bma_100, slot = "counts")), 
                          '~/Desktop/bma_100_counts.csv', 
                          sep = ',', row.names = T, col.names = T, quote = F)

bma_150_cts <- write.table(as.matrix(GetAssayData(object = bma_150, slot = "counts")), 
                           '~/Desktop/bma_150_counts.csv', 
                           sep = ',', row.names = T, col.names = T, quote = F)

bma_200_cts <- write.table(as.matrix(GetAssayData(object = bma_200, slot = "counts")), 
                           '~/Desktop/bma_200_counts.csv', 
                           sep = ',', row.names = T, col.names = T, quote = F)

bma_250_cts <- write.table(as.matrix(GetAssayData(object = bma_250, slot = "counts")), 
                           '~/Desktop/bma_250_counts.csv', 
                           sep = ',', row.names = T, col.names = T, quote = F)

bma_300_cts <- write.table(as.matrix(GetAssayData(object = bma_300, slot = "counts")), 
                           '~/Desktop/bma_300_counts.csv', 
                           sep = ',', row.names = T, col.names = T, quote = F)

bma_400_cts <- write.table(as.matrix(GetAssayData(object = bma_400, slot = "counts")), 
                           '~/Desktop/bma_400_counts.csv', 
                           sep = ',', row.names = T, col.names = T, quote = F)

bma_500_cts <- write.table(as.matrix(GetAssayData(object = bma_500, slot = "counts")), 
                           '~/Desktop/bma_500_counts.csv', 
                           sep = ',', row.names = T, col.names = T, quote = F)

# import csv files for each count matrix
bma_0_cts <- read.csv("~/Desktop/bma_0_counts.csv")
bma_50_cts <- read.csv("~/Desktop/bma_50_counts.csv")
bma_100_cts <- read.csv("~/Desktop/bma_100_counts.csv")
bma_150_cts <- read.csv("~/Desktop/bma_150_counts.csv")
bma_200_cts <- read.csv("~/Desktop/bma_200_counts.csv")
bma_250_cts <- read.csv("~/Desktop/bma_250_counts.csv")
bma_300_cts <- read.csv("~/Desktop/bma_300_counts.csv")
bma_400_cts <- read.csv("~/Desktop/bma_400_counts.csv")
bma_500_cts <- read.csv("~/Desktop/bma_500_counts.csv")

# add the first column containing the gene names to the rest of the data frame and rename the column
bma_0_cts <- rownames_to_column(bma_0_cts, var = "genes")
bma_50_cts <- rownames_to_column(bma_50_cts, var = "genes")
bma_100_cts <- rownames_to_column(bma_100_cts, var = "genes")
bma_150_cts <- rownames_to_column(bma_150_cts, var = "genes")
bma_150_cts <- rownames_to_column(bma_150_cts, var = "genes")
bma_200_cts <- rownames_to_column(bma_200_cts, var = "genes")
bma_250_cts <- rownames_to_column(bma_250_cts, var = "genes")
bma_300_cts <- rownames_to_column(bma_300_cts, var = "genes")
bma_400_cts <- rownames_to_column(bma_400_cts, var = "genes")
bma_500_cts <- rownames_to_column(bma_500_cts, var = "genes")


# create new column in each count matrix with total sum of reads for each gene across all cells

bma_0_total <- bma_0_cts %>% 
  mutate(total_0bp = rowSums(bma_0_cts[,-1]))

bma_50_total <- bma_0_cts %>% 
  mutate(total_50bp = rowSums(bma_50_cts[,-1]))

bma_100_total <- bma_100_cts %>% 
  mutate(total_100bp = rowSums(bma_100_cts[,-1]))

bma_150_total <- bma_150_cts %>% 
  mutate(total_150bp = rowSums(bma_150_cts[,-1]))

bma_200_total <- bma_200_cts %>% 
  mutate(total_200bp = rowSums(bma_200_cts[,-1]))

bma_250_total <- bma_250_cts %>% 
  mutate(total_250bp = rowSums(bma_250_cts[,-1]))

bma_300_total <- bma_300_cts %>% 
  mutate(total_300bp = rowSums(bma_300_cts[,-1]))

bma_400_total <- bma_400_cts %>% 
  mutate(total_400bp = rowSums(bma_400_cts[,-1]))

bma_500_total <- bma_500_cts %>% 
  mutate(total_500bp = rowSums(bma_500_cts[,-1]))


# create new data.frame with gene column and total columns in each extension-specific matrix
ext_all <- data.frame(bma_0_total$genes, bma_0_total$total_0bp, bma_50_total$total_50bp,bma_100_total$total_100bp, bma_150_total$total_150bp, bma_200_total$total_200bp, bma_250_total$total_250bp, bma_300_total$total_300bp, bma_400_total$total_400bp, bma_500_total$total_500bp)

# rename columns for simplicity (and sanity)
nms <- c("gene_id", "0", "50", "100", "150", "200", "250", "300", "400", "500")
setnames(ext_all, nms)


# identify the ext column with the largest count number in each row and print new column with the optimal ext length
ext_all$final <- as.numeric(colnames(ext_all[,-1])[max.col(ext_all[,-1], ties.method = "first")])


# remove the ext length columns and maintain "gene_id" and "final" columns
ext_all <- select(ext_all, -c("0", "50", "100", "150", "200", "250", "300", "400", "500"))


# read in B. malayi gtf
gtf <- read.table("~/Desktop/geneset.gtf", sep = "\t", header = FALSE, quote = "")


# append the "final" column of ext_all to the gtf based on gene_id 
gtf_new <- gtf %>% 
  mutate(gene_id = str_match(V9, "gene_id \"(.*?)\";")[,2]) # unable to pipe the next line for some reason??

gtf_new <- left_join(gtf_new, ext_all, by = "gene_id") %>% 
  mutate(final = replace_na(final, as.numeric(0)))



# extend gene, transcript, final exon, and 3' UTR by the amount in the "final" column
gtf.ext <- gtf_new %>%
  mutate(transcript_id = str_match(V9, "transcript_id \"(.*?)\";")[,2]) %>%
  mutate(exon_number = str_match(V9, "exon_number \"(\\d+)\";")[,2]) %>%
  mutate(exon_number = replace_na(exon_number, as.character(0))) %>% 
  mutate(exon_number = as.numeric(exon_number)) %>% 
  group_by(transcript_id) %>% 
  mutate(group_max = max(exon_number)) %>% 
  mutate(group_max = ifelse(exon_number == group_max & V3 == "exon", "yes", "no")) %>% 
  ungroup() %>% 
  mutate(V5 = ifelse(group_max == "yes" & V7 == "+", as.numeric(V5) + as.numeric(final), V5)) %>% 
  mutate(V4 = ifelse(group_max == "yes" & V7 == "-", as.numeric(V4) - as.numeric(final), V4)) %>% 
  mutate(V4 = ifelse(group_max == "yes" & V4 < 1, 1, V4)) %>%
  mutate(V5 = ifelse(V3 == "three_prime_utr" & V7 =="+", as.numeric(V5) + as.numeric(final), V5)) %>%
  mutate(V4 = ifelse(V3 == "three_prime_utr" & V7 =="-", as.numeric(V4) - as.numeric(final), V4)) %>%
  mutate(V4 = ifelse(V3 == "three_prime_utr" & V4 < 1, 1, V4)) %>% 
  mutate(V5 = ifelse(V3 == "transcript" & V7 =="+", as.numeric(V5) + as.numeric(final), V5)) %>%
  mutate(V4 = ifelse(V3 == "transcript" & V7 =="-", as.numeric(V4) - as.numeric(final), V4)) %>%
  mutate(V4 = ifelse(V3 == "transcript" & V4 < 1, 1, V4)) %>%
  mutate(V5 = ifelse(V3 == "gene" & V7 =="+", as.numeric(V5) + as.numeric(final), V5)) %>%
  mutate(V4 = ifelse(V3 == "gene" & V7 =="-", as.numeric(V4) - as.numeric(final), V4)) %>%
  mutate(V4 = ifelse(V3 == "gene" & V4 < 1, 1, V4))
  
colnames(gtf.ext) <- c("seqname","source","feature","start","end","score","strand","frame","attribute", "gene_id", "final", "transcript_id", "exon_number", "group_max")


#read in contig lengths
contig.len <- read.table("~/Desktop/contig_lengths.txt", sep = "\t", header = FALSE, quote = "")
colnames(contig.len) <- c("seqname","len")

#truncate if exceeds boundaries of contig 
gtf.ext <- left_join(gtf.ext,contig.len, by = "seqname") %>%
  mutate(end = ifelse(strand == "+" & end >= len, len, end)) %>%
  select(-len, -group_max, -exon_number, -transcript_id, -gene_id, -final)
  
#generate new gtf
write.table(gtf.ext,"~/Desktop/geneset.ext.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


  