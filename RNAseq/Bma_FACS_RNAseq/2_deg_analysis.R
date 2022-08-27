library(tidyverse)
library(cowplot)
library(paletteer)
library(viridis)
library(DESeq2)
library(magick)
#brew install poppler librsvg cairo
library(grConvert)
library(grImport2)

setwd("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/RNAseq/Bma_MF_ES_scRNA/")
RNAdata <- c("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/")


##########################################
### 1A. Prep raw counts for DESEQ2 analysis of Bma MF bulk RNAseq of large cells via FACS
##########################################

# read in counts (raw) and filter for samples
counts.raw <- readRDS("data/counts-es.raw.rds") 


# convert to matrix
counts.raw <- counts.raw %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")


# make metadata file
sample_list <- colnames(counts.raw)
samples <- data.frame("sample" = sample_list, stringsAsFactors = FALSE)
samples <- samples %>%
  mutate(cond = c("Round_with_wolb", "High_SCC", "High_FSC", "Round_without_wolb", "Highest_SSC", "Highest_FSC")) 
samples$cond <- as.factor(samples$cond)
saveRDS(samples, "data/deg/samples-es.rds")

# read in sample metadata df > samples
samples <- readRDS("data/deg/samples-es.rds")


##########################################
## 1B. Transform and create heatmap for quick gene expression visualization across all samples
##########################################

library(hexbin)
library(wesanderson)

gene_count <- as.matrix(counts.raw)
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ cond)

#filter 1 (require x reads total across all samples)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

#filter 2 (require x reads in each of at least y samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]
nrow(dds) #11709 > 1529


#pick and compare vst / rlog transform
vst <- vst(dds, blind = TRUE)

dds <- estimateSizeFactors(dds)
df <- bind_rows(as_tibble(log2(counts(dds, normalized = TRUE)[, 1:6]+1)) %>% 
                  mutate(transformation = "log2(x + 1)"),
                as_tibble(assay(vst)[,1:6]) %>%  mutate(transformation = "vst"))


#calculate sample distances and use hclust to cluster samples based on dist
sampleDists <- dist(t(assay(vst)), method = "euclidean") #chose vsd over rld
clust.dist <- hclust(sampleDists, method="ward.D2")

#get list order of clustered samples from hclust output
ord <- clust.dist$order

#convert original distance matrix to df, re-order, and tidy for plotting
sampleDists.df <- as.data.frame(as.matrix(sampleDists)) %>%
  rownames_to_column(var = "sample") 
sampleDists.df$sample <- factor(sampleDists.df$sample, 
                                levels = c(sampleDists.df$sample[ord]),
                                labels = c("S1", "S2", "S3", "S4", "S5", "S6"))
sampleDists.df <- sampleDists.df %>%
  pivot_longer(cols = 2:ncol(sampleDists.df),
               names_to = "sample_2",
               values_to = "dist")
sampleDists.df$sample_2 <- factor(sampleDists.df$sample_2, 
                                  levels = c(sampleDists.df$sample_2[ord]), 
                                  labels = c("S1", "S2", "S3", "S4", "S5", "S6"))


#plot
cluster.ht <- ggplot(sampleDists.df, aes(x = sample, y = sample_2, fill = dist)) +
  theme_bw() + xlab('') + ylab('') + labs(fill = "Distance") +
  theme(axis.text.x = element_text(angle=90)) +
  geom_tile() 
cluster.ht

##########################################
## 1C. PCA analyses
##########################################

#Run PCA on vst transformed data and return output for ggplot
pca <- plotPCA(vst, intgroup = c("cond"), returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

#Fig 1B
pca.ht <- ggplot(pca, aes(x = PC1, y = PC2, color = cond)) +
  geom_point(size =2.5, alpha = 0.75) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #scale_color_manual(values = c("#4A544A","#DF6589","#FFCC52")) +
  theme(legend.position = "right", legend.title = element_blank()) +
  coord_fixed()
pca.ht



##########################################
### 1F. Heatmap of genes of interest (vaccine and drug targets) - TPMs
##########################################

# read in counts (raw) and filter for samples
counts.tpm <- readRDS("data/counts-es.tpm.rds") 

# read in sample metadata df > samples
samples <- readRDS("data/deg/samples-es.rds")

# option: filter for gene list
#gene_list_dt <- read.csv("/Users/mzamanian/Box/ZamanianLab/Data/Airs_Experiments/P1_SpatialTranscriptomics_Manuscript/RNAseq/Gene_Lists/drug_targets.csv",
                        # header = TRUE, sep = ",")
#gene_list_dt <- unique(gene_list$gene_id)

#gene_list_ag <- read.csv("/Users/mzamanian/Box/ZamanianLab/Data/Airs_Experiments/P1_SpatialTranscriptomics_Manuscript/RNAseq/Gene_Lists/Antigens.csv",
 #                        header = TRUE, sep = ",")
#gene_list_ag <- unique(gene_list$gene_id)


#counts.tpm <- counts.tpm %>%
#  filter(gene_id %in% gene_list)


#join counts and sample info, widen, declare and normalize matrix
matrix_heatmap <- function (df) {
  df <- left_join(df, samples, by = "sample") %>%
    mutate(sample = paste0(cond)) %>%
    select(gene_id, sample, expression) %>%
    pivot_wider(names_from = sample, values_from = expression) %>%
    column_to_rownames(var = "gene_id")
  
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  #df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  #df.m <- log2(df.m+1)
  df.m <- t(scale(t(log2(df.m+1)),center=TRUE,scale=TRUE)) #or log2(df.m+1)
  return(df.m)
}
counts <- matrix_heatmap(counts.tpm)

#calculate gene distances and use hclust to cluster samples based on dist
geneDists <- dist(counts, method = "euclidean")
gclust.dist <- hclust(geneDists, method="ward.D2")

#get list order of clustered genes from hclust output
ord <- gclust.dist$order

#convert genecount matrix to df, re-order, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "sample", values_to = "expression")
counts$sample <- factor(counts$sample)




library("ZamanianLabThemes")
heatmap <- ggplot(counts, aes(gene_id, sample)) +
  geom_tile(aes(fill = expression), show.legend = TRUE) +
  scale_fill_viridis() +
  ylab("") + xlab("B. malayi Genes") + theme(axis.ticks = element_blank()) +
  theme_bw() + theme_zlab() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme()
#theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) + #element_text(angle = 45, hjust = 1, size = 6)
theme(legend.title=element_blank()) 

ggsave("plots/cluster.pdf", heatmap, width = 10, height = 6, units = "in")



 
