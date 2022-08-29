############
## FIG 3 Exploring antigens
############
#data wrangling/analysis
library(tidyverse)
library(Seurat)
library(DESeq2)
library(dplyr)

# plotting
library(magick)
library(pdftools)
library(cowplot)
library(ggplot2)
library(ggdendro)
library(pals)
library(ggtext)
library(ZamanianLabThemes)
library(viridis)

#other
library(here)
library(glue)


# read in seurat object
new_combined <- readRDS(here("Figures/Figure_3/10XGenomics/scmulti_integrated.RDS"))
DefaultAssay(new_combined) <- "RNA"


# pull out the normalized counts, metadata, and UMAP coordinates into dataframes for plotting in ggplot2
data <- as_tibble(new_combined@reductions$umap@cell.embeddings, rownames = 'index') # UMAP coordinates for each cell

md <- as_tibble(new_combined@meta.data, rownames = 'index') # metadata detailing ut/t identity and cluster information

counts2 <- as_tibble(new_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") 


# color palette (30 total)
dakota <- c("#d97f64", "#263946", "#bebab6", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#65838d", "#82aca7", "#a0b4ac", "#b5b9b0", "#fbc1c1", "#e89690", "#d76660", "#cac6b9", "#878787", "#cb8034", "#7f93a2", "#ac8287", "#c1d6d3", "#cd4c42", "#5c8492", "#b25757", "#fe906a", "#6f636b", "#6a9491", "#82ac92", "#a26f6a", "#184459", "#596c7f")



################
### Fig. 3a - Schematic of FACS filtration and sequencing
################
facs_scheme <- ggdraw() +
  draw_image(magick::image_read_pdf(here("Figures/Figure_3/FACS_RNAseq.pdf"))



#################
### Fig. 3b - Bma FACS RNAseq heatmap
#################

# LOAD RDS OBJECT BELOW instead of running code underneath and skip to heatmap on line 146
#counts <- readRDS("data/bm_filtered_counts.RDS")


#Prep raw counts for DESEQ2 analysis of Bma MF bulk RNAseq of large cells via FACS

#read in raw counts
counts.raw <- readRDS(here("RNAseq/Bma_FACS_RNAseq/counts-es.raw.rds"))

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
#saveRDS(samples, "data/deg/samples-es.rds")

# read in sample metadata df > samples
samples <- readRDS(here("RNAseq/Bma_FACS_RNAseq/data/deg/samples-es.rds"))


library(hexbin)
library(wesanderson)

#filter (require 1 read in at least 1 sample)
keep <- counts.raw %>% 
  mutate(keep = ifelse(rowSums(counts.raw >= 1) >1, "keep", "discard")) %>% 
  subset(keep == "keep")
keep <- rownames_to_column(keep, var = "gene_id")
keep <- keep$gene_id




# read in counts (tpm) and filter for samples that meet the filters for raw
counts.tpm <- readRDS(here("RNAseq/Bma_FACS_RNAseq/data/counts-es.tpm.rds"))
  filter(gene_id %in% keep)


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
bma_counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "sample", values_to = "expression")
bma_counts$sample <- factor(bma_counts$sample)



# rename samples
bma_counts$sample <- factor(bma_counts$sample, levels = c( "Round_with_wolb", "Round_without_wolb","High_SCC", "Highest_SSC", "High_FSC", "Highest_FSC"), labels = c("Round (*wBm* +)", "Round (*wBm* -)", "Granular", "Most Granular", "Large", "Largest"))


# plot heatmap
(bma_heatmap <- bma_counts %>% 
  ggplot(aes( y = gene_id,x = sample)) +
  geom_tile(aes(fill = expression), show.legend = TRUE) +
  scale_fill_viridis() +
  scale_x_discrete(limits = c("Largest","Large","Most Granular", "Granular","Round (*wBm* -)", "Round (*wBm* +)"), labels = c("Largest","Large","Most Granular", "Granular","Round (*wBm* -)", "Round (*wBm* +)"))+
  labs(y = "*B. malayi* genes", fill = "Z-Score") +
  theme(#text = element_text(family = "helvetica"),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_markdown(size = 8, angle = 45, hjust = 1, vjust = 1.1),
        axis.title.y = element_markdown(size = 10, angle = 90),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.title=element_markdown(size = 8),
        legend.text = element_text(size = 8),
        legend.text.align = 0.8,
        legend.position = "right",
        legend.margin = margin(0,0,0,-0.2, "cm"),
        plot.margin = margin(0.5, 0, 0, 0, "cm")))

#saveRDS(counts, "data/bm_filtered_tpm_counts.RDS")



# combine figure 3a panels
(fig3a <- plot_grid(facs_scheme, bma_heatmap, ncol = 2, rel_widths = c(1, 0.9), rel_heights = c(0.8, 0.8), scale = c(0.9, 1)))


###################
### Fig.3B - Large cell genes reverse plotted on scRNAseq UMAPS
###################
# retrieve raw counts table from scdata
raw <- as_tibble(new_combined@assays[["RNA"]]@counts, rownames = "gene_id") %>%  # gene expression matrix of raw counts
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") 

# pull out genes with z-score max from "largest" sample
genes <- bma_counts %>% 
  subset(expression >= 2.02) %>% 
  subset(sample == "Largest")


bma_facs <- as.character(genes$gene_id, stringAsFactors = FALSE)


# combine FACS-RNAseq df with the single-cell coordinate df
#filter the raw sc data for genes in "largest" with max z-score
tmp <- raw %>% 
  subset(gene_id %in% bma_facs) %>% 
  #subset(expression > 2.02) %>% 
  left_join(data) %>% 
  left_join(md) %>% 
  subset(counts > 0)

#pivot wider to sum total gene counts per cell and reduce dataframe to cell and gene count
test <- tmp %>% 
  select("gene_id", "index", "counts") %>% 
  pivot_wider(names_from = "gene_id", values_from = "counts") 


test[is.na(test)] <- 0
test$total <- rowSums(test[,-1]) 
test <- test %>% 
  select("index", "total")

# left join the total back to the dataframe
tmp <- tmp %>% 
  left_join(test)



# how many cells are expressing more than 1 of the DEG fron the FACS-RNAseq
tmp$count <- 1
test <- tmp %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 

test[is.na(test)] <- 0
test$markers <- rowSums(test[,-1]) 
test <- test %>% select("index", "markers")

tmp <- tmp %>% 
  left_join(test) 


# plot summed gene expression totals per cell
(umap <- tmp%>%
  subset(total >=3.5) %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data= data, size = 0.05, alpha = 0.5, color = "grey")+
  geom_point(data = subset(tmp, total >= 3.5 & total < 4), aes(color = total), size = 0.5)+
  geom_point(data = subset(tmp, total >= 4 & total < 5), aes(color = total), size = 0.5)+
  geom_point(data = subset(tmp, total >= 5 & total < 6), aes(color = total), size = 0.5)+
  geom_point(data = subset(tmp, total >= 6 & total < 10), aes(color = total), size = 0.5)+
  geom_point(data = subset(tmp, total >= 10), aes(color = total), size = 0.5)+
  geom_point(data = subset(tmp, total >= 20), aes(color = total), size = 0.5)+
  geom_point(data = subset(tmp, total >= 30), aes(color = total), size = 0.5)+
  geom_point(data = subset(tmp, total >= 40), aes(color = total), size = 0.5)+
  scale_color_viridis()+
  geom_ellipse( aes(x0= -2, y0=3.5, a = 4, b = 2, angle = 15), size = 0.15, inherit.aes = FALSE, color = "red")+
  labs(color = str_wrap("Read Count", width = 5))+
  theme(#text = element_text(family = "helvetica"),
        strip.background = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = ggplot2::element_text(size = 8),
        axis.text.y = ggplot2::element_text(size = 8),
        axis.title.x = ggplot2::element_text(size = 9, vjust = 1),
        axis.title.y = ggplot2::element_text(size = 9, vjust = -2),
        strip.text = ggtext::element_markdown(size = 8), 
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.title=element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.text.align = 0.8,
        legend.position = "right",
        legend.margin = margin(0, 0, 0, -0.2, "cm"),
        plot.margin = margin(0.3, 0, -0.05, -0.1, "cm"))+
  guides(alpha = "none"))



###################
### Fig.3C - # plot for marker quantification per cell from fig3b
###################
library(janitor)

# calculate read total per gene 
cells <-tmp %>% 
  subset(total >= 3.5) %>% 
  tabyl(integrated_snn_res.0.5, index) %>% 
  adorn_totals(("col")) %>% 
  select("integrated_snn_res.0.5", "Total")

# calculate number of markers expressed in each cell (co-expression)
markers <- tmp %>% 
  subset(total >= 3.5) %>% 
  tabyl(integrated_snn_res.0.5, markers) %>% 
  pivot_longer(!integrated_snn_res.0.5, names_to = "markers", values_to = "counts")

test <- cells %>% 
  left_join(markers)

# rename cluster numbers and assign annotation ID 
test$integrated_snn_res.0.5 <- factor(test$integrated_snn_res.0.5, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"))

test <- test %>% 
  mutate(ID = case_when(
    integrated_snn_res.0.5 == "1" ~ "Unannotated",
    integrated_snn_res.0.5 == "2" ~ "MS",
    integrated_snn_res.0.5 == "3" ~ "Unannotated",
    integrated_snn_res.0.5 == "4" ~ "Unannotated",
    integrated_snn_res.0.5 == "5" ~ "Unannotated",
    integrated_snn_res.0.5 == "6" ~ "C",
    integrated_snn_res.0.5 == "7" ~ "Unannotated",
    integrated_snn_res.0.5 == "8" ~ "Unannotated",
    integrated_snn_res.0.5 == "9" ~ "MD",
    integrated_snn_res.0.5 == "10" ~ "Unannotated",
    integrated_snn_res.0.5 == "11" ~ "Neuron",
    integrated_snn_res.0.5 == "12" ~ "Neuron",
    integrated_snn_res.0.5 == "13" ~ "Neuron",
    integrated_snn_res.0.5 == "15" ~ "S",
    integrated_snn_res.0.5 == "14" ~ "CA",
    integrated_snn_res.0.5 == "16" ~ "Unannotated",
    integrated_snn_res.0.5 == "17" ~ "MD",
    integrated_snn_res.0.5 == "18" ~ "Neuron",
    integrated_snn_res.0.5 == "19" ~ "MS",
    integrated_snn_res.0.5 == "20" ~ "Unannotated",
    integrated_snn_res.0.5 == "21" ~ "Unannotated",
    integrated_snn_res.0.5 == "22" ~ "IB",
    integrated_snn_res.0.5 == "23" ~ "Neuron",
    integrated_snn_res.0.5 == "24" ~ "Neuron",
    integrated_snn_res.0.5 == "25" ~ "Neuron",
    integrated_snn_res.0.5 == "26" ~ "Neuron",
    integrated_snn_res.0.5 == "27" ~ "Neuron"))

test$ID <- factor(test$ID, levels = c("MS","MD", "C", "S", "CA", "IB", "Neuron", "Unannotated"))


# plot
(bar<- test %>% 
  ggplot(aes(x =integrated_snn_res.0.5, y = counts, fill = markers))+
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE))+
  scale_fill_viridis(discrete = T, direction = -1)+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Cluster", y = "Cell Count", fill = str_wrap("Markers Per Cell", width = 8))+
  facet_grid(cols = vars(ID), space = "free", scales = "free", drop = TRUE)+
  theme(panel.background = element_blank(),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.key = element_blank(),
        legend.title=element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.text.align = 0.8,
        legend.position = c(0.8, 0.9),
        legend.direction = "horizontal",
        legend.margin = margin(0, 1, 0, -0.2, "cm"),
        legend.box.margin = margin(0, 0.01, 0, -2),
        axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
        axis.title.x = ggplot2::element_text(size = 10, vjust = 0.5),
        axis.title.y = ggplot2::element_text(size = 10), 
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 7.5),
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        plot.margin = margin(0.5, 0.1, 0, 0.3, "cm")))





####################
### Fig. 3D - Cel. BK26 FACS-RNAseq and markers
####################


# Cel. BK26 with GFP(+) cell
cel_GFP <- ggdraw() +
  draw_image(magick::image_read_pdf("Figures/Figure_3/Fig3a_cel_GFP.pdf"))

cel_GFP <- plot_grid(NULL, cel_GFP, ncol = 2, rel_widths = c(0.001, 1))



# Cel FACS-RNAseq Heatmap

# read in counts (raw) and filter for samples
counts.raw <- readRDS("RNAseq/Cel_FACS_RNAseq/data/counts.raw.rds")


# convert to matrix
counts.raw <- counts.raw %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")


# make metadata file
sample_list <- colnames(counts.raw)
samples <- data.frame("sample" = sample_list, stringsAsFactors = FALSE)
samples <- samples %>%
  mutate(cond = c("Large", "NegCells", "Small")) 
samples$cond <- as.factor(samples$cond)
#saveRDS(samples, "data/deg/samples.rds")

# read in sample metadata df > samples
#samples <- readRDS("data/deg/samples.rds")



#Transform and create heatmap for quick gene expression visualization across all samples
keep <- counts.raw %>% 
  mutate(keep = ifelse(rowSums(counts.raw >= 1) >1, "keep", "discard")) %>% 
  subset(keep == "keep")
keep <- rownames_to_column(keep, var = "gene_id")
keep <- keep$gene_id




# read in counts (tpm) and filter for samples that meet the filters for raw
counts.tpm <- readRDS(here("RNAseq/Cel_FACS_RNAseq/data/counts.tpm.rds")) %>%
  filter(gene_id %in% keep)


#join counts and sample info, widen, declare and normalize matrix
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

# for downstream venn diagram
cel_facs <- counts %>% 
  filter(sample == "Large") %>% 
  filter(expression >= 1.15)

# rename samples
counts$sample <- factor(counts$sample, levels = c("Large","Small", "NegCells"), labels = c( "GFP(+) Large", "GFP(+) Small", "GFP(-)"))


(cel_heatmap <- counts %>% 
  ggplot(aes(gene_id,sample)) +
  geom_tile(aes(fill = expression), show.legend = TRUE) +
  scale_fill_viridis() +
  scale_y_discrete(limits = c("GFP(+) Small", "GFP(+) Large","GFP(-)"), position = "left")+
  labs(x = "*C. elegans* genes", fill = "Z-Score")+
  theme(#text = element_text(family = "helvetica"),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_markdown(size = 8, angle = 45, hjust = 1, vjust = 1.1),
    axis.title.y = element_markdown(size = 10, angle = 90, vjust = -20),
    axis.title.x = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.title=element_markdown(size = 8),
    legend.text = element_text(size = 8),
    legend.text.align = 0.8,
    legend.position = "right",
    legend.margin = margin(0,0,0,-0.2, "cm"),
    plot.margin = margin(0.5, 0, 0, 0, "cm"))+
  coord_flip())


# combine all panels for Fig3D
fig3d <- plot_grid(cel_GFP, cel_heatmap, ncol = 2, rel_widths = c(1.05, 0.95), scale = c(0.99, 1), axis= "t")


#################
### Fig.3f - Antigen dotplot
#################
antigen <- read.csv(here("Auxillary/antigens.csv"))
ant_genes <- antigen$gene_id
ant_genes <- ant_genes[!duplicated(ant_genes)]

# use Seurat DotPlot function to calculate average and percet expression for genes for dotplot
dot <- DotPlot(new_combined, features = ant_genes, assay = "RNA", scale = FALSE) # not found: WBGene00220628, WBGene00249804, WBGene00230908
dot <- dot$data 
dot <- rownames_to_column(dot, "genes")
dot <- dot %>% 
  mutate(gene_id = substr(genes, 1, 14)) %>% 
  select(-"genes")

dot <- dot %>% 
  left_join(antigen)

dot[is.na(dot)] <- 0

# rename clusters and assign annotated labels
dot$id <- factor(dot$id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"))

dot <- dot %>% 
  mutate(ID = case_when(
    id == "1" ~ "Unannotated",
    id == "2" ~ "MS",
    id == "3" ~ "Unannotated",
    id == "4" ~ "Unannotated",
    id == "5" ~ "Unannotated",
    id == "6" ~ "C",
    id == "7" ~ "Unannotated",
    id == "8" ~ "Unannotated",
    id == "9" ~ "MD",
    id == "10" ~ "Unannotated",
    id == "11" ~ "Neuron",
    id == "12" ~ "Neuron",
    id == "13" ~ "Neuron",
    id == "15" ~ "S",
    id == "14" ~ "CA",
    id == "16" ~ "Unannotated",
    id == "17" ~ "MD",
    id == "18" ~ "Neuron",
    id == "19" ~ "MS",
    id == "20" ~ "Unannotated",
    id == "21" ~ "Unannotated",
    id == "22" ~ "IB",
    id == "23" ~ "Neuron",
    id == "24" ~ "Neuron",
    id == "25" ~ "Neuron",
    id == "26" ~ "Neuron",
    id == "27" ~ "Neuron"))


dot$ID <- factor(dot$ID, levels = c("MS","MD", "C", "S", "CA", "IB", "Neuron", "Unannotated"))


# create the labels 
labels <- glue_data(
  dot, "<span style='color: {if_else(Moreno_abundant == 1, 'red', 'black')}'>{gene_name}</span>"
)
names(labels) <- dot$gene_name


# dotplot of all antigens scaled with alpha by expression  
(non_zinc <- dot %>%
  filter(zinc_finger == 0) %>% 
  ggplot(aes(y = id, x = gene_name))+
  geom_point(aes(size = pct.exp, color = avg.exp.scaled), show.legend = FALSE)+
  scale_size("Proportion (%)", range = c(-1, 3), breaks=c(0, 10, 25, 50, 75))+
  #scale_size_continuous(range = c(-1, 3), nice.breaks = TRUE)+
  scale_color_viridis_c(limits = c(0, 7))+
  scale_y_discrete(labels = labels)+
  labs(x = "Genes", y = "Cluster", size = "Proportion (%)", color = "Avg. Exp.")+
  facet_grid(cols = vars(ID), space = "free", scales = "free", drop = TRUE)+
  #facet_grid(cols = vars(ID),space = "free", scales = "free", drop = TRUE)+
  theme(panel.background = element_blank(),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        axis.text.x = element_markdown(size = 8, angle = 90 + 1e-09, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7.5, face = "italic", angle=0 + 1e-09, hjust = 1),
        axis.title.x = ggplot2::element_text(size = 10, vjust = -1),
        axis.title.y = ggplot2::element_text(size = 10, vjust = -5), 
        strip.text.x = element_text(size = 7),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.title = element_text(size = 8, vjust =1),
        legend.text = element_text(size = 8),
        legend.text.align = 0.8,
        legend.position = "bottom",
        legend.margin = margin(-0.2, 0.5, 0, -0.2, "cm"),
        plot.margin = margin(0, 0, 1.3, -0.5, "cm"),
        panel.grid = element_line(color = "#e0e0e0", size = 0.05),
        panel.spacing.x = unit(0.25, "lines"))+
     coord_flip()+
     NULL)

legend <- get_legend(non_zinc)

(zinc <- dot %>%
    filter(zinc_finger == 1) %>% 
    ggplot(aes(y = id, x = gene_name))+
    geom_point(aes(size = pct.exp, color = avg.exp.scaled), show.legend = FALSE)+
    scale_size("Proportion (%)", range = c(-1, 3), breaks=c(0, 10, 25, 50, 75))+
    #scale_size_continuous(range = c(-1, 3), nice.breaks = TRUE)+
    scale_color_viridis_c(limits = c(0, 7))+
    scale_y_discrete(labels = labels)+
    labs(x = "Genes", y = "Cluster", size = "Proportion (%)", color = "Avg. Exp.")+
    facet_grid(cols = vars(ID), space = "free", scales = "free", drop = TRUE)+
    #facet_grid(cols = vars(ID),space = "free", scales = "free", drop = TRUE)+
    theme(panel.background = element_blank(),
          axis.line = element_line (colour = "black"),
          legend.background=element_blank(),
          axis.text.x = element_markdown(size = 8, angle = 90 + 1e-09, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 7.5, face = "italic", angle=0 + 1e-09),
          axis.title.x = ggplot2::element_text(size = 10, vjust = -1),
          axis.title.y = ggplot2::element_text(size = 10, vjust = -6, color = "white"), 
          strip.text.x = element_text(size = 7),
          strip.text.y = element_blank(),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.key.height = unit(0.25, "cm"),
          legend.key.width = unit(0.35, "cm"),
          legend.title = element_text(size = 8, vjust =1),
          legend.text = element_text(size = 8),
          legend.text.align = 0.8,
          legend.position = "bottom",
          legend.margin = margin(0.25, 0.5, 0, -0.2, "cm"),
          plot.margin = margin(0, 0, 1.3, -0.5, "cm"),
          panel.grid = element_line(color = "#e0e0e0", size = 0.05),
          panel.spacing.x = unit(0.25, "lines"))+
    coord_flip()+
    NULL)


## heatmap based on antigen target to indicate diagnostic, vaccine, or prominent antigens
purpose_nonzinc <- dot %>% 
  select("gene_name", "Vaccine", "DiagnosticAntigen", "Moreno_abundant") %>% 
 pivot_longer(!gene_name, names_to = "purpose", values_to = "counts") %>% 
  mutate(purpose = case_when(purpose == "DiagnosticAntigen" ~ "Diagnos.",
                   purpose == "Vaccine" ~ "Vacc.",
                   purpose == "Moreno_abundant" ~ "Abundant")) %>% 
  left_join(antigen, by = "gene_name") 
purpose_nonzinc[is.na(purpose_nonzinc)] <- 0

purpose_nonzinc <- purpose_nonzinc %>% 
  filter(zinc_finger == 0) %>% 
  select("gene_name", "purpose", "counts", "zinc_finger")  
  distinct()
 

purpose_nonzinc <- purpose_nonzinc%>% 
  mutate(zinc_finger = case_when(zinc_finger == 1 ~ "Zinc finger TFs (C2H2 type)",
                                 zinc_finger == 0 ~ "Prominent antigens"))


heatmap_nonzinc <- purpose_nonzinc %>% 
  filter(zinc_finger == "Prominent antigens")  %>% 
  ggplot()+
  geom_tile(data = subset(purpose_nonzinc, counts == 1), aes(x = purpose, y = gene_name, fill = counts))+
  geom_tile(data = subset(purpose_nonzinc, counts == 0), aes(x = purpose, y = gene_name, fill = counts), color = "white", alpha = 0.001)+
  scale_fill_gradient2(guide = "none")+
  scale_y_discrete(position= "right")+
  facet_grid(rows = vars(zinc_finger), space = "free", scales = "free", drop = TRUE)+
  theme(panel.background = element_blank(),
        axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(color = "black"),
        plot.margin = margin(0.575, 0, 0.7, 0.01, "cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.grid = element_line(color = "#e0e0e0", size = 0.05))


# combine the dot and heatmap into single panel
(prominent_antigens <- plot_grid(non_zinc, heatmap_nonzinc, nrow = 1, rel_widths = c(1, 0.1), rel_heights = c(1.6, 0.6), align = "v", axis = "tb"))
 

################
### Figure 3G - Dot plot of zinc fingers enriched in secretory cell
################
purpose_zinc <- dot %>% 
  select("gene_name", "Vaccine", "DiagnosticAntigen", "Moreno_abundant") %>% 
  pivot_longer(!gene_name, names_to = "purpose", values_to = "counts") %>% 
  mutate(purpose = case_when(purpose == "DiagnosticAntigen" ~ "Diagnos.",
                             purpose == "Vaccine" ~ "Vacc.",
                             purpose == "Moreno_abundant" ~ "Abundant")) %>% 
  left_join(antigen, by = "gene_name") 
purpose_zinc[is.na(purpose_zinc)] <- 0

purpose_zinc <- purpose_zinc %>% 
  filter(zinc_finger == 1) %>% 
  select("gene_name", "purpose", "counts", "zinc_finger")  
distinct()


(heatmap_zinc <- purpose_zinc %>% 
  ggplot()+
  geom_tile(data = subset(purpose_zinc, counts == 1), aes(x = purpose, y = gene_name, fill = counts))+
  geom_tile(data = subset(purpose_zinc, counts == 0), aes(x = purpose, y = gene_name, fill = counts), color = "white", alpha = 0.001)+
  scale_fill_gradient2(guide = "none")+
  scale_y_discrete(position= "right")+
  facet_grid(rows = vars(zinc_finger), space = "free", scales = "free", drop = TRUE)+
    theme(panel.background = element_blank(),
          axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color = "black"),
          plot.margin = margin(0.575, 0, 0.7, 0.01, "cm"),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          panel.grid = element_line(color = "#e0e0e0", size = 0.05)))


# combine the dot and heatmap into single panel
(zinc_fingers <- plot_grid(zinc, heatmap_zinc, nrow = 1, rel_widths = c(1, 0.1), rel_heights = c(1.6, 0.4)))

  


#########################
### Figure 3E - Venn Diagram of sc, bma-facs, and cel-facs datasets
#########################
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)

# compute average expression for all differentially expressed genes in secretory cluster and grab the top 100 (?) genes
BM_marker_15 <- FindConservedMarkers(new_combined, ident.1 = 14, grouping.var = "orig.ident", verbose = TRUE) %>% 
  filter(utBM_p_val_adj <= 0.05 & utBM_avg_log2FC > 0)
BM_marker_15 <- rownames_to_column(BM_marker_15, var = "gene_id")

bma_genes <- BM_marker_15$gene_id
bma_genes <- bma_genes[!duplicated(bma_genes)]

# use Seurat DotPlot function to calculate average and percent expression for genes for dotplot
bma <- DotPlot(new_combined, features = bma_genes, assay = "RNA", scale = FALSE) 
bma <- bma$data 
bma <- rownames_to_column(bma, "genes")
bma <- bma %>% 
  mutate(gene_id = substr(genes, 1, 14)) %>% 
  select(-"genes")

# filter for genes in the secretory cell and exp >= 1, sort in decending order and take the top 100 genes
bma_tmp <- bma %>% filter(id == 14)

sc <- as.character(bma_tmp$gene_id, stringAsFactors = FALSE)



# replace Cel gene IDs with bma 1:1 orthologs


orthos <- read.csv(here("Auxillary/cel_orthos.csv")) %>%
  select(-"bma_genome_project")
colnames(orthos)[2] <- "gene_id"
gene_list_cel <- cel_facs %>% 
  left_join(orthos)
gene_list_cel <- na.omit(gene_list_cel) %>%  # left with 23 orthologous genes
  select("bma_gene_ID")

colnames(gene_list_cel)[1] <- "gene_id"
cel <- as.character(gene_list_cel$gene_id, stringAsFactors = FALSE)


# create list for venn diagram -- omit the Bma-sc data
x <- list(S1 = bma_facs, S2 = cel)


venn <- Venn(x)
data = process_data(venn, shape_id == "201f")

(venn_diagram <- ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE, color = "black", size = 0.1) +
  xlim(-4, 8)+
  #geom_sf_text(aes(label = name), data = venn_setlabel(data), position = ) +
  geom_sf_text(aes(label = count), data = venn_region(data))+
  annotate(geom = "text", label = "Cel-FACS", x = 6, y = 5.5, size = 3)+
  annotate(geom = "text", label = "(Large GFP(+))", x = 6.3, y = 4.75, size = 2.5)+
  annotate(geom = "text", label = "Bma-FACS", x = 6.2, y = 0, size = 3)+
  annotate(geom = "text", label = "(Largest)", x = 5.5, y = -0.75, size = 2.5)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  theme(legend.position = "none", 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        ))


#################
### Figure 3 - Assemble complete figure
################
row1 <- plot_grid(fig3a, fig3d, ncol = 2, rel_widths = c(1, 0.7), rel_heights = c(1,1),labels = c("A", "D"),hjust = c(0.5, -0.5), label_fontface = "plain", label_fontfamily = "Helvetica", axis = "t", align = "v", label_size = 10) + theme(plot.margin = margin(0, 0, 0, 0.25, "cm"))

row2 <- plot_grid(umap,bar, venn_diagram, ncol = 3, rel_widths = c(1.05, 1.15, 0.85), rel_heights = c(0.8, 1, 1), scale= c(0.95, 1, 0.925), labels = c("B", "C", "E"), vjust = c(1.2, 1.2, 1.2), label_fontface = "plain", label_fontfamily = "Helvetica", label_size = 10, axis = "t", align = "v")

row3 <- plot_grid(prominent_antigens, zinc_fingers, ncol = 2, rel_widths = c(1, 1),labels = c("F", "G"), vjust = c(1.25, 1.25), hjust = c(-0.1, -0.5), label_fontface = "plain", label_fontfamily = "Helvetica", axis = "t", align = "v", label_size = 10)+ theme(plot.margin = margin(0, 0, 0, 0.1, "cm"))
  


# combine rows into final figure 3
Figure3 <- plot_grid(row1, row2, NULL, row3, legend, nrow = 5, rel_heights = c(0.7,0.6, 0.04, 1.5,0.03), align = "v") + theme(plot.margin = margin(0.1, 0.1, 0.1, 0, "cm"))

# export pdf file
ggsave(Figure3, filename = "Figure3.pdf", device = cairo_pdf, width = 9, height = 12, units = "in")





#####################
### Figure 3 - supplemental figure, dot plot with normalized transcript count
####################
library(Rfast)
raw <- as_tibble(new_combined@assays[["RNA"]]@counts, rownames = "gene_id") %>% 
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") # gene expression matrix of raw counts


trans <- data %>% 
  left_join(md) %>% 
  left_join(raw) %>% 
  filter(counts > 0) %>% 
  select("gene_id", "counts", "integrated_snn_res.0.5", "index") %>% 
  filter(gene_id %in% ant_genes)

tmp <- trans %>% 
  group_by(integrated_snn_res.0.5) %>% 
  pivot_wider(names_from = index, values_from = counts)

tmp[is.na(tmp)] <- 0

# calculate the total
tmp$raw_summed <- rowsums(as.matrix(tmp[,c(-1, -2)])) 
new <- tmp %>% select("gene_id", "integrated_snn_res.0.5", "raw_summed")

total <- new %>% 
  group_by(gene_id) %>% 
  summarise(total = sum(raw_summed)) %>% 
  left_join(new) %>% 
  mutate(fraction = ((raw_summed/total)*100)) %>% 
  left_join(antigen)


total$zinc_finger[is.na(total$zinc_finger)] <- 0


total$integrated_snn_res.0.5<- factor(total$integrated_snn_res.0.5, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "2", "3", "4", "5","6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"))

total <- total %>% 
  mutate(ID = case_when(
    integrated_snn_res.0.5 == "1" ~ "Unannotated",
    integrated_snn_res.0.5 == "2" ~ "MS",
    integrated_snn_res.0.5 == "3" ~ "Unannotated",
    integrated_snn_res.0.5 == "4" ~ "Unannotated",
    integrated_snn_res.0.5 == "5" ~ "Unannotated",
    integrated_snn_res.0.5 == "6" ~ "C",
    integrated_snn_res.0.5 == "7" ~ "Unannotated",
    integrated_snn_res.0.5 == "8" ~ "Unannotated",
    integrated_snn_res.0.5 == "9" ~ "MD",
    integrated_snn_res.0.5 == "10" ~ "Unannotated",
    integrated_snn_res.0.5 == "11" ~ "Neuron",
    integrated_snn_res.0.5 == "12" ~ "Neuron",
    integrated_snn_res.0.5 == "13" ~ "Neuron",
    integrated_snn_res.0.5 == "15" ~ "S",
    integrated_snn_res.0.5 == "14" ~ "CA",
    integrated_snn_res.0.5 == "16" ~ "Unannotated",
    integrated_snn_res.0.5 == "17" ~ "MD",
    integrated_snn_res.0.5 == "18" ~ "Neuron",
    integrated_snn_res.0.5 == "19" ~ "MS",
    integrated_snn_res.0.5 == "20" ~ "Unannotated",
    integrated_snn_res.0.5 == "21" ~ "Unannotated",
    integrated_snn_res.0.5 == "22" ~ "IB",
    integrated_snn_res.0.5 == "23" ~ "Neuron",
    integrated_snn_res.0.5 == "24" ~ "Neuron",
    integrated_snn_res.0.5 == "25" ~ "Neuron",
    integrated_snn_res.0.5 == "26" ~ "Neuron",
    integrated_snn_res.0.5 == "27" ~ "Neuron"))


total$ID <- factor(total$ID, levels = c("MS","MD", "C", "S", "CA", "IB", "Neuron", "Unannotated"))

# plot
(supp_nonzinc <- total %>%
    filter(zinc_finger == 0) %>% 
    filter(fraction > 1) %>% 
    ggplot(aes(y = integrated_snn_res.0.5, x = gene_name))+
    geom_point(aes(size = fraction), color = "black", show.legend = FALSE)+
    scale_size("Total reads (%)", range = c(0, 4), breaks = c(1, 5, 10, 25, 50, 75, 100))+
    #scale_color_viridis(discrete = TRUE)+
    #scale_color_manual(values = dakota[3:8])+
    labs(x = "Genes", y = "Cluster", size = "Total reads (%)", title = "Prominent Antigens")+
    facet_grid(cols = vars(ID), rows = vars(zinc_finger), space = "free", scales = "free", drop = TRUE)+
    theme(#text=element_text(family="Helvetica"),
      panel.background = element_blank(),
      axis.line = element_line (colour = "black"),
      legend.background=element_blank(),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8, vjust = 1),
      legend.key = element_blank(),
      axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 8, hjust = 1, face = "italic"),
      axis.title.x = ggplot2::element_text(size = 8, vjust = -1),
      axis.title.y = ggplot2::element_text(size = 8, vjust = -1.5), 
      strip.text.x = element_text(size = 8),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = unit(0.5, "lines"), 
      legend.key.width = unit(0.35, "cm"),
      legend.key.height = unit(0.25, "cm"),
      #legend.key.size = unit(0.25, "cm"), 
      legend.position = "bottom",
      panel.grid = element_line(color = "#ededed", size = 0.05))+
    coord_flip()+
    guides(color= "none"))


(supp_zinc <- total %>%
    filter(zinc_finger == 1) %>% 
    filter(fraction > 1) %>% 
    ggplot(aes(y = integrated_snn_res.0.5, x = gene_name,))+
    geom_point(aes(size = fraction), color = "black", show.legend = FALSE)+
    scale_size("Total reads (%)", range = c(0, 4), breaks = c(1, 5, 10, 25, 50, 75, 100))+
    #scale_color_viridis(discrete = TRUE)+
    #scale_color_manual(values = dakota[3:8])+
    labs(x = "Genes", y = "Cluster", size = "Total reads (%)", title = "Zinc Finger TFs (C2H2)")+
    facet_grid(cols = vars(ID), rows = vars(zinc_finger), space = "free", scales = "free", drop = TRUE)+
    theme(#text=element_text(family="Helvetica"),
      panel.background = element_blank(),
      axis.line = element_line (colour = "black"),
      legend.background=element_blank(),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8, vjust = 1),
      legend.key = element_blank(),
      axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 8, hjust = 1, face = "italic"),
      axis.title.x = ggplot2::element_text(size = 8, vjust = -1),
      axis.title.y = ggplot2::element_text(size = 8, vjust = -3), 
      strip.text.x = element_text(size = 8),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = unit(0.5, "lines"), 
      legend.key.width = unit(0.35, "cm"),
      legend.key.height = unit(0.25, "cm"),
      #legend.key.size = unit(0.25, "cm"), 
      legend.position = "bottom",
      panel.grid = element_line(color = "#ededed", size = 0.05))+
    coord_flip()+
    guides(color= "none"))

legend <- get_legend(supp_zinc)

supp_plots <- plot_grid(supp_nonzinc, supp_zinc, ncol = 2)
supp_plots <- plot_grid(supp_plots, legend, nrow = 2, rel_heights = c(1, 0.1))


ggsave(supp_plots, filename = "antigens_readfraction_percluster.pdf", device = cairo_pdf, width = 10, height = 7, units = "in")




