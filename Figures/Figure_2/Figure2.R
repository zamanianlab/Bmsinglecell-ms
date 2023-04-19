############
## FIG 2 B. malayi mf single-cell atlas via 10X Genomics and cluster annotation
############
#data wrangling/analysis
library(tidyverse)
library(Seurat)
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
setwd("path/to/directory")
library(here)



# read in integrated seurat object
new_combined <- readRDS(here("Figures/Figure_2/10XGenomics/scmulti_integrated.RDS"))
DefaultAssay(new_combined) <- "RNA"


# pull out the normalized counts, metadata, and UMAP coordinates into dataframes for plotting in ggplot2
data <- as_tibble(new_combined@reductions$umap@cell.embeddings, rownames = 'index') # UMAP coordinates for each cell

md <- as_tibble(new_combined@meta.data, rownames = 'index') # metadata detailing ut/t identity and cluster information

counts <- as_tibble(new_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") 


# color palette (30 total)
dakota <- c("#cd4c42", "#5c8492", "#b25757", "#fe906a", "#6f636b", "#6a9491", "#82ac92", "#a26f6a", "#184459", "#596c7f", "#d97f64", "#263946", "#bebab6", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#65838d", "#82aca7", "#a0b4ac", "#b5b9b0", "#fbc1c1", "#e89690", "#d76660", "#cac6b9", "#878787", "#cb8034", "#7f93a2", "#ac8287", "#c1d6d3" )



####################
### Fig. 2a - Total utBM UMAP with bulk vs sc RNA-seq inset
####################

# subset data to only the untreated
data2 <- data %>% 
  left_join(counts) %>% 
  left_join(md) %>% 
  subset(counts >= 2.2) %>% 
  subset(orig.ident == "utBM") %>% 
  select("index", "UMAP_1", "UMAP_2","integrated_snn_res.0.5", "orig.ident") %>% 
  distinct()


#Assign identified cluster names
data2$integrated_snn_res.0.5 <- factor(data2$integrated_snn_res.0.5, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "Muscle (2)", "3", "4", "5", "Coelomocyte (6)", "7", "8", "Mesoderm (9)", "10", "Neuron (11)", "Neuron (12)", "Neuron (13)", "Canal-assoc. (14)", "Secretory (15)", "16", "Mesoderm (17)", "Neuron (18)", "Muscle (19)", "20", "21", "Inner body (22)", "Interneuron (23)", "Neuron (24)", "Neuron (25)", "Neuron (26)", "Neuron (27)"))



#table for cluster numbers on umap
clusters <- read.csv(here("Figures/Figure_2/fig2a_labels.csv"))


# plot
global_plot <- data2 %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data2, aes(color = integrated_snn_res.0.5), size = 0.01, show.legend = FALSE)+
  geom_text(data = clusters, aes(x = x, y = y, label = str_wrap(text, width = 8)), size = 3, fontface = "plain")+
  geom_segment(data = clusters, aes(x = xline, y = yline, xend = xend, yend = yend), color = "black", size = 0.5)+
  scale_size_area(max_size = 15)+
  #scale_color_manual(values = dakota)+
  scale_color_manual(values = c("#c1d6d3", "#5c8492", "#b25757", "#6a9491", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#ac8287", "#65838d", "#82aca7", "#fe906a", "#e3e2e1", "#e89690","#cd4c42", "#6f636b", "#82ac92", "#a26f6a", "#184459", "#596c7f","#263946", "#d97f64", "#a0b4ac", "#e3e2e1", "#fbc1c1", "#7f93a2", "#d76660", "#cac6b9", "#e3e2e1", "#cb8034"), labels = function(color) str_wrap(color, width = 8))+
  labs( color = "Cell Type")+
  theme(#text=element_text(family="Helvetica"),
        axis.text.x = ggplot2::element_text(size = 8),
        axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
        axis.title.x = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 10), 
        legend.text = element_markdown(size = 8, face = "plain"),
        legend.key.size = unit(0.2, "cm"),
        panel.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0.5, "cm"),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.25, unit = "cm"))+
  guides(color = guide_legend(override.aes = list(size=3), ncol = 1))+
  NULL



##################
### Figure 2a- inset comparing the sc raw counts to the bulk RNA raw counts for mf (READ RDS AT TOP)
##################
library(Rfast)

# Load the RDS oject instead of running all the code below
rna <- readRDS(here("Figures/Figure_2/bulk_sc_rnaseq_comp.RDS"))


###########################################################
# pull out the raw counts for each gene in each cell
raw<- as_tibble(new_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") 

# create counts matrix with raw count dataframe
#matrix <- raw %>% 
  pivot_wider(names_from = index, values_from = counts)

# sum fraw counts for each gene across all cells 
matrix$total <- rowsums(as.matrix(matrix[,2:46621]))

# extract count data for each gene
sums <- matrix %>% select("gene_id", "total")

#saveRDS(sums, "sc_count_summary.RDS")
#sc_rna <- readRDS("sc_count_summary.RDS")

#For the bulk RNA seq for bma mf:
# pull in the Rda object with counts table
bulk_rna <- readRDS(here("Figures/Figure_2/bulk_rna_tpm.RDS"))


# combine the two rna-seq datasets
rna <- left_join(bulk_rna, sc_rna, by = "gene_id")
#saveRDS(rna, "bulk_sc_rnaseq_comp.RDS")
##############################################################3
  

# trip the dataset of low count/tpm values and NAs in either dataset
rna <- rna %>%
  filter(tpm >= 1) %>% 
  filter(total > 1) %>% 
  na.omit()

# calculate correlation coefficient (r)
cor.test(log10(rna$total), log10(rna$tpm), method = "pearson") # 0.85


#calculate R^2 (coefficient of determination)
summary(lm(log10(rna$tpm) ~ log10(rna$total))) #R^2= 0.72, p = < 2.2e-16


# plot the comparison
plot <- ggplot(rna, aes(log10(total), log10(tpm))) +
  geom_point(size = 0.05)+
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.5)+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(expand = c(0,0))+
  labs(x = expression("SC:" ~Log[10]*"(Total)"), y = expression("Bulk:" ~Log[10]*"(tpm)"))+
  theme(panel.background = element_blank(),
        panel.grid = element_line(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = margin(0, 0, 0, 0, "cm"))


# combined global plot with the RNA comparison inset
global_plot <- global_plot + annotation_custom(ggplotGrob(plot), xmin = -15, xmax = -6, ymin = 5, ymax = 15)





###################################################################################
## Figure 2b - Mapping neuron classes (cholinergic, amingergic, GABAergic, etc.) ##
###################################################################################
# read in csv with neuron info
csv <- read.csv(here("Auxillary/neuron_types.csv"))

genes <- csv$gene_id
genes <- genes[!duplicated(genes)]

#calculate average gene expression per cluster using seurat's DotPlot function
dot <- DotPlot(new_combined, features = genes, assay = "RNA", scale = FALSE)
dot <- dot$data
dot <- rownames_to_column(dot, "genes")
dot <- dot %>% 
  mutate(gene_id = substr(genes, 1, 14)) %>% 
  select(-"genes")

dot <- dot %>% 
  left_join(csv)


#rename clusters
dot$id<- factor(dot$id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "2", "3", "4", "5","6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"))

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

# plot
(row4 <- dot %>% 
    ggplot(aes(y = id, x = gene_name))+
    geom_point(aes(size = pct.exp, color = avg.exp.scaled))+
    scale_size("Proportion (%)", range = c(-1, 3))+
    #scale_size_continuous(range = c(-1, 3), nice.breaks = TRUE)+
    scale_color_viridis()+
    labs(x = "Genes", y = "Cluster", size = "Proportion (%)", color = "Avg. Exp.")+
    facet_grid(cols = vars(ID), rows = vars(neurotransmitter), space = "free", scales = "free", drop = TRUE)+
    theme(text=element_text(family="Helvetica"),
          panel.background = element_blank(),
          axis.line = element_line (colour = "black"),
          legend.background=element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, vjust = 1),
          legend.key = element_blank(),
          axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5),
          axis.text.y = ggplot2::element_text(size = 8, hjust = 1, face = "italic"),
          axis.title.x = ggplot2::element_text(size = 10, vjust = -1),
          axis.title.y = ggplot2::element_text(size = 10), 
          strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 8, angle = 0),
          strip.background = element_blank(),
          panel.spacing.x = unit(0.5, "lines"), 
          legend.key.width = unit(0.35, "cm"),
          legend.key.height = unit(0.25, "cm"),
          #legend.key.size = unit(0.25, "cm"), 
          legend.position = "bottom",
          panel.grid = element_line(color = "#ededed", size = 0.05))+
    coord_flip()+
    #guides(size=guide_bins(title= str_wrap("Proportion (%)", width = 13)))+
    NULL)


#### Combine Fig2A and Fig2B to create Figure 2
(Figure2 <- plot_grid(new_plot, row4, nrow = 2, rel_widths = c(1, 1), rel_heights = c(2.2, 1.7), labels = c("A", "B"), label_fontface = "plain")+theme(plot.margin = margin(0.1, 0, 0, 0, "cm")))

# save plot
ggsave(Figure2, filename = "~/Desktop/Figure2.pdf", device = cairo_pdf, width = 6.5, height = 8, units = "in")


##############################################
########### Supplemental Plots################
##############################################


##################
### Fig2 - figure supplement 1 - histogram of gene and read counts per cell for annotated clusters
##################

# genes expressed per cell, median marked in red vertical line @ 230 genes/cell
genes <- md %>% 
  ggplot()+
  geom_histogram(aes(x = nFeature_RNA), bins = 60)+
  geom_vline(aes(xintercept = median(nFeature_RNA), col = "red"), show.legend = FALSE)+
  scale_x_continuous(breaks = c(0, 50, 250, 500,750, 1000, 1250, 1500, 1750), expand = c(0, 50))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Genes per cell", y = "Count")+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))

# reads captured per cell, median marked in red vertical line @ 267 reads/cell
counts <- md %>% 
  ggplot()+
  geom_histogram(aes(x = nCount_RNA), bins = 60)+
  geom_vline(aes(xintercept = median(nCount_RNA), col = "red"), show.legend = FALSE)+
  labs(x = "Reads per cell", y = "Count")+
  scale_x_continuous(breaks = c(0, 50, 250, 500,750, 1000, 1250, 1500, 1750), expand = c(0, 50))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))

annotated <- plot_grid(genes, counts, ncol = 2)

ggsave(annotated, filename = "figure2-figuresupplement1.pdf", device = cairo_pdf, width = 6, height = 3, units = "in")



##########################
### Fig2 - figure supplement 2 - Orthologous transcription factors and UMAPs of ciliated sensory neuron markers
##########################

## dataframe of Bma-daf-19 expression
daf19 <- counts %>% 
  subset(gene_id == "WBGene00224065") %>% 
  subset(counts > 0) %>% 
  left_join(md) %>% 
  left_join(data)

daf19$cat <- "Bma-daf-19"
daf19$gene_name <- "daf-19"


#### Distribution of genes involved in cilia assembly
cilia <- read.csv(here("Auxillary/cilia_assembly_genes.csv")) %>% 
  select("bma_ortho", "bma_genename", "component") %>% 
  unique()

cilia <- cilia[-4,]
colnames(cilia) <- c("gene_id", "gene_name", "component")

list <- cilia$gene_id

data2<- data %>%
  left_join(md) %>% 
  subset(orig.ident == "utBM") %>% 
  left_join(counts) %>% 
  subset(counts >=1) %>% 
  select("index", "UMAP_1", "UMAP_2","integrated_snn_res.0.5", "orig.ident", "gene_id", "counts") %>% 
  distinct()

# create dataframe and rename the ciliogenesis functions to make shorter
umaps <- data2 %>% 
  filter(gene_id %in% list) %>% 
  subset(orig.ident == "utBM") %>% 
  left_join(cilia) %>% 
  mutate(cat = case_when(
    component == "Kinesin-II" ~ "Kinesin-II",
    component == "IFT-dynein" ~ "IFT",
    str_detect(component, "IFT") ~ "IFT",
    component == "BBS proteins" ~ "BBS",
    component == "Motor activators" ~ "Motor activator",
    component == "Various" ~ "IFT",
    TRUE ~ component)) %>% 
  select(-"component")

daf19 <- daf19 %>% select("index", "UMAP_1", "UMAP_2", "integrated_snn_res.0.5", "orig.ident", "gene_id", "gene_name", "cat", "counts")

# add daf19 data to the ciliogenesis genes
umaps <- rbind(umaps, daf19)

# factor the ciliogenesis genes by function
umaps$cat <- factor(umaps$cat, levels = c("Bma-daf-19", "IFT", "Kinesin-II", "BBS", "Motor activator"), labels = c("DAF-19", "IFT", "Kinesin-II", "BBS", "Motor activator"))



# UMAPs for each gene
figb <- umaps %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data= data2[data2$counts > 4,], size = 0.1, alpha = 0.1, color = "grey")+
  geom_point(data = subset(umaps, counts > 1 & counts <= 1.5), aes(color = counts), size = 0.5, show.legend = TRUE)+
  geom_point(data = subset(umaps, counts > 1.5 & counts <= 3), aes(color = counts), size = 0.5, show.legend = TRUE)+
  geom_point(data = subset(umaps, counts > 3), aes(color = counts), size = 0.5, show.legend = TRUE)+
  scale_color_viridis(limits = c(1,4))+
  facet_grid(rows = vars(cat))+
  labs(color = "Avg. Exp.")+
  theme(axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.placement = "outside",
        axis.line = element_line(),
        legend.background=element_blank(),
        legend.key = element_blank(), 
        legend.position = "bottom",
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(vjust =1.1))+
  NULL



### put it all together
sup_fig<- plot_grid(tf_dot, figb, ncol = 2, rel_widths = c(1.25, 0.75), scale = c(1, 0.95), labels = c("A", "B"), label_fontface = "plain")

# save plot
ggsave(sup_fig, filename = "~/Desktop/figure2_figuresupplement2.pdf", device = cairo_pdf, width = 8.5, height = 10, units = "in")



##########################
### Fig2 - figure supplement 3 - Individual expression UMAPs of markers defining cell tyes
##########################

table <- data %>% 
  left_join(md) %>% 
  left_join(counts) %>% 
  subset(counts > 0)

# Muscle 
muscles <- c("WBGene00231447","WBGene00222011", "WBGene00224604")
muscle <- table %>% 
  filter(gene_id %in% muscles) %>% 
  subset(counts >= 2.5) 


# create a count column and pivot wider to sum across all columns to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
muscle$count <- 1
tmp <- muscle %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/3) %>% 
  select("index", "fraction") %>% 
  mutate(alpha = ifelse(fraction > 0.6, 1, 0.1))


#left join the fraction of markers per cel back onto the muscle dataframe for plotting
muscle <- muscle %>% 
  left_join(tmp) %>% 
  subset(alpha == 1)

muscle$type <- "Muscle"

# plot using the fraction of markers expressed as the color gradient
muscle_plot <- muscle %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
  geom_point(data = subset(muscle, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
  geom_point(data = subset(muscle, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
  geom_point(data = subset(muscle, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
  scale_color_viridis(limits = c(2, 7.5))+
  facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color = "white", vjust = 5),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 11),
        strip.placement = "outside",
        axis.line.x = element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(), 
        legend.title = element_blank(),
        panel.spacing.y = unit(6, "cm"))+
  NULL


# coelomocytes (cluster 5 = 1664 cells)
######## Coelomocytes
coelomocytes <- c("WBGene00268467", "WBGene00227182", "WBGene00223567","WBGene00223869")
coel <- table %>% 
  filter(gene_id %in% coelomocytes) %>% 
  subset(counts >= 2.5)


# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
coel$count <- 1
tmp <- coel %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/4) %>% 
  select("index", "fraction") %>% 
  mutate(alpha = ifelse(fraction >= 0.5, 1, 0.1))


#left join the fraction of markers per cel back onto the muscle dataframe for plotting
coel <- coel %>% 
  left_join(tmp) %>% 
  subset(alphs = 1)


coel$type <- "Coelomocytes"
# plot using the fraction of markers expressed as the color gradient
coel_plot <- coel %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
    geom_point(data = subset(coel, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(coel, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(coel, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    #scale_alpha_manual(guide = FALSE)+
    scale_color_viridis(limits = c(2, 7.5))+
    facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "white", vjust = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 3),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 11),
          strip.placement = "outside",
          axis.line.x = element_blank(),
          #axis.line.y = element_line(color = "black"),
          legend.background=element_blank(),
          legend.key = element_blank(), 
          legend.title = element_blank(),
          panel.spacing.y = unit(6, "cm"))+
    NULL

# Inner Body related (cluster 21 = 503 cells)
IB_marker <- c("WBGene00229597", "WBGene00223435", "WBGene00228562", "WBGene00224494")
IB <- table %>% 
  filter(gene_id %in% IB_marker) %>% 
  subset(counts >= 2.5) 


# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
IB$count <- 1
tmp <- IB %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/4) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction >= 0.5, 1, 0.1))


#left join the fraction of markers per cel back onto the muscle dataframe for plotting
IB <- IB %>% 
  left_join(tmp) %>% 
  subset(alpha == 1)

IB$type <- "Inner body"

# plot using the fraction of markers expressed as the color gradient
IB_plot <- IB %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
    geom_point(data = subset(IB, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(IB, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(IB, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    scale_color_viridis(limits = c(2, 7.5))+
    facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "white", vjust = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 3),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 11),
          strip.placement = "outside",
          axis.line.x = element_blank(),
          legend.background=element_blank(),
          legend.key = element_blank(), 
          legend.title = element_blank(),
          panel.spacing.y = unit(6, "cm"))+
    NULL

# canal-associated (cluster 13 = 1006 cells)
canal_markers <- c("WBGene00222948", "WBGene00226559", "WBGene00226410")
canal<- table %>% 
  filter(gene_id %in% canal_markers) %>% 
  subset(counts >= 2.2)


# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
canal$count <- 1
tmp <- canal %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/3) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction >= 0.6, 1, 0.1))


#left join the fraction of markers per cel back onto the muscle dataframe for plotting
canal <- canal %>% 
  left_join(tmp) 

canal$type <- "Canal-associated"

# plot using the fraction of markers expressed as the color gradient
canal_plot <- canal %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
  geom_point(data = subset(canal, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
  geom_point(data = subset(canal, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
  geom_point(data = subset(canal, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
  scale_color_viridis(limits = c(2, 7.5))+
  facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color = "white", vjust = 5),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 11),
        strip.placement = "outside",
        axis.line.x = element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(), 
        legend.title = element_blank(),
        panel.spacing.y = unit(6, "cm"))+
  NULL

# Mesoderm (clusters 8 and 16, 1455 +  942 = 2397 cells)
mesodermal <- c("WBGene00225364", "WBGene00222162")
mes <- table %>% 
  filter(gene_id %in% mesodermal) %>% 
  subset(counts >= 3) 


# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
mes$count <- 1
tmp <- mes %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/2) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction ==1, 1, 0.1))


#left join the fraction of markers per cel back onto the muscle dataframe for plotting
mes <- mes %>% 
  left_join(tmp)

mes$type <- "Mesoderm"

# plot using the fraction of markers expressed as the color gradient
mes_plot <- mes %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
  geom_point(data = subset(mes, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
  geom_point(data = subset(mes, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
  geom_point(data = subset(mes, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
  #scale_alpha_manual(guide = FALSE)+
  scale_color_viridis(limits = c(2, 7.5))+
  facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color = "white", vjust = 5),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 11),
        strip.placement = "outside",
        axis.line.x = element_blank(),
        #axis.line.y = element_line(color = "black"),
        legend.background=element_blank(),
        legend.key = element_blank(), 
        legend.title = element_blank(),
        panel.spacing.y = unit(6, "cm"))+
  NULL




### Combine non-neuronal markers into single dataframe

nonneuronal <- rbind(muscle, IB, coel, canal, mes) 
nonneuronal$type <- factor(nonneuronal$type, levels = c("Muscle", "Coelomocytes", "Inner body", "Canal-associated", "Mesoderm"))






# pan-neuronal 
pan_markers <- c("WBGene00223147", "WBGene00223381", "WBGene00221982", "WBGene00225764", "WBGene00226594")
pan<- table %>% 
  filter(gene_id %in% pan_markers) %>% 
  subset(counts >=2) 



# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
pan$count <- 1
tmp <- pan %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/4) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction >= 0.4, 1, 0.1))

#left join the fraction of markers per cel back onto the muscle dataframe for plotting
pan <- pan %>% 
  left_join(tmp) 

pan$type <- "Pan-neuronal"

pan_plot <- pan %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
    geom_point(data = subset(pan, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(pan, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(pan, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    scale_color_viridis(limits = c(2, 7.5))+
    facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "white", vjust = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 3),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 11),
          strip.placement = "outside",
          axis.line.x = element_blank(),
          legend.background=element_blank(),
          legend.key = element_blank(), 
          legend.title = element_blank(),
          panel.spacing.y = unit(6, "cm"))+
    NULL


# motor
motor<- table %>% 
  subset(gene_id == "WBGene00223870") %>% 
  subset(counts >=2) 

# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
motor$count <- 1
tmp <- motor %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/1) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction >= 1, 1, 0.1))

#left join the fraction of markers per cel back onto the muscle dataframe for plotting
motor <- motor %>% 
  left_join(tmp)

motor$type <- "Motor"

motor_plot <- motor %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
    geom_point(data = subset(motor, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(motor, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(motor, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    scale_color_viridis(limits = c(2, 7.5))+
    facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "white", vjust = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 3),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 11),
          strip.placement = "outside",
          axis.line.x = element_blank(),
          legend.background=element_blank(),
          legend.key = element_blank(), 
          legend.title = element_blank(),
          panel.spacing.y = unit(6, "cm"))+
    NULL


### DVA (cluster 22 = 503 cells)
# tail

dva<- table %>% 
  subset(gene_id == "WBGene00225297") %>% 
  subset(counts >=4.5) 

# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
dva$count <- 1
tmp <- dva %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/1) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction == 1, 1, 0.1))

#left join the fraction of markers per cel back onto the muscle dataframe for plotting
dva <- dva %>% 
  left_join(tmp)
dva$type <- "Interneuron"

dva_plot <- dva %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
    geom_point(data = subset(dva, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(dva, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(dva, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    scale_color_viridis(limits = c(2, 7.5))+
    facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "white", vjust = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 3),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 11),
          strip.placement = "outside",
          axis.line.x = element_blank(),
          legend.background=element_blank(),
          legend.key = element_blank(), 
          legend.title = element_blank(),
          panel.spacing.y = unit(6, "cm"))+
    NULL



# Neuropeptidergic neurons 
pep_markers <- c("WBGene00223554","WBGene00223147", "WBGene00233246", "WBGene00225067" )
pep<- table %>% 
  filter(gene_id %in% pep_markers) %>% 
  subset(counts >=2.5) 


# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
pep$count <- 1
tmp <- pep %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/4) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction >= 0.5, 1, 0.1))

#left join the fraction of markers per cel back onto the muscle dataframe for plotting
pep <- pep %>% 
  left_join(tmp)
pep$type <- "Neuropeptidergic"

pep_plot <- pep %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
    geom_point(data = subset(pep, counts <= 3), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(pep, counts <= 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    geom_point(data = subset(pep, counts > 5), aes(color = counts), size = 0.5, show.legend = FALSE)+
    scale_color_viridis(limits = c(2, 7.5))+
    facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "white", vjust = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 3),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 11),
          strip.placement = "outside",
          axis.line.x = element_blank(),
          legend.background=element_blank(),
          legend.key = element_blank(), 
          legend.title = element_blank(),
          panel.spacing.y = unit(6, "cm"))+
    NULL

# Aminergic neurons 

amine<- table %>% 
  subset(gene_id == "WBGene00225236") %>% 
  subset(counts >=2) 


# create a count column and pivot wider to sum across all columsn to show how many genes are expressed per cell and then calculate the fraction of markers expressed in each cell
amine$count <- 1
tmp <- amine %>% 
  select(gene_id, index, count) %>% 
  pivot_wider(names_from = "gene_id", values_from = "count") 
tmp[is.na(tmp)] <- 0
tmp$total <- rowSums(tmp[,-1]) 

tmp <- tmp %>% mutate(fraction = total/1) %>% 
  select("index", "fraction")%>% 
  mutate(alpha = ifelse(fraction ==1, 1, 0.1))

#left join the fraction of markers per cel back onto the muscle dataframe for plotting
amine <- amine %>% 
  left_join(tmp)

amine$type <- "Aminergic"

amine_plot <- amine %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(data= data2, size = 0.1, alpha = 0.1, color = "grey")+
    geom_point(data = subset(amine, counts <= 3), aes(color = counts), size = 0.5, show.legend = TRUE)+
    geom_point(data = subset(amine, counts <= 5), aes(color = counts), size = 0.5, show.legend = TRUE)+
    geom_point(data = subset(amine, counts > 5), aes(color = counts), size = 0.5, show.legend = TRUE)+
  labs(color = "Norm. Counts")+
    scale_color_viridis(limits = c(2, 7.5))+
    facet_grid(cols = vars(gene_id), rows = vars(type), switch = "y")+
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, color = "white", vjust = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 3),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 11),
          strip.placement = "outside",
          axis.line.x = element_blank(),
          legend.background=element_blank(),
          legend.key = element_blank(), 
          legend.title = element_text(size = 9, vjust = 0.9),
          legend.position = "bottom",
          panel.spacing.y = unit(6, "cm"),
          legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.5, "cm"))+
    NULL


legend <- get_legend(amine_plot)  # then go back and remove legend from amine plot (show.legend = FALSE)


### Combining neuron markers into a single UMAP
neuronal <- rbind(motor, dva, pan, pep, amine)

neuronal$type <- factor(neuronal$type, levels = c("Pan-neuronal", "Motor", "Interneuron", "Neuropeptidergic", "Aminergic"))


# combine the plots to make complete supplemental figure
row1<- plot_grid(muscle_plot, NULL,legend, ncol = 3, rel_widths = c(0.6, 0.1, 0.3))
row2 <- plot_grid(coel_plot, NULL, ncol = 2, rel_widths = c(0.8, 0.2))
row3 <- plot_grid(IB_plot, NULL, ncol = 2, rel_widths = c(0.8, 0.2))
row4 <- plot_grid(canal_plot, mes_plot, ncol = 2, rel_widths = c(0.6, 0.4))
row5 <- plot_grid(pan_plot, NULL, ncol = 2, rel_widths = c(0.99, 0.01))
row6 <- plot_grid(motor_plot, dva_plot, NULL, ncol = 3, rel_widths = c(0.25, 0.25, 0.5))
row7 <- plot_grid(pep_plot, amine_plot, NULL, ncol = 3, rel_widths = c(0.6, 0.25, 0.15))


UMAPS <- plot_grid(row1, row2, row3, row4, row5, row6, row7, nrow = 7)

# save plot
ggsave(UMAPS, filename = "marker_umaps.pdf", device = cairo_pdf, width = 8, height = 12, units = "in")









###################
### Create Figure 2
###################
Figure2 <- plot_grid(fig2a, row2, row3, row4, nrow = 4, rel_widths = c(1.5, 1, 1, 1), rel_heights = c(2.2, 0.9, 1.1, 1.9), labels = c("A", "B","", "C"), label_fontfamily = "helvetica", label_fontface = "plain")+theme(plot.margin = margin(0.1, 0, 0, 0, "cm"))


ggsave(Figure2, filename = "Figure2.pdf", device = cairo_pdf, width = 8, height = 12, units = "in")



###################
### Supplemenetal Figure - alternate dot plot of normalized transcription per cluster for neuron classes
##################
library(Rfast)
raw <- as_tibble(new_combined@assays[["RNA"]]@counts, rownames = "gene_id") %>% 
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") # gene expression matrix of raw counts
  

trans <- data %>% 
  left_join(md) %>% 
  left_join(raw) %>% 
  filter(counts > 0) %>% 
  select("gene_id", "counts", "integrated_snn_res.0.5", "index") %>% 
  filter(gene_id %in% genes)

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
  left_join(csv)

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
(supp_plot <- total %>% 
    filter(fraction > 1) %>% 
    ggplot(aes(y = integrated_snn_res.0.5, x = gene_name))+
    geom_point(aes(size = fraction, color = neurotransmitter))+
    scale_size("Total reads (%)", range = c(0, 4), breaks = c(1, 5, 10, 25, 50, 75, 100))+
    scale_color_manual(values = dakota[3:8])+
    labs(x = "Genes", y = "Cluster", size = "Total reads (%)")+
    facet_grid(cols = vars(ID), rows = vars(neurotransmitter), space = "free", scales = "free", drop = TRUE)+
    theme(#text=element_text(family="Helvetica"),
          panel.background = element_blank(),
          axis.line = element_line (colour = "black"),
          legend.background=element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, vjust = 1),
          legend.key = element_blank(),
          axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5),
          axis.text.y = ggplot2::element_text(size = 8, hjust = 1, face = "italic"),
          axis.title.x = ggplot2::element_text(size = 10, vjust = -1),
          axis.title.y = ggplot2::element_text(size = 10), 
          strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 8, angle = 0),
          strip.background = element_blank(),
          panel.spacing.x = unit(0.5, "lines"), 
          legend.key.width = unit(0.35, "cm"),
          legend.key.height = unit(0.25, "cm"),
          #legend.key.size = unit(0.25, "cm"), 
          legend.position = "bottom",
          panel.grid = element_line(color = "#ededed", size = 0.05))+
    coord_flip()+
    guides(color= "none"))


ggsave(supp_plot, filename = "neurons_readfraction_percluster.pdf", device = cairo_pdf, width = 6, height = 5, units = "in")




##################
### Figure 2 Supplemental marker UMAPs
##################

row1<- plot_grid(muscle_plot, NULL,legend, ncol = 3, rel_widths = c(0.6, 0.1, 0.3))
row2 <- plot_grid(coel_plot, NULL, ncol = 2, rel_widths = c(0.8, 0.2))
row3 <- plot_grid(IB_plot, NULL, ncol = 2, rel_widths = c(0.8, 0.2))
row4 <- plot_grid(canal_plot, mes_plot, ncol = 2, rel_widths = c(0.6, 0.4))
row5 <- plot_grid(pan_plot, NULL, ncol = 2, rel_widths = c(0.99, 0.01))
row6 <- plot_grid(motor_plot, dva_plot, NULL, ncol = 3, rel_widths = c(0.25, 0.25, 0.5))
row7 <- plot_grid(pep_plot, amine_plot, NULL, ncol = 3, rel_widths = c(0.6, 0.25, 0.15))


UMAPS <- plot_grid(row1, row2, row3, row4, row5, row6, row7, nrow = 7)

ggsave(UMAPS, filename = "marker_umaps.pdf", device = cairo_pdf, width = 8, height = 12, units = "in")



##################
### Figure 2 Supplementals - histogram of gene and read counts per cell for annotated clusters
##################
# genes expressed per cell, median marked in red vertical line @ 230 genes/cell
genes <- md %>% 
  ggplot()+
  geom_histogram(aes(x = nFeature_RNA), bins = 60)+
  geom_vline(aes(xintercept = median(nFeature_RNA), col = "red"), show.legend = FALSE)+
  scale_x_continuous(breaks = c(0, 50, 250, 500,750, 1000, 1250, 1500, 1750), expand = c(0, 50))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Genes per cell", y = "Count")+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))

# reads captured per cell, median marked in red vertical line @ 267 reads/cell
counts <- md %>% 
  ggplot()+
  geom_histogram(aes(x = nCount_RNA), bins = 60)+
  geom_vline(aes(xintercept = median(nCount_RNA), col = "red"), show.legend = FALSE)+
  labs(x = "Reads per cell", y = "Count")+
  scale_x_continuous(breaks = c(0, 50, 250, 500,750, 1000, 1250, 1500, 1750), expand = c(0, 50))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))

annotated <- plot_grid(genes, counts, ncol = 2)

ggsave(annotated, filename = "annotated_summary.pdf", device = cairo_pdf, width = 6, height = 3, units = "in")




#########################
### Fig2 - figure supplement 4 - Pseudobulk analysis
#########################
### Pseudobulk analysis of B. malayi mf scRNAseq data
# read in processed RDS object for pheatmap plotting. RDS object represents all processing between the ## marks below
centered <- readRDS("~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/other/bma_mf_pseudobulk.RDS")

################################################################
# aggregate the raw read counts for utBM and tBM clusters
split <- SplitObject(new_combined, split.by = "orig.ident")

utBM_agg <- AggregateExpression(split$utBM, assays = "RNA", slot = "counts", verbose = TRUE, )
utBM_agg <- as.data.frame(utBM_agg$RNA)

tBM_agg <- AggregateExpression(split$tBM, assays = "RNA", slot = "counts", verbose = TRUE, )
tBM_agg <- as.data.frame(tBM_agg$RNA)


#rename columns so utBM and tBM clusters are distinguishable
colnames(utBM_agg) <- paste(colnames(utBM_agg), "utBM", sep= "_")
colnames(tBM_agg) <- paste(colnames(tBM_agg), "tBM", sep= "_")

# rename first column in order to combine the two dataframes
utBM_agg <- rownames_to_column(utBM_agg, var = "gene_id")
tBM_agg <- rownames_to_column(tBM_agg, var = "gene_id")


# leftjoin the two dataframes 
expression <- utBM_agg %>%  left_join(tBM_agg)


### TMM normalization using edgeR
#turn the dataframe back to a matrix
rownames(joined) <- joined[,1]
joined <- joined %>% select(-"gene_id")

# create a vector with the sample names
samples <- data.frame(samples = colnames(joined))
samples <-factor(levels = colnames(joined))

# condition
condition <- samples %>% 
  mutate(condition = ifelse(grepl('utBM', samples), "Control", "Treated")) %>% 
  column_to_rownames(var = "samples")
condition <- factor(condition$condition)       


# build the DGEList object
dge <- DGEList(counts = joined, group=condition, samples = samples)

#normalize the libraries using the default trimmed mean of M values (TMM)
dge <- calcNormFactors(dge)
tmm <- cpm(dge)  


# remove genes that have < 10 counts
#filter 1 (require x reads total across all samples)
keep <- rowSums(tmm) > 10
tmm <- tmm[keep,]  

# add 1 to each count for log transformation and median-centering
tmm <- tmm + 1

# log2 transform
log <- log2(tmm)

#median-center the transformed data
centered <- t(apply(log,1,function(
    x){x-median(x)
}))

# change colnmaes to match the renamed clusters
colnames(centered) <- c("1_utBM", "2_utBM",  "3_utBM",  "4_utBM",  "5_utBM",  "6_utBM" ,"7_utBM", "8_utBM", "9_utBM",  "10_utBM",  "11_utBM", "12_utBM", "13_utBM", "14_utBM", "15_utBM", "16_utBM", "17_utBM", "18_utBM", "19_utBM", "20_utBM", "21_utBM", "22_utBM", "23_utBM", "24_utBM", "25_utBM", "26_utBM", "27_utBM", "1_tBM", "2_tBM",  "3_tBM",  "4_tBM", "5_tBM", "6_tBM",  "7_tBM", "8_tBM", "9_tBM", "10_tBM", "11_tBM", "12_tBM",  "13_tBM", "14_tBM", "15_tBM", "16_tBM" , "17_tBM",  "18_tBM",  "19_tBM" , "20_tBM", "21_tBM", "22_tBM", "23_tBM", "24_tBM", "25_tBM",  "26_tBM", "27_tBM")

#saveRDS(centered, "~/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/other/bma_mf_pseudobulk.RDS")
###################################################



### Heirarchically cluster the genes and "samples"
library(ClassDiscovery)
# use the uncentered pearson correlation using the classdiscovery package
rows <- distanceMatrix(as.matrix(t(centered)), "uncentered correlation")
rowclus <- hclust(rows, method = "complete") # cluster the genes

cols <- distanceMatrix(as.matrix(centered), "uncentered correlation")
colclus <- hclust(cols, method = "complete") #cluster the samples


## generate heatmap of pseudobulk expression with axis labels
rg <- max(abs(centered))

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")

heatmap <- pheatmap(centered, cluster_rows = rowclus, cluster_cols = colclus, show_rownames = FALSE, color = colorRampPalette(c("deepskyblue", "black", "yellow"))(40), breaks = seq(-rg, rg, length.out = 40))

setHook("grid.newpage", NULL, "replace")
#add x and y axis labels
library(grid)
grid.text("Pseudobulk sample", y = -0.025, gp=gpar(fontsize=13))
grid.text("Genes", x=-0.025, rot=90, gp=gpar(fontsize=13))


# save plot
ggsave(filename = "~/Desktop/figure2_figuresupplement4.pdf", device = cairo_pdf, width = 7, height = 6.52, units = "in")



###################
### Figure 2 - figure supplement 5 - alternate dot plot of normalized transcription per cluster for neuron classes
##################
library(Rfast)
raw <- as_tibble(new_combined@assays[["RNA"]]@counts, rownames = "gene_id") %>% 
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") # gene expression matrix of raw counts


trans <- data %>% 
  left_join(md) %>% 
  left_join(raw) %>% 
  filter(counts > 0) %>% 
  select("gene_id", "counts", "integrated_snn_res.0.5", "index") %>% 
  filter(gene_id %in% genes)

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
  left_join(csv)

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
(supp_plot <- total %>% 
    filter(fraction > 1) %>% 
    ggplot(aes(y = integrated_snn_res.0.5, x = gene_name))+
    geom_point(aes(size = fraction, color = neurotransmitter))+
    scale_size("Total reads (%)", range = c(0, 4), breaks = c(1, 5, 10, 25, 50, 75, 100))+
    scale_color_manual(values = dakota[3:8])+
    labs(x = "Genes", y = "Cluster", size = "Total reads (%)")+
    facet_grid(cols = vars(ID), rows = vars(neurotransmitter), space = "free", scales = "free", drop = TRUE)+
    theme(#text=element_text(family="Helvetica"),
      panel.background = element_blank(),
      axis.line = element_line (colour = "black"),
      legend.background=element_blank(),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8, vjust = 1),
      legend.key = element_blank(),
      axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 8, hjust = 1, face = "italic"),
      axis.title.x = ggplot2::element_text(size = 10, vjust = -1),
      axis.title.y = ggplot2::element_text(size = 10), 
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8, angle = 0),
      strip.background = element_blank(),
      panel.spacing.x = unit(0.5, "lines"), 
      legend.key.width = unit(0.35, "cm"),
      legend.key.height = unit(0.25, "cm"),
      #legend.key.size = unit(0.25, "cm"), 
      legend.position = "bottom",
      panel.grid = element_line(color = "#ededed", size = 0.05))+
    coord_flip()+
    guides(color= "none"))


ggsave(supp_plot, filename = "~/Desktop/Figure2-figuresupplement5.pdf", device = cairo_pdf, width = 6, height = 5, units = "in")
