############
## FIG 4 Expression of anthelmintic targets via UMAP
########################
# this script generates many large objects that may need to be removed from the global environment prior to figure export in order avoid aborting the R session due to memory issues

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


# read in data
new_combined <- readRDS(here("10XGenomics/scmulti_integrated.RDS"))
DefaultAssay(new_combined) <- "RNA"


# pull out the normalized counts, metadata, and UMAP coordinates into dataframes for plotting in ggplot2
data <- as_tibble(new_combined@reductions$umap@cell.embeddings, rownames = 'index') # UMAP coordinates for each cell

md <- as_tibble(new_combined@meta.data, rownames = 'index') # metadata detailing ut/t identity and cluster information

counts <- as_tibble(new_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") 


# color palette (30 total)
dakota <- c("#d97f64", "#263946", "#bebab6", "#7a7f84", "#cab6b2", "#fae2af", "#f3933b","#65838d", "#82aca7", "#a0b4ac", "#b5b9b0", "#fbc1c1", "#e89690", "#d76660", "#cac6b9", "#878787", "#cb8034", "#7f93a2", "#ac8287", "#c1d6d3", "#cd4c42", "#5c8492", "#b25757", "#fe906a", "#6f636b", "#6a9491", "#82ac92", "#a26f6a", "#184459", "#596c7f")




##################
### Fig 4a - Drug target dot plot
##################
# read in csv for anthelmintic targets
drugs <- read.csv(here("Auxillary/drug_targets.csv"))

genes <- drugs$gene_id
genes <- genes[!duplicated(genes)]

# use Seurat dotplot function to calculate average gene expression per cluster
target <- DotPlot(new_combined, features = genes, assay = "RNA", scale = FALSE) 
target <- target$data 
target <- rownames_to_column(target, "genes")
target <- target %>% 
  mutate(gene_id = substr(genes, 1, 14)) %>% 
  select(-"genes")

target <- target %>% 
  left_join(drugs)

# rename clusters
target$id <- factor(target$id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"))


target <- target %>% 
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




target$ID <- factor(target$ID, levels = c("MS","MD", "C", "S", "CA", "IB", "Neuron", "Unannotated"))

target$type <- factor(target$type, levels = c("B-Tubulin", "BK", "GPCR", "LGIC", "GluCl", "ACC", "nAChR", "TRP", "CNG"), labels = c("β-tubulin","BK", "GPCR", "LGIC", "GluCl", "ACC", "nAChR", "TRP", "CNG"))


# dotplot of all targets scaled with alpha by expression  
fig4a <- target %>%
  ggplot(aes(y = id, x = name))+
  geom_point(aes(size = pct.exp, color = avg.exp.scaled ))+
  scale_size("Proportion (%)", range = c(-1, 3))+
  scale_color_viridis(limits = c(0, 3))+
  facet_grid(rows = vars(type),cols = vars(ID), space = "free", scales = "free", drop = TRUE)+
  labs(x = "Gene", y = "Cluster", color = "Avg. Exp.")+
  #scale_color_manual(values = dakota[20:26])+
  theme(text = element_text(family = "helvetica"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8, vjust = 1.05),
        panel.background = element_blank(),
        axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.3),
        axis.text.y = ggplot2::element_text(size = 8, face = "italic"),
        axis.title.x = ggplot2::element_text(size = 10, face = "plain"),
        axis.title.y = ggplot2::element_text(size = 10, vjust = -3, face = "plain"), 
        axis.line = element_line (colour = "black"),
        legend.background = element_blank(),
        strip.text.y = element_text(angle=0, size = 8),
        strip.text.x = element_text(size = 7.5),
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.key.height = unit(0.15, "cm"),
        legend.key.width = unit(0.15, "cm"),
        legend.box.background = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_blank(),
        panel.spacing.x = unit(0.2, "lines"),
        panel.spacing.y = unit(0.2, "lines"),
        panel.grid = element_line(color = "#e0e0e0", size = 0.05), 
        plot.margin = margin(-0.3, 0, 0, 0, "cm"))+
  coord_flip()




#####################
### Fig. 4b - drug target UMAP plots
#####################
table <- data %>% 
  left_join(md) %>% 
  left_join(counts) %>% 
  subset(counts >= 1.5)

# Glutamate-gated chloride channels (GluCls)
#subset data to only glucls
glucl_markers <- c("WBGene00221971", "WBGene00222703", "WBGene00223839", "WBGene00228311")
glucls <- table %>% 
  filter(gene_id %in% glucl_markers)


#plot glucls
(glucls_plot <- glucls %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(aes(color = gene_id), size = 0.25)+
  scale_color_manual(values = c("#f3933b","#65838d","#b25757", "#e60000"), labels = c("*avr-14*", "*glc-4*", "*glc-2*", "*glc-3*"))+
  #scale_size(limits = c(0.5, 2.5), name = "Expression")+
  #annotate(geom = "text", x = 15, y = 0, label = "Macrocyclic lactones", angle = 270, size = 3)+
  labs(title = "Macrocyclic Lactones")+
  theme(text = element_text(family = "helvetica"),
        plot.title = element_markdown(size = 10, hjust = 0.5),
        panel.background = element_blank(),
        axis.text.x = ggplot2::element_text(size = 8,  angle = 0, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
        axis.title.x = ggplot2::element_text(size = 10, color = "black"),
        axis.title.y.left = ggplot2::element_text(size = 10, color = "black"),
        axis.line.y = element_line (colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.background=element_blank(),
        #legend.spacing.x = unit(0.005, "cm"),
        #legend.spacing = unit(0.1, "lines"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.key = element_blank(),
        legend.position = "top",
        legend.margin = margin(-0.2, -0.5, 0, -1, "cm"),
        legend.text = element_markdown(size = 8),
        legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=2)))+
  NULL)

# grab plot legend
glucls_l <- get_legend(glucls_plot)




## Betatubulins 
#create dataframe subset to only the btubs
btub_markers <- c("WBGene00229959", "WBGene00224994", "WBGene00228922","WBGene00233027")
btubs <- table %>% 
  filter(gene_id %in% btub_markers) 

#plot
(btubs_plot <- btubs %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(data = subset(btubs, gene_id == "WBGene00224994" & counts >= 0.5), aes(color = gene_id), size = 0.25)+
  geom_point(data = subset(btubs, gene_id == list("WBGene00228922", "WBGene00233027")), aes(color = gene_id), size = 0.25)+
  geom_point(data = subset(btubs, gene_id == "WBGene00229959"), aes(color = gene_id), size = 0.25)+
  scale_color_manual(values = c("#f3933b","#b25757","#65838d", "#e60000"), labels = c("*btub-1*", "*mec-7*", "*btub-2*", "*tbb-4*"))+
  scale_size(range = c(0.25, 3), name = "Expression")+
  labs(title = "Benzimidazoles")+
    theme(text = element_text(family = "helvetica"),
          plot.title = element_markdown(size = 10, hjust = 0.5),
          panel.background = element_blank(),
          axis.text.x = ggplot2::element_text(size = 8,  angle = 0, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
          axis.title.x = ggplot2::element_text(size = 10, color = "black"),
          axis.title.y.left = ggplot2::element_text(size = 10, color = "black"),
          axis.line.y = element_line (colour = "black"),
          axis.line.x = element_line(colour = "black"),
          legend.background=element_blank(),
          #legend.spacing.x = unit(0.005, "cm"),
          #legend.spacing = unit(0.1, "lines"),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.key = element_blank(),
          legend.position = "top",
          legend.margin = margin(-0.2, -0.5, 0, -1, "cm"),
          legend.text = element_markdown(size = 8),
          legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size=2)))+
    NULL)
 
# grab plot legend
btubs_l <- get_legend(btubs_plot)



# slo-1 (Emodepside target)
slo1 <- table %>%
  filter(gene_id == "WBGene00226980") %>% 
  subset(counts > 2.5)


# plot
(slo1_plot <- slo1 %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(data = subset(slo1, gene_id == "WBGene00226980" & counts >= 0.5), aes(color = gene_id), size = 0.25)+
  scale_color_manual(values = "#65838d", labels = "*slo-1*")+
  scale_size(range = c(0.1, 3), name = "Expression")+
  labs(title = "Emodepside")+
    theme(text = element_text(family = "helvetica"),
          plot.title = element_markdown(size = 10, hjust = 0.5),
          panel.background = element_blank(),
          axis.text.x = ggplot2::element_text(size = 8,  angle = 0, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
          axis.title.x = ggplot2::element_text(size = 10, color = "black"),
          axis.title.y.left = ggplot2::element_text(size = 10, color = "black"),
          axis.line.y = element_line (colour = "black"),
          axis.line.x = element_line(colour = "black"),
          legend.background=element_blank(),
          #legend.spacing.x = unit(0.005, "cm"),
          #legend.spacing = unit(0.1, "lines"),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.key = element_blank(),
          legend.position = "top",
          legend.margin = margin(-0.2, -0.5, 0, -1, "cm"),
          legend.text = element_markdown(size = 8),
          legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(size=2)))+
    NULL)
# grab plot legend
slo1_l <- get_legend(slo1_plot)




# create faceted plot with all three drug target umaps
# dataframe with all three targets
glucls$type <- "glucl"
btubs$type <- "btubs"
slo1$type <- "slo1"

targets <- rbind(glucls, btubs, slo1) %>% 
  mutate(gene_name = case_when( 
    gene_id == "WBGene00221971" ~ "*avr-14*", 
    gene_id == "WBGene00222703" ~ "*glc-4*",
    gene_id == "WBGene00223839" ~ "*glc-2*",
    gene_id == "WBGene00228311" ~ "*glc-3*",
    gene_id == "WBGene00229959" ~ "*btub-2*",
    gene_id == "WBGene00224994" ~ "*btub-1*",
    gene_id == "WBGene00228922" ~ "*mec-7*",
    gene_id == "WBGene00233027" ~ "*tbb-4*",
    gene_id == "WBGene00226980" ~ "*slo-1*"))

targets$type  <- factor(targets$type, levels = c("btubs", "slo1", "glucl"), labels = c("Benzimidazoles", "Emodepside", "Macrocyclic lactones"))
targets$gene_name <- factor(targets$gene_name, levels = c("*avr-14*","*glc-4*", "*glc-2*","*glc-3*","*btub-1*","*mec-7*", "*btub-2*","*tbb-4*", "*slo-1*"), labels = c("*avr-14*","*glc-4*", "*glc-2*","*glc-3*","*btub-1*","*mec-7*", "*btub-2*","*tbb-4*", "*slo-1*"))




#plot faceted plot
targets_plot <- targets %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data[sample(nrow(data), 10000),], size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(data = subset(targets, gene_name == "*btub-1*" & counts > 0.5), aes(color =gene_name), size = 0.25, show.legend = FALSE)+
  geom_point(data = subset(targets, !gene_name == "*btub-2*"), aes(color = gene_name), size = 0.25, show.legend = FALSE)+
  geom_point(data = subset(targets, type == "Macrocyclic lactones"), aes(color = gene_name), size = 0.25, show.legend = FALSE)+
  geom_point(data = subset(targets, gene_id == "WBGene00224994" & counts >= 0.5), aes(color = gene_name), size = 0.25, show.legend = FALSE)+
  geom_point(data = subset(targets, gene_id == list("WBGene00228922", "WBGene00233027")), aes(color = gene_name), size = 0.25, show.legend = FALSE)+
  geom_point(data = subset(targets, gene_id == "WBGene00229959"), aes(color = gene_name), size = 0.25, show.legend = FALSE)+             
  geom_point(data = subset(targets, gene_id == "WBGene00226980" & counts >= 0.5), aes(color = gene_name), size = 0.25, show.legend = FALSE)+            
  scale_color_manual(values = c("#f3933b", "#f3933b","#65838d","#b25757", "#e60000", "#65838d","#b25757","#65838d", "#e60000"))+
  facet_grid(cols = vars(type))+
  theme(#text = element_text(family = "Helvetica"),
        #plot.title = element_markdown(size = 10, hjust = 0.75),
        legend.text = element_markdown(size = 8),
        legend.title = element_blank(),
        panel.background = element_blank(),
        axis.text.x = ggplot2::element_text(size = 8),
        axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
        axis.title.x = ggplot2::element_text(size = 10, color = "black"),
        axis.title.y = ggplot2::element_text(size = 10, color = "black"), 
        axis.line.y = element_line (colour = "black"),
        axis.line.x = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key = element_blank(),
        panel.spacing.y = unit(1.5, "lines"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "cm"))+
  guides(color = guide_legend(override.aes = list(size=1)))



# combine all plot legends
legends <- plot_grid(btubs_l, slo1_l,glucls_l,  ncol = 3) #rel_heights = c(0, 3, 0.25)) + theme(plot.margin = margin(0, 0, 0.5, -0.5, "cm"))


# combine all faceted umaps and legends for Fig 4b
(fig4b <- plot_grid(targets_plot, legends, ncol = 2, rel_widths = c(1, 0.3), rel_heights = c(1, 0), scale = c(1, 0.7))+ theme(plot.margin = margin(-0.4, 0, 3, 0, "cm")))

fig4b <- plot_grid(btubs_plot, slo1_plot, glucls_plot, ncol = 3, labels = c("B", "", ""), label_fontface = "plain", label_fontfamily = "helvetica", label_size = 12, vjust = 0.5)+theme(plot.margin = margin(0, 0.05, 0, 0.05, "cm"))






#######################
## Fig 4c. - Co-expression plots for GluCls, nicotinic subunits and btubs
######################
## building the trees
library(ggtree)
library(tidytree)
library(treeio)
library(ape)
library(dplyr)
library(stringr)


# load in Bma and Cel ids
Bma.id <- read.csv(here("Auxillary/Bma.Proteins.csv",
                   header = FALSE, sep = ","))
colnames(Bma.id) <-  c("protein_id", "gene_id", "gene_name") 
Bma.id <- Bma.id %>% 
  group_by(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)
Bma.protein.list <- unique(Bma.id$protein_id)

Cel.id <- read.csv(here("Auxillary/Cel.Proteins.csv",
                   header = FALSE, sep = ","))
colnames(Cel.id) <-  c("protein_id", "gene_id", "gene_name") 
Cel.id <- Cel.id %>% 
  group_by(gene_id,gene_name) %>%
  distinct(gene_id, .keep_all=TRUE)
Cel.gene.list <- unique(Cel.id$gene_id)
Cel.protein.list <- unique(Cel.id$protein_id)

all_id <- bind_rows(Bma.id, Cel.id) %>% 
  mutate(label = case_when(
    str_detect(protein_id, 'Bm') == TRUE ~ protein_id,
    TRUE ~ gene_id)
  )

# read in iqtree file
lgic.phylo <- read.iqtree(here("Phylogenetics/lgics/tree/LGIC_trim_final.aln.treefile"))

# convert tree (phylo object) to tibble (relabel) and generate d1 with other data
d1 <- as_tibble(lgic.phylo) %>% 
  mutate(species = case_when(
    label %in% Bma.protein.list ~ 'Bma',
    label %in% Cel.gene.list ~ 'Cel'
  )) %>% 
  left_join(., all_id) %>% 
  mutate(tiplab = case_when(
    is.na(gene_name) == TRUE ~ protein_id,
    TRUE ~ gene_name
  ))

#read lgics in from tree (first run code further down to get d1)
lgics <- d1 %>% dplyr::filter(species == "Bma")
lgic.list <- unique(lgics$gene_id)

gene_counts <- counts %>% 
  dplyr::filter(gene_id %in% lgic.list) %>% 
  dplyr::select(gene_id, index, counts) %>%
  pivot_wider(names_from = index, values_from = counts) %>%
  column_to_rownames(var = "gene_id")



# generate matrix (no scale normalization)
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  #df.m <- log2(df.m+1) 
  #df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  return(df.m)
}
gene_counts <- matrix_heatmap(gene_counts)



# correlation matrix on transposed counts (columns = genes)
library(Hmisc)
gene_counts.t <- t(gene_counts)
cor <- rcorr(gene_counts.t, type = c("pearson","spearman"))
cor.r <- as.data.frame(cor$r) %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols=2:52, names_to = "gene_id_2", values_to = "corr") %>%
  transmute(from = pmin(gene_id, gene_id_2), to = pmax(gene_id, gene_id_2), corr) %>%
  distinct()


#reroot
lgic.phylo <- ape::root(lgic.phylo, node = 214)

(node_labels <- ggtree(lgic.phylo, layout = "circular", branch.length = "none") %<+% d1 +
    geom_text2(aes(subset = !isTip, label = node), size = 2, hjust = -.3) +
    geom_text2(aes(subset = !isTip, label = UFboot, color = UFboot), size = 2, hjust = -.3, vjust = 3) +
    geom_tiplab(aes(label = tiplab)) +
    theme_tree2() +
    NULL)


# prepare correlation data
cor_tree <- cor.r %>% 
  left_join(., select(d1, gene_id, tiplab), by = c('from' = 'gene_id')) %>% 
  rename(from_gene_id = tiplab) %>% 
  left_join(., select(d1, gene_id, tiplab), by = c('to' = 'gene_id')) %>%
  rename(to_gene_id = tiplab) %>% 
  select(from = from_gene_id, to = to_gene_id, corr)


# subset tree (glc)
glc.phylo <- tree_subset(lgic.phylo, node = 161, levels_back = 0) 
glc.tibble <- as_tibble(glc.phylo) %>% 
  left_join(., d1, by = 'label') %>% 
  select(!contains('.y')) %>% 
  rename(node = node.x) %>% 
  dplyr::slice(1:10) %>%
  mutate(newlab = str_remove(tiplab, "Bma-"))

glc_cor <- cor_tree %>% 
  filter(from %in% glc.tibble$tiplab & to %in% glc.tibble$tiplab) %>%
  filter(from != to) %>% 
  filter(corr > 0) %>%
  mutate(color = cut(corr, 9, labels = LETTERS[1:9]))

glc.phylo@phylo$tip.label <- glc.tibble$tiplab
glc.phylo@data$SH_aLRT <- 0.1

(glc.tree <- ggtree(glc.phylo, layout = 'inward_circular', branch.length = 'SH_aLRT', xlim = 3) %<+% glc.tibble +
    geom_tippoint(aes(fill = species),
                  shape = 21,
                  size = 2) +
    geom_taxalink(data = glc_cor, mapping = aes(taxa1 = from, taxa2 = to,color = color),
                  size = 1,
                  ncp = 2,
                  offset = 1.1,
                  outward = FALSE,
                  alpha = 0.8) +
    geom_tiplab(aes(label = newlab),
                size = 3,
                align = TRUE,
                linesize = 0,
                linetype = 0,
                offset = -1,
                hjust = 0,
                fontface = 'italic') +
    scale_color_brewer(palette = 'Reds') +
    scale_fill_manual(values = c('black', 'white')) +
    theme(legend.position = "empty")+
    NULL)
#save_plot('glc_tree.pdf', glc.tree, base_height = 12)
#ggsave(glc.tree, filename = "~/Desktop/glc.tree.pdf", device = cairo_pdf, width = 3, height = 3)
glc <- ggdraw() +
  draw_image(magick::image_read_pdf(here("Figures/Figure_4/glc.tree.pdf"))




# subset tree (nachr)

nachr.phylo <- tree_subset(lgic.phylo, node = 241, levels_back = 0) 
nachr.tibble <- as_tibble(nachr.phylo) %>% 
  left_join(., d1, by = 'label') %>% 
  select(!contains('.y')) %>% 
  rename(node = node.x) %>% 
  dplyr::slice(1:30) %>%
  mutate(newlab = str_remove(tiplab, "Bma-"))

nachr_cor <- cor_tree %>% 
  filter(from %in% nachr.tibble$tiplab & to %in% nachr.tibble$tiplab) %>%
  filter(from != to) %>% 
  filter(corr > 0) %>% 
  # remove isoform links
  separate(from, into = c('from', 'from_isoform'), sep = '\\.') %>% 
  separate(to, into = c('to', 'to_isoform'), sep = '\\.') %>% 
  filter(from != to) %>% 
  mutate(to_isoform = case_when(
    is.na(to_isoform) == TRUE ~ '',
    TRUE ~ to_isoform
  )) %>% 
  mutate(from_isoform = case_when(
    is.na(from_isoform) == TRUE ~ '',
    TRUE ~ from_isoform
  )) %>% 
  mutate(
    from = str_c(from, from_isoform, sep = '.'),
    to = str_c(to, to_isoform, sep = '.')
  ) %>% 
  mutate(
    from = str_remove(from, '\\.$'), 
    to = str_remove(to, '\\.$')
  ) %>% 
  select(from, to, corr) %>% 
  mutate(color = cut(corr, 10, labels = LETTERS[1:10])) 

nachr.phylo@phylo$tip.label <- nachr.tibble$tiplab

(nachr.tree <- ggtree(nachr.phylo, layout = 'inward_circular', xlim = 5) %<+% nachr.tibble +
    geom_tippoint(aes(fill = species),
                  shape = 21,
                  size = 2) +
    geom_taxalink(data = nachr_cor, mapping = aes(taxa1 = from, taxa2 = to,color = color),
                  size = 1,
                  ncp = 2,
                  offset = 1.1,
                  outward = FALSE,
                  alpha = 0.8
    ) +
    geom_tiplab(aes(label = newlab),
                size = 3.5,
                align = TRUE,
                linesize = 0,
                linetype = 0,
                offset = -1,
                hjust = 0,
                fontface = 'italic') +
    scale_color_brewer(palette = 'Reds') +
    scale_fill_manual(values = c('black', 'white')) +
    theme(legend.position = "empty")+
    NULL)
#save_plot('nachr.tree.pdf', plot_grid(nachr.tree), base_height = 3)
#ggsave(nachr.tree, filename = "~/Desktop/nachr_tree.pdf", device = cairo_pdf, width = 5, height = 5)
nachr<- ggdraw() +
  draw_image(magick::image_read_pdf("nachr_tree.pdf")) +
  theme(plot.margin = margin(0, 0, 1, 0, "cm"))





### betatubulins
# read in btub tree file
btub.phylo <- read.iqtree(here("Phylogenetics/btubs/tree/btubs_trim_final.aln.treefile"))

# convert tree (phylo object) to tibble (relabel) and generate d1 with other data
d1 <- as_tibble(btub.phylo) %>% 
  mutate(species = case_when(
    label %in% Bma.protein.list ~ 'Bma',
    label %in% Cel.gene.list ~ 'Cel'
  )) %>% 
  left_join(., all_id) %>% 
  mutate(tiplab = case_when(
    is.na(gene_name) == TRUE ~ protein_id,
    TRUE ~ gene_name
  ))

#read lgics in from tree (first run code further down to get d1)
btubs <- d1 %>% dplyr::filter(species == "Bma")
btubs.list <- unique(btubs$gene_id)

gene_counts <- counts %>% 
  dplyr::filter(gene_id %in% btubs.list) %>% 
  dplyr::select(gene_id, index, counts) %>%
  pivot_wider(names_from = index, values_from = counts) %>%
  column_to_rownames(var = "gene_id")


gene_counts <- matrix_heatmap(gene_counts)



# correlation matrix on transposed counts (columns = genes)
library(Hmisc)
gene_counts.t <- t(gene_counts)
cor <- rcorr(gene_counts.t, type = c("pearson","spearman"))
cor.r <- as.data.frame(cor$r) %>%
  rownames_to_column(var = "gene_id") %>% 
  pivot_longer(cols=2:5, names_to = "gene_id_2", values_to = "corr") %>%
  transmute(from = pmin(gene_id, gene_id_2), to = pmax(gene_id, gene_id_2), corr) %>%
  distinct()



#reroot
#btub.phylo <- ape::root(btub.phylo, node = 14)

(node_labels <- ggtree(btub.phylo, layout = "circular", branch.length = "none") %<+% d1 +
    geom_text2(aes(subset = !isTip, label = node), size = 2, hjust = -.3) +
    geom_text2(aes(subset = !isTip, label = UFboot, color = UFboot), size = 2, hjust = -.3, vjust = 3) +
    geom_tiplab(aes(label = tiplab)) +
    theme_tree2() +
    NULL)


# prepare correlation data
cor_tree <- cor.r %>% 
  left_join(., select(d1, gene_id, tiplab), by = c('from' = 'gene_id')) %>% 
  rename(from_gene_id = tiplab) %>% 
  left_join(., select(d1, gene_id, tiplab), by = c('to' = 'gene_id')) %>%
  rename(to_gene_id = tiplab) %>% 
  select(from = from_gene_id, to = to_gene_id, corr)



btub.tibble <- as_tibble(btub.phylo) %>% 
  left_join(., d1, by = 'label') %>% 
  select(!contains('.y')) %>% 
  rename(node = node.x) %>% 
  dplyr::slice(1:10) %>%
  mutate(newlab = str_remove(tiplab, "Bma-")) %>% 
  mutate(newlab = case_when(newlab== "Bm4733" ~ "btub-1",
                            newlab == "Bm9698" ~ "btub-2",
                            TRUE ~ as.character(newlab)))

btub_cor <- cor_tree %>% 
  filter(from %in% btub.tibble$tiplab & to %in% btub.tibble$tiplab) %>%
  filter(from != to) %>% 
  filter(corr > 0) %>%
  mutate(color = cut(corr, 9, labels = LETTERS[1:9]))

btub.phylo@phylo$tip.label <- btub.tibble$tiplab
btub.phylo@data$SH_aLRT <- 0.1

(btub.tree <- ggtree(btub.phylo, layout = 'inward_circular', branch.length = 'SH_aLRT', xlim = 3) %<+% btub.tibble +
    geom_tippoint(aes(fill = species),
                  shape = 21,
                  size = 2) +
    geom_taxalink(data = btub_cor, mapping = aes(taxa1 = from, taxa2 = to,color = color),
                  size = 1,
                  ncp = 2,
                  offset = 1.1,
                  outward = FALSE,
                  alpha = 0.8) +
    geom_tiplab(aes(label = newlab),
                size = 3,
                align = TRUE,
                linesize = 0,
                linetype = 0,
                offset = -1,
                hjust = 0,
                fontface = 'italic') +
    scale_color_brewer(palette = 'Reds') +
    scale_fill_manual(values = c('black', 'white')) +
    theme(legend.position = "empty")+
    NULL)

#ggsave(btub.tree, filename = "~/Desktop/btub.tree.pdf", device = cairo_pdf, width = 3, height = 3)
btub <- ggdraw() +
  draw_image(magick::image_read_pdf(here("Figures/Figure_4/btub.tree.pdf"))





row1 <- plot_grid(fig4a, fig4b, ncol = 2, rel_widths = c(1.2, 0.8), labels = c("A", "B"), rel_heights = c(0.9, 0.5), scale = c(1, 0.9), vjust = -0.2, label_fontface = "plain", label_fontfamily = "helvetica", label_size = 12, align = "v", axis = "t") 

row2_labels <- plot_grid(NULL, NULL, NULL,ncol = 3, labels = c("GluCls", "β-tubulins", "nAChRs"),label_fontface = "plain", label_fontfamily = "helvetica", label_size = 10, hjust = c(-2.5, -0.9, -1.5), vjust = c(1.5, 1.5, -6))

row2 <- plot_grid(NULL, glc, btub, nachr, ncol = 4, labels = c("C","", "", ""), rel_widths = c(0.01, 0.9, 0.9, 1.3), rel_heights = c(1, 1.2, 1.2, 2), scale = c(1, 1.1, 1.1, 1.8), vjust = c(-2,0.75,0.75,0), hjust = c(-1, -2, -1.25, -1.65), label_fontface = "plain", label_fontfamily = "helvetica", label_size = 12) 


Figure_4 <- plot_grid(row1, row2_labels, row2,NULL, nrow=4, rel_heights = c(1.2, 0.05, 0.3, 0.05))+
  theme(plot.margin = unit(c(0.5, 0, 0.1, 0), "cm"))



(coexp <- plot_grid(glc,NULL, btub,NULL, nachr, nrow = 5, rel_heights = c(0.75,0.01, 0.75,0.01, 1), scale = c(0.9,1, 0.9, 1,1), labels = c("GluCl","", "β-tubulin", "","nAChR"), hjust = c(-6.25, 1, -4, 1, -5), vjust = c(1.5, 1, 1.5, 1, 0.75), label_fontface = "plain", label_size = 10))

row1b <- plot_grid(fig4a, coexp, ncol = 2, rel_widths = c(1.3, 0.7), labels = c("A", "C"), scale = c(1, 1), vjust = -0.2, label_fontface = "plain", label_size = 12, align = "v", axis = "t") 


Figure_4 <- plot_grid(row1b,NULL, fig4b, nrow=3, rel_heights = c(1.4,0.01, 0.45))+
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm"))




ggsave(Figure_4, filename =  here("Figures/Figure_4/Figure_4.pdf"), device = cairo_pdf, width = 8, height = 11.5, units = "in")





#################
### Supplemental Figure - Dot plot of total read fraction from raw transcripts for each cluster
#################
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
  left_join(drugs)

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
total$type <- factor(total$type, levels = c("B-Tubulin", "BK", "GPCR", "LGIC", "GluCl", "ACC", "nAChR", "TRP", "CNG"), labels = c("β-tubulin","BK", "GPCR", "LGIC", "GluCl", "ACC", "nAChR", "TRP", "CNG"))


# plot
(supp_plot <- total %>% 
    filter(fraction >= 1) %>% 
    ggplot(aes(y = integrated_snn_res.0.5, x = name))+
    geom_point(aes(size = fraction, color = type))+
    scale_size("Total reads (%)", range = c(0, 4), breaks = c(1, 5, 10, 25, 50, 75, 100))+
    scale_color_manual(values = dakota)+
    labs(x = "Genes", y = "Cluster", size = "Total reads (%)")+
    facet_grid(cols = vars(ID), rows = vars(type), space = "free", scales = "free", drop = TRUE)+
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


ggsave(supp_plot, filename = here("drugtargets_readfraction_percluster.pdf"), device = cairo_pdf, width = 6, height = 10, units = "in")
