#################
# Figure 4-  Signal peptide, transmembrane domain, and C2H2 ZF TF analysis in the secretory cell
#################

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
library(pheatmap)
library(ggtext)
library(ggforce)
library(scales)

#other
setwd("path/to/directory")
library(here)





#############################################
#### Signal Peptide Sequence Predictions ####
#############################################
CAN<- read.csv(here("Figures/Figure_4/signalp/CAN_SP_results.csv"), header = TRUE)
coel <- read.csv(here("Figures/Figure_4/signalp/coelomocyte_SP_results.csv"), header = TRUE)
IB <- read.csv(here("Figures/Figure_4/signalp/IB_SP_results.csv"), header = TRUE)
meso <- read.csv(here("Figures/Figure_4/signalp/mesoderm_SP_results.csv"), header = TRUE)
muscle <- read.csv(here("Figures/Figure_4/signalp/muscle_SP_results.csv"), header = TRUE)
secretory <- read.csv(here("Figures/Figure_4/signalp/secretory_SP_results.csv"), header = TRUE)
neuronal <- read.csv(here("Figures/Figure_4/signalp/neuronal_SP_results.csv"), header = TRUE)

# add column identifying the cell type
CAN$celltype <- "CA"
coel$celltype <- "C"
IB$celltype <- "IB"
meso$celltype <- "MD"
muscle$celltype <- "MS"
secretory$celltype <- "S"
neuronal$celltype <- "N"

# join all SP dataframes into a single dataframe
SPs <- rbind(CAN, coel, IB, meso, muscle, secretory, neuronal)


# pick out the isoforms that have a SP and remove all other isoforms
iso.keep <- SPs %>% 
  subset(str_detect(gene_id, "_")) %>% 
  mutate(keep = ifelse(SP == "NO_SP", "remove", "keep")) %>% 
  filter(keep == "keep") %>% 
  ungroup() %>% 
  select(-"keep") %>% 
  mutate(gene_id = substr(gene_id, 1, 14)) %>% 
  unique()

iso.list <- iso.keep$gene_id

#remove the isoforms from the original list and combine the two lists
final <- SPs %>% 
  mutate(detect = str_detect(SPs$gene_id, "_")) %>% 
  subset(detect == "FALSE") %>% 
  select(-"detect") %>% 
  rbind(., iso.keep) %>% 
  group_by(celltype) %>% 
  unique() 

tmp <- final %>% 
  filter(SP == "NO_SP") %>% 
  subset(gene_id %in% iso.list)
 
SPs <- anti_join(final, tmp) %>% 
  rbind(iso.keep) %>% 
  unique()

 


# some isoforms indicated SP(+) while another isoform indicated SP(-) and are not removed by unique. We want to keep the isoform that indicated SP(+)
SPs <- SPs %>% 
  mutate(Prediction = case_when(
    SP == "SP" ~ "Yes",
    SP == "NO_SP" ~ "No",
    SP == "LIPO" ~ "Yes",
  )) 

# quantify the SP frequency
tmp <- as.data.frame(table(SPs$celltype, SPs$Prediction), stringsAsFactors = FALSE)
colnames(tmp) <- c("celltype", "Prediction", "cnt")
SP <- tmp %>% 
  group_by(celltype) %>% 
  mutate(pct = (cnt/sum(cnt))*100) 
SP$pct <- round(SP$pct, digits = 1)

#factor the celltype labels to  maintain the same theme as the rest of the figures
SP$celltype <- factor(SP$celltype, levels = c("S", "MS","MD","C","CA","IB","N"))

#create new dataframe with n= summary for annotation
label <- SP %>% 
  mutate(table = paste(cnt," (",pct,"%)", sep = "")) %>% 
  group_by(celltype) %>% 
  mutate(total = sum(cnt)) %>% 
  ungroup() %>% 
  mutate(total = paste("n=",total, sep = "")) %>% 
  select("celltype", "total") %>% 
  unique()



#plot
(figa <- SP %>%
    ggplot(aes(x = "", y = pct, fill = Prediction, color = Prediction))+
    geom_bar(stat = "identity", width = 1, show.legend = TRUE)+
    coord_polar(theta = "y", start = 0, clip = "off")+
    scale_y_discrete(expand = c(0,0))+
    labs(fill = "Signal Peptide")+
    facet_grid(cols = vars(celltype))+
    geom_text(data = label, aes(label= total, x="", y = 2), inherit.aes = FALSE, size = 2.75, vjust = -4)+
    scale_fill_manual(values = c("#EE9B00","#005F73"))+
    scale_color_manual(values = c("#EE9B00","#005F73"))+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(margin = margin(b = 7)),
          panel.spacing = unit(-0.1, "cm"),
          legend.key.height = unit(0.25, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 9),
          legend.margin = margin(-0.7,0.2,0,0.3,"cm"),
          legend.position = "right")+
          #strip.text.y.left  = element_text(angle = 0),
          #plot.margin = margin(5.5, 6, 5.5, 5.5, "points"))+
    guides(color = "none"))





############################
# TMs for all major groups #
############################
CAN<- read.csv(here("Figures/Figure_4/hmmtop2.1/CAN_results.csv"), header = TRUE)
coel <- read.csv(here("Figures/Figure_4/hmmtop2.1/coelomocyte_results.csv"), header = TRUE)
IB <- read.csv(here("Figures/Figure_4/hmmtop2.1/IB_results.csv"), header = TRUE)
meso <- read.csv(here("Figures/Figure_4/hmmtop2.1/mesoderm_results.csv"), header = TRUE)
muscle <- read.csv(here("Figures/Figure_4/hmmtop2.1/muscle_results.csv"), header = TRUE)
secretory <- read.csv(here("Figures/Figure_4/hmmtop2.1/secretory_results.csv"), header = TRUE)
neuronal <- read.csv(here("Figures/Figure_4/hmmtop2.1/neuronal_results.csv"), header = TRUE)

# add column identifying the cell type
CAN$celltype <- "CA"
coel$celltype <- "C"
IB$celltype <- "IB"
meso$celltype <- "MD"
muscle$celltype <- "MS"
secretory$celltype <- "S"
neuronal$celltype <- "N"

# join all SP dataframes into a single dataframe
TMs <- rbind(CAN, coel, IB, meso, muscle, secretory, neuronal)


# remove the _# nomenclature for the isoforms in the data and remove duplicated rows
TMs <- TMs %>% mutate(gene_id = substr(gene_id, 1, 14)) %>% 
  unique()

# some isoforms indicated SP(+) while another isoform indicated SP(-) and are not removed by unique. We want to keep the isoform that indicated SP(+)
TMs <- TMs %>% 
  group_by(gene_id, celltype) %>% 
  mutate(domains = max(TMs)) %>% 
  select("gene_id", "domains") %>% 
  unique() %>% 
  ungroup()


#quantify the domain frequency and calculate percent
tmp <- as.data.frame(table(TMs$domains, TMs$celltype), stringsAsFactors = FALSE)
tmp$Var1 <- as.integer(tmp$Var1)
colnames(tmp) <- c("domains", "celltype", "cnt")
TMs <- tmp %>% 
  group_by(celltype) %>% 
  mutate(pct = (cnt/sum(cnt))*100) 
TMs$pct <- round(TMs$pct, digits = 1)

#factor by number of domains and facet
TMs$domains <- factor(TMs$domains, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 19, 23))

TMs$celltype <- factor(TMs$celltype, levels = c("S", "MS","MD","C","CA","IB","N"))


#create new dataframe with n= summary for annotation
label <- TMs %>% 
  mutate(table = paste(cnt," (",pct,"%)", sep = "")) %>% 
  group_by(celltype) %>% 
  mutate(total = sum(cnt)) %>% 
  ungroup() %>% 
  mutate(total = paste("n=",total, sep = "")) %>% 
  select("celltype", "total") %>% 
  unique()







#plot
(figb <- TMs %>%
    ggplot(aes(x = "", y = pct, fill = domains, color = domains))+
    geom_col(color = NA)+
    #geom_bar(stat = "identity", show.legend = TRUE)+
    #coord_polar(theta = "y", start = 0, clip = "off")+
    coord_polar(theta = "y", clip = "off")+
    scale_y_discrete(expand = c(0,0))+
    geom_text(data = label, aes(label= total, x="", y = 2), inherit.aes = FALSE, size = 2.75, vjust = -4)+
    scale_fill_manual(values = c("#EE9B00","#005F73","#E9D8A6", "#607057", "#9B2226","#94D2BD", "#701839", "#063633", "#a194c2", "#9bad86", "#84c4e6", "#7082a4", "#99c96c", "#8d52a5", "#fa4454", "#b82c43", "#000000", "grey"))+
    scale_color_manual(values = c("#EE9B00","#005F73","#E9D8A6", "#607057", "#9B2226","#94D2BD", "#701839", "#063633", "#a194c2", "#9bad86", "#84c4e6", "#7082a4", "#99c96c", "#8d52a5", "#fa4454", "#b82c43", "#000000", "grey"))+
    #geom_text(data = label, aes(label= total, x="", y = 2), inherit.aes = FALSE, size = 2.5, vjust = -2.5)+
    facet_grid(cols = vars(celltype))+
    labs(fill = "    TM Domains")+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.25, "cm"),
          panel.spacing = unit(-0.1, "cm"),
          strip.background = element_blank(),
          strip.text.y.left  = element_text(angle = 0),
          strip.text = element_text(margin = margin(b = 7)))+
          #legend.box.margin = margin(0, -0.5, 0, 0, "cm"),
          #plot.margin = margin(5.5, 6.5, 5.5, 5.5, "points"))+
    guides(color = "none", fill = guide_legend(ncol =3, title.position = "top")))


(figb <- plot_grid(plot_TMs, NULL, TM_tbl, ncol = 3, rel_widths = c(0.425,0.01, 1.575)))




######################################
# C2H2 ZF TF Dot Plot and Phylo Tree #
######################################
# load data
new_combined <- readRDS(here("Figures/Figure_3/10XGenomics/scmulti_integrated.RDS"))
DefaultAssay(new_combined) <- "RNA"


# read in list of orthologous Bma C2H2 TF hits from phylogenetic analysis
list <- read.csv(here("Auxillary/bma_c2h2_wbIDs.csv"), header = FALSE) %>% select("V1") %>% unique()
list <- list$V1

#calculate the expression data for the dot plot
bma <- DotPlot(new_combined, features = list, assay = "RNA", scale = FALSE) 
bma <- bma$data 
bma <- rownames_to_column(bma, "genes")
bma <- bma %>% 
  mutate(gene_id = substr(genes, 1, 14)) %>% 
  select(-"genes")


bma$id <- factor(bma$id, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"))

dot <- bma %>% 
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



#### continue on to phylogenetic tree then combine the two to make the complete tree and paired dot plot



#######################################
### Pylogenetic Tree of C2H2 ZF TFS ###
#######################################

# create phylogenetic tree 
## building the trees
library(ggtree)
library(tidytree)
library(treeio)
library(ape)
library(dplyr)
library(stringr)


# load in Bma and Cel ids
Bma.id <- read.csv(here("Auxillary/Bma.Proteins.csv"),
                   header = FALSE, sep = ",")
colnames(Bma.id) <-  c("protein_id", "gene_id", "gene_name") 
Bma.id <- Bma.id %>% 
  group_by(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)
Bma.protein.list <- unique(Bma.id$protein_id)

Cel.id <- read.csv("Auxillary/Cel.Proteins.csv",
                   header = FALSE, sep = ",")
colnames(Cel.id) <-  c("protein_id", "gene_id", "gene_name") 
Cel.id <- Cel.id %>% 
  group_by(gene_id,gene_name) %>%
  distinct(gene_id, .keep_all=TRUE)
Cel.gene.list <- unique(Cel.id$gene_id)
Cel.protein.list <- unique(Cel.id$protein_id)

all_id <- bind_rows(Bma.id, Cel.id) %>% 
  mutate(label = case_when(
    str_detect(protein_id, "Bm") == TRUE ~ protein_id,
    TRUE ~ gene_id
  ))



# read in iqtree file
c2h2.phylo <- read.iqtree(here("Phylogenetics/C2H2_TFs/tree/c2h2_trim_final.aln.treefile"))


# convert tree (phylo object) to tibble (relabel) and generate d1 with other data
d1 <- as_tibble(c2h2.phylo) %>%
  mutate(label = case_when(
    str_detect(label, "Bm") == TRUE ~ sub("\\..*", "", label),
    TRUE ~ label
  )) %>% 
  mutate(species = case_when(
    str_detect(label, 'Bm') == TRUE ~ "Bma",
    label %in% Cel.gene.list ~ 'Cel'
  )) %>%
  left_join(., all_id) %>%
  mutate(tiplab = case_when(
    is.na(gene_name) == TRUE ~ protein_id,
    TRUE ~ gene_name
  ))


#reroot
(node_labels <- ggtree(c2h2.phylo, layout = "circular", branch.length = "none") %<+% d1 +
    geom_text2(aes(subset = !isTip, label = node), size = 2, hjust = -.3) +
    geom_text2(aes(subset = !isTip, label = UFboot, color = UFboot), size = 2, hjust = -.3, vjust = 3) +
    geom_tiplab(aes(label = tiplab)) +
    theme_tree2() +
    NULL)

c2h2.tibble <- as_tibble(c2h2.phylo) %>% 
  mutate(label = case_when(
    str_detect(label, "Bm") == TRUE ~ sub("\\..*", "", label),
    TRUE ~ label
  )) %>% 
  left_join(., d1) %>% 
  mutate(newlab = str_remove(tiplab, "Bma-"))  %>% 
  unique()

c2h2.tibble <- c2h2.tibble[-284,]
c2h2.tibble <- c2h2.tibble[-504,]



data <- left_join(c2h2.tibble, dot, by = "gene_id", multiple = "all") 
data<- data %>% 
  mutate(avg.exp = case_when(
    is.na(avg.exp) ~ 0,
    TRUE ~ avg.exp)) %>% 
  mutate(pct.exp = case_when(
    is.na(pct.exp) ~ 0,
    TRUE ~ pct.exp)) %>% 
  mutate(features.plot = case_when(
    is.na(features.plot) ~ gene_id,
    TRUE ~ features.plot)) %>% 
  mutate(id = case_when(
    is.na(id) ~ "Cel",
    TRUE ~ id)) %>% 
  mutate(avg.exp.scaled = case_when(
    is.na(avg.exp.scaled) ~ 0,
    TRUE ~ avg.exp.scaled)) %>% 
  mutate(ID = case_when(
    is.na(ID) ~ "Cel",
    TRUE ~ ID))




# tree layout
(p <- ggtree(c2h2.phylo, layout = 'rectangular', branch.length = "none", size = 0.1) %<+% data +
    geom_tiplab(aes(label = newlab, color = species),
                family='Helvetica',
                fontface = "italic",
                geom = "text",
                size = 1) +
    xlim(0, 50)+
    scale_color_manual(values = c('red', 'black')) +
    theme(legend.position = "empty",
          plot.margin = margin(-0.075, 0.45, 1.7, 0, "cm"))+
    NULL)

# get the tree structure
pa <- ggplot_build(p) # the [[3]] dataframe contains the node labels and order

#pull the dataframe and order based on the y column then factor the data object by this node ordering
ordering <- as.data.frame(pa$data[[3]]) %>% 
  select("y", "node")
colnames(ordering)[1] <- "order"

data <- left_join(data, ordering)

data$ID <- factor(data$ID, levels = c("MS", "MD", "C", "S", "CA", "IB", "Neuron", "Unannotated", "Cel"))


(dot_plot <- data %>%   
    ggplot(aes(y = id, x = factor(order)))+
    geom_point(aes(size = pct.exp, color = avg.exp.scaled))+
    scale_size("Proportion (%)", range = c(-1, 3), breaks=c(5, 10, 25, 50, 75))+
    #scale_size_continuous(range = c(-1, 3), nice.breaks = TRUE)+
    scale_color_viridis()+
    labs(x = NULL, y = "Cluster", size = "Proportion (%)", color = "Avg. Exp.")+
    facet_grid(cols = vars(ID), space = "free", scales = "free", drop = TRUE)+
    #scale_x_reverse()+
    #scale_x_continuous(labels = c(as.character(data$newlab)))+
    theme(#text=element_text(family="Helvetica"),
          panel.background = element_blank(),
          axis.line = element_line (colour = "black"),
          legend.background=element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, vjust = 1),
          legend.key = element_blank(),
          axis.text.x = ggplot2::element_text(size = 8, angle = -90, vjust = 0.5),
          axis.text.y = element_blank(),
          axis.title.x = ggplot2::element_text(size = 10, vjust = -1, angle = -180),
          axis.title.y = ggplot2::element_text(size = 10), 
          strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 8, angle = 0),
          strip.background = element_blank(),
          panel.spacing.x = unit(0.25, "lines"), 
          legend.key.width = unit(0.35, "cm"),
          legend.key.height = unit(0.25, "cm"),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          #legend.key.size = unit(0.25, "cm"), 
          legend.position = "bottom",
          plot.margin = margin(0, -0.7, 0, -0.8, "cm"),
          #plot.margin = margin(0.1, -0.59, 0, -0.65, "cm"),
          #panel.grid = element_line(color = "black", size = 0.05))+
          panel.grid = element_line(color = "#ededed", size = 0.05))+
    coord_flip()+
    #guides(size=guide_bins(title= str_wrap("Proportion (%)", width = 13)))+
    NULL)



# combine the dot plot with the phylogenetic tree on the y axis
(tmp <- plot_grid(p, dot_plot, ncol = 2, rel_widths = c(0.2, 1.25), scale = c(0.9528, 1)))


ggsave(tmp, filename = "~/Desktop/phylo_dot_plot.pdf", device = cairo_pdf, width = 6.5, height = 11)

###dot plot and C2H2 tree exported and format refined in illustrator prior to reading it back in for total Figure 4 construction. 




## create figure 
(top <- plot_grid(figa, figb, nrow = 2, labels = c("A", "B"), label_size = 10, label_fontface = "plain"))

# read in the bottom pdf
tmp <-image_read_pdf("~/Desktop/phylo_dot_plot.pdf", density = 600)
bottom <- ggdraw() + draw_image(tmp)


tmp <- plot_grid(top, bottom, nrow = 2, labels = c("", "C"), rel_heights = c(0.65, 1.35), label_size = 10, label_fontface = "plain")

ggsave(tmp, filename = "~/Desktop/Figure4.pdf", device = cairo_pdf, width = 8.5, height = 9, units = "in")



#############################################
### Supplemental GO term analysis figure ####
#############################################

# GO enrichment -----------------------------------------------------------
# load in the text file that has GO IDs matched to gene or transcript IDs
# Note: this must be a text file (can't be a data frame already in memory)
#       with a very specific structure; you can regenerate these files with
#       the commands listed in the first section

library(tidyverse)
library(here)
library(topGO)
library(conflicted)
library(Rgraphviz)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")

#Create GO ID --> gene/transcript ID mappings ----------------------------
  # Note: this only needs to be performed if you have reason to believe that
  #       the mappings have been updated, otherwise you can use the previously
  #       generated files

### WormBase ParaSite species
library(biomaRt)
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# brugia
brugia_go <- getBM(mart = mart, 
                   filters = c("species_id_1010"),
                   value = list("brmalaprjna10729"),
                   attributes = c("wbps_gene_id", "wbps_transcript_id", "go_accession", "go_name_1006", "go_definition_1006", "go_linkage_type", "go_namespace_1003")) 

# Note: we use gene_id for WBP species
brugia_go_out <- dplyr::select(brugia_go, gene_id = wbps_gene_id, go_id = go_accession) %>%
  group_by(gene_id) %>%
  distinct() %>% 
  filter(go_id != "") %>% 
  # mutate(transcript_id = str_remove(transcript_id, '\\.[0-9]*$')) %>%
  summarise(go_ids = list(go_id)) %>% 
  mutate(go_ids = paste0(go_ids)) %>% 
  mutate(go_ids = str_remove_all(go_ids, "c\\(")) %>% 
  mutate(go_ids = str_remove_all(go_ids, '\\"')) %>%
  mutate(go_ids = str_remove_all(go_ids, '\\)'))

write.table(brugia_go_out, '~/Desktop/brugia_gene_go.txt', 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)



# load in the mappings (this must be a file, can't be a df)
gene_go <- topGO::readMappings('~/Library/CloudStorage/Box-Box/brugia_gene_go.txt')

# read the df for later perusal
genes_go <- read_tsv(here('~/Library/CloudStorage/Box-Box//brugia_gene_go.txt'), col_names = c('gene_id', 'go_ids')) %>%
  separate_rows(go_ids, sep = ', ') %>%
  rename(go_id = go_ids)

# get the list of possible gene_ids
gene_ids <- names(gene_go)

# read in your list of genes/transcripts of interest
interest_genes <- ungroup(SPs) %>%  select("gene_id")
interest_go <- left_join(interest_genes, genes_go, by = 'gene_id')

go_summary <- group_by(interest_go, go_id) %>%
  summarize(n = n()) %>%
  filter(!is.na(go_id))

# the final data.frame needs to have one column with all the transcript_ids
# and a second column denoting whether or not it is a transcript of interest
final_genes <- distinct(select(genes_go, gene_id)) %>%
  mutate(interest = case_when(
    gene_id %in% interest_genes$gene_id ~ 1,
    TRUE ~ 0
  )) 

# get the interest column as a factor
final_genes_tg <- as.factor(final_genes$interest)

# convert to a factor with names
names(final_genes_tg) <- final_genes$gene_id

# create the topGOdata objects
# MF == molecular function
# BP == biological process
# CC == cellular component
go_data_mf <- new("topGOdata", ontology = "MF", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)
go_data_bp <- new("topGOdata", ontology = "BP", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)
go_data_cc <- new("topGOdata", ontology = "CC", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)

# create the statistical test
fisher_test <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

all_go_data <- tibble(top_go_object = c(go_data_mf, go_data_bp, go_data_cc)) %>% # make a tibble of all 3 topGOobjects
  mutate(test = map(top_go_object, getSigGroups, fisher_test)) %>% # run the fisher test on each topGOobject
  mutate(result = map2(top_go_object, test, GenTable, ranksOf = "classic", topNodes = 20)) %>% # extract significant GO IDs to a df
  mutate(result = map2(top_go_object, result, ~mutate(.y, class = .x@ontology))) # add the GO class as a column for each nested df

# view the graph for a given class
showSigOfNodes(all_go_data[[1]][[2]], all_go_data[[2]][[2]]@score, firstSigNodes = 5, useInfo = 'all')

plot_data <- select(all_go_data, result) %>%
  unnest(cols = c(result)) %>%
  janitor::clean_names() %>%
  rename(result = result1) %>%
  mutate(class = case_when(
    class == 'BP' ~ 'Biological Process',
    class == 'MF' ~ 'Molecular Function',
    class == 'CC' ~ 'Cellular Component'
  ))

(plot_GO <- plot_data %>% 
    subset(result <= 0.05) %>% 
    #subset(class == "Molecular Function") %>% 
    ggplot() + 
    geom_point(aes(x = term, y = as.numeric(result)), color = "black", size = 1.5) +
    facet_grid(rows = vars(class), scales = "free", space='free') +
    labs(y = "p-value", x = "GO Term") + 
    scale_y_reverse()+
    coord_flip() +
    theme_minimal(base_size = 10, base_family = "Helvetica") +
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 10),
          strip.text = element_text(size =7.5),
          axis.title.y = element_text(vjust = -2, size = 10))
)



ggsave(plot_GO, filename = "~/Desktop/Figure4-figure_supplement1.pdf", device = cairo_pdf, width = 6, height = 6, units = "in")






############### Supplemental Figures ######################3

###############
### Figure 4 - Supplemental Table 1
###############
# table for SP summary
library(mmtable2)

# combine the cnt and pct information in same dataframe cell and create a column with the total SPs for each celltype
SP <- SP %>%
  mutate(table = paste(cnt," (",pct,"%)", sep = "")) %>%
  group_by(celltype) %>%
  mutate(total = sum(cnt)) %>%
  ungroup() %>%
  select("celltype", "Prediction", "table", "total")

# # include "Total" as a "Prediction" element so it represents a summary of the column in the table
tmp <- SP %>%   select("celltype", "total")
tmp$Prediction <- "Total"
colnames(tmp)[2] <- "table"

# #combine the two reduced tables
table <- SP %>%
  select("celltype", "Prediction", "table")
table <- table %>% rbind(table, tmp)

# # factor the "Prediction" column accordingly
table$Prediction <- factor(table$Prediction, levels = c("Yes", "No", "Total"))

# include "Total" as a "Prediction" element so it represents a summary of the column in the table
tmp <- SP %>%   select("celltype", "total")
tmp$Prediction <- "Total"
colnames(tmp)[2] <- "table"

#combine the two reduced tables
table <- SP %>% 
  select("celltype", "Prediction", "table") 
table <- table %>% rbind(table, tmp)

# factor the "Prediction" column accordingly
table$Prediction <- factor(table$Prediction, levels = c("Yes", "No", "Total"))


# create table
(tbl_SP <- table %>% 
    unique() %>% 
    mmtable(cells = table)+
    header_top(celltype)+
    header_left(Prediction)+
    header_format(header = Prediction, style = list(cell_text(weight = "bold")))+
    header_format(header = celltype, style = list(cell_text(weight = "bold", align = "center"))) +
    cells_format(cell_predicate = T, style = list(cell_text(weight= "normal", align = "center"))))



# the table is formatted in a way that cannot be used by cowplot and not easily reformated to a grob. Must be exported as a png and read in for forwatting in cowplot. 

#install.packages("webshot2")
library(webshot2)
gtsave(apply_formats(tbl_SP), filename = "~/Desktop/figure4_supplementaltable1.png")




##################
### Figure 4 - Supplemental Table 2
##################
# table for TM summary

# combine the cnt and pct information in same dataframe cell and create a column with the total TMs for each celltype
TMs <- TMs %>% 
  mutate(table = paste(cnt," (",pct,"%)", sep = "")) %>% 
  group_by(celltype) %>% 
  mutate(total = sum(cnt)) %>% 
  ungroup() %>% 
  select("celltype", "domains", "table", "total")

# include "Total" as a "domain" element so it represents a summary of the column in the table
tmp <- TMs %>%   select("celltype", "total")
tmp$domains <- "Total"
colnames(tmp)[2] <- "table"

#combine the two reduced tables
table <- TMs %>% 
  select("celltype", "domains", "table") 
table <- table %>% rbind(table, tmp)

# factor the "Prediction" column accordingly
table$domains <- factor(table$domains, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "15", "16", "19", "23","Total"))


# create table
(tbl_TM <- table %>% 
    unique() %>% 
    mutate(table = case_when(
      table == "0 (0%)" ~ "0",
      TRUE ~ table)) %>% 
    mmtable(cells = table)+
    header_left(domains)+
    header_top(celltype)+
    header_format(header = celltype, style = list(cell_text(weight = "bold", align = "center")))+
    header_format(header = domains, style = list(cell_text(weight = "bold"))) +
    cells_format(cell_predicate = T, style = list(cell_text(weight= "normal", align = "center"))))

# the table is formatted in a way that cannot be used by cowplot and not easily reformated to a grob. Must be exported as a png and read in for forwatting in cowplot. 

gtsave(apply_formats(tbl_TM), filename = "~/Desktop/figure4_supplementaltable2.png")




