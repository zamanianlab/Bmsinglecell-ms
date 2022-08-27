############
## FIG 5 Flow cytometry analysis for drug dose responses
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
library(ggforce)

#other
library(here)




#################
### Fig. 5a - Viability workflow (small graphic)
#################
Fig5a_diagram <- ggdraw() +
  draw_image(here("Figures/Figure_5/Fig5a_workflow.png"))



#################
### Fig. 5a - spectralflo contour plots of viability
#################
Fig5a_flow <- ggdraw() +
  draw_image(magick::image_read_pdf(here("Figures/Figure_5/viability_flowcytometry.pdf"))




################
### Fig. 5b - flow cytometry, cell viability over time
################

via <- read.csv(here("Flow/viability/20220203_viability_curve.csv"))

# calculate % viability
via <- via %>% 
  mutate(norm = (viable/total)) 

# plot
(p1<-via %>% 
    ggplot(aes(x = incubation, y =  norm, color = treatment))+
    geom_point(data = subset(via, date == 20220203), size = 1.25)+
    geom_line(data = subset(via, date == 20220203),size = 0.5)+
    scale_color_manual(labels = c("DMSO (0.1%)", "Methanol Fixed", "Untreated"), values = c("#DD8D29", "#a2b246", "#65838d"))+
    scale_y_continuous(labels = scales::percent, breaks = c(0,0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8), limits = c(0, 0.8), expand = c(0,0))+
    labs(x = "Timepoint (hr)", y = "Viability", color = "Treatment") +
    theme(#text = element_text(family = "helvetica"),
          axis.text.x = ggplot2::element_text(size = 8, hjust = 0.5,face = "plain"),
          axis.text.y = ggplot2::element_text(size = 8, hjust = 1,face = "plain"),
          axis.title.x = ggplot2::element_text(size = 8, hjust = 1, vjust = -1,face = "plain"),
          axis.title.y = ggplot2::element_text(size = 8,face = "plain"), 
          axis.line = ggplot2::element_line(size = 0.25, colour = "black"),
          axis.ticks = ggplot2::element_line(size = 0.25), 
          strip.text = element_text(size = 8,face = "plain"), 
          legend.text = element_text(size = 8,face = "plain"),
          legend.title = element_text(size = 8,face = "plain"),
          legend.key = element_blank(),
          legend.key.size = unit(0.35, "cm"),
          legend.margin = margin(0, 0, 0, 4, "cm"),
          panel.grid.major = element_line(color = "#ededed", size = 0.25),
          strip.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank())+
    NULL)

# plot violin plot of several DMSO controls
(p2 <- via %>% 
    ggplot()+
    geom_violin(data = subset(via, incubation == 2 & sample == "DMSO"),size = 0.5, aes(x = incubation, y = norm, color = treatment), show.legend = FALSE)+
    geom_point(data = subset(via, incubation == 2 & sample == "DMSO"), aes(x = incubation, y = norm, color = treatment), size = 1.25, show.legend = FALSE)+
    scale_color_manual(labels = "DMSO (0.1%)", values = "#DD8D29")+
    scale_y_continuous(labels = scales::percent, breaks = c(0,0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8), limits = c(0, 0.8), expand = c(0,0))+
    scale_x_discrete(limits = 2)+
    #theme_minimal(base_family = "Helvetica")+
    theme(axis.text.x = ggplot2::element_text(size = 8, hjust = 0.5,face = "plain"),
          axis.text.y = element_blank(),
          axis.title.x = ggplot2::element_text(size = 8,face = "plain", color = "white"), 
          axis.title.y = ggplot2::element_text(size = 8,face = "plain", color = "white"), 
          axis.line.x = ggplot2::element_line(size = 0.25, colour = "black"),
          axis.ticks.x = ggplot2::element_line(size = 0.25), 
          axis.ticks.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, -0.5, "cm"),
          panel.grid.major = element_line(color = "#ededed", size = 0.25),
          strip.background = element_blank(),
          panel.background = element_blank())+
    NULL)



viability<- plot_grid(p1 + theme(legend.position = "none"), p2, ncol = 2, rel_widths = c(1, 0.4), rel_heights = c(1, 1))+ theme(plot.margin=margin(0, 0, 0, 0, "cm"))

legend <- get_legend(p1 + theme(legend.margin = margin(0,0,0,4)))

(Fig5b <- plot_grid(viability, legend, ncol = 2, rel_widths = c(1, .2), scale = c(1, 0.1), axis = "t", align = "hv") + theme(plot.margin = margin(0, 1, 0, 0, "cm")))






################
### Fig. 5c - drug dose response curves with IVM, LEV, AZS
################

#read in drug dose response data files for IVM, AZS, and LEV
# rep1
ddr_1 <- read.csv(here("Flow/anthelmintics/20211111_bma_mf_ddrc_rep1.csv"))

ddr_1 <- ddr_1 %>% 
  mutate(treatment = case_when(treatment == "ABS" ~ "AZS",
                               TRUE ~ as.character(treatment))) 

ddr_1 <- ddr_1 %>% 
  mutate(treatment = case_when(sample == "DMSO_1" ~ "IVM",
                               sample == "DMSO_2" ~ "AZS",
                               sample == "DMSO_3" ~ "LEV",
                               TRUE ~ as.character(treatment))) %>% 
  mutate(sample = case_when(sample == "DMSO_1" ~ "DMSO",
                            sample == "DMSO_2" ~ "DMSO",
                            sample == "DMSO_3" ~ "DMSO",
                            TRUE ~ as.character(sample)))

ddr_1 <- ddr_1[,-6]
ddr_1 <- ddr_1[-c(5,9,14),]
ddr_1$date <- "20211111"


# rep2
ddr_2 <- read.csv(here("Flow/anthelmintics/20220203_bma_mf_ddrc_slim.csv"))

ddr_2<- ddr_2[,-5]
colnames(ddr_2)[5] <- "sample"

ddr_2 <- ddr_2[-c(2,7,12, 16),]
ddr_2$date <- "20220203"

ddr_2 <- ddr_2 %>% 
  mutate(sample = case_when(sample == "DMSO_1" ~ "DMSO",
                            sample == "DMSO_2" ~ "DMSO",
                            sample == "DMSO_3" ~ "DMSO",
                            TRUE ~ as.character(sample)))


# calculate normalization
ddr_1 <- ddr_1 %>%
  mutate(perc = (viable/total))

norm_1 <- ddr_1 %>% 
  group_by(treatment) %>% 
  filter(sample == "DMSO") %>% 
  mutate(dmso = perc)

norm_1 <- left_join(ddr_1, norm_1) %>% 
  group_by(treatment) %>% 
  fill(dmso,.direction = "downup") %>% 
  mutate(norm = (perc/dmso))




ddr_2 <- ddr_2 %>% 
  mutate(perc = (viable/total))

norm_2 <- ddr_2 %>% 
  group_by(treatment) %>% 
  filter(sample == "DMSO") %>% 
  mutate(dmso = perc)

norm_2 <- left_join(ddr_2, norm_2) %>% 
  group_by(treatment) %>% 
  fill(dmso,.direction = "downup") %>% 
  mutate(norm = (perc/dmso))


# add final rep for IVM (AZS and LEV have 2 reps, IVM has 3)
# rep 1
IVM_1 <- read.csv(here("Flow/anthelmintics/20211021_Bma_mf_ddr_IVM.csv"))

IVM_1 <- IVM_1 %>% 
  mutate(
    conc= factor(sample, levels = c("Untreated", "DMSO", "50nM", "1uM", "50uM", "100uM")),
    conc = case_when(
      conc == "DMSO" ~ 0.01,
      conc == "Untreated" ~ 0,
      conc == "50nM" ~ 0.05,
      str_detect(conc, "uM") ~ as.numeric(str_remove(conc, "uM"))),
    date = "20211021")

IVM_1 <- dplyr::select(IVM_1, -group, -vol)
colnames(IVM_1)[3] <- "stage"
IVM_1 <- IVM_1[-1,]
IVM_1 <- IVM_1 %>% 
  mutate(perc = (viable/total))

# calculate normalization
norm_3<- IVM_1 %>% 
  group_by(treatment) %>% 
  filter(sample == "DMSO") %>% 
  mutate(dmso = perc)

norm_3 <- left_join(IVM_1, norm_3 ) %>% 
  fill(dmso,.direction = "up") %>% 
  mutate(norm = (perc/dmso))


# combine the three dataframes
combined <- rbind(norm_1, norm_2, norm_3)


combined$treatment <- factor(combined$treatment, levels = c("AZS", "LEV", "IVM"), labels = c("AZS", "LEV", "IVM"))

# create dataframe for VIM EC50 geom_vline on plot
vline <- subset(combined, treatment == "IVM")
vline$x <- 50.97


# add function for creating logtick marks for a single facet
add_logticks  <- function (base = 10, sides = "bl", scaled = TRUE, 
                           short = unit(0.1, "cm"), mid = unit(0.2, "cm"),  long = unit(0.3, "cm"), 
                           colour = "black",  size = 0.5, linetype = 1, alpha = 1, color = NULL, 
                           data =data.frame(x = NA),... )   {
  if (!is.null(color)) 
    colour <- color
  layer(geom = "logticks", params = list(base = base, 
                                         sides = sides, scaled = scaled, short = short, 
                                         mid = mid, long = long, colour = colour, size = size, 
                                         linetype = linetype, alpha = alpha, ...), 
        stat = "identity", data = data , mapping = NULL, inherit.aes = FALSE,  position = "identity",
        show.legend = FALSE)
}


x <- data.frame(c("DMSO (0.1%)", "0.10", "1", "10", "100"))


#plot

(ddrc<- combined %>% 
    ggplot(aes(x = conc, y = norm))+
    geom_point(color = "grey", size = 1, alpha = 0.9, show.legend = FALSE)+
    geom_line(aes(group = date), color = "grey", size = 0.5, alpha = 0.5, show.legend = FALSE)+
    stat_summary(data = combined, aes(x = conc, y = norm),geom = "line", fun = "mean", color = "black", size = 0.75)+
    geom_vline(data = vline, aes(xintercept = x), color = "red", size = 0.5, alpha = 0.75)+
    facet_grid(rows = vars(treatment))+
    stat_summary(data = combined, aes(x = conc, y = norm),geom = "point", fun = "mean", color = "black", size = 1.5)+
    labs(x = expression(paste(~Log[10](Concentration)~ (µM))), y = "Viability (% DMSO Control)")+
    scale_y_continuous(labels = scales::percent, limits = c(0.25, 1.05), expand = c(0,0))+
    scale_x_log10(breaks = c(0.01, 0.10, 1, 10, 100), labels = c("DMSO (0.1%)", "0.10", "1", "10", "100"))+
    #theme_minimal(base_family = "Helvetica")+
    theme(axis.text.x = ggplot2::element_text(size = 8, hjust = 0.5, vjust = -1,face = "plain"),
          axis.text.y = ggplot2::element_text(size = 8, hjust = 1,face = "plain"),
          axis.title.x = element_text(size = 8, vjust = -1,face = "plain"),
          axis.title.y = ggplot2::element_text(size = 8,face = "plain"), 
          axis.line = ggplot2::element_line(size = 0.25, colour = "black"),
          axis.ticks = ggplot2::element_line(size = 0.25), 
          strip.text = element_text(size = 10,face = "plain"),
          panel.grid = element_line(color = "#ededed", size = 0.25),
          strip.background = element_blank(),
          panel.background = element_blank())+
    NULL) 


# add logtick marks to single facet on the bottom

(Fig5c <- ddrc + 
  add_logticks(side = "b", data = filter(combined, treatment == "IVM"), short = unit(-0.05, "cm"), mid = unit(-0.1, "cm"), long = unit(-0.2, "cm"), size = 0.25)+
  coord_cartesian(clip = "off"))







################
### Fig. 5d - spectralflo contour plots of IVM response
################
Fig5d <- ggdraw() +
  draw_image(magick::image_read_pdf(here("Figures/Figure_5/ivermectin_flowcytometry.pdf"))



###############
### Fig. 5e - 10X utBM/tBM UMAP and IVM DEG volcano plot
###############
new_combined <- readRDS(here("10XGenomics/scmulti_integrated.RDS"))
DefaultAssay(new_combined) <- "RNA"
# pull out the normalized counts, metadata, and UMAP coordinates into dataframes for plotting in ggplot2
data <- as_tibble(new_combined@reductions$umap@cell.embeddings, rownames = 'index') # UMAP coordinates for each cell

md <- as_tibble(new_combined@meta.data, rownames = 'index') # metadata detailing ut/t identity and cluster information

counts2 <- as_tibble(new_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") 




## UMAP
data2 <- data %>% 
  left_join(counts2) %>% 
  left_join(md) %>% 
  select("index", "UMAP_1", "UMAP_2", "orig.ident") %>% 
  distinct() %>% 
  mutate(orig.ident = case_when(
    orig.ident == "utBM" ~ "Untreated",
    orig.ident == "tBM" ~ "1 µM IVM"
  ))

umap <- data2 %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
    geom_point(aes(color = orig.ident), size = 0.005, alpha = 1/10, show.legend = FALSE)+
 scale_color_manual(values = c("#f3933b", "#65838d"), labels = c("IVM (1 µM)", "Untreated"))+
  labs( color = NULL)+
  facet_wrap(~orig.ident, ncol = 1)+
  geom_text(data = data2, mapping = aes(x = -9, y = 8, label = orig.ident, fontface = "plain"), check_overlap = TRUE, family = "Helvetica", size = 2.5)+
  theme(#text = element_text(family = "Helvetica"),
        #axis.text.x = ggplot2::element_text(size = 8, face = "plain"),
        #axis.text.y = ggplot2::element_text(size = 8, face = "plain"),
        #axis.title.x = ggplot2::element_text(size = 10, face = "plain"),
        #axis.title.y = ggplot2::element_text(size = 10, face = "plain"), 
        #axis.line = element_line (colour = "black"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_markdown(size = 10, face = "plain"),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(0, -0.2, 0, -0.2, "cm"))



# Differential gene expression analysis IVM(+) vs untreated

# LOAD RDS BELOW!  -- the following clusters did not have features that passed the logfc.threshold of 0.25.   Positive avg_logFC indicates the gene is more highly expressed in the first group
response <- readRDS(here("Figures/Figure_5/IVM_response.RDS"))


#-------------------------------------------------------------------------------------
BM_combined1 <- new_combined
# find markers between the two treatment groups to look at possible effects of IVM
BM_combined1$celltype.trt <- paste(Idents(BM_combined1), BM_combined1$orig.ident, sep= "_")   #create meta.data slot to store both the cell type and treatment information
BM_combined1$celltype <- Idents(BM_combined1)
Idents(BM_combined1) <- "celltype.trt"



IVM_response1 <- FindMarkers(BM_combined1, ident.1 = "0_tBM", ident.2 = "0_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response1 <- rownames_to_column(IVM_response1, var = "gene_id")
IVM_response1['Cluster'] <- "1"

IVM_response2 <- FindMarkers(BM_combined1, ident.1 = "1_tBM", ident.2 = "1_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response2 <- rownames_to_column(IVM_response2, var = "gene_id")
IVM_response2['Cluster'] <- "2"

IVM_response3 <- FindMarkers(BM_combined1, ident.1 = "2_tBM", ident.2 = "2_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response3 <- rownames_to_column(IVM_response3, var = "gene_id")
IVM_response3['Cluster'] <- "3"

IVM_response4 <- FindMarkers(BM_combined1, ident.1 = "3_tBM", ident.2 = "3_utBM",test.use = "t", verbose = TRUE, logfc.threshold = 0)
IVM_response4 <- rownames_to_column(IVM_response4, var = "gene_id")
IVM_response4['Cluster'] <- "4"

IVM_response5 <- FindMarkers(BM_combined1, ident.1 = "4_tBM", ident.2 = "4_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response5 <- rownames_to_column(IVM_response5, var = "gene_id")
IVM_response5['Cluster'] <- "5"

IVM_response6 <- FindMarkers(BM_combined1, ident.1 = "5_tBM", ident.2 = "5_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response6 <- rownames_to_column(IVM_response6, var = "gene_id")
IVM_response6['Cluster'] <- "6"

IVM_response7 <- FindMarkers(BM_combined1, ident.1 = "6_tBM", ident.2 = "6_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response7 <- rownames_to_column(IVM_response7, var = "gene_id")
IVM_response7['Cluster'] <- "7"

IVM_response8 <- FindMarkers(BM_combined1, ident.1 = "7_tBM", ident.2 = "7_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response8 <- rownames_to_column(IVM_response8, var = "gene_id")
IVM_response8['Cluster'] <- "8"

IVM_response9 <- FindMarkers(BM_combined1, ident.1 = "8_tBM", ident.2 = "8_utBM", verbose = TRUE,logfc.threshold = 0)
IVM_response9 <- rownames_to_column(IVM_response9, var = "gene_id")
IVM_response9['Cluster'] <- "9"

IVM_response10 <- FindMarkers(BM_combined1, ident.1 = "9_tBM", ident.2 = "9_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response10 <- rownames_to_column(IVM_response10, var = "gene_id")
IVM_response10['Cluster'] <- "10"

IVM_response11 <- FindMarkers(BM_combined1, ident.1 = "10_tBM", ident.2 = "10_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response11 <- rownames_to_column(IVM_response11, var = "gene_id")
IVM_response11['Cluster'] <- "11"

IVM_response12 <- FindMarkers(BM_combined1, ident.1 = "11_tBM", ident.2 = "11_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response12 <- rownames_to_column(IVM_response12, var = "gene_id")
IVM_response12['Cluster'] <- "12"

IVM_response13 <- FindMarkers(BM_combined1, ident.1 = "12_tBM", ident.2 = "12_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response13 <- rownames_to_column(IVM_response13, var = "gene_id")
IVM_response13['Cluster'] <- "13"


IVM_response14 <- FindMarkers(BM_combined1, ident.1 = "13_tBM", ident.2 = "13_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response14 <- rownames_to_column(IVM_response14, var = "gene_id")
IVM_response14['Cluster'] <- "14"

IVM_response15 <- FindMarkers(BM_combined1, ident.1 = "14_tBM", ident.2 = "14_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response15 <- rownames_to_column(IVM_response15, var = "gene_id")
IVM_response15['Cluster'] <- "15"

IVM_response16 <- FindMarkers(BM_combined1, ident.1 = "15_tBM", ident.2 = "15_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response16 <- rownames_to_column(IVM_response16, var = "gene_id")
IVM_response16['Cluster'] <- "16"


IVM_response17 <- FindMarkers(BM_combined1, ident.1 = "16_tBM", ident.2 = "16_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response17 <- rownames_to_column(IVM_response17, var = "gene_id")
IVM_response17['Cluster'] <- "17"


IVM_response18 <- FindMarkers(BM_combined1, ident.1 = "17_tBM", ident.2 = "17_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response18 <- rownames_to_column(IVM_response18, var = "gene_id")
IVM_response18['Cluster'] <- "18"

IVM_response19 <- FindMarkers(BM_combined1, ident.1 = "18_tBM", ident.2 = "18_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response19 <- rownames_to_column(IVM_response19, var = "gene_id")
IVM_response19['Cluster'] <- "19"

IVM_response20 <- FindMarkers(BM_combined1, ident.1 = "19_tBM", ident.2 = "19_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response20 <- rownames_to_column(IVM_response20, var = "gene_id")
IVM_response20['Cluster'] <- "20"


IVM_response21 <- FindMarkers(BM_combined1, ident.1 = "20_tBM", ident.2 = "20_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response21 <- rownames_to_column(IVM_response21, var = "gene_id")
IVM_response21['Cluster'] <- "21"

IVM_response22 <- FindMarkers(BM_combined1, ident.1 = "21_tBM", ident.2 = "21_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response22 <- rownames_to_column(IVM_response22, var = "gene_id")
IVM_response22['Cluster'] <- "22"

IVM_response23 <- FindMarkers(BM_combined1, ident.1 = "22_tBM", ident.2 = "22_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response23 <- rownames_to_column(IVM_response23, var = "gene_id")
IVM_response23['Cluster'] <- "23"

IVM_response24 <- FindMarkers(BM_combined1, ident.1 = "23_tBM", ident.2 = "23_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response24 <- rownames_to_column(IVM_response24, var = "gene_id")
IVM_response24['Cluster'] <- "24"

IVM_response25 <- FindMarkers(BM_combined1, ident.1 = "24_tBM", ident.2 = "24_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response25 <- rownames_to_column(IVM_response25, var = "gene_id")
IVM_response25['Cluster'] <- "25"


IVM_response26 <- FindMarkers(BM_combined1, ident.1 = "25_tBM", ident.2 = "25_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response26 <- rownames_to_column(IVM_response26, var = "gene_id")
IVM_response26['Cluster'] <- "26"

IVM_response27 <- FindMarkers(BM_combined1, ident.1 = "26_tBM", ident.2 = "26_utBM", verbose = TRUE, logfc.threshold = 0)
IVM_response27 <- rownames_to_column(IVM_response27, var = "gene_id")
IVM_response27['Cluster'] <- "27"


response <- rbind(IVM_response1, IVM_response2, IVM_response3,IVM_response4, IVM_response5, IVM_response6, IVM_response7, IVM_response8, IVM_response9, IVM_response10, IVM_response11, IVM_response12, IVM_response13, IVM_response14, IVM_response15, IVM_response16, IVM_response17, IVM_response18, IVM_response19, IVM_response20, IVM_response21, IVM_response22, IVM_response23, IVM_response24, IVM_response25, IVM_response26, IVM_response27)

#saveRDS(response, "~/Desktop/IVM_response.RDS")


#-----------------------------------------------------------------------------

# add the bma gene names
list <- read.csv(here("Auxillary/gene_list.csv")) %>%
  select("Gene.stable.ID", "Gene.name")
colnames(list)[1] <- "gene_id"
list <- unique(list)
response <- response %>% left_join(list)

# subset response to only those with significant adjusted p-value
tmp <- response %>% 
  subset(avg_log2FC >= 1 & -log10(p_val_adj) >= -log(0.05))

tmp <- response %>% 
  filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>% 
  filter(-log10(p_val_adj) >= -log(0.05) | -log10(p_val_adj) <= -log(0.05))


response <- response %>% 
  mutate(ID = case_when(
  Cluster == "1" ~ "Unannotated",
  Cluster == "2" ~ "MS",
  Cluster == "3" ~ "Unannotated",
  Cluster == "4" ~ "Unannotated",
  Cluster == "5" ~ "Unannotated",
  Cluster == "6" ~ "C",
  Cluster == "7" ~ "Unannotated",
  Cluster == "8" ~ "Unannotated",
  Cluster == "9" ~ "MD",
  Cluster == "10" ~ "Unannotated",
  Cluster == "11" ~ "Neuron",
  Cluster == "12" ~ "Neuron",
  Cluster == "13" ~ "Neuron",
  Cluster == "15" ~ "S",
  Cluster == "14" ~ "CA",
  Cluster == "16" ~ "Unannotated",
  Cluster == "17" ~ "MD",
  Cluster == "18" ~ "Neuron",
  Cluster == "19" ~ "MS",
  Cluster == "20" ~ "Unannotated",
  Cluster == "21" ~ "Unannotated",
  Cluster == "22" ~ "IB",
  Cluster == "23" ~ "Neuron",
  Cluster == "24" ~ "Neuron",
  Cluster == "25" ~ "Neuron",
  Cluster == "26" ~ "Neuron",
  Cluster == "27" ~ "Neuron"))


# assign cluster colors based on Fig2a golbal umap
response <- response %>% 
  mutate(color = case_when(
    Cluster == "1" ~ "#c1d6d3",
    Cluster == "2" ~ "#5c8492",
    Cluster == "3" ~ "#b25757",
    Cluster == "4" ~ "#6a9491",
    Cluster == "5" ~ "#7a7f84",
    Cluster == "6" ~ "#cab6b2",
    Cluster == "7" ~ "#fae2af",
    Cluster == "8" ~ "#f3933b",
    Cluster == "9" ~ "#ac8287",
    Cluster == "10" ~ "#65838d",
    Cluster == "11" ~ "#82aca7",
    Cluster == "12" ~ "#fe906a",
    Cluster == "13" ~ "#e3e2e1",
    Cluster == "14" ~ "#e89690",
    Cluster == "15" ~ "#cd4c42",
    Cluster == "16" ~ "#6f636b",
    Cluster == "17" ~ "#82ac92",
    Cluster == "18" ~ "#a26f6a",
    Cluster == "19" ~ "#184459",
    Cluster == "20" ~ "#596c7f",
    Cluster == "21" ~ "#263946",
    Cluster == "22" ~ "#d97f64",
    Cluster == "23" ~ "#a0b4ac",
    Cluster == "24" ~ "#e3e2e1",
    Cluster == "25" ~ "#fbc1c1",
    Cluster == "26" ~ "#7f93a2",
    Cluster == "27" ~ "#d76660"))

library(ggrepel)

#plot faceted by annotation
vol_panel <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val_adj)))+
    geom_point(data = subset(response, avg_log2FC <= -1 & -log10(p_val_adj) >= -log(0.05)), color = "#f3933b", show.legend = TRUE, size = 1)+
    geom_point(data = subset(response, avg_log2FC >= 1 & -log10(p_val_adj) >= -log(0.05)),color = "#f3933b", show.legend = TRUE, size = 1)+
    geom_point(data = subset(response, avg_log2FC <= 1 & avg_log2FC > -0.6),color = "#65838d", alpha = 0.5, size = 0.25, show.legend = FALSE)+
    geom_point(data = subset(response, avg_log2FC >= -1 & -log10(p_val_adj) <= -log(0.05)), color = "light grey", alpha = 0.5, size = 0.25, show.legend = FALSE)+
    geom_point(data = subset(response, -log10(p_val_adj) < -log(0.05)), color = "#65838d", alpha = 0.5, size = 0.5, show.legend = FALSE)+
    #geom_text_repel(data = subset(response, avg_log2FC >= 1 & -log10(p_val_adj) >= -log(0.05)), aes(label = Gene.name, fontface = 3),size = 2.5, nudge_y =0, nudge_x = 0, color = "black", max.iter = 100000, max.overlaps = 9, force  = 3)+  
    xlim(-3, 3)+
    geom_hline(yintercept = -log(0.05), linetype = "dotted", col = "black", alpha = 0.3)+
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", col = "black", alpha = 0.3)+
   # annotate(data = subset(response, ID == "C"), geom= "text", x = 2,y = 100, label = "Treated")+
   # annotate(geom= "text", x = -1,y = 100, label = "Untreated")+
    facet_wrap(~ID, ncol = 4, nrow = 2)+
    labs(x = expression(Log["2"]~'(Fold Change)'), y = expression(-Log["10"]~italic(P)), color = "Cluster")+
    geom_text(aes(label = ifelse(avg_log2FC <= -1 & -log10(p_val_adj) >= -log(0.05), as.character(gene_id), "")))+
    #scale_color_identity(guide = "legend", labels = c("1", "Muscle (2)", "3", "4", "5", "Coelomocyte (6)", "7", "8", "Mesoderm (9)", "10", "Neuron (11)", "Neuron (12)", "Neuron (13)", "Canal-assoc. (14)", "Secretory (15)", "16", "Mesoderm (17)", "Neuron (18)", "Muscle (19)", "20", "21", "Inner body (22)", "Unannotated (23)", "Neuron (24)", "Neuron (25)", "Neuron (26)", "Neuron (27)"))+
    theme(axis.text.x = ggplot2::element_text(size = 8, face = "plain"),
          axis.text.y = ggplot2::element_text(size = 8, face = "plain"),
          axis.title.x = ggplot2::element_text(size = 8, face = "plain"),
          axis.title.y = ggplot2::element_text(size = 8,face = "plain"), 
          legend.text = element_markdown(size = 8,face = "plain"),
          panel.background = element_blank(),
          axis.line = element_line (colour = "black"),
          legend.background=element_blank(),
          legend.key = element_blank(),
          legend.spacing.y = unit(0.5, "lines"),
          strip.background = element_blank())+
    guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))+
    NULL


# combined volcao plot colored by cluster and labeled points
volcano <- ggplot(data = response, aes(x = avg_log2FC, y = -log10(p_val_adj)))+
  geom_point(data = subset(response, avg_log2FC <= -1 & -log10(p_val_adj) >= -log(0.05)), aes(color = color), show.legend = TRUE, size = 1)+
  geom_point(data = subset(response, avg_log2FC >= 1 & -log10(p_val_adj) >= -log(0.05)), aes(color = color), show.legend = TRUE, size = 1)+
  geom_point(data = subset(response, avg_log2FC <= 1 & avg_log2FC > -0.6),color = "light grey", alpha = 0.5, size = 0.25, show.legend = FALSE)+
  geom_point(data = subset(response, avg_log2FC >= -1 & -log10(p_val_adj) <= -log(0.05)), color = "light grey", alpha = 0.5, size = 0.25, show.legend = FALSE)+
  geom_point(data = subset(response, -log10(p_val_adj) < -log(0.05)), color = "light grey", alpha = 0.5, size = 0.5, show.legend = FALSE)+
    geom_text_repel(data = subset(response, avg_log2FC >= 1 & -log10(p_val_adj) >= -log(0.05)), aes(label = Gene.name, fontface = 3),size = 2.5, nudge_y =0, nudge_x = 0, color = "black", max.iter = 100000, max.overlaps = 12, force  = 3)+  
  xlim(-3, 3)+
  geom_hline(yintercept = -log(0.05), linetype = "dotted", col = "black", alpha = 0.3)+
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", col = "black", alpha = 0.3)+
  annotate(geom= "text", x = 1.5,y = 115, label = "1 µM IVM", size = 2.5)+
  annotate(geom= "text", x = -1.5,y = 115, label = "Untreated", size = 2.5)+
  #facet_wrap(~ID)+
  labs(x = expression(Log["2"]~'(Fold Change)'), y = expression(-Log["10"]~italic(P)), color = "Cluster")+
  geom_text(aes(label = ifelse(avg_log2FC <= -1 & -log10(p_val_adj) >= -log(0.05), as.character(gene_id), "")))+
  scale_color_identity(guide = "legend", labels = c("1", "Muscle (2)", "3", "4", "5", "Coelomocyte (6)", "7", "8", "Mesoderm (9)", "10", "Neuron (11)", "Neuron (12)", "Neuron (13)", "Canal-assoc. (14)", "Secretory (15)", "16", "Mesoderm (17)", "Neuron (18)", "Muscle (19)", "20", "21", "Inner body (22)", "Unannotated (23)", "Neuron (24)", "Neuron (25)", "Neuron (26)", "Neuron (27)"))+
  theme(axis.text.x = ggplot2::element_text(size = 8, face = "plain"),
        axis.text.y = ggplot2::element_text(size = 8, face = "plain"),
        axis.title.x = ggplot2::element_text(size = 8, face = "plain"),
        axis.title.y = ggplot2::element_text(size = 8,face = "plain"), 
        legend.title = element_text(size = 8),
        legend.text = element_markdown(size = 7,face = "plain"),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.15, 0.5),
        #legend.position = "none",
        #legend.spacing.y = unit(0.5, "lines"),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.1, "cm"))+
        #legend.key.size = unit(0.1, "cm"))+
 guides(color = guide_legend(override.aes = list(size=1, alpha = 1), ncol = 1))+
  NULL




################
### Figure 5
################
row1 <- cowplot::plot_grid(NULL, Fig5a_diagram, Fig5a_flow, ncol = 3, 
                  rel_widths = c(0.01, 0.25, 1),
                  rel_heights = c(1, 0.25, 1),
                  labels = c("A", "", ""), 
                  label_fontface = "plain",
                  label_fontfamily = "Helvetica",
                  label_size = 12, 
                  vjust = 0.5)+ theme(plot.margin = margin(0.3, 0.01, 0, 0, "cm"))



row2 <- cowplot::plot_grid(NULL,Fig5b, Fig5c,  ncol =3 , 
                           rel_widths = c(0.01, 1, 1),
                           labels = c("B", "", "C"), 
                           scale = c(1, 1, 0.99),
                           label_fontface = "plain",
                           label_fontfamily = "Helvetica",
                           label_size = 12, 
                           vjust = c(0, 0),
                           align = "v",
                           axis = "b") + theme(plot.margin = margin(0, 0.01, 0, 0, "cm"))




row3 <- cowplot::plot_grid( Fig5d, NULL, ncol = 2,
                           rel_widths = c( 1.4, 0.01),
                           rel_heights = c( 1.4, 1),
                           scale = c( 1, 1),
                           labels = c("D", ""),
                           label_fontface = "plain",
                           label_fontfamily = "Helvetica",
                           label_size = 12, 
                           vjust = 1) + theme(plot.margin = margin(0, 0, 0, 0, "cm"))

row4 <- cowplot::plot_grid( umap, vol_panel, volcano, ncol = 3,
                            rel_widths = c( 0.4,1.25, 1.35),
                            rel_heights = c(0.5, 1, 1),
                            scale =c(0.75, 1, 1),
                            labels = c("E", ""),
                            label_fontface = "plain",
                            label_fontfamily = "Helvetica",
                            label_size = 12, 
                            vjust = 1.5) + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
  



Figure_5 <- cowplot::plot_grid(row1, row2, row3, row4, nrow = 4, rel_widths = c(1, 1, 1,1), rel_heights = c(1, 0.85, 1,1.1))

ggsave(Figure_5, filename = here("Figures/Figure_5.pdf"), device = cairo_pdf, width = 8.5, height = 10)


