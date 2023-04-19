############
## FIG 1 Optimization of single-cell dispersions
############
#data wrangling/analysis
library(tidyverse)
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

################
### Figure 1a
################

# load illustrator generated workflow scheme
tmp <-image_read_pdf(here("Figures/Figure_1/bma_sc_workflow_trimmed.pdf"), density = 600)
Fig1a <- ggdraw() + draw_image(tmp)



Fig1a <- plot_grid( NULL, Fig1a, ncol = 2, rel_widths = c(0.06,1), labels = c("A", NULL), label_size = 12, label_fontface = "plain", label_fontfamily = "helvetica")



#################
### Figure 1b - PD-10 filtration
#################

pd1 <- ggdraw()+
  draw_image(here("sc_protocol/optimization/pd10/20201119-p14-MZ_A01_w2_binary.png"))

pd2 <- ggdraw()+
  draw_image(here("optimization/pd10/20210820-p14-NJW_A01_binary.png"))

Fig1b <- plot_grid(NULL, pd1,NULL, pd2, ncol = 4, label_size = 12, rel_widths = c(0.2,1, 0.01, 1), label_fontface = "plain", label_fontfamily = "helvetica")



################
### Figure 1c - chitinase treatment
################

# load chitinase images = fig1c
tmp<-image_read(here("optimization/chitinase/20210325_chitinase_untreated_2.png"), density = 200)
chi1 <- ggdraw() + draw_image(tmp)


tmp<-image_read(here("optimization/chitinase/20210325_chitinase_2mgml_20min.png"), density = 200)
chi2 <- ggdraw() + draw_image(tmp)


Fig1c <- plot_grid(NULL,chi1,NULL, chi2, ncol = 4, label_size = 12, rel_widths = c(0.2,1, 0.01, 1), label_fontface = "plain", label_fontfamily = "helvetica")



################
### Figure 1d
################

# load pronase time course images and create single panel = fig1d

#15 min
tmp<-image_read(here("optimization/pronase/15min_pronase/20220208-p03-CRH_A07_s18_15min.png"), density = 300)
pro15.2 <- ggdraw() +
  draw_image(tmp)


#20min
tmp<-image_read(here("optimization/pronase/20min_pronase/20220208-p03-CRH_B07_s53_20min.png"), density = 300)
pro20.2 <- ggdraw() +
  draw_image(tmp)


#24min
tmp<-image_read(here("optimization/pronase/24min_pronase/20220208-p03-CRH_D07_s7_24min.png"), density = 300)
pro24.2 <- ggdraw() +
  draw_image(tmp)


#30min
tmp<-image_read(here("optimization/pronase/30min_pronase/20220208-p03-CRH_G07_s74_30min.png"), density = 300)
pro30.2 <- ggdraw() +
  draw_image(tmp)


#32min
tmp<-image_read(here("optimization/pronase/32min_pronase/20220208-p03-CRH_H07_s16_32min_50um.png"), density = 300)
pro32.2 <- ggdraw() +
  draw_image(tmp)



#Generate Fig1d 
Fig1d <- plot_grid(NULL, pro15.2,NULL, pro20.2,NULL, pro24.2, NULL, pro30.2, NULL, pro32.2, ncol = 10, labels = "D", vjust = 0, label_size = 12, rel_widths = c(0.1, 1,0.01, 1,0.01, 1,0.01, 1,0.01, 1), label_fontface = "plain", label_fontfamily = "helvetica")




######################
#### Figure 1e
######################

# load ImageStream plots and images

#imagestream gating of live cells
IS1 <- ggdraw()+
  draw_image(magick::image_read_pdf(here("Figures/Figure_1/20211110_imagestream_gate.pdf")))+ theme(plot.margin = margin(0.6, 0, 0, 0.7))



#generating aspect ratios and cell diameter plot
combined_es <- read.csv(here("Figures/Figure_1/imagestream_featurewizard.csv"))

# filter dataframe to get es cells to be highlighted and create a second legend
highlight_es <- combined_es %>% 
  filter(stain=="es")
aspect_ratios <- data.frame(combined_es$Aspect.Ratio_M01, combined_es$Aspect.Ratio_M05)
colnames(aspect_ratios) <- c("Aspect.Ratio_M01", "Aspect.Ratio_M05")
aspect_ratios$golden <- aspect_ratios$Aspect.Ratio_M01 >= 0.8 & aspect_ratios$Aspect.Ratio_M05 >=0.8

ratios_tidy <- aspect_ratios %>% 
  pivot_longer(cols = c(Aspect.Ratio_M01, Aspect.Ratio_M05), names_to = "ratio_type", values_to = "value")



# plot
(testplot_3 <- combined_es %>% 
    ggplot(aes(x = Diameter, y = Aspect.Ratio_M01))+
    geom_point(aes(colour = Aspect.Ratio_M01 > 0.8 & Aspect.Ratio_M05 > 0.8), alpha = 0.7, size = 0.8, position = position_jitter(width = 1, seed = 2))+
    geom_point(data = highlight_es, aes(x = Diameter, y = Aspect.Ratio_M01, alpha = "ES Pore-like"), colour = "black")+
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.78, ymax = 0.82), alpha = 0.025, fill = " gray93")+
    labs(x = "Cell Diameter (µm)", y = "Aspect Ratios", colour = "Aspect Ratios")+
    scale_alpha_manual(name = "Cell Type", values = c(1), breaks = c("ES Pore-like"))+
    scale_colour_manual(labels = c("Multiplet/Debris", "Single-cells"), values = c("#DD8D29", "#5BBCD6"))+
    theme(text = element_text(family = "helvetica"),
          plot.margin = margin(-0.5, 2.15, 0.5, 0.25, "cm"), 
          axis.line = element_line(color = "black"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(size = 9.5),
          legend.text = element_text(size = 10),
          legend.key.height = unit(0.05, "cm"),
          legend.key.width = unit(0.05, "cm"),
          legend.title = element_text(size = 10),
          legend.key = element_blank(),
          legend.key.size = unit(0.05, "cm"),
          legend.spacing = unit(0.05, "cm"),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.position = c(1.3, 0.22),
          panel.grid.major = element_blank())+
    guides(colour = guide_legend(override.aes = list(size = 2)))+
    guides(alpha = guide_legend(override.aes = list(size = 2)))+
    NULL)



xhist2 <- axis_canvas(testplot_3, axis = "x")+
  geom_density(data = combined_es, aes(x = Diameter, colour = golden))+
  scale_colour_manual(values = c("#DD8D29", "#5BBCD6"))

yhist2 <- axis_canvas(testplot_3, axis = "y")+
  #geom_histogram(data = ratios_tidy, aes(y = value, fill = golden), bins = 200, alpha = 0.8, position = "identity")+
  geom_density(data = ratios_tidy, aes(y = value, color = golden))+
  scale_color_manual(values = c("#DD8D29", "#5BBCD6"))

(IS2 <- ggdraw(insert_xaxis_grob(testplot_3, xhist2, height = grid::unit(0.3, "null"), position = "top") %>% 
                              insert_yaxis_grob(., yhist2, position = "right")))






# reduce dataframe for calculating average diameter and SD
stat <- combined_es %>% 
  select("X", "Object.Number", "stain", "Diameter", "golden") %>% 
  subset(golden == "TRUE")

## calculate the average cell size for cells meeting the aspect ratio requirements (>= 0.8 for brightfield and DRAQ5 channels)
mean(stat$Diameter)  #5.22831 um
sd(stat$Diameter) # 0.773079 um 



# load representative images for each cell classification (single, multiple, debris, es pore-like)

tmp<-image_read_pdf(here("Figures/Figure_1/imagestream/0.9aspectratio/0.9aspectratio_montage.pdf"), density = 200)
IS3.1 <- ggdraw() +
  draw_image(tmp)


tmp<-image_read_pdf(here("Figures/Figure_1/imagestream/grayregion_aspectratio/grayregion_montage_2.pdf"), density = 200)
IS3.2 <- ggdraw() +
  draw_image(tmp)



tmp<-image_read_pdf(here("Figures/Figure_1/imagestream/under0.8_aspectratio/under0.8aspectratio_montage.pdf"), density = 200)
IS3.3 <- ggdraw() +
  draw_image(tmp)


tmp <-image_read_pdf(here("Figures/Figure_1/imagestream/espores_montage.pdf"), density = 200)
IS3.4 <- ggdraw() +
  draw_image(tmp)


# panel all 4 groups of representative images
cells <- plot_grid(NULL, IS3.1, NULL, IS3.2, NULL, IS3.3, NULL, IS3.4, nrow = 2, ncol = 4, rel_widths = c(0.025, 1, 0.02,1))+
  theme(plot.margin = margin(0, 0.01, 0, 0.01))







# Generate figure rows

row1 <- Fig1a

null1 <- plot_grid(NULL, NULL, NULL, NULL, ncol = 4, labels = c("Unfiltered", "Filtered", "(-) Chitinase", "(+) Chitinase"), label_size = 10, label_fontfamily = "Helvetica",hjust = c(-1.3, -1.5, -1, -0.65), vjust = c(3, 3, 3.25, 3.25), label_fontface = "plain")

row2 <- plot_grid(Fig1b, NULL, Fig1c, ncol = 3, labels = c("B", "", "C"), vjust = 0, scale = c(0.95, 0.95, 1, 1), rel_widths = c(1, 0.01, 1), label_size = 12, label_fontface = "plain", label_fontfamily = "helvetica") 

null2 <- plot_grid(NULL, NULL, NULL, NULL, NULL, ncol = 5, labels = c("15 min.", "20 min.", "24 min.", "30 min.", "32 min."), label_size = 10, label_fontfamily = "Helvetica", hjust = c(-1.6, -1.5, -1.3, -1.3, -1.2), vjust = 2.3, label_fontface = "plain")
  
row3 <- Fig1d

row4 <- plot_grid(NULL, IS1,NULL, IS2, ncol = 4, labels = c("E", "","", ""),rel_widths = c(0.05, 0.75,0.03, 1.45), hjust = c(-0.1, 0,1, 0), vjust = c(0.75,1,1, 0.5), scale = c(1, 1.12, 1,1), label_fontface = "plain", label_fontfamily = "helvetica")  

row5 <- cells



# Generate complete figure
Figure1 <- plot_grid(row1, null1, row2, null2, row3, row4, row5, nrow = 7, scale = c(1,1,1, 1, 1, 1,1), rel_heights = c(1, 0.1, 1, 0.025, 1, 1.45, 1))+
  theme(plot.margin = unit(c(0, 0.25, 0, 0), "cm"))


# export pdf file
ggsave(Figure1, filename = "Figure1.pdf", device = cairo_pdf, width = 7, height = 10, units = "in")





