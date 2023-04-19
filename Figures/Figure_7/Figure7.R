# Figure 7 - DIC microscopy and cell culture viability curve
### Figure 6 was assembled in Illustrator. This document creates the viability plot in 7B.

# plotting
library(magick)
library(tidyverse)
library(pdftools)
library(cowplot)
library(ggplot2)
library(ggtext)
library(ggforce)

#data wrangling/analysis
library(tidyverse)
library(dplyr)

#other
setwd("path/to/directory")
library(here)



#################
### Fig. 7b - Viability of cells in culture over 96 hrs
################

viability <- read.csv(here("Auxillary/cellculture_viability.csv")) %>% 
  filter(!is.na(count))
viability$replicate <- as.character(viability$replicate)

viability <- viability %>% 
  mutate(via = count/total)

plot <- viability %>% 
  group_by(replicate, date) %>%
  filter(replicate != "3") %>% 
  ggplot()+
  geom_point(aes(x = hours, y = via, shape = replicate, color = replicate), size = 1.5)+
  scale_color_manual(values = c("#65838d", "#82aca7"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72, 84, 96), labels = c(0, 12, 24, 36, 48, 60, 72, 84, 96))+
  scale_shape_manual(values = c(17, 20))+
  labs(x = "Timepoint (hr)", y = "Viability", color = "Replicate", shape = "Replicate")+
  stat_summary(aes(x = hours, y = via, group = replicate),geom = "line", fun = "mean", color = "black", size = 0.75)+
  stat_summary(aes(x = hours, y = via, group = replicate), fun.data = mean_se, geom = "errorbar", width = 0.75)+
  theme(#text = element_text(family = "helvetica"),
    axis.text.x = ggplot2::element_text(size = 11,face = "plain"),
    axis.text.y = ggplot2::element_text(size = 11, face = "plain"),
    axis.title.x = ggplot2::element_text(size = 12, face = "plain"),
    axis.title.y = ggplot2::element_text(size = 12,face = "plain"), 
    axis.line = ggplot2::element_line(size = 0.25, colour = "black"),
    axis.ticks = ggplot2::element_line(size = 0.25), 
    strip.text = element_text(size = 8,face = "plain"), 
    legend.text = element_text(size = 10,face = "plain"),
    legend.title = element_text(size = 10,face = "plain"),
    legend.key = element_blank(),
    panel.grid.major = element_line(color = "#ededed", size = 0.25),
    strip.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank())+
guides(shape = guide_legend(override.aes = list(size = 2.5)))


ggsave(plot, filename = "~/Desktop/viability_plot.pdf", device = cairo_pdf, width = 6, height = 3, units = "in")




