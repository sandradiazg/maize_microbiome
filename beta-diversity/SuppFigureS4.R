## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Calculate Bray-Curtis dissimilarity 
## - Plot Supplementary Figure S4
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(cowplot)
library(ggplot2)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

#Plot's theme
main_theme_pcoa <- theme(panel.background=element_blank(),
                         panel.border = element_rect(color = "black", fill=NA, size = 0.6),
                         panel.grid=element_blank(),
                         axis.line=element_line(color="black", size=0.6),
                         axis.ticks=element_line(color="black", size=0.6),
                         legend.background=element_blank(),
                         legend.key = element_rect(fill = NA),
                         legend.title = element_text(colour="black", size=10, face="bold"),
                         legend.text = element_text(colour="black", size = 8),
                         legend.spacing.y = unit(0.5, 'mm'),
                         legend.spacing.x = unit(0.5, 'mm'),
                         legend.key.height = unit(4, "mm"),
                         axis.text.y = element_text(size = 10, color = "black"),
                         axis.text.x = element_text(size = 10, color = "black"),
                         axis.title.y = element_text(size = 12, color = "black"),
                         axis.title.x = element_text(size = 12, color = "black"),
                         text=element_text(size=12, color="black"))

cols_comp <- c("Bulk soil" = "#CC9999", 
               "Rhizosphere" = "#990033", 
               "Root" = "#33CCFF", 
               "Leaf" = "#66CC00", 
               "Grain" = "#FFCC33")

irrigation_cols <- c("OW" = "#4a80ff",
                     "LW" = "tan2")

tp_cols <- c("T01" = "#66cccc",
             "T02" = "#3399CC",
             "T03" = "#003333")

#Supplementary Figure S4A
bray_dist = phyloseq::distance(ps_filtered_norm, method="bray", weighted=T)
ordination = ordinate(ps_filtered_norm, method="PCoA", distance=bray_dist)
Sfig4A <- plot_ordination(ps_filtered_norm, ordination,color = "compartment") + 
  theme(aspect.ratio=1) + geom_point(size=1.5)+
  scale_colour_manual(values = cols_comp) + 
  guides(color=FALSE)+
  main_theme_pcoa

#Supplementary Figure S4B
bray_dist = phyloseq::distance(ps_filtered_norm, method="bray", weighted=T)
ordination = ordinate(ps_filtered_norm, method="PCoA", distance=bray_dist)
Sfig4B <- plot_ordination(ps_filtered_norm, ordination, color = "tpreal")+ 
  theme(aspect.ratio=1) + geom_point(size=1.5)+
  scale_colour_manual(values = tp_cols)+
  guides(color=FALSE)+
  main_theme_pcoa


#Supplementary Figure S4C
bray_dist = phyloseq::distance(ps_filtered_norm, method="bray", weighted=T)
ordination = ordinate(ps_filtered_norm, method="PCoA", distance=bray_dist)
Sfig4C <- plot_ordination(ps_filtered_norm, ordination, color = "irrigation") + 
  theme(aspect.ratio=1) + geom_point(size=1.5)+
  scale_colour_manual(values = irrigation_cols)+
  guides(color=FALSE)+
  main_theme_pcoa

#Compile figure
Sfig4 <- plot_grid(Sfig4A, Sfig4B, Sfig4C,
                   ncol = 3, nrow = 1,
                   align = "v", axis = "l"
                   #rel_heights = c(1, 0.95, 0.75)
                  )
#Save figures
ggsave("SFigure_S4A_16s.pdf", Sfig4A, units = "mm", width = 70, height = 70, dpi = 300)
ggsave("SFigure_S4B_16S.pdf", Sfig4B, units = "mm", width = 70, height = 70, dpi = 300)
ggsave("SFigure_S4C_16S.pdf", Sfig4C, units = "mm", width = 70, height = 70, dpi = 300)
ggsave("SFigure_S4_16S.pdf", Sfig4, units = "mm", width = 200, height = 70, dpi = 300)
