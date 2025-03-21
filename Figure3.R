## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Calculate Bray-Curtis dissimilarity 
## - Plot Figure 3
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(cowplot)
library(ggplot2)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

#Subset ps objects by time point
ps_t01 <- subset_samples(ps_filtered_norm, tpreal == "T01")
ps_t02 <- subset_samples(ps_filtered_norm, tpreal == "T02")
ps_t03 <- subset_samples(ps_filtered_norm, tpreal == "T03")

#Plots format
main_theme_pcoa3 <- theme(panel.background=element_blank(),
                          panel.border = element_rect(color = "black", fill=NA, size = 0.6),
                          panel.grid=element_blank(),
                          axis.line=element_line(color="black", size=0.6),
                          axis.ticks=element_line(color="black", size=0.6),
                          legend.background=element_blank(),
                          legend.key = element_rect(fill = NA),
                          axis.text.y = element_text(size = 10, color = "black"),
                          axis.text.x = element_text(size = 10, color = "black"),
                          axis.title.y = element_text(size = 12, color = "black"),
                          axis.title.x = element_text(size = 12, color = "black")
                          #text=element_text(size=10, color="black")
)

#Calculate Bray-Curtis dissimilarity and plot PCoAs
set.seed(1025)
bray_dist = phyloseq::distance(ps_t01, method="bray", weighted=T)
ordination = ordinate(ps_t01, method="PCoA", distance=bray_dist)
pcoa_t01 <- plot_ordination(ps_t01, ordination, color = "compartment") + 
labs(title = "T01") + geom_point(size=2)+scale_colour_manual(values = c("#CC9999", "#66CC00","#990033","#33CCFF"))+
stat_ellipse()+guides(color=FALSE)+main_theme_pcoa3

bray_dist = phyloseq::distance(ps_t02, method="bray", weighted=T)
ordination = ordinate(ps_t02, method="PCoA", distance=bray_dist)
pcoa_t02 <- plot_ordination(ps_t02, ordination, color = "compartment") + 
labs(subtitle = "T02") + geom_point(size=2)+scale_colour_manual(values = c("#CC9999", "#66CC00","#990033","#33CCFF"))+
stat_ellipse()+guides(color=FALSE)+main_theme_pcoa3

bray_dist = phyloseq::distance(ps_t03, method="bray", weighted=T)
ordination = ordinate(ps_t03, method="PCoA", distance=bray_dist)
pcoa_t03 <- plot_ordination(ps_t03, ordination, color = "compartment") + 
labs(subtitle = "T03") + geom_point(size=2)+scale_colour_manual(values = c("#CC9999", "#FFCC33", "#66CC00","#990033","#33CCFF"))+
stat_ellipse()+guides(color=FALSE)+main_theme_pcoa3

fig3 <- plot_grid(pcoa_t01, pcoa_t02, pcoa_t03,
                     ncol = 3, nrow = 1,
                     align = "v", axis = "l"
                     #rel_heights = c(1, 0.95, 0.75)
)

ggsave("Figure 3_16S.pdf", fig3, units = "mm", width = 180, height = 65, dpi = 300)
ggsave("Figure 3_T01_16S.pdf", pcoa_t01, units = "mm", width = 66.7, height = 60, dpi = 300)
ggsave("Figure 3_T02_16S.pdf", pcoa_t02, units = "mm", width = 66.7, height = 60, dpi = 300)
ggsave("Figure 3_T03_16S.pdf", pcoa_t03, units = "mm", width = 66.7, height = 60, dpi = 300)