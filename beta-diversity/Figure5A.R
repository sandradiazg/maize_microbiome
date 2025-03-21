## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Calculate Bray-Curtis dissimilarity 
## - Plot Figure 5A
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(cowplot)
library(ggplot2)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

#Subset ps objects by compartment
ps_BULK_SOIL <- subset_samples(ps_filtered_norm, compartment == "Bulk soil")
ps_RHIZOSPHERE <- subset_samples(ps_filtered_norm, compartment == "Rhizosphere")
ps_ROOT <- subset_samples(ps_filtered_norm, compartment == "Root")
ps_LEAF <- subset_samples(ps_filtered_norm, compartment == "Leaf")
ps_GRAIN <- subset_samples(ps_filtered_norm, compartment == "Grain")


#Plot's format
main_theme_pcoa2<- theme(panel.background=element_blank(),
                         panel.border = element_rect(color = "black", fill=NA, size = 0.5),
                         panel.grid=element_blank(),
                         axis.line=element_line(color="black", size=0.4),
                         axis.ticks=element_line(color="black", size=0.4),
                         legend.background=element_blank(),
                         axis.text.y = element_text(size = 7, color = "black"),
                         axis.text.x = element_text(size = 7, color = "black"),
                         axis.title.y = element_text(size = 8, color = "black"),
                         axis.title.x = element_text(size = 8, color = "black"),
                         text=element_text(size=8, color="black"))

irrigation_cols <- c("OW" = "#4a80ff",
                     "LW" = "tan2")

tp_shape <- c("T01" = 22,
              "T02" = 21,
              "T03" = 24)
#Fig5A
#Bulk soil
bray_dist = phyloseq::distance(ps_BULK_SOIL, method="bray", weighted=T)
ordination = ordinate(ps_BULK_SOIL, method="PCoA", distance=bray_dist)
fig5A_bs <- plot_ordination(ps_BULK_SOIL, ordination,  color = "irrigation", shape="tpreal") + 
  theme(aspect.ratio=1) + geom_point(size=0.7)+
  scale_colour_manual(values = irrigation_cols)+
  scale_shape_manual(values = c(15,16,17))+
  guides(color=FALSE, shape = FALSE)+
  main_theme_pcoa2

#Rhizosphere
bray_dist = phyloseq::distance(ps_RHIZOSPHERE, method="bray", weighted=T)
ordination = ordinate(ps_RHIZOSPHERE, method="PCoA", distance=bray_dist)
fig5A_rz <- plot_ordination(ps_RHIZOSPHERE, ordination,  color = "irrigation", shape="tpreal") + 
  theme(aspect.ratio=1) + geom_point(size=0.7)+
  scale_colour_manual(values = irrigation_cols)+
  scale_shape_manual(values = c(15,16,17))+
  guides(color=FALSE, shape = FALSE)+
  main_theme_pcoa2

#Root
bray_dist = phyloseq::distance(ps_ROOT, method="bray", weighted=T)
ordination = ordinate(ps_ROOT, method="PCoA", distance=bray_dist)
fig5A_rt <- plot_ordination(ps_ROOT, ordination,  color = "irrigation", shape="tpreal") + 
  theme(aspect.ratio=1) + geom_point(size=0.7)+
  scale_colour_manual(values = irrigation_cols)+
  scale_shape_manual(values = c(15,16,17))+
  guides(color=FALSE, shape = FALSE)+
  main_theme_pcoa2

#Leaf
bray_dist = phyloseq::distance(ps_LEAF, method="bray", weighted=T)
ordination = ordinate(ps_LEAF, method="PCoA", distance=bray_dist)
fig5A_lf <- plot_ordination(ps_LEAF, ordination,  color = "irrigation", shape="tpreal") + 
  theme(aspect.ratio=1) + geom_point(size=0.7)+
  scale_colour_manual(values = irrigation_cols)+
  scale_shape_manual(values = c(15,16,17))+
  guides(color=FALSE, shape = FALSE)+
  main_theme_pcoa2


#Grain
bray_dist = phyloseq::distance(ps_GRAIN, method="bray", weighted=T)
ordination = ordinate(ps_GRAIN, method="PCoA", distance=bray_dist)
fig5A_gr <- plot_ordination(ps_GRAIN, ordination,  color = "irrigation", shape="tpreal") + 
  theme(aspect.ratio=1)  + geom_point(size=0.7)+
  scale_colour_manual(values = irrigation_cols)+
  scale_shape_manual(values = 17)+
  guides(color=FALSE, shape = FALSE)+
  main_theme_pcoa2

#Compile in one figure
fig5A <- plot_grid(fig5A_bs, fig5A_rz, fig5A_rt,
                   fig5A_lf, fig5A_gr,
                   ncol = 5, nrow = 1,
                   align = "v", axis = "l"
                   #rel_heights = c(1, 0.95, 0.75)
)

ggsave("Figure 5_16S.pdf", fig5A, units = "mm", width = 180, height = 60, dpi = 300)
ggsave("Figure 5_bs_16S.pdf", fig5A_bs, units = "mm", width = 50, height = 60, dpi = 300)
ggsave("Figure 5_rz_16S.pdf", fig5A_rz, units = "mm", width = 50, height = 60, dpi = 300)
ggsave("Figure 5_rt_16S.pdf", fig5A_rt, units = "mm", width = 50, height = 60, dpi = 300)
ggsave("Figure 5_lf_16S.pdf", fig5A_lf, units = "mm", width = 50, height = 60, dpi = 300)
ggsave("Figure 5_gr_16S.pdf", fig5A_gr, units = "mm", width = 50, height = 60, dpi = 300)