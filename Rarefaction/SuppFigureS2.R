## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Create rarefaction curves for each compartment
## - Plot Supplementary Figure S2
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filter ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(ggplot2)
library(phyloseq.extended)


#Load ps object (note that this one is not normalize!)
ps_filter <- readRDS("ps_filter_16S.rds")


#Subset samples by compartment
ps_bs<- subset_samples(ps_filter, compartment == "Bulk soil" )
ps_rz <- subset_samples(ps_filter, compartment == "Rhizosphere")
ps_rt <- subset_samples(ps_filter, compartment == "Root")
ps_lf<- subset_samples(ps_filter, compartment == "Leaf")
ps_gr <- subset_samples(ps_filter, compartment == "Grain")

#Plot's format
main_theme_rare <- theme(panel.background=element_blank(),
                         panel.border = element_rect(color = "black", fill=NA, size = 0.6),
                         panel.grid=element_blank(),
                         axis.line=element_line(color="black", size=0.6),
                         axis.ticks=element_line(color="black", size=0.6),
                         legend.background=element_blank(),
                         legend.key = element_rect(fill = NA),
                         axis.text.y = element_text(size = 7, color = "black"),
                         axis.text.x = element_text(size = 7, color = "black"),
                         axis.title.y = element_text(size = 7, color = "black"),
                         axis.title.x = element_text(size = 7, color = "black"),
                         plot.subtitle = element_text(size = 7),
                         text=element_text(size=8, color="black"))

cols_comp <- c("Bulk soil" = "#CC9999", 
               "Rhizosphere" = "#990033", 
               "Root" = "#33CCFF", 
               "Leaf" = "#66CC00", 
               "Grain" = "#FFCC33")

#Calculate and plot rarefaction curves for each compartment
p1 <- ggrare(ps_bs, step = 1000, color = "compartment", se = FALSE)+labs(subtitle = "Bulk soil")+scale_colour_manual(values=cols_comp)+guides(color=FALSE)+main_theme_rare
p2 <- ggrare(ps_rz, step = 1000, color = "compartment", se = FALSE)+labs(subtitle = "Rhizosphere")+scale_colour_manual(values=cols_comp)+guides(color=FALSE)+main_theme_rare
p3 <- ggrare(ps_rt, step = 1000, color = "compartment", se = FALSE)+labs(subtitle = "Root")+scale_colour_manual(values=cols_comp)+guides(color=FALSE)+main_theme_rare
p4 <- ggrare(ps_lf, step = 1000, color = "compartment", se = FALSE)+labs(subtitle = "Leaf")+scale_colour_manual(values=cols_comp)+guides(color=FALSE)+main_theme_rare
p5 <- ggrare(ps_gr, step = 1000, color = "compartment", se = FALSE)+labs(subtitle = "Grain")+scale_colour_manual(values=cols_comp)+guides(color=FALSE)+main_theme_rare

#Compile figures
Sfig2 <- plot_grid(p1, p2, p3, p4, p5,
                   ncol = 3, nrow = 2
                   # labels = c('B', 'D', 'F','H','J'),
                   # label_size = 12
                   #rel_heights = c(1, 0.95, 0.75)
)
#Save figure
ggsave("SFigure_S2_16S.pdf", Sfig2, units = "mm", width = 200, height = 70, dpi = 300)
