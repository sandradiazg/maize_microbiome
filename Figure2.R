## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Estimate alpha-diversity indexes
## - Plot Figure 2
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(phyloseqCompanion)
library(tibble)
library(ggplot2)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

#calculate alpha-diversity indexes
rich = estimate_richness(ps_filtered_norm)
rich
#Merge sample_medata with alpha-diversity matrix
df_sample_metadata <- sample.data.frame(ps_filtered_norm)
df_sample_metadata <- as.data.frame(df_sample_metadata)
rich <- rownames_to_column(rich, var = "name")
rich

rich_metadata <- merge(df_sample_metadata, rich, by.x = "name",
                       by.y = "name")
rich_metadata

write.table(rich_metadata ,"results_alpha_div_16S.txt", sep = ";", quote = F)


#Plots format

main_theme_box <- theme(panel.background=element_blank(),
                        panel.border = element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(color="black", size=0.3),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_line(color="black", size=0.3),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 8, color = "black"),
                        axis.title.y = element_text(size = 10, color = "black"),
                        legend.background=element_blank(),
                        legend.key = element_rect(fill =),
                        legend.title = element_text(colour="black", size=10, face="bold"),
                        legend.text = element_text(colour="black", size = 8),
                        legend.spacing.y = unit(0.5, 'mm'),
                        legend.spacing.x = unit(0.5, 'mm'),
                        legend.key.height = unit(4, "mm"),
                        strip.text.x = element_blank(),
                        strip.background = element_rect(colour=NA, fill=NA))

compartment_names <- c(
  'Bulk soil'="Bulk soil",
  'Rhizosphere'="Rhizosphere",
  'Root'="Root",
  'Leaf'="Leaf",
  'Grain' = "Grain")#Just if you want to change the names of your treatments for the graphs

cols_comp <- c("Bulk soil" = "#CC9999", 
               "Rhizosphere" = "#990033", 
               "Root" = "#33CCFF", 
               "Leaf" = "#66CC00", 
               "Grain" = "#FFCC33")
compartment_order <- c("Bulk soil", "Rhizosphere", "Root", "Leaf", "Grain") #To determine a specific order for your treatments in your graphs
irrigation_order <- c("OW", "LW")
irrigation_names <- c("OW" = "OW",
                      "LW" = "LW")
irrigation_cols <- c("OW" = "#4a80ff",
                     "LW" = "tan2")
tp_shape <- c("T01" = 22,
              "T02" = 21,
              "T03" = 24)

tp_cols <- c("T01" = "#66cccc",
             "T02" = "#3399CC",
             "T03" = "#003333")

##Figure 2A
p_chao<-plot_richness(ps_filtered_norm, x= "compartment", color= "compartment", measures=c("Chao1")) 
p_chao$layers <- p_chao$layers[-c(1,2)]
fig2A<-p_chao+geom_boxplot(lwd=0.5, alpha=0.4, size =4, width=0.8, outlier.shape = NA, aes(x=factor(compartment, level= compartment_order), fill=sample_data(ps_filtered_norm)$compartment))+ 
  geom_jitter(stroke=0.4, alpha=0.5, size =1, aes(color = compartment, shape = tpreal, fill=sample_data(ps_filtered_norm)$compartment))+
  labs(x="", y="Chao1 Index") + 
  ylim(0,10000)+
  scale_x_discrete(labels= compartment_names)+
  scale_colour_manual(values=cols_comp)+
  scale_fill_manual(values= cols_comp)+
  scale_shape_manual(values = tp_shape)+
  facet_wrap(~tpreal) +
  guides(fill=FALSE, color = FALSE, shape = FALSE)+
  main_theme_box

#Figure 2B
p_shn<-plot_richness(ps_filtered_norm, x= "compartment", color= "compartment", measures=c("Shannon")) 
p_shn$layers <- p_shn$layers[-c(1,2)]
fig2B<-p_shn+geom_boxplot(lwd=0.5, alpha=0.4, size =4, width=0.8, outlier.shape = NA, aes(x=factor(compartment, level=compartment_order), fill=sample_data(ps_filtered_norm)$compartment))+ 
  geom_jitter(stroke=0.4, alpha=0.5, size =1, aes(color = compartment, shape = tpreal, fill=sample_data(ps_filtered_norm)$compartment))+
  labs(x="", y="Shannon Index") + 
  scale_x_discrete(labels=compartment_names)+
  scale_colour_manual(values=cols_comp)+
  scale_fill_manual(values= cols_comp)+
  scale_shape_manual(values = tp_shape)+
  facet_wrap(~tpreal) +
  guides(fill=FALSE, color = FALSE, shape = FALSE)+
  main_theme_box


#Figure 2C
p_chao<-plot_richness(ps_filtered_norm, x= "tpreal", color= "irrigation", measures=c("Chao1")) 
p_chao$layers <- p_chao$layers[-c(1,2)]
fig2C <- p_chao+
  geom_boxplot(lwd=0.35, alpha=0.4, size =4, width=0.7, outlier.shape = NA, 
               aes(x=factor(tpreal), fill=sample_data(ps_filtered_norm)$irrigation))+ 
  geom_jitter(stroke=0.3, alpha=0.5, size =0.7, position = position_jitterdodge(),
              aes(color = irrigation, shape = tpreal, fill=sample_data(ps_filtered_norm)$irrigation))+
  labs(x="", y="Chao1 Index") + 
  scale_colour_manual(name = "Irrigation",breaks = irrigation_order, 
                      labels = irrigation_names,values=irrigation_cols)+
  scale_fill_manual(name = "Irrigation", breaks= irrigation_order, 
                    labels = irrigation_names,values=irrigation_cols)+
  scale_shape_manual(values = tp_shape)+
  ylim(0,10000)+
  facet_wrap(~factor(compartment,levels = compartment_names), 
             ncol=5, labeller = as_labeller(compartment_names)) +
  guides(fill= FALSE, color = FALSE, shape = FALSE)+
  main_theme_box


p_shn<-plot_richness(ps_filtered_norm, x= "tpreal", color= "irrigation", measures=c("Shannon")) 
p_shn$layers <- p_shn$layers[-c(1,2)]
fig2D <- p_shn+
  geom_boxplot(lwd=0.35, alpha=0.4, size =4, width=0.7, outlier.shape = NA, 
               aes(x=factor(tpreal), fill=sample_data(ps_filtered_norm)$irrigation))+ 
  geom_jitter(stroke=0.3, alpha=0.5, size =0.7, position = position_jitterdodge(),
              aes(color = sample_data(ps_filtered_norm)$irrigation, shape = factor(tpreal), fill=sample_data(ps_filtered_norm)$irrigation))+
  labs(x="", y="Shannon Index") + 
  scale_colour_manual(name = "Irrigation", breaks= irrigation_order, 
                      labels = irrigation_names,values=irrigation_cols)+
  scale_fill_manual(name = "Irrigation", breaks= irrigation_order, 
                    labels = irrigation_names,values=irrigation_cols)+
  scale_shape_manual(values = tp_shape)+
  facet_wrap(~factor(compartment,levels = compartment_names), 
             ncol=5, labeller = as_labeller(compartment_names)) +
  guides(fill= FALSE, color = FALSE, shape = FALSE)+
  main_theme_box

#Save figures as PDF (vectorial format)
ggsave("Figure_2A_16S.pdf", fig2A, units = "mm", width = 85, height = 50, dpi = 300)
ggsave("Figure_2B_16S.pdf", fig2B, units = "mm", width = 85, height = 50, dpi = 300)
ggsave("Figure_2C_16S.pdf", fig2C, units = "mm", width = 85, height = 50, dpi = 300)
ggsave("Figure_2D_16S.pdf", fig2D, units = "mm", width = 85, height = 50, dpi = 300)
