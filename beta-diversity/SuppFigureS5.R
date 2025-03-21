## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Calculate Bray-Curtis dissimilarity 
## - Plot distances in Supplementary Figure S5
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(ggplot2)
library(dplyr)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

#Calculate Bra-Curtis distances between samples 
bray_dist = phyloseq::distance(ps_filtered_norm, method="bray", weighted=T)
bray_matrix <- as.matrix(bray_dist) # Convert to matrix format
bray_long <- as.data.frame(as.table(as.matrix(bray_dist)))# Convert to long format
colnames(bray_long) <- c("Sample1", "Sample2", "Distance")

# Remove self-comparisons (i.e., same sample comparisons)
bray_long <- bray_long[bray_long$Sample1 != bray_long$Sample2, ]

# Extract metadata for each sample (compartment, irrigation, time point)
bray_long$Group1 <- sample_data(ps_filtered_norm)$compartment_irrigation_tpreal[match(bray_long$Sample1, rownames(sample_data(ps_filtered_norm)))]
bray_long$Group2 <- sample_data(ps_filtered_norm)$compartment_irrigation_tpreal[match(bray_long$Sample2, rownames(sample_data(ps_filtered_norm)))]

# Extract time point (T01, T02, T03) and compartment from metadata
bray_long<-bray_long %>%
  mutate(tpreal1 = case_when(
    endsWith(Group1, "T01") ~ "T01",
    endsWith(Group1, "T02") ~ "T02",
    endsWith(Group1, "T03") ~ "T03"
  ))%>%
  mutate(tpreal2 = case_when(
    endsWith(Group2, "T01") ~ "T01",
    endsWith(Group2, "T02") ~ "T02",
    endsWith(Group2, "T03") ~ "T03"
  ))%>%
  mutate(compartment1 = case_when(
    startsWith(Group1, "Bulk soil") ~ "Bulk soil",
    startsWith(Group1, "Rhizosphere") ~ "Rhizosphere",
    startsWith(Group1, "Root") ~ "Root",
    startsWith(Group1, "Leaf") ~ "Leaf",
    startsWith(Group1, "Grain") ~ "Grain",
  )) %>%
  mutate(compartment2 = case_when(
    startsWith(Group2, "Bulk soil") ~ "Bulk soil",
    startsWith(Group2, "Rhizosphere") ~ "Rhizosphere",
    startsWith(Group2, "Root") ~ "Root",
    startsWith(Group2, "Leaf") ~ "Leaf",
    startsWith(Group2, "Grain") ~ "Grain",
  )) %>%
  mutate(irrigation1 = case_when(
    str_detect(Group1, "OW") ~ "OW",
    str_detect(Group1, "LW") ~ "LW",
  )) %>%
  mutate(irrigation2 = case_when(
    str_detect(Group2, "OW") ~ "OW",
    str_detect(Group2, "LW") ~ "LW",
  ))

# Define compartment order
compartment_order <- c("Bulk soil", "Rhizosphere", "Root", "Leaf", "Grain")

# Filter pairwise comparisons excluding "Grain" and keeping same time points
bray_dist_main <- bray_long %>%
  mutate(
    compartment1 = factor(compartment1, levels = compartment_order),
    compartment2 = factor(compartment2, levels = compartment_order)
  ) %>%
  filter(compartment1 != compartment2,  # Exclude self-comparisons
         !compartment1 %in% "Grain",   # Exclude "Grain" comparisons
         !compartment2 %in% "Grain",   # Exclude "Grain" comparisons
         tpreal1 == tpreal2) %>%       # Keep only same time point comparisons
  mutate(compartment_pair = ifelse(as.numeric(compartment1) < as.numeric(compartment2),
                                   paste(compartment1, compartment2, sep = " vs. "),
                                   paste(compartment2, compartment1, sep = " vs. "))) %>%
  group_by(tpreal = tpreal1, 
           compartment_pair) %>%
  dplyr::summarise(mean_distance = mean(Distance),
                   sd_distance = sd(Distance),
                   n = n(), .groups = "drop")

# Set factor levels for plotting
bray_dist_main$compartment_pair <- factor(bray_dist_main$compartment_pair, 
                                          levels = c("Bulk soil vs. Rhizosphere",
                                                     "Bulk soil vs. Root",
                                                     "Bulk soil vs. Leaf",
                                                     "Rhizosphere vs. Root",
                                                     "Rhizosphere vs. Leaf",
                                                     "Root vs. Leaf"))

#Supplementary Figure S5A
#Plot pairwise comparisons of Bray-Curtis distance between bulk soil, rhizosphere, root and leaf compartments over time.
plot_main <- ggplot(bray_dist_main, 
                    aes(x = tpreal, y = mean_distance, color = compartment_pair, group = compartment_pair)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = mean_distance - sd_distance, 
                    ymax = mean_distance + sd_distance), width = 0.2, size = 0.3) +
  facet_wrap(~ compartment_pair, scales = "free_y") +
  labs(x = "Time point", 
       y = "Bray-Curtis dissimilarity") +
  ylim(c(0.25,0.85))+
  guides(color = FALSE, fill = FALSE)+
  theme(            panel.border = element_rect(color = "black", fill=NA, size = 0.5),
                    panel.background = element_blank(),
                    panel.grid = element_line(color = "grey",size = 0.2),
                    axis.line=element_line(color="black", size=0.2),
                    axis.ticks=element_line(color="black", size=0.2),
                    axis.text.y = element_text(size = 7, color = "black"),
                    axis.text.x = element_text(size = 7, color = "black"),
                    axis.title.y = element_text(size = 7, color = "black"),
                    axis.title.x = element_text(size = 7, color = "black"),
                    strip.text.x = element_text(size = 5.8, color = "black", face = "bold"),
                    strip.background = element_rect(colour=NA, fill=NA))


#Save Figure
ggsave("SFigure_S5A_16S.pdf", plot_main, units = "mm", width = 100, height = 60, dpi = 300)

#Supplementary Figure S5B
# Filter comparisons involving "Grain" at time point T03
bray_dist_grain <- bray_long %>%
  filter((compartment1 == "Grain" | compartment2 == "Grain") &  # Al menos uno es "Grain"
           !(compartment1 == "Grain" & compartment2 == "Grain") & # No ambos "Grain"
           tpreal1 == "T03" & tpreal2 == "T03") %>%  # Mismo tiempo
  mutate(
    comparison = ifelse(compartment1 == "Grain", 
                        paste(compartment2, "vs. Grain"),
                        paste(compartment1, "vs. Grain"))
  )

# Set factor levels for plotting
bray_dist_grain$comparison <- factor(bray_dist_grain$comparison, 
                                     levels = c("Bulk soil vs. Grain",
                                                "Rhizosphere vs. Grain",
                                                "Root vs. Grain",
                                                "Leaf vs. Grain"))

#Plot pairwise comparisons of Bray-Curtis distance between grain and the rest of compartments (bulk soil, rhizosphere, root and leaf) at T03
plot_grain <- ggplot(bray_dist_grain, aes(x = comparison, y = Distance, fill = comparison)) +
  geom_boxplot(lwd=0.5, alpha=0.4, size =4, width=0.8, 
               outlier.shape=NA, 
               color = "#CD3278", 
               fill= "#CD3278")+ 
  #geom_jitter(stroke=0.4, alpha=0.5, size =1, color = "#CD3278", fill= "#CD3278")+
  labs(y = "Bray-Curtis dissimilarity")+
  guides(color = FALSE, fill = FALSE)+
  stat_summary( fun.data="mean_sdl",
                fun.args = list(mult=1), 
                geom = "point",  
                size = 1.8, shape = 16, 
                position = position_dodge(.5))+
  theme(panel.background=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black", size=0.3),
        axis.ticks.x = element_line(color="black", size=0.3),
        axis.ticks.y = element_line(color="black", size=0.3),
        axis.text.x = element_text(size = 6, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank())

#Save figure
ggsave("SFigure_S5B_16S.pdf", plot_grain, units = "mm", width = 70, height = 60, dpi = 300)