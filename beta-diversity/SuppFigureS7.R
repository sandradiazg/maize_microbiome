## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Calculate Jaccard distances between samples
## - Create dendrogram and plot Supplementary Figure S7
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(ggplot2)
library(ggdendro)  # Visualize dendrograms
library(dplyr)
library(RColorBrewer)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

# Merge samples by compartment and time point for dendrogram analysis
ps_merged = merge_samples(ps_filtered_norm, "compartment_tpreal")
# Calculate Jaccard distances (IMPORTANT: binary = TRUE, ensures presence/absence data is used)
jaccard_dist <- phyloseq::distance(ps_merged, method="jaccard", binary = T)
# #jaccard_similarity <- 1 - jaccard_dist #Jaccard similarity = 1 - Jaccard distance. Here we work with distances. It is important to take into account to interpretate the results.
# Perform hierarchical clustering using UPGMA (average linkage method)
hc <- hclust(jaccard_dist, method = "average")  # Method UPGMA (average)
plot(hc, main = "Dendrograma basado en distancias Jaccard", cex = 0.8)
# Cut dendrogram into 4 clusters
clusters <- cutree(hc, k = 4)  # Divide into 4 groups
table(clusters)   # Display the number of elements in each cluster

# Convert dendrogram data to ggplot-compatible format
dend_data <- dendro_data(hc)
# Add medata information to labels_df
labels_df <- data.frame(label = hc$labels)
labels_df <- labels_df %>%
  mutate(Time_point = case_when(
    endsWith(label, "T01") ~ "T01",
    endsWith(label, "T02") ~ "T02",
    endsWith(label, "T03") ~ "T03"
  )) %>%
  mutate(Compartment = case_when(
    startsWith(label, "Bulk soil") ~ "Bulk soil",
    startsWith(label, "Rhizosphere") ~ "Rhizosphere",
    startsWith(label, "Root") ~ "Root",
    startsWith(label, "Leaf") ~ "Leaf",
    startsWith(label, "Grain") ~ "Grain"
  ))
# Define colors and shapes for visualization
cols_comp <- c("Bulk soil" = "#CC9999", 
               "Rhizosphere" = "#990033", 
               "Root" = "#33CCFF", 
               "Leaf" = "#66CC00", 
               "Grain" = "#FFCC33")

tp_cols <- c("T01" = "#66cccc",
             "T02" = "#3399CC",
             "T03" = "#003333")

tp_shape2 <- c("T01" = 15,
               "T02" = 16,
               "T03" = 17)

# # Convert metadata columns to factors
labels_df$Compartment <- factor(labels_df$Compartment)
labels_df$Time_point <- factor(labels_df$Time_point)
labels_df$Cluster <- as.factor(clusters[match(labels_df$label, names(clusters))])

# Merge dendrogram data with metadata
labels_merged <- merge(dend_data$labels, labels_df, by.x = "label", by.y = "label")

# Add Time_point to dendrogram segments
dend_data$segments$Time_point <- labels_merged$Time_point[match(dend_data$segments$xend, labels_merged$x)]


# labels_merged <- labels_merged %>%
#   mutate(label = gsub("BULK_SOIL", "Bulk soil", label),
#          label = gsub("RIZOSPHERE", "Rhizosphere", label),
#          label = gsub("ROOT", "Root", label),
#          label = gsub("LEAF", "Leaf", label),
#          label = gsub("GRAIN", "Grain", label))


jaccard_dend = ggplot() +
  # Draw dendrogram lines, shorthening segments
  geom_segment(data = dend_data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               size = 0.5) +  
  
  # Add shapes by time point and the end of the segments
  geom_point(data = labels_merged, 
             aes(x = x, y = y, 
                 color = as.factor(Time_point), 
                 shape = as.factor(Time_point)), 
             size = 2.5, stroke = 0.5) +
  
  # Add text
  geom_text(data = labels_merged, 
            aes(x = x, y = y - 0.03 ,  
                label = label, 
                color = as.factor(Compartment)), 
            angle = 0, hjust = 1, size = 2) +  
  
  # Color and shape scales
  scale_color_manual(values = tp_cols) +  
  scale_shape_manual(values = tp_shape2) +  
  guides(color = FALSE, 
         shape = FALSE) +  
  
  # Expand graph
  expand_limits(y = min(dend_data$segments$y) - 0.38) + 
  coord_cartesian(clip = "off") +  #
  coord_flip() +  # Flips the graph
  #Plot's theme
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.1),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black", size = 0.1),
        axis.text.x = element_text(size = 5, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5, 5, 5, 5)) 

#Save figure
ggsave("SFigure_S7_16S.pdf", jaccard_dend, units = "mm", width = 85, height = 70, dpi = 300)