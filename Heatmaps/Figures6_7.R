## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Plot Heatmap Figures 6 and 7
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(MicEco)
library(MicrobiotaProcess)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(factoextra)
library(tibble)
library(RColorBrewer)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

#Obtain ASVs shared between at least two compartments with Venn lists
vennlist=ps_venn(ps_filtered_norm, group = "compartment", plot = FALSE)
vennlist=ldply(vennlist, data.frame) # Convert list to data frame
colnames(vennlist) <- c("dataset", "ASV") # Rename columns

# Filter ASVs
vennlist_heatmap<- vennlist[!vennlist$dataset %in% c("Bulk soil", "Rhizosphere", "Root", "Leaf", "Grain"),]
ASV_to_keep <- vennlist_heatmap$ASV  # Store ASVs to keep
#str(ASV_to_keep)

#Aggregate data by compartment and irrigation treatment
ps_merged = merge_samples(ps_filtered_norm, "compartment_irrigation")
ps_mg_rel_abund= phyloseq::transform_sample_counts(ps_merged, function(x){x / sum(x)}) # Normalize by relative abundance
#phyloseq::otu_table(ps_mg_rel_abund)[1:5, 1:5]

df_tax_table_mg_rel_abund <- as.data.frame(phyloseq::tax_table (ps_mg_rel_abund)) # Extract taxonomy data
df_otu_table_mg_rel_abund <- as.data.frame(phyloseq::otu_table(ps_mg_rel_abund)) # Extract OTU table
df_heatmap<-as.data.frame(t(df_otu_table_mg_rel_abund)) # Transpose OTU table
#str(df_heatmap)

#Select the ASVs common for at least two compartments in the otu_table
df_heatmap_ASV_to_keep <- subset(df_heatmap, rownames(df_heatmap) %in% ASV_to_keep)
# colnames(df_heatmap_ASV_to_keep) = c("Bulk soil_LW", "Bulk soil_OW", "Grain_LW", "Grain_OW", "Leaf_LW", "Leaf_OW", "Rhizosphere_LW", "Rhizosphere_OW", "Root_LW", "Root_OW")
# colnames(df_heatmap_ASV_to_keep)

# Convert selected ASVs into a matrix for clustering analysis
set.seed(101)
matriz <- as.matrix(df_heatmap_ASV_to_keep)
matriz_p <- matriz  %>% as_tibble(rownames = "ASV")

# Extract unique ASVs
unique_ASV <- matriz_p %>% 
  pull(ASV) %>%             # extract the ASV column as a vector
  unique()                   # retain only unique values

# Create a matrix with ASV abundance data
hclust_matrix <- matriz_p %>% dplyr::select(-ASV) %>% as.matrix()
rownames(hclust_matrix) <- matriz_p$ASV # Assign ASV names as rownames
hclust_matrix <- hclust_matrix[unique_ASV, ]

# Normalize (scale) the data for clustering
hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so ASV are as columns
  t() %>% 
  # apply scalling to each column of the matrix (ASVs)
  scale() %>% 
  # transpose back so ASVs are as rows again
  t()

# Determine the optimal number of clusters using the elbow method
#SSE:
#The first measure is using the sum of squared error (SSE). 
#SSE is defined as the sum of the squared distance between each member of a cluster and 
#its cluster centroid. We repeatedly test and increase number of clusters and evaluate the SSE. 
#As we increase the number of clusters the distance between any point and it’s centroid will be 
#smaller since the cluster itself is smaller. At a certain number of clusters number however, 
#the SSE will not significantly decrease with each new addition of a cluster. 
#This is the elbow and suggests a suitable number of clusters:

wss <- (nrow(hclust_matrix)-1)*sum(apply(hclust_matrix,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(hclust_matrix,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# Perform hierarchical k-means clustering
hkmeans_cluster <- hkmeans(x = hclust_matrix, hc.metric = "euclidean",
                           hc.method = "complete", k = 10) #K is the number of clusters (usually obtained with the elbow method)
hkmeans_cluster <- hkmeans_cluster$cluster %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(ASV = name, cluster = value)
head(hkmeans_cluster)

# Convert clustering results to a data frame
ASV_cluster_df <- as.data.frame(hkmeans_cluster)
#colnames(ASV_cluster_df) = c("ASV","Cluster")

# Merge taxonomic and clustering data
hclust_df <- as.data.frame(hclust_matrix)
hclust_df = hclust_df %>%
  rownames_to_column(var="ASV")
df_tax_table_ASV_to_keep <- subset(df_tax_table_mg_rel_abund, rownames(df_tax_table_mg_rel_abund) %in% ASV_to_keep)
df_tax_table_ASV_to_keep = df_tax_table_ASV_to_keep%>%
  rownames_to_column(var="ASV")
hclust_df <-merge(hclust_df,df_tax_table_ASV_to_keep, by="ASV")
hclust_df <-merge(hclust_df,ASV_cluster_df, by="ASV")

##Formatting tax names
# hclust_df$Kingdom = gsub(x = hclust_df$Kingdom , pattern = "k__", replacement = "")
# hclust_df$Phylum = gsub(x = hclust_df$Phylum , pattern = "p__", replacement = "")
# hclust_df$Class = gsub(x = hclust_df$Class , pattern = "c__", replacement = "")
# hclust_df$Order = gsub(x = hclust_df$Order , pattern = "o__", replacement = "")
# hclust_df$Family = gsub(x = hclust_df$Family , pattern = "f__", replacement = "")
# hclust_df$Genus = gsub(x = hclust_df$Genus , pattern = "g__", replacement = "")
# hclust_df$Species = gsub(x = hclust_df$Species , pattern = "s__", replacement = "")

# Replace missing taxonomic values with "Unknown"
hclust_df$Family <-replace(hclust_df$Family,is.na(hclust_df$Family),"Unknown")
hclust_df$Genus <- replace(hclust_df$Genus,is.na(hclust_df$Genus),"Unknown")
hclust_df$Species <- replace(hclust_df$Species,is.na(hclust_df$Species),"Unknown")
hclust_df = hclust_df%>%
  column_to_rownames(var="ASV")

# Define compartment and irrigation labels for annotations
comp = gsub("_OW", "", colnames(hclust_matrix))
comp = gsub("_LW", "", comp)

irr = gsub("Bulk soil_", "", colnames(hclust_matrix))
irr = gsub("Rhizosphere_", "", irr)
irr = gsub("Root_", "", irr)
irr = gsub("Leaf_", "", irr)
irr = gsub("Grain_", "", irr)

# Define color schemes for heatmap annotations
col = list(Comp = c("Bulk soil" = "#CC9999", 
                    "Rhizosphere" = "#990033", 
                    "Root" = "#33CCFF", 
                    "Leaf" = "#66CC00", 
                    "Grain" = "#FFCC33"),
           Irr = c("OW" = "#4a80ff",
                   "LW" = "tan2") )

# Create heatmap annotation
ha = HeatmapAnnotation(
  df = data.frame(Comp = comp, Irr = irr),
  simple_anno_size = unit(0.2, "cm"),
  annotation_name_side = "left",
  col = col,
  gp = gpar(col = "#333333"),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    Irr = list(at = c("OW", "LW")),
    Comp = list(at = c("Bulk soil", 
                       "Rhizosphere", 
                       "Root", 
                       "Leaf", 
                       "Grain"))
  ))

# Create the heatmap annotation

hm=Heatmap(matrix = hclust_matrix, name = "Z-score",
           clustering_distance_columns = "euclidean",
           clustering_distance_rows = "euclidean",
           clustering_method_columns = "complete",
           clustering_method_rows = "complete",
           split = hkmeans_cluster$cluster, 
           show_row_names = FALSE,
           show_column_names = FALSE, 
           top_annotation = ha)
#Save Figure as png
png("Figure_6.png", width = 210, height = 180, units ="mm", res = 300)  # Convert mm to inches
hm
dev.off()
