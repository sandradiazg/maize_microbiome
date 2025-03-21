## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and OW availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Calculate relative abundances under different conditions (Supplementary Tables S4 and S5)
## This script is for 16S dataset, but the same code was used for ITS dataset using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(tibble)
library(openxlsx)
library(phyloseqCompanion)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

## ----Pre-processing of phyloseq object, if needed ----

# df_tax<-as.data.frame(phyloseq::tax_table(ps_filtered_norm))

# #Remove taxonomic prefixes from different levels
# df_tax$Genus <- gsub("g__", "",df_tax$Genus)
# df_tax$Family <- gsub("f__", "",df_tax$Family)
# df_tax$Order <- gsub("o__", "",df_tax$Order)
# df_tax$Class <- gsub("c__", "",df_tax$Class)
# df_tax$Phylum <- gsub("p__", "",df_tax$Phylum)
# df_tax$Kingdom <- gsub("k__", "",df_tax$Kingdom)
# #df_tax<-replace(df_tax,is.na(df_tax_plot_taxa),"Unknown")
# view(df_tax)

# #Replace missing Family and Genus names with higher taxonomic levels
# df_tax$Family<-ifelse(is.na(df_tax$Family),str_c("Unk_Ord_", df_tax$Order),df_tax$Family)
# df_tax$Genus<-ifelse(is.na(df_tax$Genus),str_c("Unk_Fam_", df_tax$Family),df_tax$Genus)

# # Create modified Tax_table
# mat_taxa <- as.matrix (df_tax)
# TAX = phyloseq::tax_table(mat_taxa)
# #Creat modified ps object
# ps_filtered_norm <- phyloseq(TAX, otu_table(ps_filtered_norm), sample_data(ps_filtered_norm))


# Aggregate data at the Genus level
dat.aglo = tax_glom(ps_filtered_norm, taxrank = "Genus", NArm=FALSE )

## ---- Aggregation by compartment ----

# Merge samples by compartment and calculate relative abundances
ps_mg_plot = merge_samples(dat.aglo, "compartment")
ps_rel_abund_mg_plot = phyloseq::transform_sample_counts(ps_mg_plot, function(x){x / sum(x)})

# Process metadata table
df_metadata_mg_plot <- sample.data.frame(ps_rel_abund_mg_plot)
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  rownames_to_column (var="rnam")
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  mutate(compartment = rnam)
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  column_to_rownames(var= "rnam")

# Retain only relevant columns
df_metadata_mg_plot <- df_metadata_mg_plot[,names(df_metadata_mg_plot) %in% c("plant_number","compartment")]
col_order <- c("plant_number","compartment")
df_metadata_mg_plot <- df_metadata_mg_plot[, col_order]

# Create a phyloseq object with updated metadata
metadata_plot_taxa <- sample_data(df_metadata_mg_plot)
ps_plot_taxa <- phyloseq(phyloseq::tax_table(ps_rel_abund_mg_plot), otu_table(ps_rel_abund_mg_plot), metadata_plot_taxa)

# Extract abundance data and refine columns
abundance_table <- ps_plot_taxa %>%
  psmelt()
abundance_table <- abundance_table[,names(abundance_table) %in% c("Abundance","compartment", "Kingdom","Phylum","Class","Order","Family","Genus")]
col_order <- c("Kingdom","Phylum","Class","Order","Family","Genus","Abundance","compartment")
abundance_table1 <- abundance_table[, col_order]

## ---- Aggregation by compartment and time point ----

# Merge samples by compartment and time point
ps_mg_plot = merge_samples(dat.aglo, "compartment_tpreal")
ps_rel_abund_mg_plot = phyloseq::transform_sample_counts(ps_mg_plot, function(x){x / sum(x)})

# Process metadata table
df_metadata_mg_plot <- sample.data.frame(ps_rel_abund_mg_plot)
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  rownames_to_column (var="rnam")
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  mutate(compartment_tpreal = rnam)
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  column_to_rownames(var= "rnam")

# Extract relevant columns
df_metadata_mg_plot <- df_metadata_mg_plot[,names(df_metadata_mg_plot) %in% c("plant_number","compartment_tpreal")]
col_order <- c("plant_number","compartment_tpreal")
df_metadata_mg_plot <- df_metadata_mg_plot[, col_order]

# Assign time points based on naming patterns
df_metadata_mg_plot<-df_metadata_mg_plot %>%
  mutate(tpreal = case_when(
    endsWith(compartment_tpreal, "T01") ~ "T01",
    endsWith(compartment_tpreal, "T02") ~ "T02",
    endsWith(compartment_tpreal, "T03") ~ "T03"
  ))%>%
  mutate(compartment = case_when(
    startsWith(compartment_tpreal, "Bulk soil") ~ "Bulk soil",
    startsWith(compartment_tpreal, "Rhizosphere") ~ "Rhizosphere",
    startsWith(compartment_tpreal, "Root") ~ "Root",
    startsWith(compartment_tpreal, "Leaf") ~ "Leaf",
    startsWith(compartment_tpreal, "Grain") ~ "Grain",
  )) 

# Create updated phyloseq object
metadata_plot_taxa <- sample_data(df_metadata_mg_plot)
ps_plot_taxa <- phyloseq(phyloseq::tax_table(ps_rel_abund_mg_plot), otu_table(ps_rel_abund_mg_plot), metadata_plot_taxa)

# Extract abundance data and refine columns
abundance_table <- ps_plot_taxa %>%
  psmelt()
abundance_table <- abundance_table[,names(abundance_table) %in% c("Abundance","compartment","tpreal", "Kingdom","Phylum","Class","Order","Family","Genus")]
col_order <- c("Kingdom","Phylum","Class","Order","Family","Genus","Abundance","compartment", "tpreal")
abundance_table2 <- abundance_table[, col_order]

## ---- Aggregation by compartment, irrigation, and time point ----

# Merge samples by compartment, irrigation, and time point
ps_mg_plot = merge_samples(dat.aglo, "compartment_irrigation_tpreal")
ps_rel_abund_mg_plot = phyloseq::transform_sample_counts(ps_mg_plot, function(x){x / sum(x)})

# Process metadata table
df_metadata_mg_plot <- sample.data.frame(ps_rel_abund_mg_plot)
str(df_metadata_mg_plot)
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  rownames_to_column (var="rnam")
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  mutate(compartment_irrigation_tpreal = rnam)
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  column_to_rownames(var= "rnam")

# Extract relevant columns
df_metadata_mg_plot <- df_metadata_mg_plot[,names(df_metadata_mg_plot) %in% c("irrigation","compartment_irrigation_tpreal")]
col_order <- c("compartment_irrigation_tpreal","irrigation")
df_metadata_mg_plot <- df_metadata_mg_plot[, col_order]

# Assign time points, compartments, and irrigation types
df_metadata_mg_plot<-df_metadata_mg_plot %>%
  mutate(tpreal = case_when(
    endsWith(compartment_irrigation_tpreal, "T01") ~ "T01",
    endsWith(compartment_irrigation_tpreal, "T02") ~ "T02",
    endsWith(compartment_irrigation_tpreal, "T03") ~ "T03"
  ))%>%
  mutate(compartment = case_when(
    startsWith(compartment_irrigation_tpreal, "Bulk soil") ~ "Bulk soil",
    startsWith(compartment_irrigation_tpreal, "Rhizosphere") ~ "Rhizosphere",
    startsWith(compartment_irrigation_tpreal, "Root") ~ "Root",
    startsWith(compartment_irrigation_tpreal, "Leaf") ~ "Leaf",
    startsWith(compartment_irrigation_tpreal, "Grain") ~ "Grain",
  )) %>%
  mutate(irrigation = case_when(
    str_detect(compartment_irrigation_tpreal, "OW") ~ "OW",
    str_detect(compartment_irrigation_tpreal, "LW") ~ "LW"
  ))


# Create final phyloseq object
metadata_plot_taxa <- sample_data(df_metadata_mg_plot)
ps_plot_taxa <- phyloseq(phyloseq::tax_table(ps_rel_abund_mg_plot), otu_table(ps_rel_abund_mg_plot), metadata_plot_taxa)

# Extract abundance data and refine columns
abundance_table <- ps_plot_taxa %>%
  psmelt()
abundance_table <- abundance_table[,names(abundance_table) %in% c("Abundance","compartment","irrigation", "tpreal", "Kingdom","Phylum","Class","Order","Family","Genus")]
col_order <- c("Kingdom","Phylum","Class","Order","Family","Genus","Abundance","compartment","irrigation","tpreal")
abundance_table3 <- abundance_table[, col_order]

# Save data tables to an Excel file
dataset_names <- list('Table S4.1' = abundance_table1, 'Table S4.2' = abundance_table2, 'Table S4.3' = abundance_table3)

openxlsx::write.xlsx(dataset_names, file = "Supplementary Table S4_Taxonomy_16S.xlsx")
