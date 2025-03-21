## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Conduct ALDEx2 differential abundance analysis for each compartment
## - Obtain a list of differentially abundant taxa (Supplementary Table S9)


#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(openxlsx)
library(dplyr)
library(tibble)
library(stringr)

##ALDEx2 analysis on 16S dataset
#Load ps object
filtered_norm_16S <- readRDS("ps_filtered_norm_16S.rds")

#subsets by compartment and time point
ps_Bulk_soil_T01 <- subset_samples(filtered_norm_16S, compartment== "Bulk soil" & tpreal == "T01")
ps_Rhizosphere_T01 <- subset_samples(filtered_norm_16S, compartment== "Rhizosphere" & tpreal == "T01")
ps_Root_T01 <- subset_samples(filtered_norm_16S, compartment== "Root" & tpreal == "T01")
ps_Leaf_T01 <- subset_samples(filtered_norm_16S, compartment== "Leaf" & tpreal == "T01")

ps_Bulk_soil_T02 <- subset_samples(filtered_norm_16S, compartment== "Bulk soil" & tpreal == "T02")
ps_Rhizosphere_T02 <- subset_samples(filtered_norm_16S, compartment== "Rhizosphere" &  tpreal == "T02")
ps_Root_T02 <- subset_samples(filtered_norm_16S, compartment== "Root" & tpreal == "T02")
ps_Leaf_T02 <- subset_samples(filtered_norm_16S, compartment== "Leaf" & tpreal == "T02")

ps_Bulk_soil_T03 <- subset_samples(filtered_norm_16S, compartment== "Bulk soil" & tpreal == "T03")
ps_Rhizosphere_T03 <- subset_samples(filtered_norm_16S, compartment== "Rhizosphere" &  tpreal == "T03")
ps_Root_T03 <- subset_samples(filtered_norm_16S, compartment== "Root" & tpreal == "T03")
ps_Leaf_T03 <- subset_samples(filtered_norm_16S, compartment== "Leaf" & tpreal == "T03")
ps_Grain_T03 <- subset_samples(filtered_norm_16S, compartment== "Grain")
# Define values for compartment and time point
valid_compartments <- c("Bulk_soil", "Rhizosphere", "Root", "Leaf", "Grain")
valid_tpreal <- c("T01", "T02", "T03")

#Identify all phyloseq objects that start with "ps_"
phyloseq_objects <- ls(pattern = "^ps_")

# Create a list to store results from ALDEx2 analysis
all_results <- list()

# Iterate through all phyloseq objects to perform ALDEx2 analysis
for (obj_name in phyloseq_objects) {
  cat("Procesando:", obj_name, "\n")
  
  # Extract compartment and time point from object name using regex
  compartment <- ifelse(any(grepl(paste(valid_compartments, collapse = "|"), obj_name)),
                        gsub(paste0(".*(", paste(valid_compartments, collapse = "|"), ").*"), "\\1", obj_name),
                        "N/A")
  
  tpreal <- ifelse(any(grepl(paste(valid_tpreal, collapse = "|"), obj_name)),
                   gsub(paste0(".*(", paste(valid_tpreal, collapse = "|"), ").*"), "\\1", obj_name),
                   "N/A")
  
  # Convert object name to actual phyloseq object reference
  obj_ps <- get(obj_name)
  
  # Aggregate ASV counts at the Genus level
  genus_ps <- tax_glom(obj_ps, taxrank = "Genus", NArm = FALSE)
  
  # Check the number of unique irrigation levels
  irrigation_levels <- unique(as.character(phyloseq::sample_data(genus_ps)$irrigation))
  
  if (length(irrigation_levels) != 2) {
    cat("Skipping:", obj_name, "- 'irrigation' does not have exactly two unique levels\n")
    next
  }
  
  # Perform ALDEx2 differential abundance analysis
  aldex2_res <- ALDEx2::aldex(
    data.frame(phyloseq::otu_table(genus_ps)),
    as.character(phyloseq::sample_data(genus_ps)$irrigation),
    test = "t", effect = TRUE, denom = "iqlr"
  )
  
  # Extract taxonomic information
  taxa_info <- data.frame(phyloseq::tax_table(genus_ps)) %>%
    rownames_to_column(var = "OTU")
  # Filter significant taxa based on statistical thresholds
  sig_aldex2 <- aldex2_res %>%
    rownames_to_column(var = "OTU") %>%
    filter((we.eBH < 0.05 | wi.eBH < 0.05) & abs(effect) > 1) %>%
    arrange(effect, we.eBH)
  
  sig_aldex2 <- left_join(sig_aldex2, taxa_info)
  
  # Add metadata columns (phyloseq object name, compartment, and time point)
  sig_aldex2 <- sig_aldex2 %>%
    mutate(object_name = obj_name,
           compartment = compartment,
           tpreal = tpreal)
  
  # Store significant results in the list
  all_results[[obj_name]] <- sig_aldex2
}

# Combine results from all compartments into a single data frame
final_results_16S <- bind_rows(all_results, .id = "source_object")
final_results_16S <- subset(final_results_16S, select = -Species )
# Replace missing Family and Genus values with 'Unknown'
final_results_16S$Family<-ifelse(is.na(final_results_16S$Family),str_c("Unknown_", final_results_16S$OTU),final_results_16S$Family)
final_results_16S$Genus<-ifelse(is.na(final_results_16S$Genus),str_c("Unknown_", final_results_16S$OTU),final_results_16S$Genus)

##Repeat the process with the ITS dataset
#Load ps object
filtered_norm_ITS <- readRDS("ps_filtered_norm_ITS.rds")

#subsets by compartment and time point
ps_Bulk_soil_T01 <- subset_samples(filtered_norm_ITS, compartment== "Bulk soil" & tpreal == "T01")
ps_Rhizosphere_T01 <- subset_samples(filtered_norm_ITS, compartment== "Rhizosphere" & tpreal == "T01")
ps_Root_T01 <- subset_samples(filtered_norm_ITS, compartment== "Root" & tpreal == "T01")
ps_Leaf_T01 <- subset_samples(filtered_norm_ITS, compartment== "Leaf" & tpreal == "T01")

ps_Bulk_soil_T02 <- subset_samples(filtered_norm_ITS, compartment== "Bulk soil" & tpreal == "T02")
ps_Rhizosphere_T02 <- subset_samples(filtered_norm_ITS, compartment== "Rhizosphere" &  tpreal == "T02")
ps_Root_T02 <- subset_samples(filtered_norm_ITS, compartment== "Root" & tpreal == "T02")
ps_Leaf_T02 <- subset_samples(filtered_norm_ITS, compartment== "Leaf" & tpreal == "T02")

ps_Bulk_soil_T03 <- subset_samples(filtered_norm_ITS, compartment== "Bulk soil" & tpreal == "T03")
ps_Rhizosphere_T03 <- subset_samples(filtered_norm_ITS, compartment== "Rhizosphere" &  tpreal == "T03")
ps_Root_T03 <- subset_samples(filtered_norm_ITS, compartment== "Root" & tpreal == "T03")
ps_Leaf_T03 <- subset_samples(filtered_norm_ITS, compartment== "Leaf" & tpreal == "T03")
ps_Grain_T03 <- subset_samples(filtered_norm_ITS, compartment== "Grain")

# Define valid values for compartment and time point
valid_compartments <- c("Bulk_soil", "Rhizosphere", "Root", "Leaf", "Grain")
valid_tpreal <- c("T01", "T02", "T03")

# Identify all phyloseq objects that start with "ps_"
phyloseq_objects <- ls(pattern = "^ps_")

# Create a list to store results from ALDEx2 analysis
all_results <- list()

# Iterate through all phyloseq objects to perform ALDEx2 analysis
for (obj_name in phyloseq_objects) {
  cat("Procesando:", obj_name, "\n")
  
  # Extract compartment and time point from object name using regex
  compartment <- ifelse(any(grepl(paste(valid_compartments, collapse = "|"), obj_name)),
                        gsub(paste0(".*(", paste(valid_compartments, collapse = "|"), ").*"), "\\1", obj_name),
                        "N/A")
  
  tpreal <- ifelse(any(grepl(paste(valid_tpreal, collapse = "|"), obj_name)),
                   gsub(paste0(".*(", paste(valid_tpreal, collapse = "|"), ").*"), "\\1", obj_name),
                   "N/A")
  
  # Convert object name to actual phyloseq object reference
  obj_ps <- get(obj_name)
  
  # Aggregate ASV counts at the Genus level
  genus_ps <- tax_glom(obj_ps, taxrank = "Genus", NArm = FALSE)
  
  # Check the number of unique irrigation levels
  irrigation_levels <- unique(as.character(phyloseq::sample_data(genus_ps)$irrigation))
  
  if (length(irrigation_levels) != 2) {
    cat("Skipping:", obj_name, "- 'irrigation' does not have exactly two unique levels\n")
    next
  }
  
  # Perform ALDEx2 differential abundance analysis
  aldex2_res <- ALDEx2::aldex(
    data.frame(phyloseq::otu_table(genus_ps)),
    as.character(phyloseq::sample_data(genus_ps)$irrigation),
    test = "t", effect = TRUE, denom = "iqlr"
  )
  
  # Extract taxonomic information
  taxa_info <- data.frame(phyloseq::tax_table(genus_ps)) %>%
    rownames_to_column(var = "OTU")
  
  # Filter significant taxa based on statistical thresholds
  sig_aldex2 <- aldex2_res %>%
    rownames_to_column(var = "OTU") %>%
    filter((we.eBH < 0.05 | wi.eBH < 0.05) & abs(effect) > 1) %>%
    arrange(effect, we.eBH)
  
  sig_aldex2 <- left_join(sig_aldex2, taxa_info)
  
  # Add metadata columns (phyloseq object name, compartment, and time point)
  sig_aldex2 <- sig_aldex2 %>%
    mutate(object_name = obj_name,
           compartment = compartment,
           tpreal = tpreal)
  
  # Store significant results in the list
  all_results[[obj_name]] <- sig_aldex2
}

# Combine results from all compartments into a single data frame
final_results_ITS <- bind_rows(all_results, .id = "source_object")
final_results_ITS <- subset(final_results_ITS, select = -Species )

## Removing the prefix "g__", "f__", "o__", "c__", "p__", "k__" from taxonomic ranks for clarity
final_results_ITS$Genus <- gsub("g__", "",final_results_ITS$Genus)
final_results_ITS$Family <- gsub("f__", "",final_results_ITS$Family)
final_results_ITS$Order <- gsub("o__", "",final_results_ITS$Order)
final_results_ITS$Class <- gsub("c__", "",final_results_ITS$Class)
final_results_ITS$Phylum <- gsub("p__", "",final_results_ITS$Phylum)
final_results_ITS$Kingdom <- gsub("k__", "",final_results_ITS$Kingdom)

# Replace missing Family and Genus values with 'Unknown'
final_results_ITS$Family<-ifelse(is.na(final_results_ITS$Family),str_c("Unknown_", final_results_ITS$OTU),final_results_ITS$Family)
final_results_ITS$Genus<-ifelse(is.na(final_results_ITS$Genus),str_c("Unknown_", final_results_ITS$OTU),final_results_ITS$Genus)

#Combine 16S and ITS results
final_results <-rbind(final_results_16S, final_results_ITS)

# Save results to an Excel file
dataset_names <- list('ALDEx2_16S' = final_results_16S, 'ALDEx2_ITS' = final_results_ITS, 'ALDEx2_all' = final_results)
openxlsx::write.xlsx(dataset_names, file = "Supplementary Table S9_ALDEx2.xlsx")

