## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
### - Creation of ITS ps_object from ASV, taxonomy and sample metadata matrixes.
### - Count and taxonomy filtering. 
### - Normalization by sequencing depth.
### - RDS object.

#Clean working environment
rm(list= ls())
#Load required packages
library(tibble) #formatting tables. Row_to_column function
library(dplyr) # filter and reformat data frames  
library(phyloseq)
library(stringr) #text strings


#Load ASV matrix
ASV_mat <- read.csv("ITS_ALL_medium_asv_otu.csv", row.names=1)
#Fortmatting ASV matrix
ASV_mat <- ASV_mat %>% 
  rownames_to_column(var ="ASV")

names(ASV_mat) <- gsub(x = names(ASV_mat), pattern = "_ITS_S[0-9]+_L00[0-9]+_R1_trimmed_val_1.fq.gz", replacement = "_ITS") 
names(ASV_mat) <- gsub(x = names(ASV_mat), pattern = "^X[0-9]+", replacement = "")  
names(ASV_mat) <- gsub(x = names(ASV_mat), pattern = "^_", replacement = "")  

#load Taxonomy matrix
tax_mat <- read.csv("ITS_ALL_medium_asv_tax.csv", row.names=1)
#Fortmatting taxonomy matrix
tax_mat <- tax_mat %>% 
  rownames_to_column(var = "ASV")

#Load sample metadata table
samples_df  <- read.csv("sample_metadata_ITS.csv",sep = ";")
#Fortmatting taxonomy matrix
samples_df$Sample_Id  <- samples_df  %>% gsub(x = samples_df$Sample_Id, pattern = "^[0-9]+", replacement = "") %>% gsub(pattern = "^_", replacement = "") 
samples_df$Sample_Id  <- samples_df  %>% gsub(x = samples_df$Sample_Id, pattern = "-", replacement = ".") 
samples_df$name  <- samples_df  %>% gsub(x = samples_df$name, pattern = "^[0-9]+", replacement = "") %>% gsub(pattern = "^_", replacement = "") 
samples_df$name  <- samples_df  %>% gsub(x = samples_df$name, pattern = "-", replacement = ".") 
samples_df$compartment_tpreal <- str_c(samples_df$compartment, "_", samples_df$tpreal)
samples_df$compartment_irrigation <- str_c(samples_df$compartment, "_", samples_df$irrigation)
samples_df$compartment_irrigation_tpreal <- str_c(samples_df$compartment, "_", samples_df$irrigation, "_", samples_df$tpreal)
samples_df$irrigation_tpreal <- str_c(samples_df$irrigation, "_", samples_df$tpreal)
samples_df$plot_name<- str_c( samples_df$plot, "_", samples_df$name)

irrigation_order <- c("OW", "LW")
samples_df$irrigation <- factor(samples_df$irrigation,     # Reorder factor levels
                                irrigation_order)
compartment_order <- c("Bulk soil", "Rhizosphere", "Root", "Leaf", "Grain")
samples_df$compartment <- factor(samples_df$compartment,     # Reorder factor levels
                                 compartment_order)
# Create phyloseq object
#define the row names from the matrixes.
ASV_mat <- ASV_mat %>%
  tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample_Id")

#Transform into matrixes ASV and tax tables (sample table can be left as data frame)
ASV_mat <- as.matrix(ASV_mat)
tax_mat <- as.matrix(tax_mat)
class(ASV_mat) <- "numeric"

#Transform to phyloseq objects
ASV = otu_table(ASV_mat, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_mat)
samples = sample_data(samples_df)

ps <- phyloseq(ASV, TAX, samples)
#Check phyloseq object parameters
ps
ntaxa(ps)
nsamples(ps)

# Taxonomic Filtering
#In most cases, the organisms within a sample are well represented in the reference database. 
#When this is the case, it's advisable to filter out reads/ASVs that cannot be assigned a high-rank taxonomy label. 
#These are most likely contaminates/artifacts that don't exist in nature and should be removed.

psp = subset_taxa(ps, !is.na(Phylum)) #Remove taxa only assigned up to Kingdom level
psp
psc = subset_taxa(psp, !is.na(Class)) #Remove taxa only assigned up to Phylum level
psc
ps0 = subset_taxa(psc, !is.na(Order))#Remove taxa only assigned up to Class level
ps0


# Keep taxa with 3 or more counts in at least >2% of the samples.
ps_filter = filter_taxa(ps0, function(x) sum(x >=3) > (0.02*length(x)), TRUE) 
ps_filter

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(ps_filter))
standf = function(x, t=total) round(t * (x / sum(x)))
ps_filtered_norm = transform_sample_counts(ps_filter, standf)
ps_filtered_norm

#Create RDS object
saveRDS(ps_filter, "ps_filter_ITS.rds")
saveRDS(ps_filtered_norm, "ps_filtered_norm_ITS.rds")
