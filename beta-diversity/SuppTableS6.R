## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Run PERMANOVA analysis (Supplementary Table S6)
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.
## This script only show the analysis for the aggregated dataset. The same method was used for the rest of analysed subsets.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(vegan)
library(tibble)
library(openxlsx)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

# Set a seed for reproducibility
set.seed(1025)

# Convert sample metadata from phyloseq object to a standard data frame
metadata_abu <- as(sample_data(ps_filtered_norm), "data.frame") ## convert sample_data to data.frame

# Compute Bray-Curtis distance matrix (weighted) for beta diversity analysis
bray_dist = phyloseq::distance(ps_filtered_norm, method="bray", weighted=T)

# --------- Beta Dispersion and PERMANOVA Analysis ---------
## Compartment analysis

b <- betadisper(bray_dist, metadata_abu$tpreal)# Test homogeneity of dispersion
permutest(b, permutations = 1000)# Perform a permutation test for dispersion

# Run PERMANOVA to test the effect of 'compartment' on community composition
perm_compartment<-adonis2(bray_dist ~ compartment, 
              data = metadata_abu,permutations = 1000)


## Tpreal (time point) analysis

b <- betadisper(bray_dist, metadata_abu$tpreal) # Test homogeneity of dispersion
permutest(b, permutations = 1000) # Perform a permutation test for dispersion

# Run PERMANOVA to test the effect of 'tpreal' on community composition
perm_tpreal <- adonis2(bray_dist ~ tpreal, 
              data = metadata_abu,permutations = 1000)

## irrigation analysis
b <- betadisper(bray_dist, metadata_abu$irrigation)# Test homogeneity of dispersion
permutest(b, permutations = 1000) # Perform a permutation test for dispersion

# Run PERMANOVA to test the effect of 'irrigation' on community composition
perm_irrigation <- adonis2(bray_dist ~ irrigation, 
              data = metadata_abu,permutations = 1000)


# --------- Extract and Compile PERMANOVA Results ---------

# Get the names of all PERMANOVA result objects
variables <- ls(pattern = "^perm_")

# Initialize an empty list to store results
results <- list()

# Iterate over each PERMANOVA result and extract key statistics
for (var in variables) {
  # Retrieve the PERMANOVA object from the global environment
  obj <- get(var)
  
  # Extract key statistical values: p-value, R² (effect size), and F-statistic
  pr_f <- obj$`Pr(>F)`[1]
  r2 <- obj$R2[1]
  f_val <- obj$F[1]
  
  # Store the extracted values in a temporary data frame
  results[[var]] <- data.frame(P = pr_f, R2 = r2, F = f_val)
}

# Combine all results into a single table
table_permanova <- do.call(rbind, results)

# Assign dataset names as row labels
rownames(table_permanova) <- variables
table_permanova = table_permanova %>% rownames_to_column(var="Dataset")

# Display the final results table
print(table_permanova)

# Save the results as an Excel file
openxlsx::write.xlsx(table_permanova, file = "Supplementary Table S6_PERMANOVA.xlsx")
