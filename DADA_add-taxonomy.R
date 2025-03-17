#!/usr/bin/env Rscript

library(dada2)
library(tidyverse)
library(DECIPHER)

cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

seqtab.nochim <- args[1]
output <- args[2]
name <- args[3]
tax_db <- strsplit(args[4], ",")[[1]]
method <- args[5]
threshold <- as.integer(args[6])

dir.create(file.path(output, "03_taxonomy"), showWarnings = FALSE)
dir.create(file.path(output, "03_taxonomy", name), showWarnings = FALSE)

output <- paste0(output,"/03_taxonomy/",name,"/")

# Assign taxonomy (general)
set.seed(42) #random  generator necessary for reproducibility

seqtab.nochim <- readRDS(seqtab.nochim)

if (grepl('[Dd]ecipher|DECIPHER', method)){ # use decipher

cat(paste0("You are using DECIPHER version ", packageVersion('DECIPHER'),'\n\n'))

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs

if (grepl('\\.RData$', tax_db[1])){
  load(tax_db[1])
} else if (grepl('\\.rds$',tax_db[1])){
  trainingSet <- readRDS(tax_db[1])
}

if (is.na(threshold)){
    threshold <- 60
}

ids <- IdTaxa(dna,
             trainingSet,
             strand = "top",
             processors = NULL,
             verbose = FALSE,
             threshold = threshold)

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

} else { # use regular dada2 classificator
  
  if (is.na(threshold)){
   threshold <- 50
  }

  cat('# The taxonomic classification included in dada2 will be used\n')
  taxid <- assignTaxonomy(seqtab.nochim,
                        tax_db[1], 
                        multithread=TRUE,
                        minBoot=threshold,
                        tryRC=TRUE
)
  cat('# Taxonomy assigned to genus level\n')
  
  if (!is.na(tax_db[2])) { # add species level if db available
    taxid <- addSpecies(taxid, 
                        tax_db[2], 
                        verbose=TRUE, 
                        allowMultiple=3)
    cat('\n# Taxonomy assigned to species level\n')
  }
}

# Write to disk
saveRDS(taxid, paste0(output, name, "_tax_assignation.rds"))

cat('\n')
cat(paste0('# The obtained taxonomy file can be found in "', paste0(output, name, "_tax_assignation.rds"), '"\n'))
cat('\n# All done!\n\n')

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

#fasta header as ASV_
for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}

#ASV sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))

#ASV abundance
asv_otu <- t(seqtab.nochim)
row.names(asv_otu) <- sub(">", "", asv_headers)

#ASV taxonomy
asv_tax <- taxid
row.names(asv_tax) <- sub(">", "", asv_headers)
 dim(asv_tax)

#Merging abundance and tax table
OTU_TAX_table <- merge(asv_otu, asv_tax, by=0)

#Writing out output files

write(asv_fasta,paste0(output, name,"_asv_fasta.fa"))
write.table(asv_otu, paste0(output, name,"_asv_otu.csv"), sep=",", quote=F, col.names=NA)
write.table(asv_tax, paste0(output, name, "_asv_tax.csv"), sep=",", quote=F, col.names=NA)
write.table(OTU_TAX_table, paste0(output,name,"_OTU_TAX_table.csv"), sep=",", quote=F, col.names=NA)

