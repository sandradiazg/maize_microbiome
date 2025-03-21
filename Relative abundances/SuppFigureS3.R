## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Plot fungal relative abundances of the top ten most abundant genera in the leaf 
## - Supplementary Figure 3

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(phyloseqCompanion)

#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_ITS.rds")

#Subsetting leaf samples
ps_LEAF <- subset_samples(ps_filtered_norm, compartment == "Leaf")

# Extracting the taxonomy table from the phyloseq object and converting it to a data frame
df_tax<-as.data.frame(phyloseq::tax_table(ps_LEAF))
# Removing the prefix "g__", "f__", "o__", "c__", "p__", "k__" from taxonomic ranks for clarity
df_tax$Genus <- gsub("g__", "",df_tax$Genus)
df_tax$Family <- gsub("f__", "",df_tax$Family)
df_tax$Order <- gsub("o__", "",df_tax$Order)
df_tax$Class <- gsub("c__", "",df_tax$Class)
df_tax$Phylum <- gsub("p__", "",df_tax$Phylum)
df_tax$Kingdom <- gsub("k__", "",df_tax$Kingdom)

# Replacing missing family or genus names with "Unknown" identifiers based on related taxonomic levels
df_tax$Family<-ifelse(is.na(df_tax$Family),str_c("Unkn_Fam_Ord_", df_tax$Order),df_tax$Family)
df_tax$Genus<-ifelse(is.na(df_tax$Genus),str_c("Unk_Fam_", df_tax$Family),df_tax$Genus)
#view(df_tax)

## Metadata
# Converting the taxonomic data frame into a matrix and creating a new phyloseq object with this updated information
mat_taxa <- as.matrix (df_tax)
TAX = phyloseq::tax_table(mat_taxa)

# Creating a new phyloseq object including the OTU table, updated taxonomic information, and sample metadata
ps_taxa <- phyloseq(TAX, otu_table(ps_LEAF), sample_data(ps_LEAF))

# Aggregating the phyloseq object at the Genus level
dat.aglo = tax_glom(ps_taxa, taxrank = "Genus", NArm=FALSE )

# Merging samples based on the "irrigation_tpreal" variable, transforming counts to relative abundance
ps_mg_plot = merge_samples(dat.aglo, "irrigation_tpreal")
ps_rel_abund_mg_plot = phyloseq::transform_sample_counts(ps_mg_plot, function(x){x / sum(x)})

# Extracting the sample data and preparing it for manipulation
df_metadata_mg_plot <- sample.data.frame(ps_rel_abund_mg_plot)
#str(df_metadata_mg_plot)

# Adding a column for the row names and another for the irrigation treatments
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  rownames_to_column (var="rnam")
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  mutate(irrigation_tpreal = rnam)
df_metadata_mg_plot<- df_metadata_mg_plot%>%
  column_to_rownames(var= "rnam")

# Keeping only the relevant columns for irrigation metadata and reordering them
df_metadata_mg_plot <- df_metadata_mg_plot[,names(df_metadata_mg_plot) %in% c("irrigation","irrigation_tpreal")]
col_order <- c("irrigation_tpreal","irrigation")
df_metadata_mg_plot <- df_metadata_mg_plot[, col_order]

# Creating categorical variables for irrigation type and time treatment (T01, T02, T03)
df_metadata_mg_plot<-df_metadata_mg_plot %>%
  mutate(irrigation = case_when(
    startsWith(irrigation_tpreal, "OW") ~ "OW",
    startsWith(irrigation_tpreal, "LW") ~ "LW"
  )) %>%
  mutate(tpreal = case_when(
    endsWith(irrigation_tpreal, "T01") ~ "T01",
    endsWith(irrigation_tpreal, "T02") ~ "T02",
    endsWith(irrigation_tpreal, "T03") ~ "T03"
  ))

# Adding metadata to the phyloseq object
metadata_plot_taxa <- sample_data(df_metadata_mg_plot)

# Extracting the taxonomic table from the phyloseq object with relative abundances
df_tax_plot_taxa<-as.data.frame(phyloseq::tax_table(ps_rel_abund_mg_plot))

# Defining a list of od the top 10 most abundant genera
fun_gen4 <- c("Abrothallus",
              "Acremonium",
              "Alternaria",
              "Aspergillus",
              "Cladorrhinum",
              "Cladosporium",
              "Exserohilum",
              "Filobasidium",
              "Fusarium",
              "Stemphylium")

# Defining a color palette for the fungal genera
nb.cols <-length(fun_gen4)
set_fun <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

# Assigning colors to each genus
c_fun4 <- data.frame(group = fun_gen4, color = set_fun)
c_fun4$group <- factor(c_fun4$group, levels = fun_gen4, ordered = TRUE)

# Grouping the remaining genera as "Others"
idx4 <- ! (df_tax_plot_taxa$Genus %in% fun_gen4)
df_tax_plot_taxa$Genus[idx4] <- "Others"

# Converting the updated taxonomic data back into a matrix
tax_plot_mat <- as.matrix(df_tax_plot_taxa)

# Creating a new phyloseq object for plotting
ps_plot_taxa <- phyloseq(phyloseq::tax_table(tax_plot_mat), otu_table(ps_rel_abund_mg_plot), metadata_plot_taxa)

# Melting the data frame for abundance plotting
abundance_table <- ps_plot_taxa %>%
  psmelt()
#write.table(abundance_table ,"results_abu_table_ITS_leaf.csv", sep = ";", quote = F)

# Filtering out "Others" genera from the abundance table for plotting
abundance_table <- abundance_table[!abundance_table$Genus %in% "Others",]
#view(abundance_table)

# Defining the list of fungal genera for which to create customized labels
fun_gen_it4 <- c(expression(italic("Abrothallus")), 
                 expression(italic("Acremonium")),
                 expression(italic("Alternaria")),
                 expression(italic("Aspergillus")),
                 expression(italic("Cladorrhinum")),
                 expression(italic("Cladosporium")),
                 expression(italic("Exserohilum")),
                 expression(italic("Filobasidium")),
                 expression(italic("Fusarium")),
                 expression(italic("Stemphylium")))

#Theme for plot
theme_tax =  theme(axis.text.y = element_text(size = 8, color = "black"),
                   axis.title.y = element_text(size = 9, color = "black"),
                   axis.text.x = element_text(size = 8, color = "black"),
                   plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
                   strip.text = element_text(size = 10, face = "bold", angle = 0),
                   strip.background = element_rect(fill="white"),
                   strip.text.x.top = element_text(size = 10, color = "black", face = "bold"),
                   panel.background=element_blank(),
                   panel.border = element_blank(),
                   panel.grid=element_blank(),
                   axis.line=element_line(color="black", size=0.4),
                   axis.ticks=element_line(color="black", size=0.4),
                   legend.background=element_blank(),
                   legend.text.align = 0,
                   legend.title = element_text(face="bold", size = 8),
                   legend.text = element_text(size = 6),
                   legend.position = "right",
                   legend.spacing.x = unit(0.2, 'mm'),
                   legend.spacing.y = unit(0.2, 'mm'),
                   legend.box.margin=margin(5,5,5,5),
                   legend.key.size = unit(2, 'mm'),
                   text=element_text(size=10, color="black"))

irrigation_order <- c("OW", "LW")
irrigation_names <- c("OW" = "OW",
                      "LW" = "LW")

# Creating the bar plot for relative abundance by genus
Sfig3 <- ggplot(data = abundance_table, aes(x = factor(tpreal), y = Abundance, fill = Genus))+
  geom_bar(aes(fill = factor(Genus, levels=fun_gen4), color = factor(Genus, levels=fun_gen4)), stat = "identity")   +     #ylim(0, 0.4)+
  #labs(title="Most abudant Genera", subtitle = "ITS")+
  facet_grid(~factor(irrigation, levels= irrigation_order, labels = irrigation_names))+
  ylab("Relative abundance") +
  xlab(NULL)+
  scale_fill_manual(values = as.character(c_fun4$color), breaks= fun_gen4, labels = fun_gen_it4) +
  scale_color_manual(values = (c_fun4$color), breaks =fun_gen4,  labels = fun_gen_it4) +
  theme_tax+
  guides(fill=guide_legend(ncol = 1,byrow=FALSE), color = FALSE)

# Save Figure
ggsave("SFigure_S3.pdf", Sfig3, units = "mm", width = 100, height = 50, dpi = 400)
