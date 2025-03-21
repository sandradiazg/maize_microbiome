## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Create Venn diagrams and lists.
## - Plot Figures 4A/D, 4B/E and 5B
## This script is for 16S graphs, but the same code was used for ITS plots using the ITS ps_filtered_norm ps object.

#Clean working environment
rm(list= ls())
#Load required packages
library(phyloseq)
library(VennDiagram)
library(MicEco)
library(MicrobiotaProcess)
library(plyr)
#Load ps object
ps_filtered_norm <- readRDS("ps_filtered_norm_16S.rds")

#Figure 4A/D
ps_venn(ps_filtered_norm, group = "compartment", main= "General (16S)", plot = TRUE) #This function allows to double-check the numbers in the Venn diagram. You can also use directly this plot, but I preferred the graphics of the 'venn.diagram' function. 
vennlist_com_ls <- get_vennlist(ps_filtered_norm, factorNames="compartment")
vennlist_com=ldply(vennlist_com_ls, data.frame)

compartment_order <- c("Bulk soil", "Rhizosphere", "Root", "Leaf", "Grain")
vennlist_com_ls = vennlist_com_ls[compartment_order] #Determine a specific order for the treatments in the graphs
compartment_names <- c(
  'Bulk soil'="Bulk soil",
  'Rhizosphere'="Rhizosphere",
  'Root'="Root",
  'Leaf'="Leaf",
  'Grain' = "Grain") #Change the names of your treatments for the graphs  (if needed)
cols_comp <- c("Bulk soil" = "#CC9999", 
               "Rhizosphere" = "#990033", 
               "Root" = "#33CCFF", 
               "Leaf" = "#66CC00", 
               "Grain" = "#FFCC33")

venn.diagram(vennlist_com_ls,
             output=TRUE,
             filename="Fig4A.tiff", 
             category.names = compartment_names,
             disable.logging = TRUE,
             
             # Output features
             imagetype="tiff" ,
             height = 48,
             width = 48,
             resolution = 600,
             compression = "lzw",
             units = "mm",
             
             # Circles
             lwd = 0.8,
             lty = 1,
             fill = cols_comp,
             col = cols_comp,
             alpha = 0.4,
             scaled = FALSE,
             
             #Numbers
             cex = 0.35,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.col= cols_comp,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(0,-18,-145,145,-10),
             cat.dist = c(0.2, 0.23,0.2, 0.22,0.22),
             cat.fontfamily = "sans",
             margin = 0.1,
)

#Figure 4B/E
ps_venn(ps_filtered_norm, group ="tpreal", main= "General (16S)", plot = TRUE)
vennlist_tpreal_ls <- get_vennlist(ps_filtered_norm, factorNames="tpreal")
vennlist_tpreal=ldply(vennlist_tpreal_ls, data.frame)

tp_cols <- c("T01" = "#66cccc",
             "T02" = "#3399CC",
             "T03" = "#003333")

venn.diagram(vennlist_tpreal_ls,
             output=TRUE,
             filename="Fig2B.tiff", 
             disable.logging = TRUE,
             
             # Output features
             imagetype="tiff" ,
             height = 48,
             width = 48,
             rsolution = 600,
             compression = "lzw",
             units = "mm",
             
             # Circles
             lwd = 0.8,
             lty = 1,
             fill = tp_cols,
             col = tp_cols,
             alpha = 0.4,
             scaled = FALSE,
             
             #Numbers
             cex = 0.8,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.9,
             cat.col= tp_cols,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             # cat.pos = c(-50, 50),
             cat.dist = c(0.07, 0.07, 0.055),
             cat.fontfamily = "sans",
             margin = 0.1,
)

#Figure 5B
ps_venn(ps_filtered_norm, group = "irrigation", main= "General (16S)", plot = TRUE)
vennlist_irr_ls <- get_vennlist(ps_filtered_norm, factorNames="irrigation")
vennlist_irr=ldply(vennlist_irr_ls, data.frame)

vennlist_irr_ls <- vennlist_irr_ls[c("OW", "LW")]
irrigation_names <- c("OW" = "OW",
                      "LW" = "LW")
irrigation_cols <- c("OW" = "#4a80ff",
                     "LW" = "tan2")

venn.diagram(vennlist_irr_ls,
             output=TRUE,
             filename="Fig5B.tiff", 
             category.names = irrigation_names,
             disable.logging = TRUE,
             
             # Output features
             imagetype="tiff" ,
             height = 48,
             width = 48,
             resolution = 600,
             compression = "lzw",
             units = "mm",
             
             # Circles
             lwd = 0.8,
             lty = 1,
             fill = irrigation_cols,
             col = irrigation_cols,
             alpha = 0.4,
             scaled = FALSE,
             
             #Numbers
             cex = 0.6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.8,
             cat.col= irrigation_cols,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(0, 0),
             cat.dist = c(0.04, 0.04),
             cat.fontfamily = "sans",
             margin = 0.08,
             rotation.degree = 180
)             
