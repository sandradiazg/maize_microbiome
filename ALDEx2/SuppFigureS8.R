## Maize associated bacterial and fungal microbiomes show contrasting conformation patterns dependent on plant compartment and water availability
## Sandra Díaz-González, 2025 (CBGP, UPM)
## Description: 
## - Plot ALDEx2 results (Supplementary Figure S8).

#Clean working environment
rm(list= ls())
#Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

#Load ALDEx2 results (from script 'SuppTableS9')
mydata<-readxl::read_excel("Supplementary Table S9_ALDEx2.xlsx", sheet="ALDEx2_all")
#View(mydata)

#Plot's theme
main_theme <- theme(panel.background=element_blank(),
                        panel.border = element_rect(color="black", size=0.3, fill=NA),
                        panel.grid=element_blank(),
                        plot.title = element_text(face = "bold", size = 12),
                        axis.line=element_line(color="black", size=0.3),
                        axis.ticks.x = element_line(color="black", size=0.2),
                        axis.ticks.y = element_line(color="black", size=0.2),
                        axis.text.x = element_text(size = 7, color = "black"),
                        #axis.text.y.left = element_text(size = 8, color = "black", face = "italic"),
                        axis.text.y = element_blank(),
                        axis.title.y = element_text(size = 8, color = "black"),
                        axis.title.x = element_text(size = 8, color = "black"),
                        # face = "bold"),
                        legend.background=element_blank(),
                        #legend.key = element_rect(fill =),
                        legend.title = element_text(colour="black", size=10, face="bold"),
                        legend.text = element_text(colour="black", size = 10),
                        legend.spacing.y = unit(0.5, 'mm'),
                        legend.spacing.x = unit(0.5, 'mm'),
                        legend.key.height = unit(4, "mm"),
                        legend.position = "right",
                        strip.text.x = element_text(colour="black", face="bold",
                                                    size=8, margin = margin(t = 2, b = 2)),# Tamaño del texto de los títulos de facet
                        strip.background = element_rect(colour = "#EEE5DE", fill = "#EEE5DE"), #Color de fondo del título
                        strip.text.y = element_blank(),  # Ocultar nombres en la derecha # Ocultar nombres arriba
                        axis.ticks.y.right = element_blank(),  # Ocultar las marcas en la derecha
                        #strip.placement = "outside",           # Asegurar que las etiquetas del facet estén fuera
                        panel.spacing = unit(0.1, "lines") # Espaciado entre paneles)
                        #aspect.ratio = 1
)
# Define a theme for the genus labels
theme_tile =   theme(#legend.position = "right",  # Leyenda a la derecha
  axis.text.y = element_text(size = 8, face = "italic", hjust = 1, margin = margin(r = 5))  # Cursiva en nombres de género
)
# Define color schemes for Kingdom and Time points
king_cols <- c("Bacteria" = "#FCC5C0",
               "Fungi" = "#C6DBEF")

tp_cols <- c("T01" = "#66cccc",
             "T02" = "#3399CC",
             "T03" = "#003333")

#Subset by compartment
bs<-subset(mydata, compartment == "Bulk_soil")
rhizos<-subset(mydata, compartment == "Rhizosphere")
root<-subset(mydata, compartment == "Root")
leaf<-subset(mydata, compartment == "Leaf")
grain<-subset(mydata, compartment == "Grain")

#############
#Bulk soil#
#############
# Expand grid to ensure all genus-time point combinations exist
combinaciones <- expand.grid(Genus = unique(bs$Genus), 
                             tpreal = unique(bs$tpreal))

# Merge expanded grid with data, filling missing Kingdom values
bs_completo <- combinaciones %>%
  left_join(bs, by = c("Genus", "tpreal")) %>%
  group_by(Genus) %>%
  fill(Kingdom, .direction = "updown") %>%
  ungroup()

# Manually order Kingdom so Bacteria appears first
bs_completo <- bs_completo %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Fungi", "Bacteria"))) %>%
  arrange(Kingdom, Genus)

# Ensure Genus follows this new order
bs_completo$Genus <- factor(bs_completo$Genus, levels = unique(bs_completo$Genus))

# Print to check completeness
print(bs_completo)

# Create the main plot
a <- ggplot(bs_completo, aes(x = effect, y = Genus)) +
  #geom_tile(aes(x = -0.5, fill = Kingdom), width = 0.5, height = 0.8) +  # Columna lateral Kingdom
  geom_segment(aes(x = 0, xend = effect, y = Genus, yend = Genus, color = tpreal) , linewidth= 0.5,na.rm = TRUE) +
  geom_point(aes(color = tpreal), size = 1.5) +  # Puntos por tpreal
  geom_vline(xintercept = 0, size = 0.3,) +
  #xlim(c(0, 6))+
  #scale_x_continuous(expand = c(0.2, 0), limits= c(0, 5.5)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Evita expansión del eje Y
  scale_color_manual(values = tp_cols)+
  xlab("Effect size (LW vs. OW)") +
  ylab(NULL) +
  facet_wrap(~ tpreal, scales = "free_x", ncol = 3) +
  guides(color = FALSE)+# Facet con ejes compartidos en Y
  main_theme
#
# Create a separate plot for Kingdom labels
b <- ggplot(unique(bs_completo[, c("Genus", "Kingdom")]), aes(x = 1, y = Genus, fill = Kingdom)) +
  geom_tile(width = 0.5, height = 0.8) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = king_cols)+
  guides(fill = FALSE)+
  theme_void() +  # Empty theme (no axis)
  theme_tile

# Combine the Kingdom labels plot and the main plot
final_plot_bs <- b + a + plot_layout(widths = c(0.2, 5))
final_plot_bs

##Save Figure
ggsave("SFigure_S8A.pdf", final_plot_bs, units = "mm", width = 90, height = 200, dpi = 300)

#############
#Rhizosphere#
#############
# Expand grid to ensure all genus-time point combinations exist
combinaciones <- expand.grid(Genus = unique(rhizos$Genus), 
                             tpreal = unique(rhizos$tpreal))
# Merge expanded grid with data, filling missing Kingdom values
rz_completo <- combinaciones %>%
  left_join(rhizos, by = c("Genus", "tpreal")) %>%
  group_by(Genus) %>%
  fill(Kingdom, .direction = "updown") %>%
  ungroup()

# Manually order Kingdom so Bacteria appears first
rz_completo <- rz_completo %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Fungi", "Bacteria"))) %>%
  arrange(Kingdom, Genus)

# Ensure Genus follows this new order
rz_completo$Genus <- factor(rz_completo$Genus, levels = unique(rz_completo$Genus))

# Print to check completeness
print(rz_completo)


# Create the main plot
a <- ggplot(rz_completo, aes(x = effect, y = Genus)) +
  #geom_tile(aes(x = -0.5, fill = Kingdom), width = 0.5, height = 0.8) +  # Columna lateral Kingdom
  geom_segment(aes(x = 0, xend = effect, y = Genus, yend = Genus, color = tpreal) , linewidth= 0.5,na.rm = TRUE) +
  geom_point(aes(color = tpreal), size = 1.5) +  # Puntos por tpreal
  geom_vline(xintercept = 0, size = 0.3,) +
  #xlim(c(-5,0))+
  #scale_x_continuous(expand = c(0.2, 0), limits= c(0, 2.5)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Evita expansión del eje Y
  scale_color_manual(values = tp_cols)+
  xlab("Effect size (LW vs. OW)") +
  ylab(NULL) +
  facet_wrap(~ tpreal, scales = "free_x", ncol = 3) +
  guides(color = FALSE)+# Facet con ejes compartidos en Y
  main_theme

# Create a separate plot for Kingdom labels
b <- ggplot(unique(rz_completo[, c("Genus", "Kingdom")]), aes(x = 1, y = Genus, fill = Kingdom)) +
  geom_tile(width = 0.5, height = 0.8) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = king_cols)+
  guides(fill = FALSE)+
  theme_void() +  # Tema vacío (sin ejes)
  theme_tile

# Combine the Kingdom labels plot and the main plot
final_plot_rz <- b + a + plot_layout(widths = c(0.2, 5))
final_plot_rz

##Save Figure
ggsave("SFigure_S8B.pdf", final_plot_rz, units = "mm", width = 85, height = 140, dpi = 300)

######
#Root#
######

# Expand grid to ensure all genus-time point combinations exist
combinaciones <- expand.grid(Genus = unique(root$Genus), 
                             tpreal = unique(root$tpreal))
# Merge expanded grid with data, filling missing Kingdom values
root_completo <- combinaciones %>%
  left_join(root, by = c("Genus", "tpreal")) %>%
  group_by(Genus) %>%
  fill(Kingdom, .direction = "updown") %>%
  ungroup()

# Manually order Kingdom so Bacteria appears first
root_completo <- root_completo %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Fungi", "Bacteria"))) %>%
  arrange(Kingdom, Genus)

# Ensure Genus follows this new order
root_completo$Genus <- factor(root_completo$Genus, levels = unique(root_completo$Genus))

# Print to check completeness
print(root_completo)

# Create the main plot
a <- ggplot(root_completo, aes(x = effect, y = Genus)) +
  #geom_tile(aes(x = -0.5, fill = Kingdom), width = 0.5, height = 0.8) +  # Columna lateral Kingdom
  geom_segment(aes(x = 0, xend = effect, y = Genus, yend = Genus, color = tpreal) , linewidth= 0.5,na.rm = TRUE) +
  geom_point(aes(color = tpreal), size = 1.5) +  # Puntos por tpreal
  geom_vline(xintercept = 0, size = 0.3,) +
  #xlim(c(-5,0))+
  #scale_x_continuous(expand = c(0.2, 0), limits= c(0, 2)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Evita expansión del eje Y
  scale_color_manual(values = tp_cols)+
  xlab("Effect size (LW vs. OW)") +
  ylab(NULL) +
  facet_wrap(~ tpreal, scales = "free_x", ncol = 3) +
  guides(color = FALSE)+# Facet con ejes compartidos en Y
  main_theme

# Create a separate plot for Kingdom labels
b <- ggplot(unique(root_completo[, c("Genus", "Kingdom")]), aes(x = 1, y = Genus, fill = Kingdom)) +
  geom_tile(width = 0.5, height = 0.8) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = king_cols)+
  guides(fill = FALSE)+
  theme_void() +  # Tema vacío (sin ejes)
  theme_tile

# Combine the Kingdom labels plot and the main plot
final_plot_root <- b + a + plot_layout(widths = c(0.2, 5))
final_plot_root
#Save Figure
ggsave("SFigure_S8C.pdf", final_plot_root, units = "mm", width = 80, height = 25, dpi = 300)

######
#Leaf#
######

# Expand grid to ensure all genus-time point combinations exist
combinaciones <- expand.grid(Genus = unique(leaf$Genus), 
                             tpreal = unique(leaf$tpreal))
# Merge expanded grid with data, filling missing Kingdom values
leaf_completo <- combinaciones %>%
  left_join(leaf, by = c("Genus", "tpreal")) %>%
  group_by(Genus) %>%
  fill(Kingdom, .direction = "updown") %>%
  ungroup()

# Ensure Genus follows this new order
leaf_completo <- leaf_completo %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Fungi", "Bacteria"))) %>%
  arrange(Kingdom, Genus)

# Print to check completeness
leaf_completo$Genus <- factor(leaf_completo$Genus, levels = unique(leaf_completo$Genus))

# Create the main plot
a <- ggplot(leaf_completo, aes(x = effect, y = Genus)) +
  #geom_tile(aes(x = -0.5, fill = Kingdom), width = 0.5, height = 0.8) +  # Columna lateral Kingdom
  geom_segment(aes(x = 0, xend = effect, y = Genus, yend = Genus, color = tpreal) , linewidth= 0.5,na.rm = TRUE) +
  geom_point(aes(color = tpreal), size = 1.5) +  # Puntos por tpreal
  geom_vline(xintercept = 0, size = 0.3,) +
  xlim(c(-5,0))+
  #scale_x_continuous(expand = c(0.2, 0), limits= c(0, 2)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Evita expansión del eje Y
  scale_color_manual(values = tp_cols)+
  xlab("Effect size (LW vs. OW)") +
  ylab(NULL) +
  facet_wrap(~ tpreal, scales = "free_x", ncol = 3) +
  guides(color = FALSE)+# Facet con ejes compartidos en Y
  main_theme

# Create a separate plot for Kingdom labels
b <- ggplot(unique(leaf_completo[, c("Genus", "Kingdom")]), aes(x = 1, y = Genus, fill = Kingdom)) +
  geom_tile(width = 0.5, height = 0.8) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = king_cols)+
  guides(fill = FALSE)+
  theme_void() +  # Tema vacío (sin ejes)
  theme_tile

# Combine the Kingdom labels plot and the main plot
final_plot_leaf <- b + a + plot_layout(widths = c(0.2, 5))
final_plot_leaf

#Save Figure
ggsave("SFigure_S8D.pdf", final_plot_leaf, units = "mm", width = 80, height = 25, dpi = 300)


#######
#Grain#
#######

# Expand grid to ensure all genus-time point combinations exist
combinaciones <- expand.grid(Genus = unique(grain$Genus),
                             tpreal = unique(grain$tpreal))
# Merge expanded grid with data, filling missing Kingdom values
grain_completo <- combinaciones %>%
  left_join(grain, by = c("Genus", "tpreal")) %>%
  group_by(Genus) %>%
  fill(Kingdom, .direction = "updown") %>%
  ungroup()

# Manually order Kingdom so Bacteria appears first
grain_completo <- grain_completo %>%
  mutate(Kingdom = factor(Kingdom, levels = c("Fungi", "Bacteria"))) %>%
  arrange(Kingdom, Genus)

# Ensure Genus follows this new order
grain_completo$Genus <- factor(grain_completo$Genus, levels = unique(grain_completo$Genus))

# Print to check completeness
print(grain_completo)

# Create the main plot
a <- ggplot(grain_completo, aes(x = effect, y = Genus)) +
  #geom_tile(aes(x = -0.5, fill = Kingdom), width = 0.5, height = 0.8) +  # Columna lateral Kingdom
  geom_segment(aes(x = 0, xend = effect, y = Genus, yend = Genus, color = tpreal) , linewidth= 0.5,na.rm = TRUE) +
  geom_point(aes(color = tpreal), size = 1.5) +  # Puntos por tpreal
  geom_vline(xintercept = 0, size = 0.3,) +
  #xlim(c(-5,0))+
  #scale_x_continuous(expand = c(0.2, 0), limits= c(0, 6)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Evita expansión del eje Y
  scale_color_manual(values = tp_cols)+
  xlab("Effect size (LW vs. OW)") +
  ylab(NULL) +
  #facet_wrap(~ tpreal, scales = "free_x", ncol = 2) +
  guides(color = FALSE)+# Facet con ejes compartidos en Y
  main_theme

# Create a separate plot for Kingdom labels
b <- ggplot(unique(grain_completo[, c("Genus", "Kingdom")]), aes(x = 1, y = Genus, fill = Kingdom)) +
  geom_tile(width = 0.5, height = 0.8) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = king_cols)+
  guides(fill = FALSE)+
  theme_void() +  # Tema vacío (sin ejes)
  theme_tile

# Combine the Kingdom labels plot and the main plot
final_plot_grain <- b + a + plot_layout(widths = c(0.2, 5))
final_plot_grain

# Save Figure
ggsave("SFigure_S8E.pdf", final_plot_grain, units = "mm", width = 50, height = 35, dpi = 300)