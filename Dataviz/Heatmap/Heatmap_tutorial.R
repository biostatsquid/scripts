# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Heatmap tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: `Laura Twomey
# Date: April 2023
# Version: 1.0

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(pheatmap) # for our heatmap
library(RColorBrewer) # for a colourful plot

# Set input path
path <- "~/Biostatsquid/Scripts/Heatmap/"
setwd(path)

# Create sample data ===================================================
set.seed(43)
data <- matrix(rnorm(500), 50, 10)
colnames(data) <- paste0("Sample_", 1:10)
rownames(data) <- paste0("Gene_", 1:50)

head(data)

## Get heatmap ===================================================
pheatmap(data)

# Clustering ===================================================
pheatmap(data, 
         cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D')


# Customisation ===================================================

## Adding a title
pheatmap(data, 
         main = "Super cool heatmap")

## Showing rows and columns
pheatmap(data,
         main = "Super cool heatmap",
         show_colnames = T, show_rownames = T,
         number_color = "black", 
         fontsize_number = 4)


## Showing values
pheatmap(data,
         fontsize_col = 10,
         fontsize_row = 10,
         display_numbers = TRUE,
         number_color = "black", 
         fontsize_number = 6,#
         border_color = "black") # default is grey60

## Cell colours
pheatmap(data,
         border_color = "black", # default is grey60
         number_color = "black", 
         fontsize_number = 8,
         col = brewer.pal(10, 'RdYlGn')) # https://r-graph-gallery.com/38-rcolorbrewers-palettes.html


## Legend customisation
pheatmap(data, 
         legend_breaks = c(-2, 0, 2),
         legend_labels = c("Low", "Medium", "High"))


pheatmap(data, 
         legend = FALSE)


# Split heatmap clusters
pheatmap(data, 
         cutree_rows = 2, cutree_cols = 4)


pheatmap(data, 
         border_color = FALSE,      # no border to cell
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size 
         angle_col = 45,             # angle for column labels
         na_col = "black",           # color of the cell with NA values
         legend = FALSE#,            # to draw legend or not (TRUE/FALSE)
)


# Annotations

# Create a data frame for column annotation
ann_df <- data.frame(Group = rep(c("Disease", "Control"), c(5, 5)),
                          Lymphocyte_count = rnorm(10))
row.names(ann_df) <- colnames(data)
head(ann_df)

gene_functions_df <- data.frame(gene_functions = rep(c('Oxidative_phosphorylation', 
                                                       'Cell_cycle',
                                                       'Immune_regulation',
                                                       'Signal_transduction',
                                                       'Transcription'), rep(10, 5)))
row.names(gene_functions_df) <- rownames(data)

ann_colors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
              "Cell_cycle" = "#708238",
              "Immune_regulation" = "#9E0142",
              "Signal_transduction" = "beige", 
              "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
              "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)


pheatmap(data, 
         col = brewer.pal(10, 'RdYlGn'),
         annotation_row = gene_functions_df, 
         annotation_col = ann_df, 
         annotation_colors = ann_colors,
         main = "Super heatmap with annotations") 

# Save it -----------------------------------------------------------
heat_plot <- pheatmap(data, 
                      col = brewer.pal(10, 'RdYlGn'), # choose a colour scale for your data
                      cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      annotation_row = gene_functions_df, # row (gene) annotations
                      annotation_col = ann_df, # column (sample) annotations
                      annotation_colors = ann_colors, # colours for your annotations
                      annotation_names_row = F, 
                      annotation_names_col = F,
                      fontsize_row = 10,          # row label font size
                      fontsize_col = 7,          # column label font size 
                      angle_col = 45, # sample names at an angle
                      legend_breaks = c(-2, 0, 2), # legend customisation
                      legend_labels = c("Low", "Medium", "High"), # legend customisation
                      show_colnames = T, show_rownames = F, # displaying column and row names
                      main = "Super heatmap with annotations") # a title for our heatmap


pdf(paste0(path, heat_plot), height = 10, width = 8)
heat_plot
dev.off()

