# ======================== # 
# 1_DGE.R
# ======================== # 

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(Seurat)
library(tidyverse)

# Read in data ===================================================
# I'll be using the seurat object we processed using Seurat_object.R
seu <- qs::qread("C:/Users/laura/Documents/Biostatsquid/Projects/data/scRNAseq/pbmc3k.qs")
DimPlot(seu, group.by = 'seurat_annotations')


# DGE analysis ===================================================
# Find markers for every cluster compared to all remaining cells
# This can take time for large datasets
all_markers <- FindAllMarkers(seu, 
                              only.pos = FALSE,    # set to TRUE if you want only positive (upregulated) markers
                              min.pct = 0.25,     # Detected in at least 25% of cells in either population
                              logfc.threshold = 0.25)
# View the top markers
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(top_markers)

# Save the complete marker list
write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)

# Compare two specific clusters
# Replace '0' and '1' with the clusters you want to compare
cluster0_vs_1_markers <- FindMarkers(seurat_obj, 
                                     ident.1 = 0,  # Cluster of interest
                                     ident.2 = 1,  # Cluster to compare against
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)

head(cluster0_vs_1_markers, n = 20)

##########################
# 10. VISUALIZATION OF MARKER EXPRESSION
##########################

# Visualize marker expression for the top markers of each cluster
top_genes_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) %>%
  pull(gene)

DoHeatmap(seurat_obj, features = top_genes_per_cluster)

# Violin plots of marker genes
VlnPlot(seurat_obj, features = top_genes_per_cluster[1:6], pt.size = 0)

# Feature plots showing expression on UMAP
FeaturePlot(seurat_obj, 
            features = top_genes_per_cluster[1:6], 
            min.cutoff = "q9",  # Show only cells with expression > 9th quantile
            ncol = 3)



