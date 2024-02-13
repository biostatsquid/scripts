# ---------------------- #
# SingleR_tutorial.R
# ---------------------- #
# Author: Laura Twomey
# Version: 1.0

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(42)

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(celldex)
library(SingleR)
library(Seurat)

# Load the data ============================================================
# We'll use a PBMC dataset from the R package scRNAseq
sce <- scRNAseq::KotliarovPBMCData(mode = c('rna'))
seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
seu <- NormalizeData(object = seu)

# Prepare SingleR inputs ===================================================
# 1. Get normalised or raw counts
raw_counts <- LayerData(seu, assay = "RNA", layer = 'counts') #
raw_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]
norm_counts <- LayerData(seu, assay = "RNA", layer = 'data') #
norm_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]

# 2. Get reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()
unique(ref$label.main)
unique(ref$label.fine)
# Subset to include only relevant cell types (CAREFUL!)
ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|
                 Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]
unique(ref$label.main)

# 3. Run SingleR
ct_ann <- SingleR(test = norm_counts, # we could also use sce or raw_counts
                  ref = ref, 
                  labels = ref$label.main,
                  de.method = 'wilcox')

ct_ann %>% head()
unique(ct_ann$pruned.labels)
table(ct_ann$pruned.labels)

# Inspect quality of the predictions
plotScoreHeatmap(ct_ann)
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)
summary(is.na(ct_ann$pruned.labels))

# Add to seurat object
seu <- AddMetaData(seu, ct_ann$pruned.labels, col.name = 'SingleR_HCA')

# Visualise them on the UMAP
seu <- SetIdent(seu, value = "SingleR_HCA")
DimPlot(seu, label = T , repel = T, label.size = 3) + NoLegend()

  
  
  
  
  
  
  
  
  