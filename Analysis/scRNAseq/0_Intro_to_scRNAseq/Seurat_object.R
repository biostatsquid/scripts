# ======================== # 
# Seurat_object.R
# ======================== # 

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(Seurat)
library(SeuratData) # we'll be using a dataset from seurat data https://github.com/satijalab/seurat-data 
library(tidyverse)
library(pheatmap)
library(patchwork)
# Install dataset. only need to run this once..
#InstallData("pbmc3k")

# Read in data ===================================================
pbmc <- LoadData("pbmc3k")
# Counts contains the raw counts
pbmc@assays$RNA$counts[c("CD3D", "TCL1A", "MS4A1"), 1:5]
# data should contain the normalised counts (but we haven't normalised yet!)
pbmc@assays$RNA$data[c("CD3D", "TCL1A", "MS4A1"), 1:5]

# Metadata ======================================================
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc@meta.data %>% head()
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc$percent.ribo <- PercentageFeatureSet(pbmc, pattern = "^RB-")
pbmc@meta.data %>% head()

# Standard Seurat workflow ======================================================

## Normalisation and scaling -------------------
# Normalised counts are stored in the "data" layer
pbmc <- NormalizeData(object = pbmc, vars.to.regress = 'percent.mt')
# The default method is "vst" which uses raw counts to calculate the highly variable genes
# You can change it with the argument "method".
pbmc <- FindVariableFeatures(object = pbmc)
pbmc@assays$RNA@meta.data[21:31,]
table(pbmc@assays$RNA@meta.data$vf_vst_counts_variable) # How many variable genes? This can be changed in FindVariableFeatures()
pbmc@assays$RNA@meta.data$var.features %>% tail(50)
# You can also visualise them
VariableFeaturePlot(pbmc) 
# Scaled, normalised counts are stored in the "scale.data" layer
# By default ScaleData() will scale the normalised counts in the "data" layer.
# If you haven't normalised then it will use raw counts.
pbmc <- ScaleData(object = pbmc)

# layer counts contains the raw counts
pbmc@assays$RNA$counts[c("CD8A", "TCL1A", "MS4A1"), 20:25]
# layer data contains the normalised counts 
pbmc@assays$RNA$data[c("CD8A", "TCL1A", "MS4A1"), 20:25]
# layer scale.data contains the normalised counts 
pbmc@assays$RNA$scale.data[c("TCL1A", "MS4A1"), 20:25]

# Why is CD8A missing from our scaled data?
c("CD8A", "TCL1A", "MS4A1") %in% pbmc@assays$RNA@meta.data$var.features 
# Looks like it's not a highly variable gene!

# If we wanted to scale the data of all features, run:
#pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Plot the expression of 2 genes
FeatureScatter(pbmc, feature1 = "CD3D", feature2 = "MS4A1") +
  theme_bw()

## PCA -------------------
pbmc <- RunPCA(object = pbmc)
pbmc@reductions$pca %>% head()
#  In this example, we can observe an ‘elbow’ around PC9-10, 
# suggesting that the majority of true signal is captured in the first 10 PCs.
ElbowPlot(pbmc)

## Nearest neighbours --------------
# Find nearest neighbours for each cell
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
names(pbmc@graphs) # Where is the data stored?
 
# unweighted graph
pheatmap(pbmc@graphs$RNA_nn[1:200, 1:200],
         col = c("white", "black"), border_color = "grey90", main = "KNN graph",
         legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2
)

# weighted graph
pheatmap::pheatmap(pbmc@graphs$RNA_snn[1:200, 1:200],
         col = colorRampPalette(c("white", 'darkred', "black"))(100),
         border_color = "grey90", main = "SNN graph",
         legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2
)

## Clustering --------------
?FindClusters 
# Clustering with louvain (algorithm 1) and a few different resolutions
pbmc <- FindClusters(pbmc, graph.name = "RNA_snn", resolution = c(0.1, 0.25, .5, 1, 1.5, 2), algorithm = 1)

# Where are the clusters stored?
# Each time you run clustering, the data is stored in metadata columns:
# seurat_clusters - lastest results only
# RNA_snn_res.XX - for each different resolution you test.
colnames(pbmc@meta.data)
head(pbmc@meta.data)
unique(pbmc$RNA_snn_res.0.1)
unique(pbmc$RNA_snn_res.2)

## Dimensionality reduction ------------
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
pbmc <- RunTSNE(object = pbmc, dims = 1:10)
names(pbmc@reductions)
pbmc@reductions$umap@cell.embeddings %>% head()

# Plot the UMAP
DimPlot(object = pbmc, reduction = "umap") 
DimPlot(object = pbmc, reduction = "umap") + theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
DimPlot(object = pbmc, reduction = "umap", group.by = 'RNA_snn_res.0.1') + 
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
DimPlot(object = pbmc, reduction = "umap", 
        group.by = c('RNA_snn_res.0.1', 'RNA_snn_res.0.25', 'RNA_snn_res.0.5', 
                     'RNA_snn_res.1', 'RNA_snn_res.2')) & 
  theme_bw() &
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

DimPlot(object = pbmc, reduction = "umap", 
        group.by = c('RNA_snn_res.0.1', 'RNA_snn_res.0.25', 'RNA_snn_res.0.5', 
                     'RNA_snn_res.1', 'RNA_snn_res.2', 'seurat_annotations')) & 
  theme_bw() &
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# Check for the expression of different markers
Idents(pbmc) <- 'RNA_snn_res.0.5'
p1 <- DotPlot(pbmc, c("CD8A", "FCGR3A", "MS4A1", "CD14"))
p2 <- FeaturePlot(pbmc, c("CD8A", "FCGR3A", "MS4A1", "CD14"))
p1 | p2

## DGE analysis -------------------------------
# Set the active identity 
Idents(pbmc) <- 'RNA_snn_res.0.1'
DimPlot(object = pbmc, reduction = "umap")
# DGE all clusters vs rest
dge_all <- FindAllMarkers(pbmc)
colnames(dge_all)
head(dge_all)
dge_all %>% group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val) %>%
  slice_head(n = 5) %>%
  print(n = 50)

# Or compare one cluster against another
DimPlot(pbmc)
dge_0vs2 <- FindMarkers(pbmc, ident.1 = '0', ident.2 = '2')
colnames(dge_0vs2)
dge_0vs2 %>% arrange(desc(avg_log2FC), p_val) %>% head(10)

# Once you have your annotations, you can compare one cell type against another
pbmc$seurat_annotations %>% unique()
Idents(pbmc) <- 'seurat_annotations'
DimPlot(object = pbmc, reduction = "umap")
dge_CD4naivevsmem <- FindMarkers(pbmc, ident.1 = 'Naive CD4 T', ident.2 = 'Memory CD4 T')
head(dge_CD4naivevsmem)

# SCTransform ==================================
pbmc <- SCTransform(object = pbmc, vars.to.regress = 'percent.mt')
pbmc <- RunPCA(object = pbmc, dims = 1:10)
pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = c(0.1, 0.25, .5, 1, 1.5, 2))
pbmc <- RunUMAP(object = pbmc, dims = 1:30)
DimPlot(object = pbmc, reduction = "umap", group.by = "SCT_snn_res.0.1")
colnames(pbmc@meta.data)
# layer counts contains the raw counts
pbmc@assays # Now we have 2 assays!

# The raw, normalised and scaled counts are different in the SCTransform assay.
pbmc@assays$RNA$counts[c("CD8A", "TCL1A", "MS4A1"), 20:25]
pbmc@assays$SCT$counts[c("CD8A", "TCL1A", "MS4A1"), 20:25]
# layer data contains the normalised counts 
pbmc@assays$RNA$data[c("CD8A", "TCL1A", "MS4A1"), 20:25]
pbmc@assays$SCT$data[c("CD8A", "TCL1A", "MS4A1"), 20:25]
# layer scale.data contains the normalised counts 
pbmc@assays$RNA$scale.data[c("TCL1A", "MS4A1"), 20:25]
pbmc@assays$SCT$scale.data[c("TCL1A", "MS4A1"), 20:25]

qs::qsave(pbmc, "C:/Users/laura/Documents/Biostatsquid/Projects/data/scRNAseq/pbmc3k.qs")
