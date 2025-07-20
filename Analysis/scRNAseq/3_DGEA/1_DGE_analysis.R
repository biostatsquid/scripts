# ======================== # 
# 1_DGE_analysis.R
# ======================== # 
# DGE Analysis with Seurat: Wilcoxon, DESeq2, and MAST
# This tutorial demonstrates three methods of differential gene expression analysis 
# using the IFN stimulation dataset from SeuratData
# Based on Seurat's vignette: https://satijalab.org/seurat/articles/de_vignette 

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Load libraries
library(Seurat)
library(SeuratData)
library(DESeq2)
library(MAST)
library(tidyverse)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

# Install the IFN dataset if not already available
if (!"ifnb" %in% AvailableData()[,"Dataset"]) {
  InstallData("ifnb")
}

# Set relevant paths
project_path <- "C:/Users/laura/Documents/Biostatsquid/Projects"
#in_path <- file.path(project_path, 'data/scRNAseq') # DGE results path
out_path <- file.path(project_path, "DGE/results") # DGE output path
dir.create(out_path, showWarnings = F, recursive = T)

# Load data =============================================
ifnb <- LoadData("ifnb")
print(paste("ifnb dataset dimensions:", dim(ifnb)[1], "genes x", dim(ifnb)[2], "cells"))
head(ifnb)
table(ifnb$stim)

# Basic preprocessing =============================================
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 2000)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb, features = VariableFeatures(object = ifnb))
ifnb <- RunUMAP(ifnb, dims = 1:20)
ifnb <- FindNeighbors(ifnb, dims = 1:20)
ifnb <- FindClusters(ifnb, resolution = 0.5)
Idents(ifnb) <- "seurat_annotations"
table(ifnb$seurat_annotations, ifnb$stim)

# Visualize the clusters =============================================
p1 <- DimPlot(ifnb, reduction = "umap", group.by = "seurat_annotations", label = TRUE) + 
  ggtitle("Clusters")
p2 <- DimPlot(ifnb, reduction = "umap", group.by = "stim", label = FALSE) + 
  ggtitle("Stimulation")
p1 + p2

# Method 1: Wilcoxon Test (default in Seurat) ====

# Compare CD14 Mono vs CD16 Mono within CTRL samples
ctrl_cells <- subset(ifnb, subset = stim == "CTRL")
DimPlot(ctrl_cells, reduction = "umap", group.by = "seurat_annotations", label = FALSE) 
# Find markers using Wilcoxon test
wilcox_results <- FindMarkers(ctrl_cells, 
                              ident.1 = "CD16 Mono", 
                              ident.2 = "CD14 Mono",
                              test.use = "wilcox", 
                              min.pct = 0.25)
# View top markers
cat("Wilcoxon results - top 10 genes:\n")
head(wilcox_results, n = 10)

# Visualize top DEGs
top_genes <- rownames(wilcox_results)[1:6]
VlnPlot(ctrl_cells, features = top_genes, idents = c("CD16 Mono", "CD14 Mono"), 
        ncol = 3, cols = c('red4', 'dodgerblue'))
# Feature plots of top genes
FeaturePlot(ctrl_cells, features = top_genes[1:2], ncol = 2)

# Method 2: MAST with mixed-effects model ============================
# Using MAST to account for technical and biological covariates

# Set identities for MAST
ifnb$celltype.stim <- paste(Idents(ifnb), ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
cat("Available cell type-condition combinations:\n")
unique(ifnb$celltype.stim)
head(ifnb)
# Select T cells
t_cells <- subset(ifnb, idents = c("CD4 Naive T_CTRL", "CD4 Naive T_STIM"))

# Run MAST with covariates
mast_results <- FindMarkers(
  t_cells,
  ident.1 = "CD4 Naive T_STIM", # Note that ident.1 should be the group of interest
  ident.2 = "CD4 Naive T_CTRL",
  test.use = "MAST",
  latent.vars = "nCount_RNA" # Include UMI count as a technical covariate
)
write.table(mast_results, file.path(out_path, 'cd4t_stimvsctrl_dge_mast.csv'),
            sep = '\t', col.names = T, row.names = T, quote = F)
head(mast_results)
# View top MAST results
head(mast_results, n = 10)

# Visualize top DEGs from MAST
top_mast_genes <- rownames(mast_results)[1:5]
VlnPlot(t_cells, features = top_mast_genes, 
        split.by = "stim", ncol = 3)

# Feature plots of top MAST genes
FeaturePlot(t_cells, features = top_mast_genes[1:2], 
            split.by = "stim", ncol = 2)

# Method 3: Pseudobulk + DESeq2 =====================================

# Create a pseudobulk dataset by aggregating cells from each cluster per sample
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, 
                                   group.by = c("stim", "seurat_annotations"))
head(pseudo_ifnb)
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
unique(pseudo_ifnb$celltype.stim)
Idents(pseudo_ifnb) <- "celltype.stim"
bulk.cd4T.de <- FindMarkers(object = pseudo_ifnb, 
                            ident.1 = "CD4 Naive T_STIM", 
                            ident.2 = "CD4 Naive T_CTRL",
                            test.use = "DESeq2")
# Each group (cell_type, condition) must have ≥2 replicates.
# Right now, each cell_type_condition combination (e.g., "CD4T_IFN") appears only once — 
# no replicates = DESeq2 cannot estimate dispersion.

## Add metadata ----------------------------
# First, we need to retrieve the sample information for each cell.
# load the inferred sample IDs of each cell
ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
info <- rbind(ctrl, stim)
# rename the cell IDs by substituting the '-' into '.'
info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)
# only keep the cells with high-confidence sample ID
info <- info[grep(pattern = "SNG", x = info$BEST), ]
# remove cells with duplicated IDs in both ctrl and stim groups
info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]
# now add the sample IDs to ifnb 
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
ifnb <- AddMetaData(ifnb, metadata = info)
head(ifnb)
# remove cells without donor IDs
ifnb$donor_id[is.na(ifnb$donor_id)] <- "unknown"
ifnb <- subset(ifnb, subset = donor_id != "unknown")
# Reset identities
Idents(ifnb) <- "seurat_annotations"
table(ifnb$donor_id)
## Run DESeq2 ----------------------------
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, 
                                   group.by = c("stim", "donor_id", "seurat_annotations"))
head(pseudo_ifnb)
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
unique(pseudo_ifnb$celltype.stim)
Idents(pseudo_ifnb) <- "celltype.stim"
bulk.cd4T.de <- FindMarkers(object = pseudo_ifnb, 
                                ident.1 = "CD4 Naive T_STIM", 
                                ident.2 = "CD4 Naive T_CTRL",
                                test.use = "DESeq2")
head(bulk.cd4T.de, n = 15)
bulk.cd4T.de <- bulk.cd4T.de[order(bulk.cd4T.de$p_val_adj), ]

# Visualize top DEGs from DESeq2
top_deseq2_genes <- rownames(bulk.cd4T.de)[1:5]
VlnPlot(ifnb, features = top_deseq2_genes, 
        idents = "CD4 Naive T", split.by = "stim", ncol = 3)

## Plot DGE results =======================================================
head(mast_results)
# Add a column for -log10(p_val_adj), handling p_val_adj = 0 by replacing with a small value
mast_results <- mast_results %>% mutate(log10_pval_adj = -log10(p_val_adj + 1e-300),
                                        significant = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5,
                                        diff_expressed = case_when(
                                          p_val_adj < 0.05 & avg_log2FC > 0.5 ~ 'UP',
                                          p_val_adj < 0.05 & avg_log2FC < -0.5 ~ 'DOWN',
                                          TRUE ~ 'n.s.'
                                        ))
head(mast_results)
table(mast_results$diff_expressed)

# Optionally select top genes for labeling
top_genes <- mast_results[order(mast_results$p_val_adj), ]
top_genes <- head(top_genes[mast_results$significant, ], 20)  # Label top 10 significant genes
top_genes$gene_labels <- rownames(top_genes)

# Create the volcano plot
ggplot(mast_results, aes(x = avg_log2FC, y = log10_pval_adj)) +
  geom_point(aes(color = diff_expressed), size = 3) +
  scale_color_manual(values = list('n.s.' = "grey", 'DOWN' = "dodgerblue", 'UP' = "red4")) +
  ggrepel::geom_text_repel(data = top_genes, aes(label = gene_labels),
                  size = 6, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3) +
  theme_bw() +
  labs(title = "CD4 Naive T cells STIM vs CONTROL",
       x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"Adjusted P-value"),
       color = '') +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "bottom"
        )

# Test for differences in cell proportions ========================== #

# Check if IFN stimulation changes cell type proportions
table(ifnb$seurat_annotations, ifnb$stim) %>% 
  prop.table(margin = 2) %>% 
  round(3)

# Simple proportion comparison
prop_summary <- ifnb@meta.data %>%
  group_by(stim, seurat_annotations) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(stim) %>%
  mutate(total = sum(n),
         proportion = n / total) %>%
  select(stim, seurat_annotations, proportion) %>%
  pivot_wider(names_from = stim, values_from = proportion, values_fill = 0) %>%
  mutate(fold_change = STIM / CTRL,
         log2_fc = log2(fold_change))
ggplot(prop_summary, aes(x = seurat_annotations, y = log2_fc)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportion Changes (STIM vs CTRL)",
       y = "Log2 Fold Change",
       x = "Cell Type") +
  geom_hline(yintercept = 0, linetype = "dashed")

# Run propeller with custom design
propeller_results <- speckle::propeller(ifnb,
                                        clusters = ifnb$seurat_annotations,
                                        sample = ifnb$donor_id,
                                        group = ifnb$stim)
speckle::plotCellTypeProps(ifnb, 
                           clusters = ifnb$seurat_annotations,
                           sample = ifnb$donor_id)

# ============================================================================ #
# Extra with functions:
# DESeq2 Pseudobulk Analysis 

# Function to create pseudobulk data with proper replicates
create_pseudobulk_seurat <- function(seurat_obj, group_by_vars = c("seurat_annotations", "stim", "donor_id")) {
  
  # Check if all grouping variables exist
  missing_vars <- group_by_vars[!group_by_vars %in% colnames(seurat_obj@meta.data)]
  if (length(missing_vars) > 0) {
    stop("Variables not found in metadata: ", paste(missing_vars, collapse = ", "))
  }
  
  # Use Seurat's AggregateExpression function
  pseudobulk_seurat <- AggregateExpression(
    seurat_obj, 
    assays = "RNA", 
    return.seurat = TRUE,
    group.by = group_by_vars,
    verbose = FALSE
  )
  
  # Extract counts matrix
  counts_matrix <- GetAssayData(pseudobulk_seurat, assay = "RNA", layer = "counts")
  
  # Create metadata from column names
  sample_names <- colnames(counts_matrix)
  group_parts <- strsplit(sample_names, "_")
  
  sample_metadata <- data.frame(
    row.names = sample_names,
    sample_id = sample_names,
    cell_type = sapply(group_parts, function(x) x[1]),
    condition = sapply(group_parts, function(x) x[2]),
    donor = sapply(group_parts, function(x) x[3])
  )
  
  cat("Seurat pseudobulk data:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "samples\n")
  cat("Sample breakdown:\n")
  print(table(sample_metadata$cell_type, sample_metadata$condition))
  
  return(list(counts = counts_matrix, metadata = sample_metadata))
}

# Main analysis function
run_deseq2_analysis <- function(pseudobulk_data, target_cell_type = "CD4 Naive T", 
                                condition_var = "condition", 
                                contrast_levels = c("STIM", "CTRL")) {
  
  # Filter for target cell type
  target_idx <- pseudobulk_data$metadata[[gsub("_.*", "", names(pseudobulk_data$metadata)[2])]] == target_cell_type
  
  if (sum(target_idx) < 4) { # Need at least 2 replicates per condition
    stop("Insufficient replicates for DESeq2. Need at least 2 samples per condition.")
  }
  
  # Check replicates per condition
  condition_counts <- table(pseudobulk_data$metadata[target_idx, condition_var])
  cat("Samples per condition for", target_cell_type, ":\n")
  print(condition_counts)
  
  if (any(condition_counts < 2)) {
    stop("Each condition must have at least 2 replicates for DESeq2")
  }
  
  # Create DESeq dataset
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_data$counts[, target_idx],
    colData = pseudobulk_data$metadata[target_idx, ],
    design = as.formula(paste("~", condition_var))
  )
  
  # Filter low count genes (recommended for single-cell derived pseudobulk)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  
  cat("Filtered to", nrow(dds), "genes for analysis\n")
  
  # Run DESeq2 with parameters optimized for single-cell pseudobulk
  # Based on DESeq2 vignette recommendations for single-cell analysis
  dds <- DESeq(dds, 
               test = "Wald",  # Use Wald test instead of LRT for two-group comparison
               fitType = "parametric",  # Usually works better than glmGamPoi for pseudobulk
               sfType = "ratio",
               minReplicatesForReplace = Inf,  # Don't replace outliers
               useT = TRUE,  # Use t-distribution for small sample sizes
               minmu = 1e-6)
  
  # Get results
  res <- results(dds, 
                 contrast = c(condition_var, contrast_levels[1], contrast_levels[2]),
                 alpha = 0.05)
  
  # Convert to data frame and sort by adjusted p-value
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
  
  # Remove genes with NA padj values
  res_df <- res_df[!is.na(res_df$padj), ]
  
  cat("DESeq2 analysis complete. Found", sum(res_df$padj < 0.05), "significant genes\n")
  
  return(list(results = res_df, dds = dds))
}

# Check donor distribution
print(table(ifnb$donor_id, ifnb$stim))
# Create pseudobulk data using Seurat method (recommended)
pseudobulk_data <- create_pseudobulk_seurat(ifnb, group_by_vars = c("seurat_annotations", "stim", "donor_id"))
print(pseudobulk_data$counts[1:5, 1:5])
print(head(pseudobulk_data$metadata))

# Run DESeq2 analysis 
tryCatch({
  deseq2_analysis <- run_deseq2_analysis(pseudobulk_data, 
                                         target_cell_type = "CD4 Naive T",
                                         condition_var = "condition",
                                         contrast_levels = c("STIM", "CTRL"))
  
  # Display results
  cat("\nTop 10 DESeq2 results:\n")
  print(head(deseq2_analysis$results, n = 10))
  
  # Create summary
  deseq2_results <- deseq2_analysis$results
  n_significant <- sum(deseq2_results$padj < 0.05, na.rm = TRUE)
  n_upregulated <- sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0, na.rm = TRUE)
  n_downregulated <- sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < 0, na.rm = TRUE)
  
  cat("\nDESeq2 Summary:\n")
  cat("Total genes tested:", nrow(deseq2_results), "\n")
  cat("Significant genes (padj < 0.05):", n_significant, "\n")
  cat("Upregulated in STIM:", n_upregulated, "\n")
  cat("Downregulated in STIM:", n_downregulated, "\n")
  
  # Visualization
  if (n_significant > 0) {
    # Get top significant genes
    top_genes <- head(rownames(deseq2_results[deseq2_results$padj < 0.05, ]), 6)
    
    if (length(top_genes) > 0) {
      cat("\nCreating visualization for top genes:", paste(top_genes, collapse = ", "), "\n")
      
      # Set correct identity for visualization
      Idents(ifnb) <- "seurat_annotations"
      
      # Create violin plot
      vlnplot_deseq2 <- VlnPlot(ifnb, 
                                features = top_genes[1:min(3, length(top_genes))], 
                                idents = "CD4 Naive T", 
                                split.by = "stim", 
                                ncol = 3,
                                cols = c("CTRL" = "blue", "STIM" = "red"))
      print(vlnplot_deseq2)
      
      # Create volcano plot
      volcano_data <- deseq2_results %>%
        mutate(
          log10_padj = -log10(pmax(padj, 1e-300)),
          significant = padj < 0.05 & abs(log2FoldChange) > 1,
          direction = case_when(
            padj < 0.05 & log2FoldChange > 1 ~ "Up",
            padj < 0.05 & log2FoldChange < -1 ~ "Down",
            TRUE ~ "NS"
          ),
          gene_name = rownames(.)
        )
      
      # Select genes for labeling
      label_genes <- volcano_data %>%
        filter(significant) %>%
        arrange(padj) %>%
        head(10)
      
      volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = log10_padj)) +
        geom_point(aes(color = direction), alpha = 0.7) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
        geom_text_repel(data = label_genes, 
                        aes(label = gene_name),
                        size = 3, max.overlaps = 10) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.7) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.7) +
        theme_bw() +
        labs(title = "DESeq2: CD4 Naive T cells (STIM vs CTRL)",
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value",
             color = "Direction") +
        theme(plot.title = element_text(hjust = 0.5))
      
      print(volcano_plot)
    }
  }
  
}, error = function(e) {
  cat("Error in DESeq2 analysis:", e$message, "\n")
  cat("This might be due to insufficient replicates or other data issues.\n")
})
