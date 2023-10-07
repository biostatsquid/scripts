# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Volcano_plots
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Date: October 2022
# Version: 1.0

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
rm(list = ls(all.names = TRUE)) 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer)
library(Seurat)
library(scRNAseq)
library(glmGamPoi) # for DGE
library(scran) # for scRNAseq steps
library(ggrepel) # for nice volcano plot

# Set input path
path = "C:/Users/laura/Documents/Biostatsquid/Scripts/Volcano_plots/"
glmgampoi_path <- path
DEG_results_path <- path
setwd(path)

# Theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black')
            ))

# Functions =====================================================================
# Seurat pipeline on expression data 
seurat_pipeline <- function(x){
  #for debug
  #x <- data
  
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 5000)
  x <- ScaleData(x)
  # dim reduction
  x <- RunPCA(x)
  # clustering
  x <- FindNeighbors(x, dims = 1:15)
  #x <- FindClusters(x, resolution = 0.8, algorithm=2)
  
  #x <- RunTSNE(x, dims = 1:15)
  x <- RunUMAP(x, dims = 1:15)
}

# Run glmgampoi for deg
run_glmgampoi <- function(sce_subset, design_formula, padj_method = 'BH', sce_subset_name = 'scesubsetname', reference = NULL) {
  # for debug
  # sce_subset <- sce_epivsdermis_withinterface
  # design_formula <- ~ patient + epivsdermis
  # sce_subset_name <- name_of_subset
  # reference <- reference_group
  
  if (!is.character(sce_subset_name)){
    stop('Provide a valid sce_subset_name for the fit RDS object')
  } 
  
  if (!class(sce_subset)[1] == 'SingleCellExperiment'){
    stop('sce_subset must be of class SingleCellExperiment')
  }
  
  if (class(design_formula) != 'formula'){
    stop('design_formula must be of class formula')
  }
  
  if (length(design_formula[[2]]) == 1){
    print(paste('Your predictor variable is: ', design_formula[[2]]))
    predictor_variable <- design_formula[[2]]
  } else {
    print(paste('Your predictor variable is: ', design_formula[[2]][[length(design_formula[[2]])]]))  
    predictor_variable <- design_formula[[2]][[length(design_formula[[2]])]]
  }
  
  # Checking if there are more than 2 unique values in the predictor variable
  if(length(unique(sce_subset@colData@listData[[predictor_variable]])) > 2){
    print(unique(sce_subset@colData@listData[[predictor_variable]]))
    stop('There are more than 2 unique values in the predictor variable. \n Subset the sce_subset and drop unused levels, change the design formula or change the contrast in test_de')
  }
  
  # the level that is not the reference is usually the 'interesting' category you want to study up or down regulation of
  interesting_category <- unique(c(reference, levels(sce_subset@colData@listData[[predictor_variable]])))[2]
  interesting_category <- stringr::str_to_title(interesting_category)
  
  # Convert counts to dense matrix
  counts(sce_subset) <- as.matrix(counts(sce_subset))
  
  # Create a 'glmGamPoi' object: fit one GLM model for each gene 
  print('Fitting GLM model')
  fit <- glm_gp(sce_subset, 
                design = design_formula,  # design formula: factor of interest to test for during differential expression testing (add last!) and the known sources of variation to control for. -1 to remove the intercept
                on_disk = FALSE, # not to load on memory, it goes faster
                reference_level = reference,
                size_factors = sce_subset$sizeFactor) # to correct for the fact that each sample is typically of different size
  saveRDS(fit, file = paste0(glmgampoi_path, sce_subset_name, '_', gsub(' ', '', Reduce(paste, deparse(design_formula[[2]]))), '_fit.RDS'))
  #summary(fit)
  #colnames(fit$Beta)
  
  # Test for Differential Expression
  # Conduct a quasi-likelihood ratio test for a Gamma-Poisson fit to find differentially expressed genes between conditions
  # The contrast argument specifies what we want to compare, e.g. treatment placebo vs treatment drug. It is set as 'conditionreference', e.g., treatmentPLACEBO
  
  print('Testing for differential expression')
  de_res <- glmGamPoi::test_de(fit, 
                               contrast = tail(colnames(fit$model_matrix), n=1), # Note! For non-binary condition you might need to change this
                               full_design = fit$model_matrix,
                               pseudobulk_by = sce_subset$unique_indexes, # If the data has multiple samples, it is a good idea to aggregate the cell counts by samples. This is called "pseudobulk". In this case we will assume each spot is an individual ST sample
                               pval_adjust_method = padj_method)
  print(paste('DEG successfully computed, condition reference = ', tail(colnames(fit$model_matrix), n=1)))
  
  # Store result in data frame with pre-determined columns
  df <- data.frame(gene_symbol = de_res$name,
                   pval = de_res$pval,
                   padj = de_res$adj_pval,
                   log2fc = de_res$lfc)
  
  # The large `lfc` values come from groups were nearly all counts are 0. Setting them to Inf makes the plots look nicer
  df$log2fc <- ifelse(abs(df$log2fc) > 20, sign(df$log2fc) * Inf, df$log2fc)
  
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
  # add a column of NAs
  df$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  df$diffexpressed[df$log2fc > 0.6 & df$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  df$diffexpressed[df$log2fc < -0.6 & df$padj < 0.05] <- "DOWN"
  df$Note <- NA
  df[1, 'Note'] <- paste('Comparison = ', predictor_variable)
  df[2, 'Note'] <- paste('Upregulation and downregulation in ', interesting_category)
  
  write.csv(df, file = paste0(glmgampoi_path, sce_subset_name, '_', gsub(' ', '', Reduce(paste, deparse(design_formula[[2]]))), '_degresults.csv'))
  
  return(df)
}

# Get DEG summary
# Input = df from run_glmgampoi. 
# Returns a dataframe with the total gene count, # upregulated and downregulated genes for conditions: 'All', 'pval_0.05', 'pval_0.01', 'padj_0.05', 'padj_0.01'
get_dge_summary <- function(dataframe){
  all <- dataframe %>% 
    filter(pval>0) %>%
    summarise('total_gene_count' = sum(pval>0), 
              'log2fc>0_up' = sum(pval>0 & log2fc > 0),
              'log2fc<0_down' = sum(pval>0 & log2fc < 0))
  pval_0.05 <- dataframe %>% 
    filter(pval<0.05) %>%
    summarise('total_gene_count' = sum(pval<0.05), 
              'log2fc>0_up' = sum(pval<0.05 & log2fc > 0),
              'log2fc<0_down' = sum(pval<0.05 & log2fc < 0))
  pval_0.01 <- dataframe %>% 
    filter(pval<0.01) %>%
    summarise('total_gene_count' = sum(pval<0.01), 
              'log2fc>0_up' = sum(pval<0.01 & log2fc > 0),
              'log2fc<0_down' = sum(pval<0.01 & log2fc < 0))
  padj_0.05 <- dataframe %>% 
    filter(padj<0.05) %>%
    summarise('total_gene_count' = sum(padj<0.05), 
              'log2fc>0_up' = sum(padj<0.05 & log2fc > 0),
              'log2fc<0_down' = sum(padj<0.05 & log2fc < 0))
  padj_0.01 <- dataframe %>% 
    filter(padj<0.01) %>%
    summarise('total_gene_count' = sum(padj<0.01), 
              'log2fc>0_up' = sum(padj<0.01 & log2fc > 0),
              'log2fc<0_down' = sum(padj<0.01 & log2fc < 0))
  summary_df <- rbind(all, pval_0.05, pval_0.01, padj_0.05, padj_0.01)
  rownames(summary_df) <- c('All', 'pval_0.05', 'pval_0.01', 'padj_0.05', 'padj_0.01')
  return(summary_df)
}

# Import data =====================================================================
#data <- BacherTCellData(filtered = TRUE, ensembl = FALSE, location = TRUE)

# # Obtain expression matrix
# exprsMat <- data@assays@data@listData[["counts"]] # matrix of expression counts, where each row is a gene and each column a cell
# dim(exprsMat) #33538  rows (genes) x 104417 columns (cells)
# 
# # Remove all genes that are expressed in <5% of all cells
# exprsMat <- exprsMat[apply(exprsMat == 0, 1, sum) <= 0.97 * ncol(exprsMat), ]
# dim(exprsMat) #12216 rows (genes) x 1977 columns (cells)
# 
# #Obtain cell types (column names)
# data_ct <- data@colData@listData[["new_cluster_names"]] # cell types
# data_cids <- 1:length(data_ct) # cell ids
# colnames(exprsMat) <- paste(data_ct, data_cids, sep = "..")
# head(rownames(exprsMat), 10) # genes
# head(colnames(exprsMat), 10) # CellID_Celltype
# 
# ct_levels <- unique(data_ct) # unique cell types of the dataset (6)
# ct_levels
# 
# #Exploring data: see how many cells of each type:
# summary(as.factor(gsub("\\.\\..*$", "", colnames(exprsMat))))
# # Central memory              Cycling      Cytotoxic / Th1             Tfh-like  Transitional memory Type-1 IFN signature 
# # 18391                  606                14262                45732                24165                 1261 
# 
# #Convert to Seurat
# data_seu <- CreateSeuratObject(counts = exprsMat) #Note: colnames should be cell ids, we're working with celltype_cellids
# data_seu  <- seurat_pipeline(data_seu)
# 
# saveRDS(data_seu, file = 'data_seu.RDS')

# Import Seurat object directly
#data_seu <- readRDS('data_seu.RDS')

# DGE --------------------------------------------------------------

## Tfh-like cells: severe vs healthy ####

### 1. Sub-setting data ####

name_of_subset <- 'severevshealthy'
reference_group <- 'healthy'
sce_tfhlike <- data[, which(data$new_cluster_names %in% c('Tfh-like'))]

sce_tfhlike <- sce_tfhlike[, which(sce_tfhlike$severity %in% c("healthy", "severe"))]
# Remove empty levels with droplevels if necessary because glm_gp() will complain otherwise
sce_tfhlike$severity <- as.factor(sce_tfhlike$severity)

#Add size_factors
sce_tfhlike <- computeSumFactors(sce_tfhlike)

### 2. Exploratory analysis - skip for now
### 3. Running DGE: glmGamPoi ####

#sce_tfhlike <- run_glmgampoi(sce_tfhlike, design_formula = ~ batch + severity, sce_subset_name = name_of_subset, reference = reference_group)
list.files(path)
dge_thflike <- read.csv(paste0(glmgampoi_path, 'severevshealthy_batch+severity_degresults.csv'), row.names = 1)
fit_thflike <- readRDS(paste0(glmgampoi_path, 'severevshealthy_batch+severity_fit.RDS'))

get_dge_summary(dge_thflike)

# Most different genes (padj)
head(dge_thflike[order(dge_thflike$padj), ])



dataframe <- dge_thflike
interesting_category <- gsub('Upregulation and downregulation in  ', '', unique(dataframe$Note)[2])
title_of_comparison <- 'Tfh-like in severe Covid versus healthy patients'

# VOLCANO PLOT
# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dataframe$delabel <- ifelse(dataframe$gene_symbol %in% head(dataframe[order(dataframe$padj), "gene_symbol"], 30), dataframe$gene_symbol, NA)

# Calculate the -log10(padj) to be plotted
dataframe$minuslog10padj <- -log10(dataframe$padj)
# The max_y_scale has to be adjusted because it cannot be infinity
max_y_scale <- max(dataframe[!is.infinite(dataframe$log2fc),]$minuslog10padj) * 1.5
infinity_line <- FALSE
genes_with_inf_minuslogpadj <- NULL
if(is.infinite(max_y_scale)){
  max_y_scale <- max(dataframe[!is.infinite(dataframe$log2fc),]$minuslog10padj[!is.infinite(dataframe[!is.infinite(dataframe$log2fc),]$minuslog10padj)])
  max_y_scale <- max_y_scale * 2
  
  genes_with_inf_minuslogpadj <- dataframe[is.infinite(dataframe$minuslog10padj),'gene_symbol']
  print('Genes with really high -log10(FDR): ')
  print(genes_with_inf_minuslogpadj)
  dataframe[is.infinite(dataframe$minuslog10padj),'minuslog10padj'] <- max_y_scale * 1.9/2 # set a fake value for minuslog10 infinity
  infinity_line <- FALSE # this adds a horizontal line marking limits of genes with -log20FDR of around infinity
}

legend_title <- interesting_category

# Change the title of the legend depending on the category
switch (interesting_category,
        'severe' = legend_title <-'Tfh-like in severe covid'
)

max_y_scale <- 100
# plot adding up all layers we have seen so far
volcano_plot <- ggplot(data = dataframe[!is.infinite(dataframe$log2fc),], 
                       aes(x = log2fc, y = -log10(pval), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  coord_cartesian(ylim = c(0, max_y_scale), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set this limits
  #theme_classic() +
  #geom_text_repel(max.overlaps = Inf) + # To show all labels
  #scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated")) +
  #ggtitle(paste(title_of_comparison)) + 
  labs(color = 'Condition', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), labels = c("Downregulated", "Not significant", "Upregulated"))

# # Let's fix the legend separately to deal with exceptions
# if(interesting_category == 'Ga_border'|interesting_category == 'Sa_border'){
#   volcano_plot <- volcano_plot + 
#     scale_color_manual(values = c("red", "black", "#008000"), labels = c("Upregulated in core", "Unchanged", "Upregulated in border")) 
# } else {
#   volcano_plot <- volcano_plot + 
#     scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated"))
# }

volcano_plot <- volcano_plot + 
  scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Not significant", "Upregulated"))


# Add infinity line and readjust genes with minuslog10padj of infinity so they show on the plot
if(infinity_line == TRUE){
  volcano_plot <- volcano_plot + 
    geom_hline(aes(yintercept = max_y_scale * 1.8/2), color = "gray", linetype = "dashed") + 
    geom_text(x = 5, aes(0, max_y_scale * 1.8/2), label = 'Infinity threshold', vjust = 1.5, color = "gray")
}


#### raw p values
ggplot(data = dataframe[!is.infinite(dataframe$log2fc),], 
       aes(x = log2fc, y = pval, col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  coord_cartesian(ylim = c(0, 1), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set this limits
  #theme_classic() +
  #geom_text_repel(max.overlaps = Inf) + # To show all labels
  #scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated")) +
  #ggtitle(paste(title_of_comparison)) + 
  labs(color = 'Condition', #legend_title, 
       x = expression("log"[2]*"FC"), y = "p-value") +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), labels = c("Downregulated", "Not significant", "Upregulated"))






# Plot DEG summary statistics
summarise_deg_results <- function(dataframe, title_of_comparison, sce_subset = NULL){
  # For debug
  # dataframe <- dge_epivsdermis_withinterface
  # title_of_comparison <- comparison_title
  # sce_subset <- sce_epivsdermis_withinterface
  
  interesting_category <- gsub('Upregulation and downregulation in  ', '', unique(dataframe$Note)[2])
  
  # Distribution of p values values
  pval_plot <- ggplot(dataframe, aes(x = pval)) + 
    geom_density() + 
    theme_bw() + 
    ggtitle(paste(title_of_comparison, "\nDistribution of p values"))
  
  # Distribution of p adj values
  padj_plot <- ggplot(dataframe, aes(x = padj)) + 
    geom_density() + 
    theme_bw() + 
    ggtitle(paste(title_of_comparison, "\nDistribution of p adj values"))
  
  # Distribution of p adj values
  log2fc_plot <- ggplot(dataframe, aes(x = log2fc)) + 
    geom_density() + 
    theme_bw() + 
    ggtitle(paste(title_of_comparison, "\nDistribution of log2fc values"))
  
  # VOLCANO PLOT
  # Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
  dataframe$delabel <- ifelse(dataframe$gene_symbol %in% head(dataframe[order(dataframe$padj), "gene_symbol"], 30), dataframe$gene_symbol, NA)
  
  # Calculate the -log10(padj) to be plotted
  dataframe$minuslog10padj <- -log10(dataframe$padj)
  # The max_y_scale has to be adjusted because it cannot be infinity
  max_y_scale <- max(dataframe[!is.infinite(dataframe$log2fc),]$minuslog10padj) * 1.5
  infinity_line <- FALSE
  genes_with_inf_minuslogpadj <- NULL
  if(is.infinite(max_y_scale)){
    max_y_scale <- max(dataframe[!is.infinite(dataframe$log2fc),]$minuslog10padj[!is.infinite(dataframe[!is.infinite(dataframe$log2fc),]$minuslog10padj)])
    max_y_scale <- max_y_scale * 2
    
    genes_with_inf_minuslogpadj <- dataframe[is.infinite(dataframe$minuslog10padj),'gene_symbol']
    print('Genes with really high -log10(FDR): ')
    print(genes_with_inf_minuslogpadj)
    dataframe[is.infinite(dataframe$minuslog10padj),'minuslog10padj'] <- max_y_scale * 1.9/2 # set a fake value for minuslog10 infinity
    infinity_line <- FALSE # this adds a horizontal line marking limits of genes with -log20FDR of around infinity
  }
  
  legend_title <- interesting_category
  
  # Change the title of the legend depending on the category
  switch (interesting_category,
          'severe' = legend_title <-'Tfh-like in severe covid'
  )
  
  # plot adding up all layers we have seen so far
  volcano_plot <- ggplot(data = dataframe[!is.infinite(dataframe$log2fc),], 
                         aes(x = log2fc, y = minuslog10padj, col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point() + 
    coord_cartesian(ylim = c(0, max_y_scale)) + # since some genes can have minuslog10padj of inf, we set this limits
    theme_classic() +
    geom_text_repel(max.overlaps = Inf) + # To show all labels
    #scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated")) +
    ggtitle(paste(title_of_comparison)) + 
    labs(color = legend_title, x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR"))
  
  # # Let's fix the legend separately to deal with exceptions
  # if(interesting_category == 'Ga_border'|interesting_category == 'Sa_border'){
  #   volcano_plot <- volcano_plot + 
  #     scale_color_manual(values = c("red", "black", "#008000"), labels = c("Upregulated in core", "Unchanged", "Upregulated in border")) 
  # } else {
  #   volcano_plot <- volcano_plot + 
  #     scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated"))
  # }
  
  volcano_plot <- volcano_plot + 
    scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Not significant", "Upregulated"))
  
  
  # Add infinity line and readjust genes with minuslog10padj of infinity so they show on the plot
  if(infinity_line == TRUE){
    volcano_plot <- volcano_plot + 
      geom_hline(aes(yintercept = max_y_scale * 1.8/2), color = "gray", linetype = "dashed") + 
      geom_text(x = 5, aes(0, max_y_scale * 1.8/2), label = 'Infinity threshold', vjust = 1.5, color = "gray")
  }
  
  # Save volcano plot separately
  pdf(file = paste0(DEG_results_path, gsub(' ', '', tolower(title_of_comparison)), '_volcanoplot.pdf'))
  print(volcano_plot)
  dev.off()
  
  # Additionally we will plot the distribution of the genes that had a -log10padj of infinity, to make sure nothing funny is going on
  # we will only plot 10 of them for now
  if(!is.null(sce_subset) & !is.null(genes_with_inf_minuslogpadj)){
    comparison_column <- gsub('Comparison =  ', '', dataframe$Note[1])
    genes_distribution_plots <- scater::plotExpression(sce_subset,
                                                       features = genes_with_inf_minuslogpadj[1:10][!is.na(genes_with_inf_minuslogpadj[1:10])], #rownames(sce_cluster0vshealthydermis)[1:6], # the genes to plot
                                                       x = comparison_column, #rownames(sce_cluster0vshealthydermis)[10], #metadata to show in x axis
                                                       exprs_values = 'counts',
                                                       colour_by = comparison_column)
    pdf(file = paste0(DEG_results_path, comparison_column, '_genedistributionplots_top10inflogfdr.pdf'), height = 15, width = 10)
    print(genes_distribution_plots)
    dev.off()
  }
  
  summary_plot <- ggarrange(pval_plot,
                            padj_plot,
                            log2fc_plot,
                            ggtexttable(get_dge_summary(dataframe)),
                            volcano_plot,
                            nrow = 3, ncol = 2, heights = c(1,1,1.7), widths = c(1,1))
  
  return(summary_plot)
}


summarise_deg_results(dge_thflike, 'Tfh-like in severe Covid versus healthy patients')


