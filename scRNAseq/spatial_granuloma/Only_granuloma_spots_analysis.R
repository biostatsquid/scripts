# Name: Only_granuloma_spots_analysis.R
# Author: Laura Twomey
# Date of creation: 19 July 2022
# Clustering and GSEA analysis of only granuloma spots

# Can we further separate border and core? 

# Clean Environment ----------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects including hidden objects
gc() #free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F) #avoid truncated output in R console and scientific notation

# Set working directory ----------------------------------------------------
path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Granulomas_only/"
setwd(path)
exploratory_analysis_path = "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Granulomas_only/DEG_exploratory/"
DEG_results_path = "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Granulomas_only/DEG_results/"
glmgampoi_path = "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Granulomas_only/glmGamPoi_results/"
pea_path = "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Granulomas_only/PEA/"

# Libraries and functions ------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(anndata)
library(reticulate)
library(hdf5r)
library(pheatmap)
library(fgsea)
library(umap)
library(ggpubr) # for ggarrange
library(SingleCellExperiment)
library(glmGamPoi)
library(ggrepel)

## Function: Patterns to python format 
to_python_format <- function(list_of_patterns){
  list_in_python_format <- cat(paste(shQuote(gsub('\\.', '_pattern', list_of_patterns), type="cmd"), collapse=", "))
  return(list_in_python_format)
}

to_python_format2 <- function(list_of_patterns){
  list_in_python_format <- cat(paste(shQuote(gsub('\\.', '-', list_of_patterns), type="cmd"), collapse=", "))
  return(list_in_python_format)
}

# Run glmgampoi for deg
run_glmgampoi <- function(sce_subset, design_formula, padj_method = 'BH', sce_subset_name = 'subsetname', reference = NULL) {
  # for debug
  # sce_subset <- sce_cluster0vshealthydermis
  # design_formula <- ~ patient + cluster0vshealthydermis
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
  
  # the level that is not the reference is usually the interesting category you want to study up or down regulation of
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
                size_factors = sce_subset$size_factors) # to correct for the fact that each sample is typically of different size
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

# Plot DEG summary statistics
summarise_deg_results <- function(dataframe, title_of_comparison, sce_subset = NULL){
  # For debug
  # dataframe <- dge_cluster2vshealthydermis
  # title_of_comparison <- comparison_title
  # sce_subset <- sce_cluster2vshealthydermis
  interesting_category <- gsub('Upregulation and downregulation in  ', '', unique(dataframe$Note)[2])
  
  # Custom names for the legend title
  interesting_category <- gsub('Cluster3', 'Border-1', gsub('Cluster2', 'Inner core', gsub('Cluster1', 'Border-2', gsub('Cluster0', 'Outer core', interesting_category))))

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
  
  # plot adding up all layers we have seen so far
  volcano_plot <- ggplot(data = dataframe[!is.infinite(dataframe$log2fc),], 
                         aes(x = log2fc, y = minuslog10padj, col = diffexpressed, label = delabel)) +
    geom_point() + 
    coord_cartesian(ylim = c(0, max_y_scale)) + # since some genes can have minuslog10padj of inf, we set this limits
    theme_classic() +
    geom_text_repel(max.overlaps = Inf) + # To show all labels
    scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated")) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") + 
    ggtitle(paste(title_of_comparison)) + 
    labs(color = interesting_category, x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) 
  
  # Add infinity line and readjust genes with minuslog10padj of infinity so they show on the plot
  if(infinity_line == TRUE){
    volcano_plot <- volcano_plot + 
      geom_hline(aes(yintercept = max_y_scale * 1.8/2), color = "gray", linetype = "dashed") + 
      geom_text(x = 5, aes(0, max_y_scale * 1.8/2), label = 'Infinity threshold', vjust = 1.5, color = "gray")
  }
  
  # Save volcano plot separately
  pdf(file = paste0(DEG_results_path, gsub(' ', '', title_of_comparison), '_volcanoplot.pdf'))
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

# This function takes the normalised counts, the column name from the normalised_counts headers dataframe 
# with the values for the comparison column, and the names of the two values 
# (groups) we want to compare, and returns a heatmap
# dge_df and number_top_genes are needed to plot only the top genes (e.g., top 50 genes sorted by padj from the DEG results)
get_deg_heatmap <- function(norm_counts, comparison_column, group1, group2, dge_df, number_top_genes){
  #For debug
  # norm_counts <- normalised_counts
  # comparison_column <- 'dge_cluster0vshealthydermis'
  # group1 <- 'dermis_lesional'
  # group2 <- 'dermis_nonlesional'
  # dge_df <- dge_dermislesionalvsnonlesional
  # number_top_genes <- 50
  
  # Some checks to avoid weird errors
  if(sum(c(group1, group2) %in% unique(normalised_counts_headers[[comparison_column]])) < 2){
    stop('Recheck your comparison groups, they must be contained in the comparison column')
  }
  
  # First we must filter the genes for the two groups
  heatmapdf <- norm_counts # get the normalised counts
  colnames(heatmapdf) <- normalised_counts_headers[[comparison_column]] # set as column names the values of the comparison variable
  heatmapdf <- as.data.frame(do.call(cbind, # merge columns with the same name and calculate the mean expression value for all spots of each group
                                     by(t(heatmapdf), INDICES = names(heatmapdf), FUN = colMeans)))
  # Subset to only the two groups we want to compare
  heatmapdf <- heatmapdf[, grepl(paste0(group1, '|', group2), colnames(heatmapdf))]
  # Subset to only the top genes
  top_genes <- head(dge_df[order(dge_df$padj), "gene_symbol"], number_top_genes) # get the 50 most different genes (padj)
  heatmapdf <- heatmapdf[top_genes, ]
  
  # Plot heatmap
  heatmap <- pheatmap(heatmapdf,
                      show_colnames = T, show_rownames = T,
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D')
  return(heatmap)
}


## Function: get_gsea_analysis_spatialde
# This function runs gsea analysis for a set of background genes (for now KEGG, Reactome or GO), given a list of deg spatialde_results from 
# the spatialDE pipeline (list of genes are in $genes column). Pathways can be filtered by padj and gene count.
# We will use signed p value (sign(deg_results$log2fc)*(-log10(deg_results$pval))) as rankings
run_gsea_analysis_deg <- function(background_genes, deg_results, filename = 'example', number_of_top_pathways_up = 25, number_of_top_pathways_down = 5){
  # for debug
  # background_genes <- 'GO'
  # deg_results <- dge_cluster0vshealthydermis
  # filename <- name_of_subset
  
  # Load in background genes
  if(background_genes == 'KEGG'){
    pwl1 <- readRDS('Background_genes/kegg.RData')
    background_genes <- '_kegg'
  } else if(background_genes == 'Reactome'){
    pwl1 <- readRDS('Background_genes/reactome.RData')
    background_genes <- '_reactome'
  } else if(background_genes == 'GO'){
    pwl1 <- readRDS('Background_genes/go.bp.RData')
    background_genes <- '_gobp'
  } else {
    stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
  }
  pwl1$term <- as.character(pwl1$term) # gsea doesn't like it when the entries are factors
  names(pwl1) <- gsub('GOBP_', '', gsub('KEGG_', '', gsub('REACTOME_', '', names(pwl1))))
  
  ## Prepare input: FGSEA needs a named vector with the rankings and the genes as names
  ## We will use  signed p value as ranking: sign(deg_results$log2fc)*(-log10(deg_results$pval)))
  rankings <- sign(deg_results$log2fc)*(-log10(deg_results$pval)) # we will use the signed p values from spatial DGE as ranking
  names(rankings) <- deg_results$gene_symbol # genes as names
  # Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
  max_ranking <- max(rankings[is.finite(rankings)])
  min_ranking <- min(rankings[is.finite(rankings)])
  rankings <- replace(rankings, rankings > max_ranking, max_ranking * 5)
  rankings <- replace(rankings, rankings < min_ranking, min_ranking * 5)
  rankings <- sort(rankings, decreasing = TRUE)
  
  # Easy peasy! Run fgsea with the pathways
  GSEAres <- fgsea(pathways = pwl1, # List of gene sets to check
                   stats = rankings,
                   scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                   minSize = 4,
                   maxSize = 500)
  saveRDS(GSEAres, file = paste0('GSEA/', paste0(filename, background_genes, '_gsea_results.RDS')))
  data.table::fwrite(GSEAres, file = paste0('GSEA/', paste0(filename, background_genes, '_gsea_results.tsv')), sep = "\t", sep2 = c("", " ", ""))
  
  # Check top pathways
  topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
  topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_top30pathways.pdf')), width = 20, height = 15)
  plotGseaTable(pwl1[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
  dev.off()
  
  # Select only independent pathways, removing redundancies/similar pathways
  collapsedPathways <- collapsePathways(GSEAres[order(pval)][padj < 0.01], pwl1, rankings)
  mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
  plotGseaTable(pwl1[mainPathways], rankings, GSEAres, gseaParam = 0.5)
  dev.off()
  
  return(GSEAres)
}

## Function:  GSEA heatmaps 
# This function takes the results from gsea pipeline  returns a heatmap
gsea_get_heatmap <- function(list_of_gsea_results, file_name = 'filename', pathway_cutoff = number_of_pathways_cutoff){
  # for debug
  # list_of_gsea_results <- all_gsea_results
  # pathway_cutoff <- 50
  # i <- 1

  if(!is.character(file_name)){
    stop('file_name must be a character')
  }
  
  list_of_pathways <- c()
  # First let's clean the gsea results
  for(i in 1:length(list_of_gsea_results)){
    list_of_gsea_results[[i]]$group <- names(list_of_gsea_results)[[i]]
    # list_of_gsea_results[[i]]$signedpval <- sign(list_of_gsea_results[[i]]$ES)*(-log10(list_of_gsea_results[[i]]$pval)) 
    # list_of_gsea_results[[i]] <- list_of_gsea_results[[i]][,c("pathway", 'pval', 'signedpval', 'group')]
    
    list_of_gsea_results[[i]]$signedpadj <- sign(list_of_gsea_results[[i]]$ES)*(-log10(list_of_gsea_results[[i]]$padj)) 
    list_of_gsea_results[[i]] <- list_of_gsea_results[[i]][,c("pathway", 'padj', 'signedpadj', 'group')]
  }
  
  df <- as.data.frame(do.call(rbind, list_of_gsea_results))
  
  # We will also select the top pathways for each group
  for(i in 1:length(list_of_gsea_results)){
    # pathways_temp <- head(unique(df[order(-abs(list_of_gsea_results[[i]]$signedpval)),'pathway']), pathway_cutoff) # top xx significant by pval for each group
    # list_of_pathways <- c(list_of_pathways, pathways_temp)
    
    pathways_temp <- head(unique(df[order(-abs(list_of_gsea_results[[i]]$signedpadj)),'pathway']), pathway_cutoff) # top xx significant by pval for each group
    list_of_pathways <- c(list_of_pathways, pathways_temp)
  }
  list_of_pathways <- unique(list_of_pathways) # remove the redundant pathways
  
  # Get heatmap
  #heatmap_df <- tidyr::pivot_wider(df[,c('signedpval', 'group', 'pathway')], values_from = signedpval, names_from = group, values_fill = 0)
  heatmap_df <- tidyr::pivot_wider(df[,c('signedpadj', 'group', 'pathway')], values_from = signedpadj, names_from = group, values_fill = 0)
  heatmap_df <- as.data.frame(heatmap_df)
  
  rownames(heatmap_df) <- heatmap_df$pathway # to be able to set them as row names
  heatmap_df <- heatmap_df[!grepl('pathway', colnames(heatmap_df))] # now remove the pathways column
  head(heatmap_df)
  
  # # We will just present the top 50 (pathway_cutoff), sorted by the absolute value of the signed pval) pathways
  if(nrow(heatmap_df) > pathway_cutoff){
    #selected_pathways <- head(unique(df[order(-abs(df$signedpval)),'pathway']), pathway_cutoff) # order by minuslogpadj and get ones with larger minuslogpadj
    selected_pathways <- list_of_pathways
    heatmap_df <- subset(heatmap_df, rownames(heatmap_df) %in% selected_pathways) # subset to less pathways for visualisation purposes if there are too many
  }
  
  # Make the rownames for the heatmap a bit nicer
  rownames(heatmap_df) <- gsub('t cell', 'T cell',
                               gsub('b cell', 'B cell',
                                    gsub('built from .*', ' (...)',
                                         gsub('mhc', 'MHC',
                                              gsub('mhc class i', 'MHC I', 
                                                   gsub('mhc class ii', 'MHC II', 
                                                        stringr::str_to_sentence(gsub('_', ' ', rownames(heatmap_df)))))))))
  heatmap <- pheatmap(heatmap_df,
                      show_colnames = T, show_rownames = T,
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      #annotation_col = col_annotations_df,
                      fontsize_col = 10,
                      fontsize_row = 10#,
                      #annotation_colors = ann_colors
  )
  
  height_of_heatmap <- nrow(heatmap_df) * 0.15 # to make the heatmap longer in case there are many rows
  
  pdf(file = paste0('GSEA/', file_name, '_heatmap.pdf'), width = 15, height = height_of_heatmap)
  grid::grid.draw(heatmap$gtable)
  dev.off()
  
  return(heatmap)
}

# Data import  ------------------------------------------------

# The frequency table has genes x sample_SPECIMENs and 1 or 0 depending if the gene is present or not
freqtable_all <- readRDS(file = paste(path, 'freqtable_all.csv', sep = '/'))
col_annotations_df <- read.csv('col_annotations_df.csv', row.names = 1)
rownames(col_annotations_df) <- col_annotations_df$sample_SPECIMEN_pattern
col_annotations_df$KMeans5 <- factor(col_annotations_df$KMeans5)

# Select here the columns you want to plot as annotations in the heatmap
heatmap_ann <- col_annotations_df %>% dplyr::select(Disease, matches('signedFDR|KMeans5'), -matches('leiden_epidermis|leiden_granuloma'))
# Rename the columns so they appear nicely in the heatmap
heatmap_ann <- heatmap_ann %>% dplyr::rename(`Core granuloma (Leiden)` = signedFDR.leiden_core_granuloma,
                                             `Border granuloma (Leiden)` = signedFDR.leiden_border_granuloma,
                                             #`Epidermis (Leiden)` = signedFDR.leiden_epidermis,
                                             #`Granuloma (Leiden)` = signedFDR.leiden_granuloma,
                                             `Granuloma (annotation)` = signedFDR.manual_granuloma, 
                                             `Epidermis (annotation)` = signedFDR.epidermis_interface, 
                                             `K-means clustering` = KMeans5)
# Set the order of the annotations in the heatmap
col_order <- c("Core granuloma (Leiden)", "Border granuloma (Leiden)", "Granuloma (annotation)", 
               "Epidermis (annotation)", "Disease", "K-means clustering")
heatmap_ann <- heatmap_ann[, col_order]

# Load in genes for each pattern & sample_SPECIMEN 
patterns_genes_df <- read.csv(paste(path, 'all_df.csv', sep = '/'), row.names = 1)
patterns_genes_df <- tidyr::separate(patterns_genes_df, slide_pattern, into = c('sample_SPECIMEN', 'pattern'), sep = '_pattern', remove = FALSE)
head(patterns_genes_df)

# Split the dataframe into a list of sub-dataframes per sample_SPECIMEN and pattern
patterns_genes_df_list <- split(patterns_genes_df, list(patterns_genes_df$sample_SPECIMEN, patterns_genes_df$pattern))
names(patterns_genes_df_list) <- gsub('-', '.', gsub('\\.', '_pattern', names(patterns_genes_df_list))) # rename so sample_SPECIMEN_pattern names match R version

# Set colours  ------------------------------------------------
col.set <- list(spots = c("DERMIS" = '#E0EEE0', 
                          "INTERFACE" = 'deepskyblue',
                          "VESSEL" = 'darkgreen',
                          "EPIDERMIS" = 'blue',
                          "HAIR FOLLICLE" = "#543005",
                          "MUSCLE" = 'darkcyan',
                          "SEBACEOUS GLAND" = 'mistyrose',
                          "SWEAT GLAND" = 'yellow',
                          "GNL" = 'orchid',
                          "GSS" = 'blueviolet',
                          "GSC" = 'mediumvioletred',
                          "GA (1)" = '#F46D43', 
                          "GA (2)" = '#D53E4F', 
                          "GA (3)" = '#9E0142',
                          "UNDETERMINED" = 'black'),
                spots_general = c("DERMIS" = '#E0EEE0', 
                                  "INTERFACE" = 'deepskyblue',
                                  "VESSEL" = 'darkgreen',
                                  "EPIDERMIS" = 'blue',
                                  "HAIR FOLLICLE" = "#543005",
                                  "MUSCLE" = 'darkcyan',
                                  "SEBACEOUS GLAND" = 'mistyrose',
                                  "SWEAT GLAND" = 'yellow',
                                  "GA" = 'firebrick',  
                                  "GNL" = 'orchid',
                                  "GSS" = 'blueviolet',
                                  "GSC" = 'mediumvioletred',
                                  "UNDETERMINED" = 'black'),
                dermis = c("UNDETERMINED" = 'black',
                           "upper EPIDERMIS" = 'blue',
                           "middle EPIDERMIS" = 'dodgerblue',
                           "basal EPIDERMIS" = 'skyblue',
                           "DERdepth1" = '#006837',
                           "DERdepth2" = '#238443',
                           "DERdepth3" = '#41AB5D',
                           "DERdepth4" = '#78C679',
                           "DERdepth5" = '#ADDD8E',
                           "DERdepth6" = '#D9F0A3',
                           "DERdepth7" = '#F7FCB9'),
                patient = c("91253" = "#F46D43",#'red',#granuloma annulare
                            "45703" = "#D53E4F",#'indianred',#granuloma annulare
                            "50107" = "#9E0142",#'firebrick',#granuloma annulare
                            "95096" = "orchid", # necrobiosis lipoidica
                            "82301" = "mediumvioletred", #sarcoidosis-cutaneous
                            "72859" = "blueviolet"), #sarcoidosis-suspected
                DISEASE = c("granuloma annulare" = "#9E0142",#'firebrick'
                            "necrobiosis lipoidica" = "orchid", # necrobiosis lipoidica
                            "sarcoidosis" = "blueviolet"), #sarcoidosis
                Disease = c("Granuloma annulare" = "#9E0142",#'firebrick'
                            "Necrobiosis lipoidica" = "orchid", # necrobiosis lipoidica
                            "Sarcoidosis" = "blueviolet"), #sarcoidosis
                lesions = c("LESIONAL" = 'tomato',
                            "NON LESIONAL" = 'darkolivegreen',
                            "UNDETERMINED" = 'lightgrey'),
                samples = c("P17851_1001" = '#F46D43',
                            "P17851_1002" = '#FDAE61',
                            "P17851_1003" = '#D53E4F',
                            "P17851_1004" = 'coral3',
                            "P18554_1001" = '#9E0142',
                            "P18554_1002" = 'firebrick3',
                            "P18554_1003" = 'orchid',
                            "P18554_1004" = 'magenta',
                            "P18554_1005" = 'blueviolet',
                            "P18554_1006" = 'mediumpurple',
                            "P18554_1007" = 'mediumvioletred',
                            "P18554_1008" = 'violet'),
                leiden_r1.3 = c("0"  = 'darkolivegreen', #mainly contains non-lesional spots
                                "1"  = '#D9F0A3', # dermdepth6
                                "2"  = '#238443', # derdepth2
                                "3"  = 'firebrick', # granuloma
                                "4"  = '#78C679', # derdepth4
                                "5"  = '#78C679', # derdepth4
                                "6"  = "#41AB5D", #dermis depth 3
                                "7"  = "#006837", #derm depth1
                                "8"  = '#ADDD8E', # dermdepth5
                                "9"  = "#238443", #derm depth2
                                "10" = '#78C679', # derdepth4
                                "11" = 'blue', #epidermis
                                "12" = 'orchid', #purple for very low % granuloma --> this is mostly ga
                                "13" = '#F46D43', #orange for very low % granuloma --> this is mostly gnl
                                "14" = 'dodgerblue', #mostly interface, has mid epidermis colour
                                "15" = "deepskyblue", #mostly interface, has interface colour
                                "16" = '#cfafaf', #sebaceous gland. 
                                "17" = 'yellow', # sweat gland
                                "18" = 'darkcyan', # muscle
                                "19" = '#006837'), #derdepth1
                leiden_r0.3_g_bc = c('0' = '#F46D43',
                                     '1' = '#78C679',
                                     '2' = '#9E0142',
                                     '3' = '#006837'),
                leiden_r0.5_g_bc = c('0' = '#006837',
                                     '1' = '#D9F0A3',
                                     '2' = '#F46D43',
                                     '3' = '#D53E4F',
                                     '4' = '#9E0142',
                                     '5' = 'mediumvioletred'),
                `K-means clustering` = c("1"="#88CFA4", "2"="#F88D52", "3"="skyblue",
                                         "4"="#9E0142", "5"="#5E4FA2"))
BuRd <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
GrBu <- colorRampPalette(colors = c(rev(brewer.pal(11, "PiYG"))[1:4],
                                    brewer.pal(11, "RdYlBu")[c(6, 8:11)]))

# Importing adata ####
# Set to T if you want to save the clustering .h5 files created in the jupyter notebook as RDS. 
# This code will convert spot type, skin layers, patients, leiden clusters..., to factors and order them
if(F){
  ad <- read_h5ad("/Users/mendenlab/work/spatial_granuloma/results/current/final/Granuloma_spots_only.h5")
  saveRDS(ad, file = "/Users/mendenlab/work/spatial_granuloma/results/current/final/Granuloma_spots_only.RDS")
  names(ad$obsm)
  names(ad$obs)
  
  df <- cbind(ad$obs,
              X_g_uc_umap1 = ad$obsm$X_umap_granuloma_uc[, 1], X_g_uc_umap2 = ad$obsm$X_umap_granuloma_uc[, 2], # no batch correction! save the UMAPs computed for the clusters 3, 12 and 13 (granuloma)
              X_g_bc_umap1 = ad$obsm$X_umap_granuloma_bc[, 1], X_g_bc_umap2 = ad$obsm$X_umap_granuloma_bc[, 2], # save the batch-corrected UMAPs computed for the clusters 3, 12 and 13 (granuloma)
              spatial1 = ad$obsm$spatial[, 1], spatial2 = ad$obsm$spatial[, 2])
  
  df$spot_type <- factor(as.character(df$spot_type),
                         levels = c("DERMIS", "INTERFACE",
                                    "EPIDERMIS",
                                    "GA", "GA (1)", "GA (2)", "GA (3)",
                                    "GSS",
                                    "GSC",
                                    "GNL",
                                    "VESSEL", "MUSCLE", "HAIR FOLLICLE", 
                                    "SEBACEOUS GLAND", "SWEAT GLAND", "UNDETERMINED"),
                         ordered = T)
  df$skin_layer <- factor(as.character(df$skin_layer), levels = names(col.set$dermis), ordered = T)
  
  df[df$sample=="P17851_1001", c("spatial1", "spatial2")] <- 30000-df[df$sample == "P17851_1001", c("spatial1", "spatial2")]
  df[df$sample=="P18554_1002", c("spatial1", "spatial2")] <- 30000-df[df$sample == "P18554_1002", c("spatial1", "spatial2")]
  df[df$sample=="P18554_1006", c("spatial1", "spatial2")] <- 30000-df[df$sample == "P18554_1006", c("spatial1", "spatial2")]
  
  df$patient <- factor(df$patient, levels = c("91253", "45703", "50107", "95096", "82301", "72859"), ordered = T)
  
  df$spot_type[df$patient == "91253" & df$spot_type == "GA"] <- "GA (1)"
  df$spot_type[df$patient == "45703" & df$spot_type == "GA"] <- "GA (2)"
  df$spot_type[df$patient == "50107" & df$spot_type == "GA"] <- "GA (3)"
  
  df <- df[order(df$spot_type), ]
  if(!is.null(df)){
    saveRDS(df, file = "/Users/mendenlab/work/spatial_granuloma/results/current/final/Only_granulomas_df.RDS")
  }
} else {
  df <- readRDS("/Users/mendenlab/work/spatial_granuloma/results/current/final/Only_granulomas_df.RDS")
}

# Heatmap of only granuloma patterns from clusters 2 and 4------------------------------------------------
# clusters_2_4_patterns <- col_annotations_df[col_annotations_df$KMeans5 == 2|col_annotations_df$KMeans5 == 4, 'sample_SPECIMEN_pattern']
# freqtable_g <- freqtable_all[, clusters_2_4_patterns] # only keep patterns belonging to K means cluster 1 and 2
# freqtable_g <- freqtable_g[rowSums(freqtable_g[])>0,] # remove rows with only 0s
# heatmap_ann_g <- heatmap_ann[clusters_2_4_patterns,]
# heatmap_g <- pheatmap(freqtable_g,
#                       show_colnames = F, show_rownames = F,
#                       clustering_method = 'ward.D',
#                       annotation_col = heatmap_ann_g,
#                       treeheight_row = 0,
#                       treeheight_col = 10,
#                       annotation_colors = col.set
# )
# pdf(width = 8, height = 10, file = paste0(path, "/Heatmaps/onlygranulomapatterns_heatmap.pdf"))
# heatmap_g
# dev.off()

# Clustering granulomas ------------------------------------------------
r_list <- c('0.3', '0.5', '0.8', '1.0', '1.3', '1.5', '3', '5')

## No batch correction ##### 
# (using uncorrected data - uc)
# Barplots of spot type, donor and disease count per cluster
list_of_barplots <- list()
for(leiden_r in r_list){
  list_of_barplots[[(length(list_of_barplots) + 1)]] <- ggarrange(
    ggarrange(
      ggplot(data=df,
             aes(x= .data[[paste0('leiden_r', leiden_r, '_g_uc')]], fill = spot_type)) + # this is like saying x = df$leiden_r0.5
        geom_bar(stat="count", position = "fill") +
        scale_fill_manual("", values = col.set$spots) +
        theme_bw(),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    ggarrange(
      ggplot(data=df,
             aes(x=.data[[paste0('leiden_r', leiden_r, '_g_uc')]], fill = patient)) +
        geom_bar(stat="count", position = "stack") +
        scale_fill_manual("", values = col.set$patient) +
        theme_bw(),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    ggarrange(
      ggplot(data=df,
             aes(x=.data[[paste0('leiden_r', leiden_r, '_g_uc')]], fill = DISEASE)) +
        geom_bar(stat="count", position = "stack") +
        scale_fill_manual("", values = col.set$DISEASE) +
        theme_bw(),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    ggarrange(
      ggplot(data=df,
             aes(x=.data[[paste0('leiden_r', leiden_r, '_g_uc')]], fill = sample_SPECIMEN)) +
        geom_bar(stat="count", position = "stack") +
        theme_bw() + 
        theme(legend.position = 'none') +
        labs(x = 'coloured by sample_SPECIMEN'),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    nrow = 1, ncol = 4, widths = c(1,1,1,1), heights = c(1,1,1,1),
    common.legend = F)
}

pdf(file = paste(path, 'Clustering/Barplots_leiden_uc.pdf', sep = ''), width = 18, height = 30)
ggarrange(list_of_barplots[[1]], list_of_barplots[[2]], list_of_barplots[[3]], list_of_barplots[[4]], 
          list_of_barplots[[5]], list_of_barplots[[6]], list_of_barplots[[7]], list_of_barplots[[8]],
          nrow = 8, ncol = 1, common.legend = F)
dev.off()

# UMAPS
UMAP_spottype_uc <- ggplot(df, aes(x = X_g_uc_umap1, y = X_g_uc_umap2, col = spot_type)) +
  geom_point(size = 0.5) +
  scale_color_manual("", values = col.set$spots) +
  labs(title = '', x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4))) 

list_of_UMAPs <- list()
for(leiden_r in r_list){
  #print(leiden_r)
  list_of_UMAPs[[(length(list_of_UMAPs) + 1)]] <- ggplot(df, aes(x = X_g_uc_umap1, y = X_g_uc_umap2, col = .data[[paste0('leiden_r', leiden_r, '_g_uc')]])) +
    geom_point(size = 0.5) +
    labs(title = paste0("Leiden clustering granuloma spots (r = ", leiden_r,')'),
         x = "UMAP 1",
         y = "UMAP 2") +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    #guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme(legend.position = 'none')
}

pdf(file = paste(path, 'Clustering/UMAPs_Leiden_uc.pdf', sep = ''), width = 18, height = 10)
ggarrange(list_of_UMAPs[[1]], list_of_UMAPs[[2]], list_of_UMAPs[[3]], list_of_UMAPs[[4]], 
          list_of_UMAPs[[5]], list_of_UMAPs[[6]], list_of_UMAPs[[7]], UMAP_spottype,
          nrow = 2, ncol = 4, common.legend = F)
dev.off()

### Leiden spatial
spatial_spottype <- ggplot(data=df,
                           aes(x=spatial1, y=spatial2, col=spot_type)) +
  geom_point(size=0.2) +
  scale_color_manual('', values = col.set$spots) +
  scale_y_reverse() +
  labs(title = 'Spatial representation granuloma spots') +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
        #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  facet_wrap(.~sample,
             scales = "fixed",
             nrow = 2)

list_of_spatialplots <- list()
for(leiden_r in r_list){
  list_of_spatialplots[[(length(list_of_spatialplots) + 1)]] <- ggplot(data=df,
                                                                       aes(x=spatial1, y=spatial2, col=.data[[paste0('leiden_r', leiden_r, '_g_uc')]])) +
    geom_point(size=0.2) +
    #scale_color_manual('', values = col.set$leiden_r1.3) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", leiden_r, ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
          #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample,
               scales = "fixed",
               nrow = 2)
}

pdf(file = paste(path, 'Clustering/Spatial_leiden_uc.pdf', sep = ''), width = 10, height = 30)
ggarrange(spatial_spottype, list_of_spatialplots[[1]], list_of_spatialplots[[2]], list_of_spatialplots[[3]], 
          list_of_spatialplots[[4]], list_of_spatialplots[[5]], list_of_spatialplots[[6]], list_of_spatialplots[[7]], 
          nrow = 8, ncol = 1, common.legend = F)
dev.off()

## Batch corrected data #####
## Leiden clusters 

# Barplots of spot type, donor and disease count per cluster
list_of_barplots <- list()
for(leiden_r in r_list){
  list_of_barplots[[(length(list_of_barplots) + 1)]] <- ggarrange(
    ggarrange(
      ggplot(data=df,
             aes(x= .data[[paste0('leiden_r', leiden_r, '_g_bc')]], fill = spot_type)) + # this is like saying x = df$leiden_r0.5
        geom_bar(stat="count", position = "fill") +
        scale_fill_manual("", values = col.set$spots) +
        theme_bw(),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    ggarrange(
      ggplot(data=df,
             aes(x=.data[[paste0('leiden_r', leiden_r, '_g_bc')]], fill = patient)) +
        geom_bar(stat="count", position = "stack") +
        scale_fill_manual("", values = col.set$patient) +
        theme_bw(),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    ggarrange(
      ggplot(data=df,
             aes(x=.data[[paste0('leiden_r', leiden_r, '_g_bc')]], fill = DISEASE)) +
        geom_bar(stat="count", position = "stack") +
        scale_fill_manual("", values = col.set$DISEASE) +
        theme_bw(),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    ggarrange(
      ggplot(data=df,
             aes(x=.data[[paste0('leiden_r', leiden_r, '_g_bc')]], fill = sample_SPECIMEN)) +
        geom_bar(stat="count", position = "stack") +
        theme_bw() + 
        theme(legend.position = 'none') +
        labs(x = 'coloured by sample_SPECIMEN'),
      NULL,
      nrow = 1, ncol = 2,
      widths = c(1, 0.03)),
    nrow = 1, ncol = 4, widths = c(1,1,1,1), heights = c(1,1,1,1),
    common.legend = F)
}

pdf(file = paste(path, 'Clustering/Barplots_leiden_bc.pdf', sep = ''), width = 18, height = 30)
ggarrange(list_of_barplots[[1]], list_of_barplots[[2]], list_of_barplots[[3]], list_of_barplots[[4]], 
          list_of_barplots[[5]], list_of_barplots[[6]], list_of_barplots[[7]], list_of_barplots[[8]],
          nrow = 8, ncol = 1, common.legend = F)
dev.off()

# UMAPS
UMAP_spottype <- ggplot(df, aes(x = X_g_bc_umap1, y = X_g_bc_umap2, col = spot_type)) +
  geom_point(size = 0.5) +
  #scale_x_reverse() +
  #scale_y_reverse() +
  scale_color_manual("", values = col.set$spots) +
  labs(title = paste0("Leiden clustering granuloma spots (spot type)"), x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4)))


list_of_UMAPs <- list()
for(leiden_r in r_list){
  #print(leiden_r)
  umap_plot <- ggplot(df, aes(x = X_g_bc_umap1, y = X_g_bc_umap2, col = .data[[paste0('leiden_r', leiden_r, '_g_bc')]])) +
    geom_point(size = 0.5) +
    labs(title = paste0("Leiden clustering granuloma spots (r = ", leiden_r,')'),
         x = "UMAP 1",
         y = "UMAP 2") +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    #guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme(legend.position = 'none')
  if(leiden_r == '0.3'){
    umap_plot <- umap_plot + scale_color_manual("", values = col.set$leiden_r0.3_g_bc)
  } else if (leiden_r == '0.5'){
    umap_plot <- umap_plot + scale_color_manual("", values = col.set$leiden_r0.5_g_bc)
  }
  list_of_UMAPs[[(length(list_of_UMAPs) + 1)]] <- umap_plot
}

pdf(file = paste(path, 'Clustering/UMAPs_Leiden_bc.pdf', sep = ''), width = 18, height = 10)
ggarrange(list_of_UMAPs[[1]], list_of_UMAPs[[2]], list_of_UMAPs[[3]], list_of_UMAPs[[4]], 
          list_of_UMAPs[[5]], list_of_UMAPs[[6]], list_of_UMAPs[[7]], UMAP_spottype,
          nrow = 2, ncol = 4, common.legend = F)
dev.off()

### Leiden spatial
spatial_spottype <- ggplot(data=df,
                           aes(x=spatial1, y=spatial2, col=spot_type)) +
  geom_point(size=0.1) +
  scale_color_manual('', values = col.set$spots) +
  scale_y_reverse() +
  labs(title = 'Spatial representation granuloma spots (coloured by spot type)') +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
        #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  facet_wrap(.~sample,
             scales = "fixed",
             nrow = 2)

spatial_leiden13 <- ggplot(data=df,
                           aes(x=spatial1, y=spatial2, col=leiden_r1.3_patient)) +
  geom_point(size=0.1) +
  scale_color_manual('', values = col.set$leiden_r1.3) +
  scale_y_reverse() +
  labs(title = 'Spatial representation granuloma spots (coloured by original Leiden clusters r1.3)') +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
        #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  facet_wrap(.~sample,
             scales = "fixed",
             nrow = 2)

list_of_spatialplots <- list()
for(leiden_r in r_list){
  spatial_plots <- ggplot(data = df, aes(x = spatial1, y = spatial2, col = .data[[paste0('leiden_r', leiden_r, '_g_bc')]])) +
    geom_point(size=0.1) +
    #scale_color_manual('', values = col.set$leiden_r1.3) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", leiden_r, ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
          #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample,
               scales = "fixed",
               nrow = 2)
  if(leiden_r == '0.3'){
    spatial_plots <- spatial_plots + scale_color_manual("", values = col.set$leiden_r0.3_g_bc)
  } else if (leiden_r == '0.5'){
    spatial_plots <- spatial_plots + scale_color_manual("", values = col.set$leiden_r0.5_g_bc)
  }
  list_of_spatialplots[[(length(list_of_spatialplots) + 1)]] <- spatial_plots
}

pdf(file = paste(path, 'Clustering/Spatial_leiden_bc.pdf', sep = ''), width = 10, height = 33)
ggarrange(spatial_spottype, spatial_leiden13, list_of_spatialplots[[1]], list_of_spatialplots[[2]], list_of_spatialplots[[3]], 
          list_of_spatialplots[[4]], list_of_spatialplots[[5]], list_of_spatialplots[[6]], list_of_spatialplots[[7]], 
          nrow = 9, ncol = 1, common.legend = F)
dev.off()


## UC vs BC figure ####
UMAP_leiden03_uc <- ggplot(df, aes(x = X_g_uc_umap1, y = X_g_uc_umap2, col = .data[[paste0('leiden_r', '0.3', '_g_uc')]])) +
  geom_point(size = 0.5) +
  labs(title = '', x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = 'right', legend.title = element_blank()) 

UMAP_leiden03_bc <- ggplot(df, aes(x = X_g_bc_umap1, y = X_g_bc_umap2, col = .data[['leiden_r0.3_g_bc']])) +
  geom_point(size = 0.5) +
  labs(title = '', x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = 'right', legend.title = element_blank()) + 
  scale_color_manual("", values = col.set$leiden_r0.3_g_bc)

barplots_uc_bc <- ggarrange(
  ggplot(data=df, aes(x=.data[[paste0('leiden_r', '0.3', '_g_uc')]], fill = patient)) +
    geom_bar(stat="count", position = "stack") +
    scale_fill_manual("", values = col.set$patient) +
    labs(title = 'Patient', x = "Leiden cluster (r = 0.3)", y = 'Count') +
    theme_bw() +
    theme(#axis.text.y = element_text(size = 12), 
      plot.title = element_text(vjust = - 10, hjust = 0.9),
          axis.line.y = element_blank(), legend.position = 'none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  ggplot(data=df, aes(x=.data[[paste0('leiden_r', '0.3', '_g_uc')]], fill = DISEASE)) +
    geom_bar(stat="count", position = "stack") +
    scale_fill_manual("", values = col.set$DISEASE) +
    labs(title = 'Disease', x = "Leiden cluster (r = 0.3)", y = 'Count') +
    theme_bw() +
    theme(#axis.text.y = element_text(size = 13),
      plot.title = element_text(vjust = - 10, hjust = 0.9),
          axis.line.y = element_blank(), legend.position = 'none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  as_ggplot(get_legend(ggplot(data=df, aes(x=.data[[paste0('leiden_r', '0.3', '_g_uc')]], fill = patient)) +
                         geom_bar(stat="count", position = "stack") + 
                         scale_fill_manual("Patient", values = col.set$patient) +
                         theme(legend.position = 'right'))),
  ggplot(data=df, aes(x=.data[[paste0('leiden_r', '0.3', '_g_bc')]], fill = patient)) +
    geom_bar(stat="count", position = "stack") +
    scale_fill_manual("", values = col.set$patient) +
    labs(title = 'Patient', x = "Leiden cluster (r = 0.3)", y = 'Count') +
    theme_bw() +
    theme(#axis.text.y = element_text(size = 13),
      plot.title = element_text(vjust = - 10, hjust = 0.9),
          axis.line.y = element_blank(), legend.position = 'none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  ggplot(data=df, aes(x=.data[[paste0('leiden_r', '0.3', '_g_bc')]], fill = DISEASE)) +
    geom_bar(stat="count", position = "stack") +
    scale_fill_manual("", values = col.set$DISEASE) +
    labs(title = 'Disease', x = "Leiden cluster (r = 0.3)", y = 'Count') +
    theme_bw() +
    theme(#axis.text.y = element_text(size = 13),
      plot.title = element_text(vjust = - 10, hjust = 0.9),
          axis.line.y = element_blank(), legend.position = 'none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  as_ggplot(get_legend(ggplot(data=df, aes(x=.data[[paste0('leiden_r', '0.3', '_g_uc')]], fill = DISEASE)) +
                         geom_bar(stat="count", position = "stack") + 
                         scale_fill_manual("Disease", values = col.set$DISEASE, labels = c('GA', 'SC', 'NL')) +
                         theme(legend.position = 'right'))),
  nrow = 2, ncol = 3, widths = c(1, 1, 0.5))

UMAPs_barplots_uc_bc <- ggarrange(
  ggarrange(
    as_ggplot(get_legend(UMAP_leiden03_uc)),
    UMAP_leiden03_uc + theme(legend.position = 'none', legend.title = element_blank()) + ggtitle(''),
    UMAP_spottype_uc + theme(legend.position = 'none') + ggtitle(''),  
    as_ggplot(get_legend(UMAP_leiden03_bc)),
    UMAP_leiden03_bc + theme(legend.position = 'none') + ggtitle(''),
    UMAP_spottype + theme(legend.position = 'none') + ggtitle(''), 
    ncol = 3, nrow = 2, widths = c(0.3, 1, 1)),
  barplots_uc_bc, widths = c(5,4))

pdf(file = paste(path, 'Clustering/UMAPs_and_barplots_uc_bc.pdf', sep = ''), width = 12, height = 8)
UMAPs_barplots_uc_bc
dev.off()

# For now let's go with r = 0.3 since it seems to distinguish border/core spots without overclustering

# DGE Analysis - clusters vs healthy dermis #####

## Data import ####
reexport_dge_df = F # only set as T if you want to re-export everything, else only read in the data
if(reexport_dge_df == T){
  df <- readRDS("/Users/mendenlab/work/spatial_granuloma/results/current/final/Only_granulomas_df.RDS") # this contains the clusters for only granulomas
  df_all <- readRDS("/Users/mendenlab/work/spatial_granuloma/results/current/final/Granuloma_DGE_df.RDS") # this contains all spots
  df$sample <- as.double(df$sample)
  merged_df <- left_join(df_all, 
                         df[,c('sample', 'array_row', 'array_col', 'leiden_r0.3_g_bc', 'leiden_r0.5_g_bc', 'X_g_bc_umap1', 'X_g_bc_umap2')], 
                         by = c('sample', 'array_row', 'array_col')) # common columns on which to merge 
  
  merged_df$cluster0vshealthydermis <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 0 & !is.na(merged_df$leiden_r0.3_g_bc == 0)), 'CLUSTER0',
                                              ifelse((merged_df$dermis_nonlesional_noborder == 1), 'HEALTHYDERMIS', 'other')))
  merged_df$cluster1vshealthydermis <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 1 & !is.na(merged_df$leiden_r0.3_g_bc == 1)), 'CLUSTER2',
                                              ifelse(merged_df$dermis_nonlesional_noborder == 1, 'HEALTHYDERMIS', 'other')))
  merged_df$cluster2vshealthydermis <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 2 & !is.na(merged_df$leiden_r0.3_g_bc == 2)), 'CLUSTER2',
                                              ifelse(merged_df$dermis_nonlesional_noborder == 1, 'HEALTHYDERMIS', 'other')))
  merged_df$cluster3vshealthydermis <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 3 & !is.na(merged_df$leiden_r0.3_g_bc == 3)), 'CLUSTER3',
                                              ifelse(merged_df$dermis_nonlesional_noborder == 1, 'HEALTHYDERMIS', 'other')))
  
  merged_df$cluster0vsrest <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 0 & !is.na(merged_df$leiden_r0.3_g_bc == 0)), 'CLUSTER0',
                                            ifelse((merged_df$leiden_r0.3_g_bc != 0 & !is.na(merged_df$leiden_r0.3_g_bc)), 'OTHERCLUSTERS',
                                                   'other')))
  merged_df$cluster1vsrest <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 1 & !is.na(merged_df$leiden_r0.3_g_bc == 1)), 'CLUSTER2',
                                            ifelse((merged_df$leiden_r0.3_g_bc != 1 & !is.na(merged_df$leiden_r0.3_g_bc)), 'OTHERCLUSTERS',
                                                   'other')))
  merged_df$cluster2vsrest <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 2 & !is.na(merged_df$leiden_r0.3_g_bc == 2)), 'CLUSTER2',
                                            ifelse((merged_df$leiden_r0.3_g_bc != 2 & !is.na(merged_df$leiden_r0.3_g_bc)), 'OTHERCLUSTERS',
                                                   'other')))
  merged_df$cluster3vsrest <- factor(ifelse((merged_df$leiden_r0.3_g_bc == 3 & !is.na(merged_df$leiden_r0.3_g_bc == 3)), 'CLUSTER3',
                                            ifelse((merged_df$leiden_r0.3_g_bc != 3 & !is.na(merged_df$leiden_r0.3_g_bc)), 'OTHERCLUSTERS',
                                                   'other')))
  
  # for dge we had to remove many indexes, so to flip the spatial coordinates of samples P17851_1001, P18554_1002 and P18554_1006 we will use the new values (see python new_adata_for_dge script for more)
  merged_df[merged_df$sample=="1", c("spatial1", "spatial2")] <- 30000-merged_df[merged_df$sample == "1", c("spatial1", "spatial2")]
  merged_df[merged_df$sample=="6", c("spatial1", "spatial2")] <- 30000-merged_df[merged_df$sample == "6", c("spatial1", "spatial2")]
  merged_df[merged_df$sample=="10", c("spatial1", "spatial2")] <- 30000-merged_df[merged_df$sample == "10", c("spatial1", "spatial2")]
  
  # We only need to save this once in the deg folder and the final results folder
  if(!is.null(merged_df)){
    saveRDS(merged_df,
            file = "/Users/mendenlab/work/spatial_granuloma/results/current/final/Granuloma_DGE_df_withonlygranuloma.RDS")
    saveRDS(merged_df,
            file = "/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_df_withonlygranuloma.RDS")
  }
  
  sce <- readRDS("/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_sce.RDS")
  sce$leiden_r0.3_g_bc <- merged_df$leiden_r0.3_g_bc
  sce$leiden_r0.5_g_bc <- merged_df$leiden_r0.5_g_bc
  sce$X_g_bc_umap1 <- merged_df$X_g_bc_umap1
  sce$X_g_bc_umap2 <- merged_df$X_g_bc_umap2
  sce$cluster0vshealthydermis <- merged_df$cluster0vshealthydermis
  sce$cluster1vshealthydermis <- merged_df$cluster1vshealthydermis
  sce$cluster2vshealthydermis <- merged_df$cluster2vshealthydermis
  sce$cluster3vshealthydermis <- merged_df$cluster3vshealthydermis
  sce$cluster0vsrest <- merged_df$cluster0vsrest
  sce$cluster1vsrest <- merged_df$cluster1vsrest
  sce$cluster2vsrest <- merged_df$cluster2vsrest
  sce$cluster3vsrest <- merged_df$cluster3vsrest
  sce$spatial1 <- merged_df$spatial1
  sce$spatial2 <- merged_df$spatial2
  
  # Check if they merged correctly - yes:)
  # sce_df <- as.data.frame(sce@colData)
  # tail(merged_df[merged_df$leiden_r0.3_g_bc == "2" & !is.na(merged_df$leiden_r0.3_g_bc),])
  # tail(sce_df[sce_df$leiden_r0.3_g_bc == "2" & !is.na(sce_df$leiden_r0.3_g_bc),])
  # ggplot(data = sce_df, aes(x = spatial1, y = spatial2, col = cluster2vsrest)) +
  #   geom_point(size=0.1) +
  #   #scale_color_manual('', values = col.set$leiden_r1.3) +
  #   scale_y_reverse() +
  #   labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
  #   theme_bw() +
  #   theme(aspect.ratio = 1,
  #         axis.title = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank(),
  #         plot.background = element_rect(fill = NA, colour=NA),
  #         #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
  #         #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
  #         axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         legend.position = "none") +
  #   facet_wrap(.~sample,
  #              scales = "fixed",
  #              nrow = 2)
  saveRDS(sce,
          file = "/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_sce_withonlygranuloma.RDS")

} else {
  df <- readRDS("/Users/mendenlab/work/spatial_granuloma/results/current/final/Granuloma_DGE_df_withonlygranuloma.RDS")
  sce <- readRDS("/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_sce_withonlygranuloma.RDS")
  normalised_counts <- readRDS("/Volumes/Drive/spatial_granuloma/output/DEG/normalised_counts.RDS")
  normalised_counts_headers <- sce@colData@listData
}

# Plot the clusters for only granuloma spots in the spatial slides

spatial_spottype <- ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = spot_type)) +
  geom_point(size=0.1) +
  #scale_color_manual('', values = col.set$leiden_r1.3) +
  scale_y_reverse() +
  labs(title = paste0('Spatial representation granuloma spots (coloured by spot type)')) +
  theme_bw() +
  scale_color_manual('Spot type   ', values = col.set$spots_general) +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
        #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(.~sample,
             scales = "fixed",
             nrow = 2)

spatial_leidenonlyg <- ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = .data[[paste0('leiden_r', '0.3', '_g_bc')]])) +
  geom_point(size=0.07) +
  #scale_color_manual('', values = col.set$leiden_r1.3) +
  scale_y_reverse() +
  labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')'), colour = 'Cluster') +
  theme_bw() +
  scale_color_manual('Leiden clusters (r = 0.3)   ', values = col.set$leiden_r0.3_g_bc, breaks = c(2, 0, 3, 1), labels = c('Inner core', 'Outer core', 'Border-1', 'Border-2')) +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        #legend.key.size = unit(3,"line"),
        #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
        #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  #guides(color = guide_legend(override.aes = list(size = 1))) +
  facet_wrap(.~sample,
             scales = "fixed",
             nrow = 2)

spatial_manual_leiden_onlyg <- ggarrange(spatial_spottype + theme(legend.position = 'none'),  
          as_ggplot(get_legend(spatial_spottype + geom_point(size = 2) + theme(legend.position = 'bottom'))),
          spatial_leidenonlyg + theme(legend.position = 'none'),  
          as_ggplot(get_legend(spatial_leidenonlyg + geom_point(size = 2) + theme(legend.position = 'bottom'))),
          ncol = 1, nrow = 4, heights = c(1, 0.3, 1, 0.3), widths = c(1)) + theme(plot.margin = margin(0,0,0,0))

pdf(file = paste0(DEG_results_path, 'spatial_spottypeandleiden_onlyg.pdf'), height = 10, width = 10)
spatial_manual_leiden_onlyg
dev.off()

## Cluster 0 vs healthy dermis ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster0vshealthydermis'
reference_group <- 'HEALTHYDERMIS'
sce_cluster0vshealthydermis <- sce[, which(sce$cluster0vshealthydermis %in% c("CLUSTER0", "HEALTHYDERMIS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster0vshealthydermis$cluster0vshealthydermis <- droplevels(sce_cluster0vshealthydermis$cluster0vshealthydermis)
#sce_cluster0vshealthydermis$cluster0vshealthydermis <- factor(sce_cluster0vshealthydermis$cluster0vshealthydermis, levels = c("healthydermis", "epidermis"), ordered = TRUE)

### 2. Exploratory analysis ####

comparison_title <- 'Outer core vs healthy dermis' #'Cluster 0 vs healthy dermis' 
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster0vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster0vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster0vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster0vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster0vshealthydermis)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER0' = 'mediumvioletred', 'HEALTHYDERMIS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)

pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster0vshealthydermis <- run_glmgampoi(sce_cluster0vshealthydermis, design_formula = ~ patient + cluster0vshealthydermis, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster0vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster0vshealthydermis_patient+cluster0vshealthydermis_degresults.csv'), row.names = 1)
fit_cluster0vshealthydermis <- readRDS(paste0(glmgampoi_path, 'sce_cluster0vshealthydermis_patient+cluster0vshealthydermis_fit.RDS'))

get_dge_summary(dge_cluster0vshealthydermis)

# Most different genes (padj)
head(dge_cluster0vshealthydermis[order(dge_cluster0vshealthydermis$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster0vshealthydermis, comparison_title, sce_cluster0vshealthydermis)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster0vshealthydermis', 
                           group1 = 'CLUSTER0', group2 = 'HEALTHYDERMIS', 
                           dge_df = dge_cluster0vshealthydermis, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()

## Cluster 1 vs healthy dermis ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster1vshealthydermis'
reference_group <- 'HEALTHYDERMIS'
sce_cluster1vshealthydermis <- sce[, which(sce$cluster1vshealthydermis %in% c("CLUSTER2", "HEALTHYDERMIS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster1vshealthydermis$cluster1vshealthydermis <- droplevels(sce_cluster1vshealthydermis$cluster1vshealthydermis)
#sce_cluster1vshealthydermis$cluster1vshealthydermis <- factor(sce_cluster1vshealthydermis$cluster1vshealthydermis, levels = c("healthydermis", "epidermis"), ordered = TRUE)

### 2. Exploratory analysis ####
comparison_title <- 'Border-2 vs healthy dermis' #'Cluster1 vs healthy dermis'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster1vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster1vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster1vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster1vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster1vshealthydermis)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER2' = 'mediumvioletred', 'HEALTHYDERMIS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)
pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster1vshealthydermis <- run_glmgampoi(sce_cluster1vshealthydermis, design_formula = ~ patient + cluster1vshealthydermis, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster1vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster1vshealthydermis_patient+cluster1vshealthydermis_degresults.csv'), row.names = 1)
fit_cluster1vshealthydermis <- readRDS(paste0(glmgampoi_path, 'sce_cluster1vshealthydermis_patient+cluster1vshealthydermis_fit.RDS'))

get_dge_summary(dge_cluster1vshealthydermis)

# Most different genes (padj)
head(dge_cluster1vshealthydermis[order(dge_cluster1vshealthydermis$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster1vshealthydermis, comparison_title, sce_cluster1vshealthydermis)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster1vshealthydermis', 
                           group1 = 'CLUSTER2', group2 = 'HEALTHYDERMIS', 
                           dge_df = dge_cluster1vshealthydermis, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()

## Cluster 2 vs healthy dermis ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster2vshealthydermis'
reference_group <- 'HEALTHYDERMIS'
sce_cluster2vshealthydermis <- sce[, which(sce$cluster2vshealthydermis %in% c("CLUSTER2", "HEALTHYDERMIS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster2vshealthydermis$cluster2vshealthydermis <- droplevels(sce_cluster2vshealthydermis$cluster2vshealthydermis)
#sce_cluster2vshealthydermis$cluster2vshealthydermis <- factor(sce_cluster2vshealthydermis$cluster2vshealthydermis, levels = c("healthydermis", "epidermis"), ordered = TRUE)

### 2. Exploratory analysis ####

comparison_title <- 'Inner core vs healthy dermis' #'Cluster 2 vs healthy dermis'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster2vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster2vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster2vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster2vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster2vshealthydermis)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER2' = 'mediumvioletred', 'HEALTHYDERMIS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)
pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster2vshealthydermis <- run_glmgampoi(sce_cluster2vshealthydermis, design_formula = ~ patient + cluster2vshealthydermis, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster2vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster2vshealthydermis_patient+cluster2vshealthydermis_degresults.csv'), row.names = 1)
fit_cluster2vshealthydermis <- readRDS(paste0(glmgampoi_path, 'sce_cluster2vshealthydermis_patient+cluster2vshealthydermis_fit.RDS'))

get_dge_summary(dge_cluster2vshealthydermis)

# Most different genes (padj)
head(dge_cluster2vshealthydermis[order(dge_cluster2vshealthydermis$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster2vshealthydermis, comparison_title, sce_cluster2vshealthydermis)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster2vshealthydermis', 
                           group1 = 'CLUSTER2', group2 = 'HEALTHYDERMIS', 
                           dge_df = dge_cluster2vshealthydermis, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()

## Cluster 3 vs healthy dermis ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster3vshealthydermis'
reference_group <- 'HEALTHYDERMIS'
sce_cluster3vshealthydermis <- sce[, which(sce$cluster3vshealthydermis %in% c("CLUSTER3", "HEALTHYDERMIS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster3vshealthydermis$cluster3vshealthydermis <- droplevels(sce_cluster3vshealthydermis$cluster3vshealthydermis)
#sce_cluster3vshealthydermis$cluster3vshealthydermis <- factor(sce_cluster3vshealthydermis$cluster3vshealthydermis, levels = c("healthydermis", "epidermis"), ordered = TRUE)

### 2. Exploratory analysis ####

comparison_title <- 'Border-1 vs healthy dermis' #'Cluster 3 vs healthy dermis'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster3vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(x = factor(patient), fill = factor(cluster3vshealthydermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(DISEASE), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(spot_type), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster3vshealthydermis@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster3vshealthydermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster3vshealthydermis)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER3' = 'mediumvioletred', 'HEALTHYDERMIS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)
pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster3vshealthydermis <- run_glmgampoi(sce_cluster3vshealthydermis, design_formula = ~ patient + cluster3vshealthydermis, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster3vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster3vshealthydermis_patient+cluster3vshealthydermis_degresults.csv'), row.names = 1)
fit_cluster3vshealthydermis <- readRDS(paste0(glmgampoi_path, 'sce_cluster3vshealthydermis_patient+cluster3vshealthydermis_fit.RDS'))

get_dge_summary(dge_cluster3vshealthydermis)

# Most different genes (padj)
head(dge_cluster3vshealthydermis[order(dge_cluster3vshealthydermis$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster3vshealthydermis, comparison_title, sce_cluster3vshealthydermis)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster3vshealthydermis', 
                           group1 = 'CLUSTER3', group2 = 'HEALTHYDERMIS', 
                           dge_df = dge_cluster3vshealthydermis, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()

# GSEA  ------------------------------------------------
setwd(pea_path)

## Cluster0 ####
name_of_subset <- 'cluster0vshealthydermis'
gsea_cluster0_go <- run_gsea_analysis_deg('GO', dge_cluster0vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster0_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster0vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster0_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster0vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Cluster1 ####
name_of_subset <- 'cluster1vshealthydermis'
gsea_cluster1_go <- run_gsea_analysis_deg('GO', dge_cluster1vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster1_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster1vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster1_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster1vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Cluster2 ####
name_of_subset <- 'cluster2vshealthydermis'
gsea_cluster2_go <- run_gsea_analysis_deg('GO', dge_cluster2vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster2_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster2vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster2_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster2vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Cluster3 ####
name_of_subset <- 'cluster3vshealthydermis'
gsea_cluster3_go <- run_gsea_analysis_deg('GO', dge_cluster3vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster3_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster3vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster3_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster3vshealthydermis, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Heatmaps ####
gsea_cluster0_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster0vshealthydermis_gobp_gsea_results.RDS'))
gsea_cluster1_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster1vshealthydermis_gobp_gsea_results.RDS'))
gsea_cluster2_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster2vshealthydermis_gobp_gsea_results.RDS'))
gsea_cluster3_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster3vshealthydermis_gobp_gsea_results.RDS'))

all_gsea_results <- list(gsea_cluster0_go, gsea_cluster1_go, gsea_cluster2_go, gsea_cluster3_go)
names(all_gsea_results) <- c('gsea_cluster0_go', 'gsea_cluster1_go', 'gsea_cluster2_go', 'gsea_cluster3_go')
gsea_clusters_heatmap <- gsea_get_heatmap(all_gsea_results, file_name = 'all_granuloma_clusters_go', pathway_cutoff = 50)

all_gsea_results <- list(gsea_cluster0_go, gsea_cluster2_go, gsea_cluster3_go)
names(all_gsea_results) <- c('gsea_cluster0_go', 'gsea_cluster2_go', 'gsea_cluster3_go')
gsea_clusters_heatmap <- gsea_get_heatmap(all_gsea_results, file_name = 'only3clusters_clusters_go', pathway_cutoff = 50)


gsea_cluster0_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster0vshealthydermis_reactome_gsea_results.RDS'))
gsea_cluster1_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster1vshealthydermis_reactome_gsea_results.RDS'))
gsea_cluster2_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster2vshealthydermis_reactome_gsea_results.RDS'))
gsea_cluster3_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster3vshealthydermis_reactome_gsea_results.RDS'))

all_gsea_results <- list(gsea_cluster0_reactome, gsea_cluster1_reactome, gsea_cluster2_reactome, gsea_cluster3_reactome)
names(all_gsea_results) <- c('gsea_cluster0_reactome', 'gsea_cluster1_reactome', 'gsea_cluster2_reactome', 'gsea_cluster3_reactome')
gsea_clusters_heatmap <- gsea_get_heatmap(all_gsea_results, file_name = 'all_granuloma_clusters_reactome', pathway_cutoff = 50)

gsea_cluster0_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster0vshealthydermis_kegg_gsea_results.RDS'))
gsea_cluster1_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster1vshealthydermis_kegg_gsea_results.RDS'))
gsea_cluster2_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster2vshealthydermis_kegg_gsea_results.RDS'))
gsea_cluster3_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster3vshealthydermis_kegg_gsea_results.RDS'))

all_gsea_results <- list(gsea_cluster0_kegg, gsea_cluster1_kegg, gsea_cluster2_kegg, gsea_cluster3_kegg)
names(all_gsea_results) <- c('gsea_cluster0_kegg', 'gsea_cluster1_kegg', 'gsea_cluster2_kegg', 'gsea_cluster3_kegg')
gsea_clusters_heatmap <- gsea_get_heatmap(all_gsea_results, file_name = 'all_granuloma_clusters_kegg', pathway_cutoff = 50)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Clusters vs the rest ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Cluster 0 vs rest ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster0vsrest'
reference_group <- 'OTHERCLUSTERS'
sce_cluster0vsrest <- sce[, which(sce$cluster0vsrest %in% c("CLUSTER0", "OTHERCLUSTERS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster0vsrest$cluster0vsrest <- droplevels(sce_cluster0vsrest$cluster0vsrest)

### 2. Exploratory analysis ####

comparison_title <- 'Outer core' #'Cluster 0 vs rest' 
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(x = factor(patient), fill = factor(cluster0vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(x = factor(patient), fill = factor(cluster0vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(DISEASE), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(DISEASE), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(spot_type), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(spot_type), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster0vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster0vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster0vsrest)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER0' = 'mediumvioletred', 'OTHERCLUSTERS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)

pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster0vsrest <- run_glmgampoi(sce_cluster0vsrest, design_formula = ~ patient + cluster0vsrest, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster0vsrest <- read.csv(paste0(glmgampoi_path, 'sce_cluster0vsrest_patient+cluster0vsrest_degresults.csv'), row.names = 1)
fit_cluster0vsrest <- readRDS(paste0(glmgampoi_path, 'sce_cluster0vsrest_patient+cluster0vsrest_fit.RDS'))

get_dge_summary(dge_cluster0vsrest)

# Most different genes (padj)
head(dge_cluster0vsrest[order(dge_cluster0vsrest$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster0vsrest, comparison_title, sce_cluster0vsrest)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster0vsrest', 
                           group1 = 'CLUSTER0', group2 = 'OTHERCLUSTERS', 
                           dge_df = dge_cluster0vsrest, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()

## Cluster 1 vs rest ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster1vsrest'
reference_group <- 'OTHERCLUSTERS'
sce_cluster1vsrest <- sce[, which(sce$cluster1vsrest %in% c("CLUSTER1", "OTHERCLUSTERS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster1vsrest$cluster1vsrest <- droplevels(sce_cluster1vsrest$cluster1vsrest)

### 2. Exploratory analysis ####

comparison_title <- 'Border-2' #'Cluster 1 vs rest'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(x = factor(patient), fill = factor(cluster1vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(x = factor(patient), fill = factor(cluster1vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(DISEASE), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(DISEASE), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(spot_type), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(spot_type), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster1vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster1vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster1vsrest)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER1' = 'mediumvioletred', 'OTHERCLUSTERS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)

pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster1vsrest <- run_glmgampoi(sce_cluster1vsrest, design_formula = ~ patient + cluster1vsrest, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster1vsrest <- read.csv(paste0(glmgampoi_path, 'sce_cluster1vsrest_patient+cluster1vsrest_degresults.csv'), row.names = 1)
fit_cluster1vsrest <- readRDS(paste0(glmgampoi_path, 'sce_cluster1vsrest_patient+cluster1vsrest_fit.RDS'))

get_dge_summary(dge_cluster1vsrest)

# Most different genes (padj)
head(dge_cluster1vsrest[order(dge_cluster1vsrest$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster1vsrest, comparison_title, sce_cluster1vsrest)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster1vsrest', 
                           group1 = 'CLUSTER1', group2 = 'OTHERCLUSTERS', 
                           dge_df = dge_cluster1vsrest, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()

## Cluster 2 vs rest ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster2vsrest'
reference_group <- 'OTHERCLUSTERS'
sce_cluster2vsrest <- sce[, which(sce$cluster2vsrest %in% c("CLUSTER2", "OTHERCLUSTERS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster2vsrest$cluster2vsrest <- droplevels(sce_cluster2vsrest$cluster2vsrest)

### 2. Exploratory analysis ####

comparison_title <- 'Inner core' #'Cluster 2 vs rest'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(x = factor(patient), fill = factor(cluster2vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(x = factor(patient), fill = factor(cluster2vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(DISEASE), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(DISEASE), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(spot_type), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(spot_type), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster2vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster2vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster2vsrest)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER2' = 'mediumvioletred', 'OTHERCLUSTERS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)

pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster2vsrest <- run_glmgampoi(sce_cluster2vsrest, design_formula = ~ patient + cluster2vsrest, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster2vsrest <- read.csv(paste0(glmgampoi_path, 'sce_cluster2vsrest_patient+cluster2vsrest_degresults.csv'), row.names = 1)
fit_cluster2vsrest <- readRDS(paste0(glmgampoi_path, 'sce_cluster2vsrest_patient+cluster2vsrest_fit.RDS'))

get_dge_summary(dge_cluster2vsrest)

# Most different genes (padj)
head(dge_cluster2vsrest[order(dge_cluster2vsrest$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster2vsrest, comparison_title, sce_cluster2vsrest)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster2vsrest', 
                           group1 = 'CLUSTER2', group2 = 'OTHERCLUSTERS', 
                           dge_df = dge_cluster2vsrest, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()


## Cluster 3 vs rest ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_cluster3vsrest'
reference_group <- 'OTHERCLUSTERS'
sce_cluster3vsrest <- sce[, which(sce$cluster3vsrest %in% c("CLUSTER3", "OTHERCLUSTERS"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_cluster3vsrest$cluster3vsrest <- droplevels(sce_cluster3vsrest$cluster3vsrest)

### 2. Exploratory analysis ####
comparison_title <- 'Border-1' #'Cluster 3 vs rest'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(x = factor(patient), fill = factor(cluster3vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(x = factor(patient), fill = factor(cluster3vsrest))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(sample_SPECIMEN), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(DISEASE), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(DISEASE), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(spot_type), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(spot_type), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_cluster3vsrest@colData), aes(factor(leiden_r1.3_patient), fill = factor(cluster3vsrest))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

all_plot_and_spatial <- ggarrange(
  ggplot(data = as.data.frame(sce@colData), aes(x = spatial1, y = spatial2, col = cluster3vsrest)) +
    geom_point(size=0.5) +
    scale_color_manual('', values = c('CLUSTER3' = 'mediumvioletred', 'OTHERCLUSTERS' = '#ADDD8E', 'other' = 'grey')) +
    scale_y_reverse() +
    labs(title = paste0("Spatial representation of Leiden clusters (r = ", '0.3', ')')) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(fill = NA, colour=NA),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none") +
    facet_wrap(.~sample, scales = "fixed", nrow = 2) +
    theme(legend.position = 'bottom'),
  all_plots, nrow = 2, heights = c(1, 4)
)

pdf(file = paste0(exploratory_analysis_path, name_of_subset, '.pdf'), height = 30, width = 20)
all_plot_and_spatial
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_cluster3vsrest <- run_glmgampoi(sce_cluster3vsrest, design_formula = ~ patient + cluster3vsrest, sce_subset_name = name_of_subset, reference = reference_group)
dge_cluster3vsrest <- read.csv(paste0(glmgampoi_path, 'sce_cluster3vsrest_patient+cluster3vsrest_degresults.csv'), row.names = 1)
fit_cluster3vsrest <- readRDS(paste0(glmgampoi_path, 'sce_cluster3vsrest_patient+cluster3vsrest_fit.RDS'))

get_dge_summary(dge_cluster3vsrest)

# Most different genes (padj)
head(dge_cluster3vsrest[order(dge_cluster3vsrest$padj), ])

### 4. Plotting DGE results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_cluster3vsrest, comparison_title, sce_cluster3vsrest)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# Heatmap
heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'cluster3vsrest', 
                           group1 = 'CLUSTER3', group2 = 'OTHERCLUSTERS', 
                           dge_df = dge_cluster3vsrest, number_top_genes = 50)

pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
heatmap
dev.off()


# GSEA  ------------------------------------------------
setwd(pea_path)

## Cluster0 ####
name_of_subset <- 'cluster0vsrest'
gsea_cluster0_go <- run_gsea_analysis_deg('GO', dge_cluster0vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster0_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster0vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster0_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster0vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Cluster1 ####
name_of_subset <- 'cluster1vsrest'
gsea_cluster1_go <- run_gsea_analysis_deg('GO', dge_cluster1vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster1_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster1vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster1_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster1vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Cluster2 ####
name_of_subset <- 'cluster2vsrest'
gsea_cluster2_go <- run_gsea_analysis_deg('GO', dge_cluster2vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster2_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster2vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster2_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster2vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Cluster3 ####
name_of_subset <- 'cluster3vsrest'
gsea_cluster3_go <- run_gsea_analysis_deg('GO', dge_cluster3vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster3_kegg <- run_gsea_analysis_deg('KEGG', dge_cluster3vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_cluster3_reactome <- run_gsea_analysis_deg('Reactome', dge_cluster3vsrest, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Heatmaps ####
gsea_cluster0_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster0vsrest_gobp_gsea_results.RDS'))
gsea_cluster1_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster1vsrest_gobp_gsea_results.RDS'))
gsea_cluster2_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster2vsrest_gobp_gsea_results.RDS'))
gsea_cluster3_go <- readRDS(paste0(pea_path, '/GSEA/', 'cluster3vsrest_gobp_gsea_results.RDS'))

all_gsea_results <- list(gsea_cluster0_go, gsea_cluster1_go, gsea_cluster2_go, gsea_cluster3_go)
names(all_gsea_results) <- c('Outer core', 'Border-2', 'Inner core', 'Border-1') #c('Cluster 0', 'Cluster 1', 'Cluster 2', 'Cluster 3')
gsea_clusters_heatmap <- gsea_get_heatmap(all_gsea_results, file_name = 'all_granuloma_clustersvsrest_go', pathway_cutoff = 50)

gsea_cluster0_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster0vsrest_reactome_gsea_results.RDS'))
gsea_cluster1_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster1vsrest_reactome_gsea_results.RDS'))
gsea_cluster2_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster2vsrest_reactome_gsea_results.RDS'))
gsea_cluster3_reactome <- readRDS(paste0(pea_path, '/GSEA/', 'cluster3vsrest_reactome_gsea_results.RDS'))

all_gsea_results <- list(gsea_cluster0_reactome, gsea_cluster1_reactome, gsea_cluster2_reactome, gsea_cluster3_reactome)
names(all_gsea_results) <- c('Outer core', 'Border-2', 'Inner core', 'Border-1') #c('Cluster 0', 'Cluster 1', 'Cluster 2', 'Cluster 3')
gsea_clusters_heatmap <- gsea_get_heatmap(all_gsea_results, file_name = 'all_granuloma_clustersvsrest_reactome', pathway_cutoff = 50)

gsea_cluster0_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster0vsrest_kegg_gsea_results.RDS'))
gsea_cluster1_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster1vsrest_kegg_gsea_results.RDS'))
gsea_cluster2_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster2vsrest_kegg_gsea_results.RDS'))
gsea_cluster3_kegg <- readRDS(paste0(pea_path, '/GSEA/', 'cluster3vsrest_kegg_gsea_results.RDS'))

all_gsea_results <- list(gsea_cluster0_kegg, gsea_cluster1_kegg, gsea_cluster2_kegg, gsea_cluster3_kegg)
names(all_gsea_results) <- c('Outer core', 'Border-2', 'Inner core', 'Border-1') #c('Cluster 0', 'Cluster 1', 'Cluster 2', 'Cluster 3')
gsea_clusters_heatmap <- gsea_get_heatmap(all_gsea_results, file_name = 'all_granuloma_clustersvsrest_kegg', pathway_cutoff = 50)



# 1vs1 scatterplots ####
dge_cluster0vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster0vshealthydermis_patient+cluster0vshealthydermis_degresults.csv'), row.names = 1)
dge_cluster1vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster1vshealthydermis_patient+cluster1vshealthydermis_degresults.csv'), row.names = 1)
dge_cluster2vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster2vshealthydermis_patient+cluster2vshealthydermis_degresults.csv'), row.names = 1)
dge_cluster3vshealthydermis <- read.csv(paste0(glmgampoi_path, 'sce_cluster3vshealthydermis_patient+cluster3vshealthydermis_degresults.csv'), row.names = 1)

colours_scatterplot <- list('Both' = c("#D53E4F", "#3288BD"),
                            `Outer core` = c('#F46D43', '#f4c643'),  
                            'Border-2' = c('#78C679', '#aeddae'),
                            `Inner core` = c('#9E0142', '#fe79b0'), 
                            'Border-1' = c('#006837', '#00b660'))

get_scatterplot <- function(dge_x, dge_y, name_x, name_y){
  # For debug
  # dge_x <- dge_cluster0vshealthydermis
  # dge_y <- dge_cluster2vshealthydermis
  # name_x <- 'Outer core'
  # name_y <- 'Inner core'
  
  # 1. Get the signed -log10(pval)
  dge_x$signedpval <- -log10(dge_x$pval) * sign(dge_x$log2fc)
  dge_y$signedpval <- -log10(dge_y$pval) * sign(dge_y$log2fc)
  
  # 2. Merge both dataframes
  df <- merge(dge_x, dge_y, by = 'gene_symbol')
  
  # Remove genes that are not differentially expressed in either one
  # df <- df[df$diffexpressed.x != 'NO' & df$diffexpressed.y != 'NO',]
  
  # 3. Plot
  labNames <- c(name_x, name_y)
  
  # Set max and min for x and y axis (don't consider infinity values)
  max_y_scale_x <- max(df[!is.infinite(df$signedpval.x),]$signedpval.x) * 1.1
  max_y_scale_y <- max(df[!is.infinite(df$signedpval.y),]$signedpval.y) * 1.1
  min_y_scale_x <- min(df[!is.infinite(df$signedpval.x),]$signedpval.x) * 1.1
  min_y_scale_y <- min(df[!is.infinite(df$signedpval.y),]$signedpval.y) * 1.1
  
  genes_with_inf_minuslogpadj <- df[is.infinite(df$signedpval.x)|is.infinite(df$signedpval.y),'gene_symbol'] # get the names of genes which have inf signedpval
  
  # Change infinity values to the max so they will be plotted nicely
  # df[is.infinite(df$signedpval.x) & df$signedpval.x > 0,'signedpval.x'] <- max_y_scale_x * 0.9
  # df[is.infinite(df$signedpval.x) & df$signedpval.x < 0,'signedpval.x'] <- min_y_scale_x * 0.9
  # df[is.infinite(df$signedpval.y) & df$signedpval.y > 0,'signedpval.y'] <- max_y_scale_y * 0.9
  # df[is.infinite(df$signedpval.y) & df$signedpval.y > 0,'signedpval.y'] <- min_y_scale_y * 0.9
  
  # Set colours
  df <- df %>% mutate(colours = case_when(diffexpressed.x == 'UP' & diffexpressed.y != 'UP' ~ paste0('UP in ', name_x), 
                                          diffexpressed.y == 'UP' & diffexpressed.x != 'UP' ~ paste0('UP in ', name_y),
                                          diffexpressed.x == 'DOWN' & diffexpressed.y != 'DOWN' ~ paste0('DOWN in ', name_x), 
                                          diffexpressed.y == 'DOWN' & diffexpressed.x != 'DOWN' ~ paste0('DOWN in ', name_y),
                                          diffexpressed.x == 'UP' & diffexpressed.y == 'UP' ~ paste0('UP in ', name_x, ' & ', name_y), 
                                          diffexpressed.y == 'DOWN' & diffexpressed.x == 'DOWN' ~ paste0('DOWN in ', name_x, ' & ', name_y)))
  df[is.na(df$colours), ]$colours <- 'Other'
  df$colours <- factor(df$colours, levels = c(paste0('UP in ', name_x, ' & ', name_y), paste0('UP in ', name_x), paste0('UP in ', name_y),
                                              paste0('DOWN in ', name_x), paste0('DOWN in ', name_y), paste0('DOWN in ', name_x, ' & ', name_y), 'Other'), 
                       ordered = TRUE)
  
  only_diffexpressed_x <- df[df$colours == paste0('UP in ', name_x) | df$colours == paste0('DOWN in ', name_x),]
  only_diffexpressed_y <- df[df$colours == paste0('UP in ', name_y) | df$colours == paste0('DOWN in ', name_y),]
  
  # Create a new column "delabel" to de, that will contain the name of the some genes (NA in case they are not)
  
  # This will show annotations for top 30 upregulated and top 30 downregulated by p adj + 10 genes that are up or downregulated in only one of the categories 
  # df$genelabel <- ifelse(df$gene_symbol %in% head(df[order(df$padj.x), "gene_symbol"], 30) |
  #                          df$gene_symbol %in% head(df[order(df$padj.y), "gene_symbol"], 30) |
  #                          df$gene_symbol %in% head(only_diffexpressed_x[order(only_diffexpressed_x$padj.x),'gene_symbol'], 10) |
  #                          df$gene_symbol %in% head(only_diffexpressed_y[order(only_diffexpressed_y$padj.y),'gene_symbol'], 10), 
  #                        df$gene_symbol, NA)
  
  df$genelabel <- ifelse(df$gene_symbol %in% head(only_diffexpressed_x[order(only_diffexpressed_x$padj.x),'gene_symbol'], 15) |
                           df$gene_symbol %in% head(only_diffexpressed_y[order(only_diffexpressed_y$padj.y),'gene_symbol'], 15), #|
                         #df$gene_symbol %in% genes_with_inf_minuslogpadj,
                         df$gene_symbol, NA)
  
  # # Add some of the genes that are up or downregulated for only one disease (previous line selects only up and down in both genes, in most cases) 
  # df$genelabel <- ifelse(df$gene_symbol %in% head(df[order(df[df$diffexpressed.x != 'NO',]$padj.x), "gene_symbol"], 10) |
  #                          df$gene_symbol %in% head(df[order(df[df$diffexpressed.y != 'NO',]$padj.y), "gene_symbol"], 10), df$gene_symbol, NA)
  
  # Colours for each of the diseases
  plot_colours <- c(colours_scatterplot$Both[1], colours_scatterplot[[name_x]][1], colours_scatterplot[[name_y]][1], 
                    colours_scatterplot[[name_x]][2], colours_scatterplot[[name_y]][2], colours_scatterplot$Both[2], 'grey')
  
  # Scatterplot
  scatterplot <- ggplot(data = df, aes(x = signedpval.x, y = signedpval.y, col = colours, label = genelabel)) +
    geom_point(alpha = 0.7) + 
    coord_cartesian(ylim = c(min_y_scale_y, max_y_scale_y), xlim = c(min_y_scale_x, max_y_scale_x), clip = 'off') + # since some genes can have minuslog10padj of inf, we set this limits
    theme_classic() +
    theme(axis.title=element_text(size = 25, face = "bold"), axis.text = element_text(size = 20), legend.text = element_text(size = 13),
          plot.margin = margin(2.5, 2.5, 2.5, 2.5, "cm")) +
    geom_label_repel(max.overlaps = Inf, force = 0.8, box.padding = 0.5, show.legend = FALSE, size = 6) + # To show all labels
    scale_color_manual(values = plot_colours) + #values = c(brewer.pal(6, 'Spectral'), 'grey')
    #geom_label_repel(fill = 'white') +
    # scale_x_continuous(breaks = seq(min_y_scale_x, max_y_scale_x, (max_y_scale_x - min_y_scale_x)/6)) + 
    # scale_y_continuous(breaks = seq(min_y_scale_y, max_y_scale_y, (max_y_scale_y - min_y_scale_y)/6)) +
    #ggtitle(paste(title_of_comparison)) + 
    labs(color = ' ', x = bquote(.(labNames[1]) *": signed p-value"), y = bquote(.(labNames[2]) *": signed p-value")) +
    # labs(color = ' ', x = bquote(.(labNames[1]) *":"~ "-log"[10]*"pval"), y = bquote(.(labNames[2]) *":"~ "-log"[10]*"pval")) +
    guides(fill = guide_legend(title = "", override.aes = aes(label = "")))
  #geom_text(colour = 'black', fill = 'white') 
  #theme(panel.background = element_rect(fill = 'lightgray', color = 'black'))
  
  return(scatterplot)
}

innervsoutercore <- get_scatterplot(dge_cluster2vshealthydermis, dge_cluster0vshealthydermis, 'Inner core', 'Outer core')
outercorevsborder1 <- get_scatterplot(dge_cluster0vshealthydermis, dge_cluster3vshealthydermis, 'Outer core', 'Border-1')
border2vsinnercore <- get_scatterplot(dge_cluster1vshealthydermis, dge_cluster2vshealthydermis, 'Border-2', 'Inner core')
border2vsoutercore <- get_scatterplot(dge_cluster1vshealthydermis, dge_cluster0vshealthydermis, 'Border-2', 'Outer core')
border2vsborder1 <- get_scatterplot(dge_cluster1vshealthydermis, dge_cluster3vshealthydermis, 'Border-2', 'Border-1') 


spotdata <- as.data.frame(sce@colData)
example_areas <- ggplot(data = spotdata[spotdata$sample_SPECIMEN == 'P18554_1002_50107-A'|spotdata$sample_SPECIMEN == 'P18554_1003_95096-A', ], aes(x = spatial1, y = spatial2, col = leiden_r0.3_g_bc)) +
  geom_point(size = 4) +
  #scale_color_manual('', values = col.set$leiden_r1.3) +
  scale_y_reverse() +
  theme_bw() +
  scale_color_manual('', values = col.set$leiden_r0.3_g_bc, breaks = c(2, 0, 3, 1), labels = c('Inner core', 'Outer core', 'Border-1', 'Border-2')) +
  guides(col = guide_legend(nrow = 2)) +
  theme(aspect.ratio = 1,
        axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), line = element_blank(),
        legend.position = 'bottom', legend.text = element_text(size = 30),
        panel.background = element_blank(),
        plot.margin = margin(2.5, 2.5, 2.5, 2.5, "cm")) 
  


pdf(file = paste0(DEG_results_path, 'comparison_scatterplots_4regions.pdf'), height = 30, width = 40, onefile=F)
ggarrange(
  innervsoutercore + theme(legend.position = 'None'),
  outercorevsborder1 + theme(legend.position = 'None'),
  example_areas,
  as_ggplot(get_legend(innervsoutercore + geom_point(size = 5) + theme(legend.position = 'bottom'))),
  as_ggplot(get_legend(outercorevsborder1 + geom_point(size = 5) + theme(legend.position = 'bottom'))),
  NULL,
  border2vsinnercore + theme(legend.position = 'None') + coord_cartesian(ylim = c(-300,300)),
  border2vsoutercore + theme(legend.position = 'None') + coord_cartesian(ylim = c(-300,300)),
  border2vsborder1 + theme(legend.position = 'None') + coord_cartesian(ylim = c(-300,300)),
  as_ggplot(get_legend(border2vsinnercore + geom_point(size = 5) + theme(legend.position = 'bottom'))),
  as_ggplot(get_legend(border2vsoutercore + geom_point(size = 5) + theme(legend.position = 'bottom'))),
  as_ggplot(get_legend(border2vsborder1 + geom_point(size = 5) + theme(legend.position = 'bottom'))),
  ncol = 3, nrow = 4, heights = c(1, 0.1, 1, 0.1),
  labels = c('A', 'B', '', '', '', '', 'C', 'D', 'E', '', '', ''), font.label = list(size = 30, color = "black")
)
dev.off()



