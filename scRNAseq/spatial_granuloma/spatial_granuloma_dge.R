# Name: spatial_granuloma_dge.R
# Author: Laura Twomey
# Date of creation: 09 May 2022
# Differential gene expression analysis on spatial granuloma samples

# Clean Environment ----------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects including hidden objects
gc() #free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F) #avoid truncated output in R console and scientific notation

# Set working directory and other relevant paths ----------------------------------------------------
setwd("/Users/mendenlab/work/spatial_granuloma/scripts/summary")
getwd()
exploratory_analysis_path = "/Volumes/Drive/spatial_granuloma/output/DEG/DEG_exploratory/"
DEG_results_path = "/Volumes/Drive/spatial_granuloma/output/DEG/DEG_results/"
glmgampoi_path = "/Volumes/Drive/spatial_granuloma/output/DEG/glmGamPoi_results/"
deg_path <- "/Volumes/Drive/spatial_granuloma/output/DEG/"

# Libraries and functions ------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(anndata)
library(reticulate)
library(hdf5r)
library(SingleCellExperiment)
library(glmGamPoi) # dor dge
library(ggpubr) # for ggarange
library(ggrepel) # for nice volcano plot
library(pheatmap) # for heatmap

# Colours --------------------------------------------------------------------
col.set.ylgn <- brewer.pal(n=9, "YlGn")
col.set.spec <- brewer.pal(n=11, "Spectral")

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
                lesions = c("LESIONAL" = 'tomato',
                            "NON LESIONAL" = 'darkolivegreen',
                            "UNDETERMINED" = 'lightgrey'),
                patients = c("91253" = "#F46D43",#'red',#granuloma annulare
                             "45703" = "#D53E4F",#'indianred',#granuloma annulare
                             "50107" = "#9E0142",#'firebrick',#granuloma annulare
                             "95096" = "orchid", # necrobiosis lipoidica
                             "82301" = "mediumvioletred", #sarcoidosis-cutaneous
                             "72859" = "blueviolet"), #sarcoidosis-suspected
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
                leiden_r0.5 = c("0" = "#E0EEE0",
                                "1" = "#41AB5D",
                                "2" = "#ADDD8E",
                                "3" = "yellow",
                                "4" = "orchid",
                                "5" = "dodgerblue",
                                "6" = "blueviolet",
                                "7" = "#D53E4F",
                                "8" = "#F46D43",
                                "9" = "mistyrose",
                                "10" = "blue",
                                "11" = "#9E0142",
                                "12" = "darkcyan"),
                leiden_r0.8 = c("0"  = 'yellow',
                                "1"  = '#ADDD8E',
                                "2"  = '#E0EEE0',
                                "3"  = '#D9F0A3',
                                "4"  = "#EAADE7",
                                "5"  = 'deepskyblue',
                                "6"  = "orchid",
                                "7"  = "#D53E4F",
                                "8"  = "darkgreen",
                                "9"  = "#F46D43",
                                "10" = "blueviolet",
                                "11" = 'skyblue',
                                "12" = 'blue',
                                "13" = "mediumvioletred",
                                "14" = 'mistyrose',
                                "15" = "#9E0142",
                                "16" = 'darkcyan'),
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
                                "19" = '#006837')) #derdepth1


# Functions  --------------------------------------------------------------------

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
    'Dermis_lesional' = legend_title <-'Lesional dermis',
    'Ga_granuloma' = legend_title <-'G.annulare',
    'Sa_granuloma' = legend_title <-'Sarcoidosis',
    'Nl_granuloma' = legend_title <-'Necrobiosis lipoidica',
    'Ga_border' = legend_title <- 'G.annulare', # for the border and core comparison we will show the upregulated for both border and core
    'Sa_border' = legend_title <- 'Sarcoidosis',  # for the border and core comparison we will show the upregulated for both border and core
    'Leiden_granuloma' = legend_title <- 'Granuloma (Leiden)',
    'Granuloma' = legend_title <- 'Granuloma (manual)'
  )
  
  # plot adding up all layers we have seen so far
  volcano_plot <- ggplot(data = dataframe[!is.infinite(dataframe$log2fc),], 
                         aes(x = log2fc, y = minuslog10padj, col = diffexpressed, label = delabel)) +
    geom_point() + 
    coord_cartesian(ylim = c(0, max_y_scale)) + # since some genes can have minuslog10padj of inf, we set this limits
    theme_classic() +
    geom_text_repel(max.overlaps = Inf) + # To show all labels
    #scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated")) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") + 
    ggtitle(paste(title_of_comparison)) + 
    labs(color = legend_title, x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR"))
  
  # Let's fix the legend separately to deal with exceptions
  if(interesting_category == 'Ga_border'|interesting_category == 'Sa_border'){
    volcano_plot <- volcano_plot + 
      scale_color_manual(values = c("red", "black", "#008000"), labels = c("Upregulated in core", "Unchanged", "Upregulated in border")) 
  } else {
    volcano_plot <- volcano_plot + 
      scale_color_manual(values = c("blue", "black", "red"), labels = c("Downregulated", "Unchanged", "Upregulated"))
  }
  
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

# This function takes the normalised counts, the column name from the normalised_counts headers dataframe 
# with the values for the comparison column, and the names of the two values 
# (groups) we want to compare, and returns a heatmap
# dge_df and number_top_genes are needed to plot only the top genes (e.g., top 50 genes sorted by padj from the DEG results)
# get_deg_heatmap <- function(norm_counts, comparison_column, group1, group2, dge_df, number_top_genes){
#   #For debug
#   # norm_counts <- normalised_counts
#   # comparison_column <- 'dermis_lesvsnonles'
#   # group1 <- 'dermis_lesional'
#   # group2 <- 'dermis_nonlesional'
#   # dge_df <- dge_dermislesionalvsnonlesional
#   # number_top_genes <- 50
#   
#   # Some checks to avoid weird errors
#   if(sum(c(group1, group2) %in% unique(normalised_counts_headers[[comparison_column]])) < 2){
#     stop('Recheck your comparison groups, they must be contained in the comparison column')
#   }
#   
#   # First we must filter the genes for the two groups
#   heatmapdf <- norm_counts # get the normalised counts
#   colnames(heatmapdf) <- normalised_counts_headers[[comparison_column]] # set as column names the values of the comparison variable
#   heatmapdf <- as.data.frame(do.call(cbind, # merge columns with the same name and calculate the mean expression value for all spots of each group
#                                      by(t(heatmapdf), INDICES = names(heatmapdf), FUN = colMeans)))
#   # Subset to only the two groups we want to compare
#   heatmapdf <- heatmapdf[, grepl(paste0(group1, '|', group2), colnames(heatmapdf))]
#   # Subset to only the top genes
#   top_genes <- head(dge_df[order(dge_df$padj), "gene_symbol"], number_top_genes) # get the 50 most different genes (padj)
#   heatmapdf <- heatmapdf[top_genes, ]
#   
#   # Plot heatmap
#   heatmap <- pheatmap(heatmapdf,
#                       show_colnames = T, show_rownames = T,
#                       cluster_rows = TRUE, cluster_cols = TRUE,
#                       clustering_distance_cols = 'euclidean',
#                       clustering_distance_rows = 'euclidean',
#                       clustering_method = 'ward.D')
#   return(heatmap)
# }


# Import data --------------------------------------------------------------
reexport_dge_df = F # only set as T if you want to re-export everything, else only read in the data
if(reexport_dge_df == T){
  ad <- read_h5ad("../../results/current/final/adata_deg.h5")
  names(ad$obsm)
  names(ad$obs)
  # Create a dataframe with the columns we need
  df <- as.data.frame(cbind(sample = ad$obs$sample,
              slide = ad$obs$slide,
              array_row = ad$obs$array_row,
              array_col = ad$obs$array_col,
              DISEASE = ad$obs$DISEASE,
              SAMPLE = ad$obs$SAMPLE,
              patient = ad$obs$patient,
              skin_layer = ad$obs$skin_layer,
              spot_type = ad$obs$spot_type,
              GRANULOMA = ad$obs$GRANULOMA,
              size_factors = ad$obs$size_factors,
              leiden_r1.3_patient = ad$obs$leiden_r1.3_patient,
              manual_border_2 = ad$obs$manual_border_2,
              leiden_core_granuloma = ad$obs$leiden_core_granuloma,
              leiden_border_granuloma = ad$obs$leiden_border_granuloma,
              leiden_epidermis = ad$obs$leiden_epidermis,
              leiden_nongranuloma = ad$obs$leiden_nongranuloma,
              dermis_lesional = ad$obs$dermis_lesional,
              dermis_nonlesional = ad$obs$dermis_nonlesional,
              dermis_lesional_noborder = ad$obs$dermis_lesional_noborder,
              dermis_nonlesional_noborder = ad$obs$dermis_nonlesional_noborder,
              epidermis_interface = ad$obs$epidermis_interface,
              dermis_lesvsnonles = ad$obs$dermis_lesvsnonles,
              ga_corevsborder = ad$obs$ga_corevsborder,
              sa_corevsborder = ad$obs$sa_corevsborder,
              ga_gvsd = ad$obs$ga_gvsd,
              sa_gvsd = ad$obs$sa_gvsd,
              nl_gvsd = ad$obs$nl_gvsd,
              manual_gvsd = ad$obs$manual_gvsd,
              leiden_gvsd = ad$obs$leiden_gvsd,
              epivsdermis = ad$obs$epivsdermis,
              indexes = ad$obs$indexes,
              spatial1 = ad$obsm$spatial[, 1], 
              spatial2 = ad$obsm$spatial[, 2]))
  # Fix up the dataframe columns
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
  df$skin_layer <- factor(as.character(df$skin_layer),
                          levels = names(col.set$dermis),
                          ordered = T)
  
  df[df$sample=="P17851_1001",
     c("spatial1", "spatial2")] <- 30000-df[df$sample=="P17851_1001",
                                            c("spatial1", "spatial2")]
  df[df$sample=="P18554_1002",
     c("spatial1", "spatial2")] <- 30000-df[df$sample=="P18554_1002",
                                            c("spatial1", "spatial2")]
  df[df$sample=="P18554_1006",
     c("spatial1", "spatial2")] <- 30000-df[df$sample=="P18554_1006",
                                            c("spatial1", "spatial2")]
  
  df$patient <- factor(df$patient,
                       levels = c("91253", "45703", "50107",
                                  "95096", "82301", "72859"),
                       ordered = T)
  
  df$spot_type <- df$spot_type
  df$spot_type[df$patient=="91253" & df$spot_type=="GA"] <- "GA (1)"
  df$spot_type[df$patient=="45703" & df$spot_type=="GA"] <- "GA (2)"
  df$spot_type[df$patient=="50107" & df$spot_type=="GA"] <- "GA (3)"
  
  # # Create columns for each deg comparison we want to do
  # # Dermis lesional vs non-lesional
  # df$dermis_lesvsnonles <- factor(ifelse(df$dermis_lesional_noborder == 1, 'dermis_lesional', 
  #                                 ifelse(df$dermis_nonlesional_noborder == 1, 'dermis_nonlesional', 'other')))
  # 
  # # GA: border vs core
  # df$ga_corevsborder <- factor(ifelse((df$DISEASE == "granuloma annulare" & df$leiden_core_granuloma == 1), 'ga_core', 
  #                                        ifelse((df$DISEASE == "granuloma annulare" & df$leiden_border_granuloma == 1), 'sa_border', 'other')))
  # 
  # # SA: border vs core
  # df$sa_corevsborder <- factor(ifelse((df$DISEASE == "sarcoidosis" & df$leiden_core_granuloma == 1), 'sa_core', 
  #                                     ifelse((df$DISEASE == "sarcoidosis" & df$leiden_border_granuloma == 1), 'sa_border', 'other')))
  # 
  # # GA: granuloma vs dermis
  # df$ga_gvsd <- factor(ifelse((df$DISEASE == "granuloma annulare" & df$GRANULOMA == 1), 'ga_granuloma', 
  #                                     ifelse(df$dermis_nonlesional_noborder == 1, 'healthydermis', 'other')))
  # 
  # # SA: granuloma vs dermis
  # df$sa_gvsd <- factor(ifelse((df$DISEASE == "sarcoidosis" & df$GRANULOMA == 1), 'sa_granuloma', 
  #                             ifelse(df$dermis_nonlesional_noborder == 1, 'healthydermis', 'other')))
  # 
  # # NL: granuloma vs dermis
  # df$nl_gvsd <- factor(ifelse((df$DISEASE == "necrobiosis lipoidica" & df$GRANULOMA == 1), 'nl_granuloma', 
  #                             ifelse(df$dermis_nonlesional_noborder == 1, 'healthydermis', 'other')))
  # 
  # # Granuloma manual: granuloma vs dermis
  # df$manual_gvsd <- factor(ifelse(df$GRANULOMA == 1, 'granuloma', 
  #                             ifelse(df$dermis_nonlesional_noborder == 1, 'healthydermis', 'other')))
  # 
  # # Granuloma Leiden: granuloma vs dermis
  # df$leiden_gvsd <- factor(ifelse((df$leiden_core_granuloma == 1 | df$leiden_border_granuloma == 1), 'leiden_granuloma', 
  #                                 ifelse(df$dermis_nonlesional_noborder == 1, 'healthydermis', 'other')))
  # 
  # # Epidermis vs dermis
  # df$epivsdermis <- factor(ifelse(df$epidermis_interface == 1, 'epidermis', 
  #                                 ifelse(df$dermis_nonlesional_noborder == 1, 'dermis', 'other')))
 
  # We only need to save this once in the deg folder and the final results folder
  if(!is.null(df)){
    saveRDS(df,
            file = "/Users/mendenlab/work/spatial_granuloma/results/current/final/Granuloma_DGE_df.RDS")
    saveRDS(df,
            file = "/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_df.RDS")
  }
  
  # Convert adata to SingleCellExperiment from raw counts
  sce <- SingleCellExperiment(
    assays      = list(counts = t(ad$layers[["counts"]])),
    colData     = ad$obs,
    rowData     = ad$var,
    reducedDims = list(#diffmap = ad$obsm[["X_diffmap"]],
      #graph_fr = ad$obsm[["X_draw_graph_fr"]],
      PCA = ad$obsm[["X_pca"]],
      #tSNE = ad$obsm[["X_tsne"]],
      umap = ad$obsm[["X_umap"]],
      spatial = ad$obsm[["spatial"]])
  )
  saveRDS(sce,
          file = "/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_sce.RDS")
  
  # Convert adata to SingleCellExperiment from normalised counts
  sce_normalised <- SingleCellExperiment(
    assays      = list(counts = t(ad$X)),
    colData     = ad$obs,
    rowData     = ad$var,
    reducedDims = list(#diffmap = ad$obsm[["X_diffmap"]],
      #graph_fr = ad$obsm[["X_draw_graph_fr"]],
      PCA = ad$obsm[["X_pca"]],
      #tSNE = ad$obsm[["X_tsne"]],
      umap = ad$obsm[["X_umap"]],
      spatial = ad$obsm[["spatial"]])
  )
  saveRDS(sce_normalised,
          file = "/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_sce_normalised.RDS")
  
  normalised_counts <- data.frame(as.matrix(counts(sce_normalised)))
  saveRDS(normalised_counts,
          file = "/Volumes/Drive/spatial_granuloma/output/DEG/normalised_counts.RDS")
  normalised_counts_headers <- sce_normalised@colData@listData
  saveRDS(normalised_counts_headers,
          file = "/Volumes/Drive/spatial_granuloma/output/DEG/normalised_counts_headers.RDS")
  
} else {
  df <- readRDS("/Users/mendenlab/work/spatial_granuloma/results/current/final/Granuloma_DGE_df.RDS")
  sce <- readRDS("/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_sce.RDS")
  #sce_normalised <- readRDS("/Volumes/Drive/spatial_granuloma/output/DEG/Granuloma_DGE_sce_normalised.RDS")
  normalised_counts <- readRDS("/Volumes/Drive/spatial_granuloma/output/DEG/normalised_counts.RDS")
  normalised_counts_headers <- readRDS("/Volumes/Drive/spatial_granuloma/output/DEG/normalised_counts_headers.RDS")
}


# DGE --------------------------------------------------------------

## Dermis vs Epidermis (with interface) ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_epivsdermis_withinterface'
reference_group <- 'healthydermis'
sce_epivsdermis_withinterface <- sce[, which(sce$epivsdermis %in% c("epidermis", "healthydermis"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_epivsdermis_withinterface$epivsdermis <- droplevels(sce_epivsdermis_withinterface$epivsdermis)
#sce_epivsdermis_withinterface$epivsdermis <- factor(sce_epivsdermis_withinterface$epivsdermis, levels = c("healthydermis", "epidermis"), ordered = TRUE)

### 2. Exploratory analysis ####

comparison_title <- 'Dermis vs Epidermis (with interface)'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(epivsdermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withinterface@colData), aes(x = factor(patient), fill = factor(epivsdermis))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withinterface@colData), aes(factor(sample_SPECIMEN), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withinterface@colData), aes(factor(DISEASE), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withinterface@colData), aes(factor(spot_type), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withinterface@colData), aes(factor(leiden_r1.3_patient), fill = factor(epivsdermis))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_epivsdermis_withinterface <- run_glmgampoi(sce_epivsdermis_withinterface, design_formula = ~ patient + epivsdermis, sce_subset_name = name_of_subset, reference = reference_group)
dge_epivsdermis_withinterface <- read.csv(paste0(glmgampoi_path, 'sce_epivsdermis_withinterface_patient+epivsdermis_degresults.csv'), row.names = 1)
fit_epivsdermis_withinterface <- readRDS(paste0(glmgampoi_path, 'sce_epivsdermis_withinterface_patient+epivsdermis_fit.RDS'))

get_dge_summary(dge_epivsdermis_withinterface)

# Most different genes (padj)
head(dge_epivsdermis_withinterface[order(dge_epivsdermis_withinterface$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_epivsdermis_withinterface, comparison_title, sce_subset = sce_epivsdermis_withinterface)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'epivsdermis', 
#                 group1 = 'epidermis', group2 = 'healthydermis', 
#                 dge_df = dge_epivsdermis_withinterface, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 15, width = 10)
# heatmap
# dev.off()

## Dermis vs Epidermis (without interface) ####

### 1. Sub-setting data ####
name_of_subset <- 'sce_epivsdermis_withoutinterface'
reference_group <- 'healthydermis'
sce_epivsdermis_withoutinterface <- sce[, which(sce$epivsdermis_withoutinterface %in% c("epidermis", "healthydermis"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_epivsdermis_withoutinterface$epivsdermis_withoutinterface <- droplevels(sce_epivsdermis_withoutinterface$epivsdermis_withoutinterface)

### 2. Exploratory analysis ####

comparison_title <- 'Dermis vs Epidermis (without interface)'
n_col = 2
n_row = 5
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(epivsdermis_withoutinterface))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withoutinterface@colData), aes(x = factor(patient), fill = factor(epivsdermis_withoutinterface))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withoutinterface@colData), aes(factor(sample_SPECIMEN), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withoutinterface@colData), aes(factor(DISEASE), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withoutinterface@colData), aes(factor(spot_type), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_epivsdermis_withoutinterface@colData), aes(factor(leiden_r1.3_patient), fill = factor(epivsdermis_withoutinterface))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_epivsdermis_withoutinterface <- run_glmgampoi(sce_epivsdermis_withoutinterface, design_formula = ~ patient + epivsdermis_withoutinterface, sce_subset_name = name_of_subset, reference = reference_group)
dge_epivsdermis_withoutinterface <- read.csv(paste0(glmgampoi_path, 'sce_epivsdermis_withoutinterface_patient+epivsdermis_withoutinterface_degresults.csv'), row.names = 1)
fit_epivsdermis_withoutinterface <- readRDS(paste0(glmgampoi_path, 'sce_epivsdermis_withoutinterface_patient+epivsdermis_withoutinterface_fit.RDS'))

get_dge_summary(dge_epivsdermis_withoutinterface)

# Most different genes (padj)
head(dge_epivsdermis_withoutinterface[order(dge_epivsdermis_withoutinterface$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_epivsdermis_withoutinterface, comparison_title, sce_subset = sce_epivsdermis_withoutinterface)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'epivsdermis_withoutinterface', 
#                            group1 = 'epidermis', group2 = 'healthydermis', 
#                            dge_df = dge_epivsdermis_withoutinterface, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## Dermis lesional vs non-lesional ####
# Covariates: dermis depth (optional: disease, patient)

### 1. Sub-setting data ####
name_of_subset <- 'sce_dermislesionalvsnonlesional'
reference_group <- 'dermis_nonlesional'
sce_dermislesionalvsnonlesional <- sce[, which(sce$dermis_lesvsnonles %in% c("dermis_lesional", "dermis_nonlesional"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_dermislesionalvsnonlesional$dermis_lesvsnonles <- droplevels(sce_dermislesionalvsnonlesional$dermis_lesvsnonles)

### 2. Exploratory analysis ####

comparison_title <- 'Lesional dermis vs healthy dermis'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(dermis_lesvsnonles))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_dermislesionalvsnonlesional@colData), aes(x = factor(patient), fill = factor(dermis_lesvsnonles))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_dermislesionalvsnonlesional@colData), aes(factor(sample_SPECIMEN), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_dermislesionalvsnonlesional@colData), aes(factor(DISEASE), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_dermislesionalvsnonlesional@colData), aes(factor(spot_type), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_dermislesionalvsnonlesional@colData), aes(factor(leiden_r1.3_patient), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_dermislesionalvsnonlesional@colData), aes(factor(skin_layer), fill = factor(dermis_lesvsnonles))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_dermislesionalvsnonlesional <- run_glmgampoi(sce_dermislesionalvsnonlesional, design_formula = ~ patient + skin_layer + dermis_lesvsnonles, sce_subset_name = name_of_subset, reference = reference_group)
dge_dermislesionalvsnonlesional <- read.csv(paste0(glmgampoi_path, 'sce_dermislesionalvsnonlesional_patient+skin_layer+dermis_lesvsnonles_degresults.csv'), row.names = 1)
fit_dermislesionalvsnonlesional <- readRDS(paste0(glmgampoi_path, 'sce_dermislesionalvsnonlesional_patient+skin_layer+dermis_lesvsnonles_fit.RDS'))

get_dge_summary(dge_dermislesionalvsnonlesional)

# Most different genes (padj)
head(dge_dermislesionalvsnonlesional[order(dge_dermislesionalvsnonlesional$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_dermislesionalvsnonlesional, comparison_title, sce_subset = sce_dermislesionalvsnonlesional)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'dermis_lesvsnonles', 
#                            group1 = 'dermis_lesional', group2 = 'dermis_nonlesional', 
#                            dge_df = dge_dermislesionalvsnonlesional, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## GA: core vs border ####
# Covariates: patient

### 1. Sub-setting data ####
name_of_subset <- 'sce_gacorevsborder'
reference_group <- 'ga_core'
sce_ga_corevsborder <- sce[, which(sce$ga_corevsborder %in% c("ga_core", "ga_border"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_ga_corevsborder$ga_corevsborder <- droplevels(sce_ga_corevsborder$ga_corevsborder)

### 2. Exploratory analysis ####

comparison_title <- 'Granuloma annulare: core vs border'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(ga_corevsborder))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_ga_corevsborder@colData), aes(x = factor(patient), fill = factor(ga_corevsborder))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_ga_corevsborder@colData), aes(factor(sample_SPECIMEN), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_ga_corevsborder@colData), aes(factor(DISEASE), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_ga_corevsborder@colData), aes(factor(spot_type), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_ga_corevsborder@colData), aes(factor(leiden_r1.3_patient), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_ga_corevsborder@colData), aes(factor(skin_layer), fill = factor(ga_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_ga_corevsborder <- run_glmgampoi(sce_ga_corevsborder, design_formula = ~ patient + ga_corevsborder, sce_subset_name = name_of_subset, reference = reference_group)
dge_ga_corevsborder <- read.csv(paste0(glmgampoi_path, 'sce_gacorevsborder_patient+ga_corevsborder_degresults.csv'), row.names = 1)
fit_ga_corevsborder <- readRDS(paste0(glmgampoi_path, 'sce_gacorevsborder_patient+ga_corevsborder_fit.RDS'))

get_dge_summary(dge_ga_corevsborder)

# Most different genes (padj)
head(dge_ga_corevsborder[order(dge_ga_corevsborder$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_ga_corevsborder, comparison_title, sce_subset = sce_ga_corevsborder)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'ga_corevsborder', 
#                            group1 = 'ga_core', group2 = 'ga_border', 
#                            dge_df = dge_ga_corevsborder, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## SA: core vs border ####
# Covariates: patient

### 1. Sub-setting data ####
name_of_subset <- 'sce_sacorevsborder'
reference_group <- 'sa_core'
sce_sa_corevsborder <- sce[, which(sce$sa_corevsborder %in% c("sa_core", "sa_border"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_sa_corevsborder$sa_corevsborder <- droplevels(sce_sa_corevsborder$sa_corevsborder)

### 2. Exploratory analysis ####

comparison_title <- 'Sarcoidosis: core vs border'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(sa_corevsborder))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_sa_corevsborder@colData), aes(x = factor(patient), fill = factor(sa_corevsborder))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_sa_corevsborder@colData), aes(factor(sample_SPECIMEN), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_sa_corevsborder@colData), aes(factor(DISEASE), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_sa_corevsborder@colData), aes(factor(spot_type), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_sa_corevsborder@colData), aes(factor(leiden_r1.3_patient), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_sa_corevsborder@colData), aes(factor(skin_layer), fill = factor(sa_corevsborder))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_sa_corevsborder <- run_glmgampoi(sce_sa_corevsborder, design_formula = ~ patient + sa_corevsborder, sce_subset_name = name_of_subset, reference = reference_group)
dge_sa_corevsborder <- read.csv(paste0(glmgampoi_path, 'sce_sacorevsborder_patient+sa_corevsborder_degresults.csv'), row.names = 1)
fit_sa_corevsborder <- readRDS(paste0(glmgampoi_path, 'sce_sacorevsborder_patient+sa_corevsborder_fit.RDS'))

get_dge_summary(dge_sa_corevsborder)

# Most different genes (padj)
head(dge_sa_corevsborder[order(dge_sa_corevsborder$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_sa_corevsborder, comparison_title, sce_subset = sce_sa_corevsborder)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'sa_corevsborder', 
#                            group1 = 'sa_core', group2 = 'sa_border', 
#                            dge_df = dge_sa_corevsborder, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## GA: granuloma vs dermis ####
# Covariates: patient

### 1. Sub-setting data ####
name_of_subset <- 'sce_ga_gvsd'
reference_group <- 'healthydermis'
sce_ga_gvsd<- sce[, which(sce$ga_gvsd %in% c("ga_granuloma", "healthydermis"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_ga_gvsd$ga_gvsd<- droplevels(sce_ga_gvsd$ga_gvsd)

### 2. Exploratory analysis ####

comparison_title <- 'Granuloma annulare: granuloma vs healthy dermis spots'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(ga_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_ga_gvsd@colData), aes(x = factor(patient), fill = factor(ga_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_ga_gvsd@colData), aes(factor(sample_SPECIMEN), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_ga_gvsd@colData), aes(factor(DISEASE), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_ga_gvsd@colData), aes(factor(spot_type), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_ga_gvsd@colData), aes(factor(leiden_r1.3_patient), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_ga_gvsd@colData), aes(factor(skin_layer), fill = factor(ga_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_ga_gvsd<- run_glmgampoi(sce_ga_gvsd, design_formula = ~ patient + ga_gvsd, sce_subset_name = name_of_subset, reference = reference_group)
dge_ga_gvsd<- read.csv(paste0(glmgampoi_path, 'sce_ga_gvsd_patient+ga_gvsd_degresults.csv'), row.names = 1)
fit_ga_gvsd<- readRDS(paste0(glmgampoi_path, 'sce_ga_gvsd_patient+ga_gvsd_fit.RDS'))

get_dge_summary(dge_ga_gvsd)

# Most different genes (padj)
head(dge_ga_gvsd[order(dge_ga_gvsd$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_ga_gvsd, comparison_title, sce_subset = sce_ga_gvsd)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'ga_gvsd', 
#                            group1 = 'ga_granuloma', group2 = 'healthydermis', 
#                            dge_df = dge_ga_gvsd, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## SA: granuloma vs dermis ####
# Covariates: patient

### 1. Sub-setting data ####
name_of_subset <- 'sce_sa_gvsd'
reference_group <- 'healthydermis'
sce_sa_gvsd<- sce[, which(sce$sa_gvsd %in% c("sa_granuloma", "healthydermis"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_sa_gvsd$sa_gvsd<- droplevels(sce_sa_gvsd$sa_gvsd)

### 2. Exploratory analysis ####

comparison_title <- 'Sarcoidosis: granuloma vs healthy dermis spots'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(sa_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_sa_gvsd@colData), aes(x = factor(patient), fill = factor(sa_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_sa_gvsd@colData), aes(factor(sample_SPECIMEN), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_sa_gvsd@colData), aes(factor(DISEASE), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_sa_gvsd@colData), aes(factor(spot_type), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_sa_gvsd@colData), aes(factor(leiden_r1.3_patient), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_sa_gvsd@colData), aes(factor(skin_layer), fill = factor(sa_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_sa_gvsd<- run_glmgampoi(sce_sa_gvsd, design_formula = ~ patient + sa_gvsd, sce_subset_name = name_of_subset, reference = reference_group)
dge_sa_gvsd<- read.csv(paste0(glmgampoi_path, 'sce_sa_gvsd_patient+sa_gvsd_degresults.csv'), row.names = 1)
fit_sa_gvsd<- readRDS(paste0(glmgampoi_path, 'sce_sa_gvsd_patient+sa_gvsd_fit.RDS'))

get_dge_summary(dge_sa_gvsd)

# Most different genes (padj)
head(dge_sa_gvsd[order(dge_sa_gvsd$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_sa_gvsd, comparison_title, sce_subset = sce_sa_gvsd)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'sa_gvsd', 
#                            group1 = 'sa_granuloma', group2 = 'healthydermis', 
#                            dge_df = dge_sa_gvsd, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## NL: granuloma vs dermis ####
# Covariates: None

### 1. Sub-setting data ####
name_of_subset <- 'sce_nl_gvsd'
reference_group <- 'healthydermis'
sce_nl_gvsd<- sce[, which(sce$nl_gvsd %in% c("nl_granuloma", "healthydermis"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_nl_gvsd$nl_gvsd<- droplevels(sce_nl_gvsd$nl_gvsd)

### 2. Exploratory analysis ####

comparison_title <- 'Necrobiosis: granuloma vs healthy dermis spots'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(nl_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_nl_gvsd@colData), aes(x = factor(patient), fill = factor(nl_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_nl_gvsd@colData), aes(factor(sample_SPECIMEN), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_nl_gvsd@colData), aes(factor(DISEASE), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_nl_gvsd@colData), aes(factor(spot_type), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_nl_gvsd@colData), aes(factor(leiden_r1.3_patient), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_nl_gvsd@colData), aes(factor(skin_layer), fill = factor(nl_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()

### 3. Running DGE: glmGamPoi ####
#dge_nl_gvsd<- run_glmgampoi(sce_nl_gvsd, design_formula = ~ nl_gvsd, sce_subset_name = name_of_subset, reference = reference_group)
dge_nl_gvsd<- read.csv(paste0(glmgampoi_path, 'sce_nl_gvsd_nl_gvsd_degresults.csv'), row.names = 1)
fit_nl_gvsd<- readRDS(paste0(glmgampoi_path, 'sce_nl_gvsd_nl_gvsd_fit.RDS'))

get_dge_summary(dge_nl_gvsd)

# Most different genes (padj)
head(dge_nl_gvsd[order(dge_nl_gvsd$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_nl_gvsd, comparison_title, sce_subset = sce_nl_gvsd)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'nl_gvsd', 
#                            group1 = 'nl_granuloma', group2 = 'healthydermis', 
#                            dge_df = dge_nl_gvsd, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

# ## Granuloma Manual: granuloma vs dermis ####
# # Covariates: patient
# 
# ### 1. Sub-setting data ####
# name_of_subset <- 'sce_manual_gvsd'
# reference_group <- 'healthydermis'
# sce_manual_gvsd<- sce[, which(sce$manual_gvsd %in% c("granuloma", "healthydermis"))]
# # Remove empty levels because glm_gp() will complain otherwise
# sce_manual_gvsd$manual_gvsd<- droplevels(sce_manual_gvsd$manual_gvsd)
# 
# ### 2. Exploratory analysis ####
# 
# comparison_title <- 'Manual annotations: granuloma vs healthy dermis spots'
# n_col = 2
# n_row = 6
# all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(manual_gvsd))) +
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
#                        ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(x = factor(patient), fill = factor(manual_gvsd))) +
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
#                        ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          coord_flip() +
#                          ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
#                        ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(sample_SPECIMEN), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') +
#                          theme_bw() + 
#                          coord_flip() +
#                          ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
#                        ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
#                        ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(DISEASE), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
#                        ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
#                        ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(spot_type), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
#                        ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
#                        ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(leiden_r1.3_patient), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
#                        ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
#                        ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(skin_layer), fill = factor(manual_gvsd))) + 
#                          geom_bar(position = 'dodge') +
#                          geom_text(aes(label = ..count..), stat = 'count') + 
#                          theme_bw() + 
#                          ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
#                        ncol = n_col, nrow = n_row
# )
# 
# png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
# all_plots
# dev.off()
# 
# 
# ### 3. Running DGE: glmGamPoi ####
# #dge_manual_gvsd<- run_glmgampoi(sce_manual_gvsd, design_formula = ~ patient + manual_gvsd, sce_subset_name = name_of_subset, reference = reference_group)
# dge_manual_gvsd<- read.csv(paste0(glmgampoi_path, 'sce_manual_gvsd_patient+manual_gvsd_degresults.csv'), row.names = 1)
# fit_manual_gvsd<- readRDS(paste0(glmgampoi_path, 'sce_manual_gvsd_patient+manual_gvsd_fit.RDS'))
# 
# get_dge_summary(dge_manual_gvsd)
# 
# # Most different genes (padj)
# head(dge_manual_gvsd[order(dge_manual_gvsd$padj), ])
# 
# ### 4. Plotting results ####
# 
# # MA plots For now we will skip this
# # MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 
# 
# # Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
# deg_results_summary <- summarise_deg_results(dge_manual_gvsd, comparison_title, condition_reference = 'Healthy dermis')
# pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
# deg_results_summary
# dev.off()
# 
# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'manual_gvsd', 
#                            group1 = 'granuloma', group2 = 'healthydermis', 
#                            dge_df = dge_manual_gvsd, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## Granuloma Manual 2: granuloma vs dermis ####
# Covariates: patient

### 1. Sub-setting data ####
name_of_subset <- 'sce_manual_gvsd_morecovariates'
reference_group <- 'healthydermis'
sce_manual_gvsd <- sce[, which(sce$manual_gvsd %in% c("granuloma", "healthydermis"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_manual_gvsd$manual_gvsd <- droplevels(sce_manual_gvsd$manual_gvsd)

### 2. Exploratory analysis ####

comparison_title <- 'Manual annotations: granuloma vs healthy dermis spots'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(manual_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(x = factor(patient), fill = factor(manual_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(sample_SPECIMEN), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(DISEASE), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(spot_type), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(leiden_r1.3_patient), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_manual_gvsd@colData), aes(factor(skin_layer), fill = factor(manual_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_manual_gvsd <- run_glmgampoi(sce_manual_gvsd, design_formula = ~ patient + skin_layer + manual_gvsd, sce_subset_name = name_of_subset, reference = reference_group)
dge_manual_gvsd <- read.csv(paste0(glmgampoi_path, 'sce_manual_gvsd_morecovariates_patient+skin_layer+manual_gvsd_degresults.csv'), row.names = 1)
fit_manual_gvsd <- readRDS(paste0(glmgampoi_path, 'sce_manual_gvsd_morecovariates_patient+skin_layer+manual_gvsd_fit.RDS'))

get_dge_summary(dge_manual_gvsd)

# Most different genes (padj)
head(dge_manual_gvsd[order(dge_manual_gvsd$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_manual_gvsd, comparison_title, sce_subset = sce_manual_gvsd)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'manual_gvsd', 
#                            group1 = 'granuloma', group2 = 'healthydermis', 
#                            dge_df = dge_manual_gvsd, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

## Granuloma Leiden: granuloma vs dermis ####
# Covariates: patient

### 1. Sub-setting data ####
name_of_subset <- 'sce_leiden_gvsd'
reference_group <- 'healthydermis'
sce_leiden_gvsd <- sce[, which(sce$leiden_gvsd %in% c("leiden_granuloma", "healthydermis"))]
# Remove empty levels because glm_gp() will complain otherwise
sce_leiden_gvsd$leiden_gvsd <- droplevels(sce_leiden_gvsd$leiden_gvsd)

### 2. Exploratory analysis ####

comparison_title <- 'Leiden clusters: granuloma vs dermis'
n_col = 2
n_row = 6
all_plots <- ggarrange(ggplot(data = as.data.frame(sce@colData), aes(x = factor(patient), fill = factor(leiden_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across patients")), 
                       ggplot(data = as.data.frame(sce_leiden_gvsd@colData), aes(x = factor(patient), fill = factor(leiden_gvsd))) +
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across patients")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(sample_SPECIMEN), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across specimens")), 
                       ggplot(data = as.data.frame(sce_leiden_gvsd@colData), aes(factor(sample_SPECIMEN), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') +
                         theme_bw() + 
                         coord_flip() +
                         ggtitle(paste(comparison_title, "\nSubset, distribution across specimens")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(DISEASE), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across diseases")), 
                       ggplot(data = as.data.frame(sce_leiden_gvsd@colData), aes(factor(DISEASE), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across diseases")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(spot_type), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across spot types")), 
                       ggplot(data = as.data.frame(sce_leiden_gvsd@colData), aes(factor(spot_type), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset,distribution across spot types")), 
                       ggplot(data = as.data.frame(sce@colData), aes(factor(leiden_r1.3_patient), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across leiden clusters")), 
                       ggplot(data = as.data.frame(sce_leiden_gvsd@colData), aes(factor(leiden_r1.3_patient), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across leiden clusters")),
                       ggplot(data = as.data.frame(sce@colData), aes(factor(skin_layer), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nAll spots, distribution across skin layers")), 
                       ggplot(data = as.data.frame(sce_leiden_gvsd@colData), aes(factor(skin_layer), fill = factor(leiden_gvsd))) + 
                         geom_bar(position = 'dodge') +
                         geom_text(aes(label = ..count..), stat = 'count') + 
                         theme_bw() + 
                         ggtitle(paste(comparison_title, "\nSubset, distribution across skin layers")),
                       ncol = n_col, nrow = n_row
)

png(filename = paste0(exploratory_analysis_path, name_of_subset, '.png'), height = n_row * 800, width = n_col * 800)
all_plots
dev.off()


### 3. Running DGE: glmGamPoi ####
#dge_leiden_gvsd <- run_glmgampoi(sce_leiden_gvsd, design_formula = ~ patient + leiden_gvsd, sce_subset_name = name_of_subset, reference = reference_group)
dge_leiden_gvsd<- read.csv(paste0(glmgampoi_path, 'sce_leiden_gvsd_patient+leiden_gvsd_degresults.csv'), row.names = 1)
fit_leiden_gvsd<- readRDS(paste0(glmgampoi_path, 'sce_leiden_gvsd_patient+leiden_gvsd_fit.RDS'))

get_dge_summary(dge_leiden_gvsd)

# Most different genes (padj)
head(dge_leiden_gvsd[order(dge_leiden_gvsd$padj), ])

### 4. Plotting results ####

# MA plots For now we will skip this
# MA plots are commonly used to represent log fold-change versus mean expression between two treatments (Figure 4). This is visually displayed as a scatter plot with base-2 log fold-change along the y-axis and normalized mean expression along the x-axis. Data points with extreme values along the y-axis represent the genes that have highly differential expression levels (although, not necessarily differentially expressed). Typically, lower mean expression values will have more variability in log fold-change than the higher expression value. This results in a fanning effect of the data points as the graph moves from right to left. 

# Plot density distr of pvalue, padj, log2fc, summary table and volcano plot
deg_results_summary <- summarise_deg_results(dge_leiden_gvsd, comparison_title, sce_subset = sce_leiden_gvsd)
pdf(file = paste0(DEG_results_path, name_of_subset, '_deg_results_summary.pdf'), height = 15, width = 10)
deg_results_summary
dev.off()

# # Heatmap
# heatmap <- get_deg_heatmap(norm_counts = normalised_counts, comparison_column = 'leiden_gvsd', 
#                            group1 = 'leiden_granuloma', group2 = 'healthydermis', 
#                            dge_df = dge_leiden_gvsd, number_top_genes = 50)
# 
# pdf(file = paste0(DEG_results_path, name_of_subset, '_heatmap.pdf'), height = 12, width = 8)
# heatmap
# dev.off()

# 1vs1 scatterplots ####

dge_ga_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_ga_gvsd_patient+ga_gvsd_degresults.csv', sep = '/'), row.names = 1)
dge_sa_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_sa_gvsd_patient+sa_gvsd_degresults.csv', sep = '/'), row.names = 1)
dge_nl_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_nl_gvsd_nl_gvsd_degresults.csv', sep = '/'), row.names = 1)

colours_scatterplot <- list('Both' = c("#D53E4F", "#3288BD"),
                            'G.annulare' = c("#9E0142", "#66C2A5"),
                            'Sarcoidosis' = c("#FC8D59", "#99D594"),
                            'N.lipoidica' = c('#ff33cc', "#5E4FA2"))
                
get_scatterplot <- function(dge_x, dge_y, name_x, name_y){
  # For debug
  # dge_x <- dge_ga_gvsd
  # dge_y <- dge_sa_gvsd
  # name_x <- 'G.annulare'
  # name_y <- 'Sarcoidosis'
  
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
    theme(axis.title=element_text(size = 25, face = "bold"), axis.text = element_text(size = 20), legend.text = element_text(size = 13)) +
    geom_text_repel(max.overlaps = Inf, force = 0.8, box.padding = 0.5, show.legend = FALSE, size = 6) + # To show all labels
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

pdf(file = paste0(DEG_results_path, 'comparison_scatterplots_3diseases.pdf'), height = 30, width = 12)
ggarrange(
  get_scatterplot(dge_ga_gvsd, dge_sa_gvsd, 'G.annulare', 'Sarcoidosis'),
  get_scatterplot(dge_ga_gvsd, dge_nl_gvsd, 'G.annulare', 'N.lipoidica'),
  get_scatterplot(dge_nl_gvsd, dge_sa_gvsd, 'N.lipoidica', 'Sarcoidosis'),
  ncol = 1, nrow = 3, labels = c('A', 'B', 'C'), font.label = list(size = 30, color = "black")
)
dev.off()

pdf(file = paste0(DEG_results_path, 'comparison_scatterplots_3diseases_horizontal.pdf'), height = 12, width = 35)
ggarrange(
  get_scatterplot(dge_ga_gvsd, dge_sa_gvsd, 'G.annulare', 'Sarcoidosis') + theme(legend.position = 'bottom'),
  get_scatterplot(dge_ga_gvsd, dge_nl_gvsd, 'G.annulare', 'N.lipoidica') + theme(legend.position = 'bottom'),
  get_scatterplot(dge_nl_gvsd, dge_sa_gvsd, 'N.lipoidica', 'Sarcoidosis') + theme(legend.position = 'bottom'),
  ncol = 3, nrow = 1, labels = c('A', 'B', 'C'), font.label = list(size = 30, color = "black")
)
dev.off()

as_ggplot(get_legend(get_scatterplot(dge_ga_gvsd, dge_sa_gvsd, 'G.annulare', 'Sarcoidosis')))
# Manual vs leiden s plots ####
dge_manual_gvsd <- read.csv(paste0(glmgampoi_path, 'sce_manual_gvsd_morecovariates_patient+skin_layer+manual_gvsd_degresults.csv'), row.names = 1)
dge_leiden_gvsd <- read.csv(paste0(glmgampoi_path, 'sce_leiden_gvsd_patient+leiden_gvsd_degresults.csv'), row.names = 1)

# Get signed p val
dge_manual_gvsd$signedpval <- -log10(dge_manual_gvsd$pval) * sign(dge_manual_gvsd$log2fc)
dge_leiden_gvsd$signedpval <- -log10(dge_leiden_gvsd$pval) * sign(dge_leiden_gvsd$log2fc)

# Merge both
dge_gvsd <- merge(dge_manual_gvsd, dge_leiden_gvsd, by = 'gene_symbol', suffixes = c('_manual', '_leiden'))
head(dge_gvsd)

# Sort by signedpval_manual and plot Manual 
dge_gvsd <- dge_gvsd[order(dge_gvsd$signedpval_manual),]
dge_gvsd$signedpval_order <- 1:nrow(dge_gvsd)

manual_scatterplot <- ggplot(data = dge_gvsd, aes(x = signedpval_manual, y = signedpval_order)) +
  geom_point() + 
  coord_cartesian(ylim = c(0, 15000), clip = 'off') + 
  scale_y_continuous(breaks = seq(0,15000, 5000)) +
  theme_classic() +
  ggtitle('Manual annotations') +
  labs(x = 'signed(p-value)', y = 'Rankings') #+

leiden_scatterplot <- ggplot(data = dge_gvsd, aes(x = signedpval_leiden, y = signedpval_order)) +
  geom_point() + 
  coord_cartesian(ylim = c(0, 15000), clip = 'off') + 
  scale_y_continuous(breaks = seq(0,15000, 5000)) +
  theme_classic() +
  ggtitle('Leiden clustering') +
  labs(x = 'signed(p-value)', y = 'Rankings')

pdf(file = paste0(DEG_results_path, 'manualvsleiden_rankings_pval_scatterplots.pdf'), height = 8, width = 14)
ggarrange(
  manual_scatterplot, leiden_scatterplot,
  ncol = 2, nrow = 1
)
dev.off()
