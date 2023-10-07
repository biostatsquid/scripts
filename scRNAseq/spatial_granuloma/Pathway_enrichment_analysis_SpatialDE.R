# Name: Pathway_enrichment_analysis_SpatialDE.R
# Author: Laura Twomey
# Date of creation: 02 June 2022
# Pathway enrichment analysis for SpatialDE results

# Clean Environment ----------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects including hidden objects
gc() #free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F) #avoid truncated output in R console and scientific notation

# Set working directory ----------------------------------------------------
path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/PEA"
setwd(path)
wilcox_path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Intensity_ratios"
clusterp_path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/PEA/clusterProfiler/"
mgsa_path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/PEA/MGSA/"
  
# Libraries and functions ------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(anndata)
library(reticulate)
library(hdf5r)
library(pheatmap)
library(mgsa)
library(fgsea)
library(clusterProfiler)
library(umap)
library(ggpubr) # for ggarrange

## Function: Patterns to python format ####
to_python_format <- function(list_of_patterns){
  list_in_python_format <- cat(paste(shQuote(gsub('\\.', '_pattern', list_of_patterns), type="cmd"), collapse=", "))
  return(list_in_python_format)
}

to_python_format2 <- function(list_of_patterns){
  list_in_python_format <- cat(paste(shQuote(gsub('\\.', '-', list_of_patterns), type="cmd"), collapse=", "))
  return(list_in_python_format)
}

## Function: get_clusterp_analysis_spatialde ####
# This function runs clusterp analysis for a set of background genes (for now KEGG, Reactome or GO), given a list of spatialde_results from 
# the spatialDE pipeline (list of genes are in $genes column). Pathways can be filtered by padj and gene count.
# If get_heatmap is set to TRUE (default) it also calls the heatmap from clusterp_get_heatmap_spatialde

get_clusterp_analysis_spatialde <- function(background_genes, spatialde_results_list, padj_cutoff = 0.05, genecount_cutoff = 5, filename = 'example', get_heatmap = TRUE, number_of_pathways_cutoff = 25){
  # for debug
  # background_genes <- 'KEGG'
  # spatialde_results_list <- cluster3_patterns
  # padj_cutoff <- 0.05
  # genecount_cutoff <- 5
  
  # Load the background genes
  if(background_genes == 'KEGG'){
    pwl2 <- readRDS('Background_genes/kegg.RDS')
    background_genes <- '_kegg'
  } else if(background_genes == 'Reactome'){
    pwl2 <- readRDS('Background_genes/reactome.RDS')
    background_genes <- '_reactome'
  } else if(background_genes == 'GO'){
    pwl2 <- readRDS('Background_genes/go.bp.RDS')
    background_genes <- '_gobp'
  } else {
    stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
  }
  
  # Run clusterProfiler on each sub-dataframe
  clusterp_all <- lapply(names(spatialde_results_list),
                         function(x) enricher(gene = spatialde_results_list[[x]]$genes,
                                              TERM2GENE = pwl2))
  names(clusterp_all) <- names(spatialde_results_list)
  
  #Convert the clusterp_all list of enrichResults for each sample_pattern to a dataframe with the pathways
  clusterp_all_df<- lapply(names(clusterp_all), function(x) rbind(clusterp_all[[x]]@result))
  names(clusterp_all_df) <- names(clusterp_all)
  clusterp_all_df = do.call(rbind, clusterp_all_df)
  #head(clusterp_all_df)
  colnames(clusterp_all_df) <- gsub('ID', 'pathways', colnames(clusterp_all_df)) # rename the pathways column
  clusterp_all_df <- clusterp_all_df %>% mutate(minuslog10padj = -log10(p.adjust),
                                                sample_SPECIMEN_pattern = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(clusterp_all_df)))
  rownames(clusterp_all_df) <- NULL
  
  # Subset to those pathways that have p adj < cutoff and gene count > cutoff
  target_pws <- unique(clusterp_all_df$pathways[clusterp_all_df$p.adjust < padj_cutoff & clusterp_all_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
  clusterp_all_df <- clusterp_all_df[clusterp_all_df$pathways %in% target_pws, ]
  
  print('Saving clusterprofiler results')
  write.csv(clusterp_all_df, paste0(clusterp_path, filename, background_genes, '.csv'), row.names = FALSE)
  
  # # Check distribution of minus log 10 p adj values
  # print('Saving minuslog(padj) distribution plot')
  # png(filename = paste0(path, 'clusterProfiler/', filename, background_genes, '_minuslog10padj_distribution.png'), height = 600, width = 600)
  # ggplot(data = clusterp_all_df, aes(x = minuslog10padj)) +
  #   geom_histogram(col = 'darkgreen', fill = 'aquamarine', bins = 100) + 
  #   scale_x_continuous(breaks = seq(0, 155, 10)) + 
  #   theme_bw() +
  #   ggtitle(paste(filename, "Distribution of -log10(padj) values")) 
  # dev.off()
  
  if(get_heatmap == TRUE){
    
    # Get the cluster number
    clusternumber <- gsub('cluster', '', filename)
    
    # Set the cluster in the annotations dataframe 
    col_annotations_df_deg$`K-means clustering` <- clusternumber
      
    clusterp_get_heatmap_spatialde(dataframe = clusterp_all_df, file_name = paste0(filename, background_genes, '_clusterp_heatmap'), pathway_cutoff = number_of_pathways_cutoff, annotations_df = col_annotations_df_deg)
  }
  return(clusterp_all_df)
}

## Function: Get heatmap_df ####
clusterp_get_heatmap_spatialde <- function(dataframe, file_name = 'filename',  pathway_cutoff = number_of_pathways_cutoff, annotations_df){
  # for debug
  # dataframe <- clusterp_all_df
  # pathway_cutoff <- 50
  if(!is.data.frame(dataframe)){
    stop('results_filtered must be a dataframe')
  }
  if(!is.character(file_name)){
    stop('file_name must be a character')
  }
  
  # First convert the dataframe into the heatmap dataframe we need to build the heatmap
  heatmap_df <- tidyr::pivot_wider(dataframe[,c('pathways', 'minuslog10padj', 'sample_SPECIMEN_pattern')], values_from = minuslog10padj, names_from = sample_SPECIMEN_pattern, values_fill = 0)
  heatmap_df <- as.data.frame(heatmap_df)
  rownames(heatmap_df) <- heatmap_df$pathways # to be able to set them as row names
  heatmap_df <- heatmap_df[!grepl('pathways', colnames(heatmap_df))] # now remove the pathways column
  head(heatmap_df)
  
  if(nrow(heatmap_df) > pathway_cutoff){
    selected_pathways <- head(unique(dataframe[order(-dataframe$minuslog10padj),'pathways']), pathway_cutoff) # order by minuslogpadj and get ones with larger minuslogpadj
    heatmap_df <- subset(heatmap_df, rownames(heatmap_df) %in% selected_pathways) # subset to less pathways for visualisation purposes if there are too many
  }

  #write.csv(heatmap_df, paste(path, 'clusterProfiler/heatmap_df.csv', sep = '/'), row.names = TRUE)
  
  # Make the rownames for the heatmap a bit nicer
  rownames(heatmap_df) <- gsub('GOBP_', '', gsub('KEGG_', '', gsub('REACTOME_', '', rownames(heatmap_df))))
  rownames(heatmap_df) <- gsub('(H|h)iv', 'HIV', 
                               gsub('pd 1', 'PD-1',
                                    gsub('ecm', 'ECM', 
                                         gsub('(I|i)nterleukin', 'IL', 
                                              gsub('(R|r)na', 'RNA', 
                                                   gsub('(D|d)na', 'DNA',
                                                        gsub(' i ', ' I ', 
                                                             gsub('(N|n)adh ', 'NADH ', 
                                                                  gsub('(N|n)ad ', 'NAD ',
                                                                       gsub('t cell', 'T cell',
                                                                            gsub('b cell', 'B cell',
                                                                                 gsub('built from .*', ' (...)',
                                                                                      gsub('mhc', 'MHC',
                                                                                           gsub('mhc class i', 'MHC I', 
                                                                                                gsub('mhc class ii', 'MHC II', 
                                                                                                     stringr::str_to_sentence(gsub('_', ' ', rownames(heatmap_df))))))))))))))))))
  
  # Remove patient names from column names
  colnames(heatmap_df) <- gsub('P17851_1001', 'P1_1001', 
                               gsub('P17851_1002', 'P1_1002',
                                    gsub('P17851_1003', 'P2_1003',
                                         gsub('P17851_1004', 'P2_1004',
                                              gsub('P18554_1001', 'P3_1001', 
                                                   gsub('P18554_1002', 'P3_1002', 
                                                        gsub('P18554_1003', 'P4_1003', 
                                                             gsub('P18554_1004', 'P4_1004',         
                                                                  gsub('P18554_1005', 'P1_1005', 
                                                                       gsub('P18554_1006', 'P1_1006', 
                                                                            gsub('P18554_1007', 'P1_1007', 
                                                                                 gsub('P18554_1008', 'P1_1008', colnames(heatmap_df)))))))))))))
    
  ## Get heatmap 
  heatmap <- pheatmap(heatmap_df,
                      show_colnames = T, show_rownames = T,
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      annotation_col = annotations_df,
                      fontsize_col = 10,
                      fontsize_row = 12,
                      annotation_colors = ann_colors
  )
  
  height_of_heatmap <- nrow(heatmap_df) * 0.5 # to make the heatmap longer in case there are many rows
  width_of_heatmap <- ncol(heatmap_df) * 0.8 # to make the heatmap longer in case there are many columns
  
  pdf(file = paste(clusterp_path, file_name, '.pdf', sep = ''), width = width_of_heatmap, height = height_of_heatmap)
  grid::grid.draw(heatmap$gtable)
  dev.off()
  
  return(heatmap)
}

## Function: get_mgsa_analysis_spatialde ####
# This function runs mgsa analysis for a set of background genes (for now KEGG, Reactome or GO), given a list of spatialde_results from 
# the spatialDE pipeline (list of genes are in $genes column). Pathways can be filtered by padj and gene count.

get_mgsa_analysis_spatialde <- function(background_genes, spatialde_results_list, filename = 'example', get_heatmap = TRUE, number_of_pathways_cutoff = 25){
  # for debug
  # background_genes <- 'KEGG'
  #spatialde_results_list <- cluster3_patterns
  
  # Load the background genes
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

  ##  Run mgsa for subset of genes in spatialde_results_list
  mgsa_fit_list <- lapply(names(spatialde_results_list), function(x) mgsa(spatialde_results_list[[x]]$genes, pwl1))
  names(mgsa_fit_list) <- names(spatialde_results_list)
  saveRDS(mgsa_fit_list, file = paste0(mgsa_path, 'Fit_results/', paste0(filename, background_genes, '_mgsa_fit_list.RData')))
  print('MGSA fit finished!')
  
  ## Get the enriched pathways 
  res_list <- lapply(names(mgsa_fit_list), function(x) setsResults(mgsa_fit_list[[x]]))
  names(res_list) <- names(spatialde_results_list)
  saveRDS(res_list, file = paste0(mgsa_path, 'Fit_results/', paste0(filename, background_genes, '_mgsa_res_list.RData'))) 
  
  res_list <- Map(cbind, res_list, sample_SPECIMEN_pattern = names(res_list))
  
  res_all_df = do.call(rbind, res_list)
  res_all_df$pathways <- gsub('^.*GOBP_|^.*KEGG_|^.*REACTOME_', '', rownames(res_all_df))
  write.csv(res_all_df, paste0(mgsa_path, 'Fit_results/', paste0(filename, background_genes, 'MGSA_res_all_.csv')))
  
  target_pathways <- unique(res_all_df$pathways[res_all_df$estimate > 0.5]) # target pathways have an estimate higher than 0.5
  res_all_df_filtered <- res_all_df[res_all_df$pathways %in% target_pathways, ]
  write.csv(res_all_df_filtered, paste0(mgsa_path, paste0(filename, background_genes, '_mgsa_res_filtered.csv')))
  
  if(get_heatmap == TRUE){
    # Get the cluster number
    clusternumber <- gsub('cluster', '', filename)
    
    # Set the cluster in the annotations dataframe 
    col_annotations_df_deg$`K-means clustering` <- clusternumber
    
    mgsa_get_heatmap_spatialde(dataframe = res_all_df_filtered, file_name = paste0(filename, background_genes, '_mgsa_heatmap'),  pathway_cutoff = number_of_pathways_cutoff, annotations_df = col_annotations_df_deg)
  }
  
  return(res_all_df_filtered)
}

## Function: Get mgsa_get_heatmap_spatialde ####
mgsa_get_heatmap_spatialde <- function(dataframe, file_name = 'filename', pathway_cutoff = number_of_pathways_cutoff, annotations_df){
  # for debug
  # dataframe <- res_all_df_filtered 
  if(!is.data.frame(dataframe)){
    stop('results_filtered must be a dataframe')
  }
  if(!is.character(file_name)){
    stop('file_name must be a character')
  }
  
  # First convert the dataframe into the heatmap dataframe we need to build the heatmap
  heatmap_df <- tidyr::pivot_wider(dataframe[,c('pathways', 'estimate', 'sample_SPECIMEN_pattern')], values_from = estimate, names_from = sample_SPECIMEN_pattern, values_fill = 0)
  heatmap_df <- as.data.frame(heatmap_df)
  rownames(heatmap_df) <- heatmap_df$pathways # to be able to set them as row names
  heatmap_df <- heatmap_df[!grepl('pathways', colnames(heatmap_df))] # now remove the pathways column
  head(heatmap_df)
  #write.csv(heatmap_df, paste(path, 'clusterProfiler/heatmap_df.csv', sep = '/'), row.names = TRUE)

  if(nrow(heatmap_df) > pathway_cutoff){
    selected_pathways <- head(unique(dataframe[order(-dataframe$estimate),'pathways']), pathway_cutoff) # order by minuslogpadj and get ones with larger minuslogpadj
    heatmap_df <- subset(heatmap_df, rownames(heatmap_df) %in% selected_pathways) # subset to less pathways for visualisation purposes if there are too many
  }
  
  # Make the rownames for the heatmap a bit nicer
  rownames(heatmap_df) <- gsub('GOBP_', '', gsub('KEGG_', '', gsub('REACTOME_', '', rownames(heatmap_df))))
  rownames(heatmap_df) <- gsub('(H|h)iv', 'HIV', 
                                    gsub('pd 1', 'PD-1',
                                         gsub('ecm', 'ECM', 
                                              gsub('(I|i)nterleukin', 'IL', 
                                                   gsub('(R|r|)na', 'RNA', 
                                                        gsub('(D|d|)na', 'DNA',
                                                             gsub(' i ', ' I ', 
                                                                  gsub('(N|n)adh ', 'NADH ', 
                                                                       gsub('(N|n)ad ', 'NAD ',
                                                                            gsub('t cell', 'T cell',
                                                                                 gsub('b cell', 'B cell',
                                                                                      gsub('built from .*', ' (...)',
                                                                                           gsub('mhc', 'MHC',
                                                                                                gsub('mhc class i', 'MHC I', 
                                                                                                     gsub('mhc class ii', 'MHC II', 
                                                                                                          stringr::str_to_sentence(gsub('_', ' ', rownames(heatmap_df))))))))))))))))))
  
  # Remove patient names from column names
  colnames(heatmap_df) <- gsub('P17851_1001', 'P1_1001', 
                               gsub('P17851_1002', 'P1_1002',
                                    gsub('P17851_1003', 'P2_1003',
                                         gsub('P17851_1004', 'P2_1004',
                                              gsub('P18554_1001', 'P3_1001', 
                                                   gsub('P18554_1002', 'P3_1002', 
                                                        gsub('P18554_1003', 'P4_1003', 
                                                             gsub('P18554_1004', 'P4_1004',         
                                                                  gsub('P18554_1005', 'P1_1005', 
                                                                       gsub('P18554_1006', 'P1_1006', 
                                                                            gsub('P18554_1007', 'P1_1007', 
                                                                                 gsub('P18554_1008', 'P1_1008', colnames(heatmap_df)))))))))))))
  heatmap <- pheatmap(heatmap_df,
                      show_colnames = T, show_rownames = T,
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      annotation_col = annotations_df,
                      fontsize_col = 10,
                      fontsize_row = 12,
                      annotation_colors = ann_colors
  )
  
  height_of_heatmap <- nrow(heatmap_df) * 0.5 # to make the heatmap longer in case there are many rows
  width_of_heatmap <- ncol(heatmap_df) * 0.8 # to make the heatmap longer in case there are many columns
  
  pdf(file = paste(mgsa_path, file_name, '.pdf', sep = ''), width = width_of_heatmap, height = height_of_heatmap)
  grid::grid.draw(heatmap$gtable)
  dev.off()
  
  return(heatmap)
}


# Read in files  ------------------------------------------------

# The frequency table has genes x sample_SPECIMENs and 1 or 0 depending if the gene is present or not
freqtable_all <- readRDS(file = paste(path, 'freqtable_all.csv', sep = '/'))

## Reading in Wilcoxon results, annotations dataframes and merging into one  
# # Read in Wilcoxon results
# wilcox_df <- read.csv(file = paste(wilcox_path, 'wilcox_df.csv', sep = '/'), row.names = 1)
# wilcox_df$sample_SPECIMEN_pattern <- gsub('-', '.', wilcox_df$sample_SPECIMEN_pattern)
# head(wilcox_df)
# # Calculate the FDR and the signedFDR from the p_values
# wilcox_df <- wilcox_df %>%
#   group_by(annotation) %>%
#   mutate(FDR = p.adjust(pvalue, method = 'BH'), signedFDR = -log10(FDR) * sign(statistic))
# head(wilcox_df)
# 
# wilcox_df$annotation <- factor(wilcox_df$annotation,
#                                levels = c("manual_granuloma", "leiden_granuloma",
#                                           "leiden_border_granuloma", "leiden_core_granuloma",
#                                           "leiden_epidermis", 'epidermis_interface'),
#                                ordered = T)
# # Mark significant pvalues (FDR<0.05 and statistic > 0)
# wilcox_df$Significant <- as.numeric(wilcox_df$FDR < 0.05 & wilcox_df$statistic > 0)
# colnames(wilcox_df)
# wilcox_plot_df <- tidyr::pivot_wider(wilcox_df, id_cols = sample_SPECIMEN_pattern, names_sep = '.', names_from = annotation,
#                               values_from = c(statistic, pvalue, FDR, signedFDR, Significant))
# head(wilcox_plot_df)
# 
# # Get the patient/sample dataframe
# ann_df <- readRDS(file = paste(path, 'ann_df.csv', sep = '/'))
# ann_df <- unique(ann_df)
# ann_df$sample_SPECIMEN <- paste(ann_df$sample, '_', ann_df$patient, '.', ann_df$SAMPLE, sep = '')
# ann_df <- ann_df %>% dplyr::slice(rep(1:n(), 8)) # we need to tweak it a bit to create 'fake' sample_specimen_patterns so we can match
# ann_df$pattern <- c(rep(0, 43), rep(1, 43), rep(2, 43), rep(3, 43), rep(4, 43), rep(5, 43), rep(6, 43), rep(7, 43))
# ann_df$sample_SPECIMEN_pattern <- paste(ann_df$sample_SPECIMEN, '_pattern', ann_df$pattern, sep = '')
# rownames(ann_df) <- ann_df$sample_SPECIMEN_pattern
# head(ann_df)
# 
# # Merge both and create an annotations dataframe
# col_annotations_df <- ann_df[grepl(paste(colnames(freqtable_all), collapse = '|'), ann_df$sample_SPECIMEN_pattern), c('patient', 'sample_SPECIMEN_pattern', 'DISEASE')]
# col_annotations_df <- merge(col_annotations_df, wilcox_plot_df, by = 'sample_SPECIMEN_pattern')
# rownames(col_annotations_df) <- col_annotations_df$sample_SPECIMEN_pattern
# # Save it
# saveRDS(col_annotations_df, 'wilcox_results_with_annotations.RDS')


# Set annotations ------------------------------------------------
# col_annotations_df <- readRDS('wilcox_results_with_annotations.RDS')
# col_annotations_df <- col_annotations_df[colnames(freqtable_all), ]
# ## KMeans clustering #####
# set.seed(1234)
# for (x in 3:10){
#   kmeans <- kmeans(t(freqtable_all), centers = x)
#   col_annotations_df[, paste0("KMeans", x)] <- as.factor(kmeans$cluster)
# }
#col_annotations_df$cluster <- kmeans$cluster

## UMAP #####
# this implementation automatically preserves the seed for reproducibility
# custom.config = umap.defaults
# custom.config$random_state <- 42
# umap.defaults$random_state <- 42
# umap <- umap::umap(t(freqtable_all), config = umap.defaults)
# colnames(umap$layout) <- c("UMAP1", "UMAP2")
# col_annotations_df <- cbind(col_annotations_df, umap$layout)

# ## PCA
# pca <- prcomp(t(freqtable_all))
# col_annotations_df <- data.frame(col_annotations_df, pca$x[, 1:10])

# # Since UMAPs were created by Ines and the seed doesn't seem to work to reproduce them, we will just copy paste those columns onto our dataframe
# anno_cols2 <- readRDS('anno_cols2.RDS') # this one also has PCA and UMAPs information (Ines did the calculations in a separate R script)
# anno_cols2$sample_SPECIMEN_pattern <- row.names(anno_cols2)
# col_annotations_df <- merge(col_annotations_df, anno_cols2[,grepl('PC|UMAP|sample_SPECIMEN_pattern', colnames(anno_cols2))], by = 'sample_SPECIMEN_pattern')

## anno_cols2 from Ines ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Another option is to read in the UMAPs, KMeans, PCAs that Ines did in a separate script (spatial_granuloma_visualizations.R) 
# anno_cols2 <- readRDS('anno_cols2.RDS') # this one also has PCA and UMAPs information (Ines did the calculations in a separate R script)
# anno_cols2$sample_SPECIMEN_pattern <- row.names(anno_cols2)
# col_annotations_df <- merge(col_annotations_df, anno_cols2[,grepl('KMeans|PC|UMAP|sample_SPECIMEN_pattern', colnames(anno_cols2))], by = 'sample_SPECIMEN_pattern')
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#col_annotations_df$Disease <- stringr::str_to_sentence(col_annotations_df$DISEASE) # First word to capital letter
#rownames(col_annotations_df) <- col_annotations_df$sample_SPECIMEN_pattern

#write.csv(col_annotations_df, 'col_annotations_df.csv') # Save df
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

# Annotations for heatmap if DEG results (only Patient, Disease and Cluster) with rownames changed to remove patient ID
col_annotations_df_deg <- col_annotations_df
rownames(col_annotations_df_deg) <- gsub('P17851_1001', 'P1_1001', 
                                         gsub('P17851_1002', 'P1_1002',
                                              gsub('P17851_1003', 'P2_1003',
                                                   gsub('P17851_1004', 'P2_1004',
                                                        gsub('P18554_1001', 'P3_1001', 
                                                             gsub('P18554_1002', 'P3_1002', 
                                                                  gsub('P18554_1003', 'P4_1003', 
                                                                       gsub('P18554_1004', 'P4_1004',         
                                                                            gsub('P18554_1005', 'P1_1005', 
                                                                                 gsub('P18554_1006', 'P1_1006', 
                                                                                      gsub('P18554_1007', 'P1_1007', 
                                                                                           gsub('P18554_1008', 'P1_1008', rownames(col_annotations_df_deg)))))))))))))
col_annotations_df_deg <- col_annotations_df_deg[,c('patient', 'DISEASE')]
col_annotations_df_deg$patient <- factor(col_annotations_df_deg$patient)
names(col_annotations_df_deg) <- c('Patient', 'Disease')
col_annotations_df_deg$Disease <- stringr::str_to_sentence(col_annotations_df_deg$Disease)
col_annotations_df_deg$`K-means clustering` <- '1'

# Set colours  ------------------------------------------------

# Set patients for subsetting freqtable
ann_colors = list(
  patient = c("91253" = "#F46D43",#'red',#granuloma annulare
              "45703" = "#D53E4F",#'indianred',#granuloma annulare
              "50107" = "#9E0142",#'firebrick',#granuloma annulare
              "95096" = "orchid", # necrobiosis lipoidica
              "82301" = "mediumvioletred", #sarcoidosis-cutaneous
              "72859" = "blueviolet"), #sarcoidosis-suspected
  Patient = c("91253" = "#F46D43",#'red',#granuloma annulare
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
  leiden_core_granuloma = rev(brewer.pal(11, "RdBu")),
  leiden_border_granuloma = rev(brewer.pal(11, "RdBu")),
  leiden_epidermis = rev(brewer.pal(11, "RdBu")),
  leiden_granuloma = rev(brewer.pal(11, "RdBu")),
  manual_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_core_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_border_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_epidermis = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.manual_granuloma = rev(brewer.pal(11, "RdBu")),
  `Core granuloma (Leiden)` = rev(brewer.pal(11, "RdBu")),
  `Border granuloma (Leiden)` = rev(brewer.pal(11, "RdBu")),
  `Epidermis (Leiden)` = rev(brewer.pal(11, "RdBu")),
  `Epidermis (annotation)` = rev(brewer.pal(11, "RdBu")),
  `Granuloma (annotation)` = rev(brewer.pal(11, "RdBu")),
  KMeans5 = c("1"="#88CFA4", "2"="#F88D52", "3"="skyblue",
              "4"="#9E0142", "5"="#5E4FA2"),
  `K-means clustering` = c("1"="#88CFA4", "2"="#F88D52", "3"="skyblue",
                           "4"="#9E0142", "5"="#5E4FA2")
  #sample_SPECIMEN_SAMPLE = sample(color, 20)
)

BuRd <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
GrBu <- colorRampPalette(colors = c(rev(brewer.pal(11, "PiYG"))[1:4],
                                    brewer.pal(11, "RdYlBu")[c(6, 8:11)]))

# Heatmaps genes ------------------------------------------------
## All  ------------------------------------------------
heatmap <- pheatmap(freqtable_all,
                    show_colnames = F, show_rownames = F,
                    clustering_method = 'ward.D',
                    annotation_col = heatmap_ann,
                    treeheight_row = 0,
                    treeheight_col = 10,
                    annotation_colors = ann_colors)

# pdf(width = 10, height = 14,
#     file = "Figures/Figure3_heat_spatial_patterns.pdf")
# heatmap
# dev.off()

## Only granuloma (clusters 2, 4)  ------------------------------------------------
# Heatmap of only granuloma patterns from clusters 2 and 4
clusters_2_4_patterns <- col_annotations_df[col_annotations_df$KMeans5 == 2|col_annotations_df$KMeans5 == 4, 'sample_SPECIMEN_pattern']
freqtable_g <- freqtable_all[, clusters_2_4_patterns] # only keep patterns belonging to K means cluster 1 and 2
freqtable_g <- freqtable_g[rowSums(freqtable_g[])>0,] # remove rows with only 0s
heatmap_ann_g <- heatmap_ann[clusters_2_4_patterns,]
heatmap_g <- pheatmap(freqtable_g,
                      show_colnames = F, show_rownames = F,
                      clustering_method = 'ward.D',
                      annotation_col = heatmap_ann_g,
                      treeheight_row = 0,
                      treeheight_col = 10,
                      annotation_colors = ann_colors
)

# pdf(width = 10, height = 14,
#     file = "Figures/onlygranulomapatterns_heatmap.pdf")
# heatmap
# dev.off()

# UMAPS: K-Means clustering ------------------------------------------------
# Getting patterns from a particular cluster:
to_python_format2(col_annotations_df$sample_SPECIMEN_pattern[col_annotations_df$KMeans5 == '3'])

kmeans <-  ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = KMeans5)) +
  geom_point() +
  scale_color_manual("K-means (5)",
                     values = c("1"="#88CFA4", "2"="#F88D52", "3"="dodgerblue",
                                "4"="#9E0142", "5"="#5E4FA2"),
                     labels = c("1", "Granuloma-2", "Interface", "Granuloma-1", "Epidermis")) +
  theme_bw() +
  ggtitle("K-Means clustering (K=5)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size=rel(1.06)),
    legend.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) 

disease <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = Disease)) +
  geom_point() +
  scale_color_manual("Disease",
                     values = ann_colors$Disease,
                     labels = c("Granuloma\nannulare", "Necrobiosis\nlipoidica", "Sarcoidosis")) +
  theme_bw() +
  ggtitle("Disease") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size=rel(1.06)),
    legend.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) 

manualg <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.manual_granuloma)) +
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = BuRd(100)) +
  #scale_colour_gradientn(colours = BuRd(100), limits=c(-100, 100)) +
  ggtitle("Granuloma (annotation)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

epidermis <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.epidermis_interface)) +
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = GrBu(100)) +
  #scale_colour_gradientn(colours = GrBu(100), limits=c(-100, 100)) +
  ggtitle("Epidermis (annotation)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

core <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.leiden_core_granuloma)) + 
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = BuRd(100)) +
  #scale_colour_gradientn(colours = BuRd(100), limits=c(-100, 100)) +
  ggtitle("Granuloma core (Leiden cluster)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

border <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.leiden_border_granuloma)) + 
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = BuRd(100)) +
  #scale_colour_gradientn(colours = BuRd(100), limits=c(-100, 100)) +
  ggtitle("Granuloma border (Leiden cluster)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

all_UMAPs <- ggarrange(
  ggarrange(kmeans + theme(legend.position = "none"),
            as_ggplot(get_legend(kmeans)),
            nrow = 1, ncol = 2,
            widths = c(1, 0.5)),
  ggarrange(disease + theme(legend.position = "none"),
            as_ggplot(get_legend(disease)),
            nrow = 1, ncol = 2,
            widths = c(1, 0.5)),
  ggarrange(manualg + theme(legend.position = "none"),
            as_ggplot(get_legend(manualg)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ggarrange(epidermis + theme(legend.position = "none"),
            as_ggplot(get_legend(epidermis)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ggarrange(core + theme(legend.position = "none"),
            as_ggplot(get_legend(core)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ggarrange(border + theme(legend.position = "none"),
            as_ggplot(get_legend(border)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ncol = 2, nrow = 3) 

# pdf(width = 10, height = 12,
#     file = "Figures/Figure3_UMAP_spatial_patterns.pdf")
# all_UMAPs
# dev.off()

# Background genes ------------------------------------------------
# We will use the same as for normal DEG - see Pathway_enrichment_analysis_DEG.R

# All pattern genes ------------------------------------------------
# Get a list of sample_SPECIMEN_patterns, each element is a dataframe with a column 'genes' belonging to that pattern

# Load in genes for each pattern & sample_SPECIMEN 
patterns_genes_df <- read.csv(paste(path, 'all_df.csv', sep = '/'), row.names = 1)
patterns_genes_df <- tidyr::separate(patterns_genes_df, slide_pattern, into = c('sample_SPECIMEN', 'pattern'), sep = '_pattern', remove = FALSE)
head(patterns_genes_df)

# Split the dataframe into a list of sub-dataframes per sample_SPECIMEN and pattern
patterns_genes_df_list <- split(patterns_genes_df, list(patterns_genes_df$sample_SPECIMEN, patterns_genes_df$pattern))
names(patterns_genes_df_list) <- gsub('-', '.', gsub('\\.', '_pattern', names(patterns_genes_df_list))) # rename so sample_SPECIMEN_pattern names match R version

# PEA ------------------------------------------------

## Cluster 1 ------------------------------------------------
# ### 1. Getting patterns ####
# cluster1 <- col_annotations_df$sample_SPECIMEN_pattern[col_annotations_df$KMeans5 == '1']
# cluster1_patterns <- patterns_genes_df_list[grepl(paste(cluster1, collapse = '|'), names(patterns_genes_df_list))]
# name_of_subset <- 'cluster1'
# 
# ### 2. clusterProfiler ####
# kegg_results_clusterp <- get_clusterp_analysis_spatialde('KEGG', cluster1_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
# reactome_results_clusterp <- get_clusterp_analysis_spatialde('Reactome', cluster1_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
# gobp_results_clusterp <- get_clusterp_analysis_spatialde('GO', cluster1_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
# 
# ### 3. MGSA ####
# kegg_results_mgsa <- get_mgsa_analysis_spatialde('KEGG', cluster1_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
# reactome_results_mgsa <- get_mgsa_analysis_spatialde('Reactome', cluster1_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
# gobp_results_mgsa <- get_mgsa_analysis_spatialde('GO', cluster1_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)

## Cluster 2 ------------------------------------------------
### 1. Getting patterns ####
cluster2 <- col_annotations_df$sample_SPECIMEN_pattern[col_annotations_df$KMeans5 == '2']
cluster2_patterns <- patterns_genes_df_list[grepl(paste(cluster2, collapse = '|'), names(patterns_genes_df_list))]
name_of_subset <- 'cluster2'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_spatialde('KEGG', cluster2_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_clusterp <- get_clusterp_analysis_spatialde('Reactome', cluster2_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_clusterp <- get_clusterp_analysis_spatialde('GO', cluster2_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)

### 3. MGSA ####
kegg_results_mgsa <- get_mgsa_analysis_spatialde('KEGG', cluster2_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_mgsa <- get_mgsa_analysis_spatialde('Reactome', cluster2_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_mgsa <- get_mgsa_analysis_spatialde('GO', cluster2_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)

## Cluster 3 ------------------------------------------------
### 1. Getting patterns ####
cluster3 <- col_annotations_df$sample_SPECIMEN_pattern[col_annotations_df$KMeans5 == '3']
cluster3_patterns <- patterns_genes_df_list[grepl(paste(cluster3, collapse = '|'), names(patterns_genes_df_list))]
name_of_subset <- 'cluster3'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_spatialde('KEGG', cluster3_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_clusterp <- get_clusterp_analysis_spatialde('Reactome', cluster3_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_clusterp <- get_clusterp_analysis_spatialde('GO', cluster3_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)

### 3. MGSA ####
kegg_results_mgsa <- get_mgsa_analysis_spatialde('KEGG', cluster3_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_mgsa <- get_mgsa_analysis_spatialde('Reactome', cluster3_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_mgsa <- get_mgsa_analysis_spatialde('GO', cluster3_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)

# Cluster 4 ------------------------------------------------
### 1. Getting patterns ####
cluster4 <- col_annotations_df$sample_SPECIMEN_pattern[col_annotations_df$KMeans5 == '4']
cluster4_patterns <- patterns_genes_df_list[grepl(paste(cluster4, collapse = '|'), names(patterns_genes_df_list))]
name_of_subset <- 'cluster4'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_spatialde('KEGG', cluster4_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_clusterp <- get_clusterp_analysis_spatialde('Reactome', cluster4_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_clusterp <- get_clusterp_analysis_spatialde('GO', cluster4_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)

### 3. MGSA ####
kegg_results_mgsa <- get_mgsa_analysis_spatialde('KEGG', cluster4_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_mgsa <- get_mgsa_analysis_spatialde('Reactome', cluster4_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_mgsa <- get_mgsa_analysis_spatialde('GO', cluster4_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)

## Cluster 5 ------------------------------------------------
### 1. Getting patterns ####
cluster5 <- col_annotations_df$sample_SPECIMEN_pattern[col_annotations_df$KMeans5 == '5']
cluster5_patterns <- patterns_genes_df_list[grepl(paste(cluster5, collapse = '|'), names(patterns_genes_df_list))]
name_of_subset <- 'cluster5'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_spatialde('KEGG', cluster5_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_clusterp <- get_clusterp_analysis_spatialde('Reactome', cluster5_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_clusterp <- get_clusterp_analysis_spatialde('GO', cluster5_patterns, padj_cutoff = 0.05, genecount_cutoff = 5, filename = name_of_subset, number_of_pathways_cutoff = 25)

### 3. MGSA ####
kegg_results_mgsa <- get_mgsa_analysis_spatialde('KEGG', cluster5_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
reactome_results_mgsa <- get_mgsa_analysis_spatialde('Reactome', cluster5_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)
gobp_results_mgsa <- get_mgsa_analysis_spatialde('GO', cluster5_patterns, filename = name_of_subset, number_of_pathways_cutoff = 25)

