# Name: spatial_granuloma_dge.R
# Author: Laura Twomey
# Date of creation: 16 May 2022
# Heatmaps and potential marker genes

# Clean Environment ----------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects including hidden objects
gc() #free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F) #avoid truncated output in R console and scientific notation

# Set working directory ----------------------------------------------------
path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Heatmaps_all"
setwd(path)
# path_ga <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/GA_slides/AEH/"
# path_gn_gs <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/GN_GS_slides/AEH/"

# Libraries and functions ------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(anndata)
library(reticulate)
library(hdf5r)
library(SingleCellExperiment)
library(pheatmap)
library(mgsa)
library(fgsea)
library(xlsx)

# For distinctive colours
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# col=sample(color, 17)
# n = 17
# pie(rep(1,n), col=sample(color, n))

## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: Patterns to python format ####
to_python_format <- function(list_of_patterns){
  list_in_python_format <- cat(paste(shQuote(gsub('\\.', '_pattern', list_of_patterns), type="cmd"), collapse=", "))
  return(list_in_python_format)
}


# Read in files  ------------------------------------------------

# Merge both frequency tables into 1, and save
# list.files(path = path_ga, pattern = 'freqtable')
# freqtable_GA <- read.csv(paste0(path_ga, "freqtable_GA_patterns.csv"), row.names = 1)
# list.files(path = path_gn_gs, pattern = 'freqtable')
# freqtable_GN_GS <- read.csv(paste0(path_gn_gs, "freqtable_GN_GS_patterns.csv"), row.names = 1)
# freqtable_all <- merge(freqtable_GA, freqtable_GN_GS, by=0, all = T)
# freqtable_all[is.na(freqtable_all)] <- 0
# rownames(freqtable_all) <- freqtable_all$Row.names
# freqtable_all$Row.names <- NULL
# saveRDS(freqtable_all, file = paste(path, 'freqtable_all.csv', sep = '/'))

# Annotations dataframe
# setwd("/Users/mendenlab/work/spatial_granuloma/scripts/summary")
# ad <- anndata::read_h5ad("../../results/current/final/Granuloma_QC_clustering.h5")
# names(ad$obsm)
# names(ad$obs)
# ann_df <- cbind(ad$obs)
# ann_df <- ann_df[grepl('sample|slide|disease|DISEASE|SAMPLE|^patient$', colnames(ann_df))]
# saveRDS(ann_df, file = paste(path, 'ann_df.csv', sep = '/'))

freqtable_all <- readRDS(file = paste(path, 'freqtable_all.csv', sep = '/'))
ann_df <- readRDS(file = paste(path, 'ann_df.csv', sep = '/'))
ann_df <- unique(ann_df)
ann_df$sample_SPECIMEN <- paste(ann_df$sample, ann_df$patient, sep = '_')
ann_df$sample_SPECIMEN_SAMPLE <- paste(ann_df$sample_SPECIMEN, ann_df$SAMPLE, sep = '.')
ann_df <- ann_df %>% dplyr::slice(rep(1:n(), 8))
ann_df$pattern <- c(rep(0, 43), rep(1, 43), rep(2, 43), rep(3, 43), rep(4, 43), rep(5, 43), rep(6, 43), rep(7, 43))
ann_df$sample_SPECIMEN_pattern <- paste(ann_df$sample_SPECIMEN_SAMPLE, '_pattern', ann_df$pattern, sep = '')
rownames(ann_df) <- ann_df$sample_SPECIMEN_pattern

# !!!! Check that you are importing the most recent version of the intensities dataframe
# Use this for the normal intensity ratio (mean, corrected by the min)
#intensity_ratio_plot_df <- read.csv(file = paste(path, 'intensity_ratio_plot_df.csv', sep = '/'), row.names = 1)
# Use this for the softmax corrected intensities
intensity_ratio_path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/Intensity_ratios"
intensity_ratio_plot_df <- read.csv(file = paste(intensity_ratio_path, 'intensity_ratio_soft.csv', sep = '/'))
head(intensity_ratio_plot_df)
colnames(intensity_ratio_plot_df) <- gsub('mse_value', 'intensity_ratio', colnames(intensity_ratio_plot_df)) # rename if it is MSE dataframe to avoid changing the next few lines all the time
intensity_ratio_plot_df <- intensity_ratio_plot_df[,c('annotation', 'sample_SPECIMEN_pattern', 'intensity_ratio')] # mse_Value or intensity_ratio
intensity_ratio_plot_df <- as.data.frame(tidyr::pivot_wider(intensity_ratio_plot_df, names_from = 'annotation', values_from = 'intensity_ratio'))
intensity_ratio_plot_df$sample_SPECIMEN_pattern<- gsub('-', '.', intensity_ratio_plot_df$sample_SPECIMEN_pattern)
intensity_ratio_plot_df[,grepl('leiden_', colnames(intensity_ratio_plot_df))] <- round(intensity_ratio_plot_df[,grepl('leiden_', colnames(intensity_ratio_plot_df))], 2)
rownames(intensity_ratio_plot_df) <- intensity_ratio_plot_df$sample_SPECIMEN_pattern
intensity_ratio_plot_df <- intensity_ratio_plot_df[,grepl('leiden_|sample_SPECIMEN_pattern', colnames(intensity_ratio_plot_df))]
head(intensity_ratio_plot_df)

# Set colours and add annotations ------------------------------------------------

# Set patients for subsetting freqtable
gannulare_patients <- c('91253', '45703', '50107')
necrobiosis_patients <- c('95096')
sarcoidosis_patients <- c('72859')

# interesting granuloma patterns
g_patterns <- c('P18554_1001_50107-A_pattern1',
                'P17851_1002_91253-A_pattern5',
                'P17851_1003_45703-A_pattern6',
                'P17851_1004_45703-A_pattern7',
                'P18554_1002_50107-A_pattern5',
                'P17851_1001_91253-A_pattern6',
                'P17851_1002_91253-B_pattern1',
                'P17851_1002_91253-A_pattern7',
                'P17851_1002_91253-A_pattern6',
                'P17851_1003_45703-A_pattern7',
                'P17851_1004_45703-A_pattern6',
                'P18554_1001_50107-A_pattern3',
                'P18554_1002_50107-A_pattern1', 
                'P18554_1007_72859-A_pattern3',
                'P18554_1006_82301-B_pattern7',
                'P18554_1004_95096-A_pattern3',
                'P18554_1005_82301-B_pattern1',
                'P18554_1003_95096-A_pattern6',
                'P18554_1006_82301-A_pattern5',
                'P18554_1004_95096-A_pattern7',
                'P18554_1005_82301-A_pattern4',
                'P18554_1007_72859-B_pattern5',
                'P18554_1008_72859-A_pattern4',
                'P18554_1008_72859-B_pattern3',
                'P18554_1004_95096-A_pattern2',
                'P18554_1008_72859-A_pattern3',
                'P18554_1005_82301-A_pattern0',
                'P18554_1008_72859-B_pattern1',
                'P18554_1006_82301-A_pattern4',
                'P18554_1007_72859-B_pattern6')
g_patterns_list <- paste(g_patterns, collapse = '|')
g_patterns_list <- gsub('-', '.', g_patterns_list)

all_g_patterns <- c("P18554_1002_50107-A_pattern5", "P18554_1003_95096-A_pattern6", "P18554_1006_82301-A_pattern5",
                  "P17851_1002_91253-B_pattern1","P18554_1007_72859-A_pattern3","P17851_1003_45703-A_pattern6",
                  "P17851_1004_45703-A_pattern7", "P17851_1001_91253-A_pattern6", "P18554_1006_82301-B_pattern7",
                  "P18554_1004_95096-A_pattern3", "P18554_1005_82301-B_pattern1", "P18554_1007_72859-B_pattern5",
                  "P18554_1005_82301-A_pattern4", "P18554_1008_72859-A_pattern4", "P17851_1002_91253-A_pattern5",
                  "P18554_1008_72859-B_pattern3", "P18554_1001_50107-A_pattern1", "P18554_1004_95096-A_pattern7")

all_g_patterns2_1 <- c('P17851_1002_91253-A_pattern1', 'P18554_1005_82301-A_pattern0', 'P18554_1008_72859-A_pattern3', 
                     'P18554_1007_72859-B_pattern6', 'P18554_1008_72859-B_pattern1', 'P18554_1001_50107-A_pattern7', 
                     'P18554_1002_50107-A_pattern4')

# Subset to only the sample_SPECIMEN_patterns we are interested in (e.g., no sample = 0)
col_annotations_df <- ann_df[grepl(paste(colnames(freqtable_all), collapse = '|'), ann_df$sample_SPECIMEN_pattern), c('patient', 'sample_SPECIMEN_pattern', 'DISEASE')]
col_annotations_df$g_pattern <- ifelse(grepl(g_patterns_list, row.names(col_annotations_df)), 'YES', 'NO')
col_annotations_df_new <- merge(col_annotations_df, intensity_ratio_plot_df, by = 'sample_SPECIMEN_pattern')
rownames(col_annotations_df_new) <- col_annotations_df_new$sample_SPECIMEN_pattern

ann_colors = list(
  patient = c("91253" = "#F46D43",#'red',#granuloma annulare
              "45703" = "#D53E4F",#'indianred',#granuloma annulare
              "50107" = "#9E0142",#'firebrick',#granuloma annulare
              "95096" = "orchid", # necrobiosis lipoidica
              "82301" = "mediumvioletred", #sarcoidosis-cutaneous
              "72859" = "blueviolet"), #sarcoidosis-suspected
  DISEASE = c("granuloma annulare" = "#9E0142",#'firebrick'
              "necrobiosis lipoidica" = "orchid", # necrobiosis lipoidica
              "sarcoidosis" = "blueviolet"), #sarcoidosis
  g_pattern = c('YES' = 'blue', 'NO' = 'grey')
  #sample_SPECIMEN_SAMPLE = sample(color, 20)
)

# Heatmaps  ------------------------------------------------
## All  ------------------------------------------------
heatmap <- pheatmap(freqtable_all,
         show_colnames = F, show_rownames = F,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = col_annotations_df_new,
         # annotation_row = data.frame(Ascore = sc.results.batch[["May22singleB"]]$est.table$`A score`,
         #                             phenotypes[phenotypes$dataset=="May22singleB", c("drug", "growth")],
         #                             row.names = rownames(featureMatrix[phenotypes$dataset=="May22singleB",])),
         annotation_colors = ann_colors
         )
heatmap
# png(file = paste(path, 'heatmap_all.png', sep = '/'), width = 800, height = 900)
# heatmap
# dev.off()
# 
# colnames(freqtable_all[,heatmap$tree_col$order]) %in% gsub('-', '.', ga_patterns)
# colnames(freqtable_all[,heatmap$tree_col$order])[119:136]
# paste(colnames(freqtable_all[,heatmap$tree_col$order])[1:15], collapse = '', '')
# paste(shQuote(gsub('\\.', '-', colnames(freqtable_all[,heatmap$tree_col$order])[1:15])), collapse=", ")

# freqtable_g <- freqtable_all[,colnames(freqtable_all[,heatmap$tree_col$order])[119:136]]
# freqtable_g$Total <- rowSums(freqtable_g)
# freqtable_g$genes <- rownames(freqtable_g)
# freqtable_g <- freqtable_g[freqtable_g$Total > 4,]
# 
# list_potential_granuloma_markers <- freqtable_g[with(freqtable_g, order(Total, decreasing = TRUE)),]$genes
# write.table(list_potential_granuloma_markers, file = paste(path, 'list_potential_granuloma_markers.txt', sep = '/'))
# 
# 
# gsub('\\.', '-', colnames(freqtable_all[,heatmap$tree_col$order])[119:136])

