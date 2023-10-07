# Name: PCA_UMAP.R
# Author: Laura Twomey
# Date of creation: 20 June 2022
# PCA, UMAPs, tSNEs

# Clean Environment ----------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects including hidden objects
gc() #free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F) #avoid truncated output in R console and scientific notation

# Set working directory ----------------------------------------------------
path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/PEA"
setwd(path) 
pca_path <- "/Volumes/Drive/spatial_granuloma/output/SpatialDE/PCAs"

# Libraries and functions ------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(anndata)
library(reticulate)
library(hdf5r)
library(pheatmap)

# Read in files

# freqtable_all is a binary matrix of all genes (rownames) showing if they are present or not in each of the sample_SPECIMEN_patterns (columns)
freqtable_all <- readRDS(file = paste(path, 'freqtable_all.csv', sep = '/'))
# df_all has all the data per spot for each pattern (sample_SPECIMEN_pattern)
df_all <- read.csv(file = '../Intensity_ratios/patterns_samples_merged.csv', sep = ',', row.names = 1)
# wilcox_df with the results from wilcoxon test
wilcox_plot_df <- read.csv(file = '../Intensity_ratios/wilcox_df.csv', sep = ',', row.names = 1)
head(wilcox_plot_df)

