# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# scRNAseq doublet detection tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Version: 1.0

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(Seurat)
library(scDblFinder)
library(DoubletFinder)

set.seed(42)

# Functions ===================================================================
#----------------------------------------------------------#
# run_doubletfinder_custom
#----------------------------------------------------------#
# run_doubletfinder_custom runs Doublet_Finder() and returns a dataframe with the cell IDs and a column with either 'Singlet' or 'Doublet'
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[['SampleID']]), '...........')) 
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_var <- (stdv^2/sum(stdv^2)) * 100
  cumulative_var <- cumsum(percent_var)
  co1 <- which(cumulative_var > 90)[1]
  co2 <- which(diff(percent_var) < 0.1)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  # Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  # Get the multiplet rate if not provided
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # Subset and save
  # head(sample@meta.data['doublet_finder'])
  # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
  # singlets$ident
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}


# Relevant paths ======================================================================
project_path <- "~/Biostatsquid/Projects/scRNAseq_tbr"
in_path <- file.path(project_path, 'data')
out_path <- file.path(project_path, 'results') # for all scRNAseq preprocessing results
doublet_folder <- file.path(out_path, 'Doublet_detection') # subfolder for Doublet detection results
#dir.create(doublet_folder, recursive = TRUE, showWarnings = FALSE) 

# Inputs ======================================================================
# Read in your Seurat object
seu <- readRDS(file.path(in_path, "seu_filt.rds"))
head(seu@meta.data)
table(seu$SampleID)

# Doublet detection =================================================================

## scDblFinder -------------------------------
# Run scDbltFinder
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce, samples = "SampleID") #dbr = multiplet_rate)
table(sce$scDblFinder.class)
sce@colData@listData %>% as.data.frame() %>% head()

# Explore results and add to seurat object
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder)
rownames(meta_scdblfinder) <- sce@colData@rownames
head(meta_scdblfinder)
seu <- AddMetaData(object = seu, metadata = meta_scdblfinder %>% dplyr::select('scDblFinder.class'))
head(seu@meta.data)
table(seu$scDblFinder.class)
rm(list = c('meta_scdblfinder', 'sce'))

# Doublet stats
# Check how doublets singlets differ in QC measures per sample.
VlnPlot(seu, group.by = 'SampleID', split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')

doublets_summary <- seu@meta.data %>% 
  group_by(SampleID, scDblFinder.class) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  group_by(SampleID) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(scDblFinder.class, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)
doublets_summary
write.table(doublets_summary, file = file.path(doublet_folder, paste0('scDblFinder_doublets_summary.txt')), quote = FALSE, row.names = FALSE, sep = '\t')

## DoubletFinder -------------------------------
# DoubletFinder should be run on a per sample basis
samp_split <- SplitObject(seu, split.by = "SampleID") 
# Get Doublet/Singlet IDs by DoubletFinder()
samp_split <- lapply(samp_split, run_doubletfinder_custom) # get singlet/doublet assigned to each of the cell IDs (each element of the list is a different sample)
head(samp_split)
sglt_dblt_metadata <- data.frame(bind_rows(samp_split)) # merge to a single dataframe
rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names # assign cell IDs to row names to ensure match
sglt_dblt_metadata$row_names <- NULL
head(sglt_dblt_metadata)
seu <- AddMetaData(seu, sglt_dblt_metadata, col.name = 'doublet_finder')

# Check how doublets singlets differ in QC measures per sample.
VlnPlot(seu, group.by = 'SampleID', split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')

# Get doublets per sample
doublets_summary <- seu@meta.data %>% 
  group_by(SampleID, doublet_finder) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  group_by(SampleID) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(doublet_finder, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)
write.table(doublets_summary, file = file.path(doublet_folder, paste0(proj, '_doubletfinder_doublets_summary.txt')), quote = FALSE, row.names = FALSE, sep = '\t')

# Doublet QC: scDblFinder vs DoubletFinder ==========================
# Percentage of doublets per sample
seu@meta.data %>% group_by(SampleID) %>% summarise(percent_doublet_doubletfinder = mean(doublet_finder == 'Doublet') * 100,
                                                   percent_doublet_scdblfinder = mean(scDblFinder.class == 'doublet') * 100)

table('scDblFinder' = seu$scDblFinder.class, 'DoubletFinder' = seu$doublet_finder)
meta <- seu@meta.data
meta <- meta %>% mutate(doublet_consensus = 
                          ifelse(doublet_finder == 'Singlet' & scDblFinder.class == 'singlet', 'singlet', 
                                 ifelse(doublet_finder == 'Doublet' & scDblFinder.class == 'doublet', 'doublet', 
                                        ifelse(doublet_finder == 'Singlet' & scDblFinder.class == 'doublet', 'DFsinglet_SCDBLdoublet', 'DFdoublet_SCDBLsinglet'))))

table(meta$doublet_consensus, meta$SampleID)

ggplot(meta, aes(x = doublet_consensus, y = nFeature_RNA, fill = doublet_consensus)) +
  geom_violin() +
  ylab('Number of genes/cell') +
  facet_wrap(~SampleID) +
  theme_bw()

ggplot(meta, aes(x = doublet_consensus, y = nCount_RNA, fill = doublet_consensus)) +
  geom_violin() +
  ylab('Number of molecules/cell') +
  facet_wrap(~SampleID)

ggplot(meta, aes(x = nFeature_RNA, fill = doublet_consensus)) +
  geom_density(alpha = 0.1) +
  xlab('Number of genes/cell') +
  facet_wrap(~SampleID)

ggplot(meta, aes(x = nCount_RNA, fill = doublet_consensus)) +
  geom_density(alpha = 0.1) +
  xlab('Number of molecules/cell') +
  facet_wrap(~SampleID)

# Save seurat with doublet annotation
saveRDS(seu, file = file.path(in_path, 'seu_dblt.rds'))

# Removing doublets ============================================================
# If all is ok subset and remove doublets
# seu_dblt <- readRDS(file.path(in_path, 'seu_dblt.rds'))
seu_dblt <- subset(seu,
                   subset = doublet_finder == 'Singlet' | scDblFinder.class == 'singlet')
table('scDblFinder' = seu_dblt$scDblFinder.class, 'DoubletFinder' = seu_dblt$doublet_finder)

stopifnot('Something went wrong, genes have been removed!' = dim(seu)[[1]] - dim(seu_dblt)[[1]] == 0)
cat('Total doublets removed: ', dim(seu)[[2]] - dim(seu_dblt)[[2]]) # removed 3500 doublets in total

saveRDS(seu_dblt, file = file.path(in_path, 'seu_dblt_filt.rds'))

cat('\n\n----------- Doublet detection & filtering complete!------------\n\n')
rm(seu_dblt)
rm(seu)
