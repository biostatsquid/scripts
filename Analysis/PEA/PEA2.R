# Name: Pathway_enrichment_analysis_DEG.R
# Author: Laura Twomey
# Date of creation: 03 July 2022
# Pathway enrichment analysis for DEG results

# Clean Environment ----------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects including hidden objects
gc() #free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F) #avoid truncated output in R console and scientific notation

# Set working directory ----------------------------------------------------
path <- "/Volumes/Drive/spatial_granuloma/output/DEG/PEA/"
setwd(path)
deg_path <- "/Volumes/Drive/spatial_granuloma/output/DEG/"
mgsa_path <- "/Volumes/Drive/spatial_granuloma/output/DEG/PEA/MGSA/"
clusterp_path <- "/Volumes/Drive/spatial_granuloma/output/DEG/PEA/clusterProfiler/"
gsea_path <-  "/Volumes/Drive/spatial_granuloma/output/DEG/PEA/GSEA/"
glmgampoi_path <- "/Volumes/Drive/spatial_granuloma/output/DEG/glmGamPoi_results/"

# Libraries  ------------------------------------------------
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
# library('org.Hs.eg.db') # to get ensembl entrez ids to match genes
# columns(org.Hs.eg.db)

# Functions ------------------------------------------------
## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: get_clusterp_analysis_deg ####
# This function runs clusterp analysis for a set of background genes (for now KEGG, Reactome or GO), given a list of deg_results from 
# the DEG pipeline (list of genes are in $gene_symbol column). Pathways can be filtered by padj and gene count.
# If get_heatmap is set to TRUE (default) it also calls the heatmap from clusterp_get_heatmap_DEG
get_clusterp_analysis_deg <- function(background_genes, deg_results_list, padj_cutoff = 0.05, genecount_cutoff = 5, filename = 'example', get_heatmap = TRUE, number_of_pathways_cutoff = 30){
  # for debug
  # background_genes <- 'Reactome'
  # deg_results_list <- genes_df_list
  # padj_cutoff <- 0.05
  # genecount_cutoff <- 5
  
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
  
  # Run clusterProfiler on each sub-dataframe
  go_all <- lapply(names(deg_results_list),
                   function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
                                        TERM2GENE = pwl2))
  names(go_all) <- names(deg_results_list)
  
  #Convert the go_all list of enrichResults for each sample_pattern to a dataframe with the pathways
  go_all_df<- lapply(names(go_all), function(x) rbind(go_all[[x]]@result))
  names(go_all_df) <- names(go_all)
  go_all_df = do.call(rbind, go_all_df)
  #head(go_all_df)
  colnames(go_all_df) <- gsub('ID', 'pathways', colnames(go_all_df)) # rename the pathways column
  go_all_df <- go_all_df %>% mutate(minuslog10padj = -log10(p.adjust),
                                    diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(go_all_df)))
  rownames(go_all_df) <- NULL
  
  # Subset to those pathways that have p adj < cutoff and gene count > cutoff
  target_pws <- unique(go_all_df$pathways[go_all_df$p.adjust < padj_cutoff & go_all_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
  go_all_df <- go_all_df[go_all_df$pathways %in% target_pws, ]
  
  print('Saving clusterprofiler results')
  write.csv(go_all_df, paste0(clusterp_path, filename, background_genes, '.csv'), row.names = FALSE)
  
  # # Check distribution of minus log 10 p adj values
  # print('Saving minuslog(padj) distribution plot')
  # png(filename = paste0(path, 'clusterProfiler/', filename, background_genes, '_minuslog10padj_distribution.png'), height = 600, width = 600)
  # ggplot(data = go_all_df, aes(x = minuslog10padj)) +
  #   geom_histogram(col = 'darkgreen', fill = 'aquamarine', bins = 100) + 
  #   scale_x_continuous(breaks = seq(0, 155, 10)) + 
  #   theme_bw() +
  #   ggtitle(paste(filename, "Distribution of -log10(padj) values")) 
  # dev.off()
  
  if(get_heatmap == TRUE){
    clusterp_get_heatmap_DEG(dataframe = go_all_df, file_name = paste0(filename, background_genes, '_clusterp_heatmap'), pathway_cutoff = number_of_pathways_cutoff)
  }
  return(go_all_df)
}

## Function: Get heatmap_df ####
clusterp_get_heatmap_DEG <- function(dataframe, file_name = 'filename',  pathway_cutoff = number_of_pathways_cutoff){
  # for debug
  # dataframe <- go_all_df 
  if(!is.data.frame(dataframe)){
    stop('results_filtered must be a dataframe')
  }
  if(!is.character(file_name)){
    stop('file_name must be a character')
  }
  
  # First convert the dataframe into the heatmap dataframe we need to build the heatmap
  heatmap_df <- tidyr::pivot_wider(dataframe[,c('pathways', 'minuslog10padj', 'diffexpressed')], values_from = minuslog10padj, names_from = diffexpressed, values_fill = 0)
  heatmap_df <- as.data.frame(heatmap_df)
  rownames(heatmap_df) <- heatmap_df$pathways # to be able to set them as row names
  heatmap_df <- heatmap_df[!grepl('pathways', colnames(heatmap_df))] # now remove the pathways column
  head(heatmap_df)
  #write.csv(heatmap_df, paste(path, 'clusterProfiler/heatmap_df.csv', sep = '/'), row.names = TRUE)
  if(nrow(heatmap_df) > pathway_cutoff){
    selected_pathways <- head(unique(dataframe[order(-dataframe$minuslog10padj),'pathways']), pathway_cutoff) # order by minuslogpadj and get ones with larger minuslogpadj
    heatmap_df <- subset(heatmap_df, rownames(heatmap_df) %in% selected_pathways) # subset to less pathways for visualisation purposes if there are too many
  }
  
  ## Get heatmap 
  heatmap <- pheatmap(heatmap_df,
                      show_colnames = T, show_rownames = T,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      #annotation_col = col_annotations_df_new,
                      fontsize_col = 10,
                      fontsize_row = 10,
                      #annotation_colors = ann_colors
  )
  
  pdf(file = paste(clusterp_path, file_name, '.pdf', sep = ''), width = 15, height = 15)
  grid::grid.draw(heatmap$gtable)
  dev.off()
  
  return(heatmap)
}

## Function: Run MGSA analysis and get heatmap ####
# This function runs mgsa analysis for a set of background genes (for now KEGG, Reactome or GO), given a list of deg_results from 
# the DEG pipeline (list of genes are in $gene_symbol column). Pathways can be filtered by padj and gene count.
run_mgsa_analysis_deg <- function(background_genes, deg_results_list, filename = 'example', get_heatmap = TRUE, number_of_pathways_cutoff = 30){
  # for debug
  #background_genes <- 'Reactome'
  #deg_results_list <- genes_df_list
  
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
  pwl1$term <- as.character(pwl1$term) # mgsa doesn't like it when the entries are factors
  
  ##  Run mgsa for subset of genes in deg_results_list
  mgsa_fit_list <- lapply(names(deg_results_list), function(x) mgsa(deg_results_list[[x]]$gene_symbol, pwl1))
  names(mgsa_fit_list) <- names(deg_results_list)
  saveRDS(mgsa_fit_list, file = paste0(mgsa_path, 'Fit_results/', paste0(filename, background_genes, '_mgsa_fit_list.RData')))
  print('MGSA fit finished!')
  
  ## Get the enriched pathways 
  res_list <- lapply(names(mgsa_fit_list), function(x) setsResults(mgsa_fit_list[[x]]))
  names(res_list) <- names(deg_results_list)
  saveRDS(res_list, file = paste0(mgsa_path, 'Fit_results/', paste0(filename, background_genes, '_mgsa_res_list.RData'))) 
  
  res_list <- Map(cbind, res_list, group = names(res_list))
  
  res_all_df = do.call(rbind, res_list)
  res_all_df$pathways <- gsub('^.*GOBP_|^.*KEGG_|^.*REACTOME_', '', rownames(res_all_df))
  write.csv(res_all_df, paste0(mgsa_path, 'Fit_results/', paste0(filename, background_genes, 'MGSA_res_all_.csv')))
  
  target_pathways <- unique(res_all_df$pathways[res_all_df$estimate > 0.5]) # target pathways have an estimate higher than 0.5
  res_all_df_filtered <- res_all_df[res_all_df$pathways %in% target_pathways, ]
  write.csv(res_all_df_filtered, paste0(mgsa_path, paste0(filename, background_genes, '_mgsa_res_filtered.csv')))
  
  if(get_heatmap == TRUE){
    mgsa_get_heatmap(results_filtered = res_all_df_filtered, file_name = paste0(filename, background_genes, '_mgsa_heatmap'), pathway_cutoff = number_of_pathways_cutoff)
  }
  return(res_all_df_filtered)
}

## Function: MGSA heatmaps ####
# This function takes the filtered results from mgsa pipeline already subset to estimate > 0.5 and returns a heatmap
mgsa_get_heatmap <- function(results_filtered, file_name = 'filename',  pathway_cutoff = number_of_pathways_cutoff){
  # for debug
  #results_filtered <- res_all_df_filtered 
  if(!is.data.frame(results_filtered)){
    stop('results_filtered must be a dataframe')
  }
  if(!is.character(file_name)){
    stop('file_name must be a character')
  }
  # Get heatmap
  heatmap_df <- tidyr::pivot_wider(results_filtered[,c('estimate', 'group', 'pathways')], values_from = estimate, names_from = group, values_fill = 0)
  heatmap_df <- as.data.frame(heatmap_df)
  
  rownames(heatmap_df) <- heatmap_df$pathways # to be able to set them as row names
  heatmap_df <- heatmap_df[!grepl('pathways', colnames(heatmap_df))] # now remove the pathways column
  head(heatmap_df)
  
  if(nrow(heatmap_df) > pathway_cutoff){
    selected_pathways <- head(unique(results_filtered[order(-results_filtered$estimate),'pathways']), pathway_cutoff) # order by minuslogpadj and get ones with larger minuslogpadj
    heatmap_df <- subset(heatmap_df, rownames(heatmap_df) %in% selected_pathways) # subset to less pathways for visualisation purposes if there are too many
  }
  
  heatmap <- pheatmap(heatmap_df,
                      show_colnames = T, show_rownames = T,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      #annotation_col = col_annotations_df,
                      fontsize_col = 10,
                      fontsize_row = 10#,
                      #annotation_colors = ann_colors
  )
  
  pdf(file = paste0(mgsa_path, file_name, '.pdf'), width = 15, height = 15)
  grid::grid.draw(heatmap$gtable)
  dev.off()
  
  return(heatmap)
}

## Function: get_gsea_analysis_spatialde ####
# This function runs gsea analysis for a set of background genes (for now KEGG, Reactome or GO), given a list of deg spatialde_results from 
# the spatialDE pipeline (list of genes are in $genes column). Pathways can be filtered by padj and gene count.
# We will use signed p value (sign(deg_results$log2fc)*(-log10(deg_results$pval))) as rankings
run_gsea_analysis_deg <- function(background_genes, deg_results, filename = 'example', number_of_top_pathways_up = 25, number_of_top_pathways_down = 5){
  # for debug
  # background_genes <- 'GO'
  # deg_results <- dge_epivsdermis_withinterface
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

## Function: GSEA heatmaps ####
# This function takes the filtered results from gsea pipeline already subset to estimate > 0.5 and returns a heatmap
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
  
  height_of_heatmap <- nrow(heatmap_df) * 0.2 # to make the heatmap longer in case there are many rows
  width_of_heatmap <- ncol(heatmap_df) * 3 # to make the heatmap longer in case there are many columns
  
  pdf(file = paste0('GSEA/Heatmaps/', file_name, '_heatmap.pdf'), width = width_of_heatmap, height = height_of_heatmap)
  grid::grid.draw(heatmap$gtable)
  dev.off()
  
  return(heatmap)
}



# Background genes ------------------------------------------------
## 15107 genes in data
# genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_epivsdermis_withinterface_patient+epivsdermis_degresults.csv', sep = '/'), row.names = 1)
# genes_in_data <- genes_df$gene_symbol

# For MGSA and GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
# gmt_files <- list.files(path = path, pattern = '.gmt') # cluster_profiler needs a special parsing of the .gmt file, we cannot use pwl1 we created for MGSA analysis
# for (file in gmt_files){
#   # Read in gmt file
#   gmt <- gmtPathways(file)
#   hidden <- unique(unlist(gmt))
#   # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
#   Ass <- matrix(NA, dimnames = list(hidden, names(gmt)),
#                 nrow = length(hidden), ncol = length(gmt))
#   for (i in 1:dim(Ass)[2]){
#     Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
#   }
#   
#   #Subset to the genes that are present in our data to avoid bias
#   hidden1 <- intersect(genes_in_data, hidden)
#   Ass1 <- Ass[hidden1, colnames(Ass)[which(colSums(Ass[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
#   # And get the list again
#   pwl1 <- matrix_to_list(Ass1) # for this we use the function we previously defined
#   saveRDS(pwl1, file = paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', as.character(file))), '.RData', sep = ''))
# }
# 
# pwl1 <- readRDS('pwl1.RData') # background gene sets (pathways) already filtered for those genes that are present in our data to avoid bias + gene sets with more than 5 genes annotated

# For clusterProfiler
# Filter out the gmt files for KEGG, Reactome and GOBP
# gmt_files <- list.files(path = path, pattern = '.gmt') # cluster_profiler needs a special parsing of the .gmt file, we cannot use pwl1 we created for MGSA analysis
# for (file in gmt_files){
#   pwl2 <- read.gmt(file) # all genes from db (GOBP, Reactome, KEGG...)
#   pwl2 <- pwl2[pwl2$gene %in% genes_in_data,] # background gene sets (pathways) already filtered for those genes (15107) that are present in our data to avoid bias
#   saveRDS(pwl2, paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', as.character(file))), '.RDS', sep = ''))
# }


# PEA analysis  ------------------------------------------------
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = '.csv') # results from DEG

## Epidermis vs dermis (with interface) ####
### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_epivsdermis_withinterface_patient+epivsdermis_degresults.csv', sep = '/'), row.names = 1)
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df$diffexpressed <- gsub('DOWN', 'Epidermis (+ interface)', gsub('UP', 'Healthy dermis', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'epivsdermis_withinterface'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## Epidermis vs dermis (without interface) ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'without.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_epivsdermis_withoutinterface_patient+epivsdermis_withoutinterface_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Epidermis', gsub('UP', 'Healthy dermis', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'epivsdermis_withoutinterface'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## Dermis lesional vs non-lesional ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'dermis.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_dermislesionalvsnonlesional_patient+skin_layer+dermis_lesvsnonles_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Lesional dermis', gsub('UP', 'Non-lesional dermis', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'dermis_lesvsnonles'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## GA: core vs border ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'ga_.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_gacorevsborder_patient+ga_corevsborder_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Border', gsub('UP', 'Core', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'ga_corevsborder'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## SA: core vs border ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'sa_.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_sacorevsborder_patient+sa_corevsborder_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Border', gsub('UP', 'Core', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'sa_corevsborder'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)


## GA: granuloma vs dermis ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'ga_gvsd.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_ga_gvsd_patient+ga_gvsd_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Granuloma (G. ann)', gsub('UP', 'Healthy dermis', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'ga_gvsd'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## SA: granuloma vs dermis ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'sa_gvsd.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_sa_gvsd_patient+sa_gvsd_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Healthy dermis', gsub('UP', 'Granuloma (sarcoidosis)', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'sa_gvsd'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## NL: granuloma vs dermis ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'nl_gvsd.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_nl_gvsd_nl_gvsd_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Healthy dermis', gsub('UP', 'Granuloma (Necrobiosis l.)', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'nl_gvsd'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)











## Granuloma Manual: granuloma vs dermis ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'manual.*gvsd.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_manual_gvsd_patient+manual_gvsd_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Granuloma (manual)', gsub('UP', 'Healthy dermis', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'manual_gvsd'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## Granuloma Manual 2: granuloma vs dermis ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'manual.*gvsd.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_manual_gvsd_morecovariates_patient+skin_layer+manual_gvsd_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Granuloma (manual)', gsub('UP', 'Healthy dermis', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'manual2_gvsd'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)

## Granuloma Leiden: granuloma vs dermis ####
list.files(paste0(deg_path, '/glmGamPoi_results'), pattern = 'leiden.*gvsd.*csv') # results from DEG

### 1. Reading in data ####
genes_df <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_leiden_gvsd_patient+leiden_gvsd_degresults.csv', sep = '/'), row.names = 1)
unique(genes_df$conditionreference)[1] # get the name of the reference of the condition
genes_df <- genes_df[genes_df$diffexpressed != 'NO', ] # remove non-significant genes
genes_df$diffexpressed <- gsub('DOWN', 'Healthy dermis', gsub('UP', 'Granuloma (Leiden)', genes_df$diffexpressed))
unique(genes_df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
genes_df_list <- split(genes_df, genes_df$diffexpressed)
name_of_comparison <- 'leiden_gvsd'

### 2. clusterProfiler ####
kegg_results_clusterp <- get_clusterp_analysis_deg('KEGG', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_clusterp <- get_clusterp_analysis_deg('Reactome', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_clusterp <- get_clusterp_analysis_deg('GO', genes_df_list, padj_cutoff = 0.05, genecount_cutoff = 5, name_of_comparison, number_of_pathways_cutoff = 30)

### 3. MGSA ####
kegg_results_mgsa <- run_mgsa_analysis_deg('KEGG', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
reactome_results_mgsa <- run_mgsa_analysis_deg('Reactome', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)
gobp_results_mgsa <- run_mgsa_analysis_deg('GO', genes_df_list, name_of_comparison, number_of_pathways_cutoff = 30)


# GSEA  ------------------------------------------------
setwd(path)

## Epidermis vs dermis (with interface) ####
name_of_subset <- 'epivsdermis_withinterface'
dge_epivsdermis_withinterface <- read.csv(paste0(glmgampoi_path, 'sce_epivsdermis_withinterface_patient+epivsdermis_degresults.csv'), row.names = 1)

gsea_epivsdermis_withinterface_go <- run_gsea_analysis_deg('GO', dge_epivsdermis_withinterface, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_epivsdermis_withinterface_kegg <- run_gsea_analysis_deg('KEGG', dge_epivsdermis_withinterface, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_epivsdermis_withinterface_reactome <- run_gsea_analysis_deg('Reactome', dge_epivsdermis_withinterface, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Epidermis vs dermis (without interface) ####
name_of_subset <- 'epivsdermis_withoutinterface'
dge_epivsdermis_withoutinterface <- read.csv(paste0(glmgampoi_path, 'sce_epivsdermis_withoutinterface_patient+epivsdermis_withoutinterface_degresults.csv'), row.names = 1)

gsea_epivsdermis_withoutinterface_go <- run_gsea_analysis_deg('GO', dge_epivsdermis_withoutinterface, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_epivsdermis_withoutinterface_kegg <- run_gsea_analysis_deg('KEGG', dge_epivsdermis_withoutinterface, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_epivsdermis_withoutinterface_reactome <- run_gsea_analysis_deg('Reactome', dge_epivsdermis_withoutinterface, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Dermis lesional vs nonlesional ####
name_of_subset <- 'dermis_lesvsnonles'
dge_dermislesionalvsnonlesional <- read.csv(paste0(glmgampoi_path, 'sce_dermislesionalvsnonlesional_patient+skin_layer+dermis_lesvsnonles_degresults.csv'), row.names = 1)

gsea_dermislesionalvsnonlesional_go <- run_gsea_analysis_deg('GO', dge_dermislesionalvsnonlesional, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_dermislesionalvsnonlesional_kegg <- run_gsea_analysis_deg('KEGG', dge_dermislesionalvsnonlesional, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_dermislesionalvsnonlesional_reactome <- run_gsea_analysis_deg('Reactome', dge_dermislesionalvsnonlesional, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Dermis lesional vs nonlesional ####
name_of_subset <- 'dermis_lesvsnonles'
dge_dermislesionalvsnonlesional <- read.csv(paste0(glmgampoi_path, 'sce_dermislesionalvsnonlesional_patient+skin_layer+dermis_lesvsnonles_degresults.csv'), row.names = 1)

gsea_dermislesionalvsnonlesional_go <- run_gsea_analysis_deg('GO', dge_dermislesionalvsnonlesional, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_dermislesionalvsnonlesional_kegg <- run_gsea_analysis_deg('KEGG', dge_dermislesionalvsnonlesional, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_dermislesionalvsnonlesional_reactome <- run_gsea_analysis_deg('Reactome', dge_dermislesionalvsnonlesional, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## GA: core vs border ####
name_of_subset <- 'gacorevsborder'
dge_gacorevsborder <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_gacorevsborder_patient+ga_corevsborder_degresults.csv', sep = '/'), row.names = 1)

gsea_gacorevsborder_go <- run_gsea_analysis_deg('GO', dge_gacorevsborder, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_gacorevsborder_kegg <- run_gsea_analysis_deg('KEGG', dge_gacorevsborder, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_gacorevsborder_reactome <- run_gsea_analysis_deg('Reactome', dge_gacorevsborder, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## SA: core vs border ####
name_of_subset <- 'sacorevsborder'
dge_sacorevsborder <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_sacorevsborder_patient+sa_corevsborder_degresults.csv', sep = '/'), row.names = 1)

gsea_sacorevsborder_go <- run_gsea_analysis_deg('GO', dge_sacorevsborder, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_sacorevsborder_kegg <- run_gsea_analysis_deg('KEGG', dge_sacorevsborder, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_sacorevsborder_reactome <- run_gsea_analysis_deg('Reactome', dge_sacorevsborder, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## GA: granuloma vs dermis ####
name_of_subset <- 'ga_gvsd'
dge_ga_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_ga_gvsd_patient+ga_gvsd_degresults.csv', sep = '/'), row.names = 1)

gsea_ga_gvsd_go <- run_gsea_analysis_deg('GO', dge_ga_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_ga_gvsd_kegg <- run_gsea_analysis_deg('KEGG', dge_ga_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_ga_gvsd_reactome <- run_gsea_analysis_deg('Reactome', dge_ga_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## SA: granuloma vs dermis ####
name_of_subset <- 'sa_gvsd'
dge_sa_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_sa_gvsd_patient+sa_gvsd_degresults.csv', sep = '/'), row.names = 1)

gsea_sa_gvsd_go <- run_gsea_analysis_deg('GO', dge_sa_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_sa_gvsd_kegg <- run_gsea_analysis_deg('KEGG', dge_sa_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_sa_gvsd_reactome <- run_gsea_analysis_deg('Reactome', dge_sa_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## NL: granuloma vs dermis ####
name_of_subset <- 'nl_gvsd'
dge_nl_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_nl_gvsd_nl_gvsd_degresults.csv', sep = '/'), row.names = 1)

gsea_nl_gvsd_go <- run_gsea_analysis_deg('GO', dge_nl_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_nl_gvsd_kegg <- run_gsea_analysis_deg('KEGG', dge_nl_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_nl_gvsd_reactome <- run_gsea_analysis_deg('Reactome', dge_nl_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Granuloma Manual2: Granuloma vs dermis ####
name_of_subset <- 'manual_gvsd_morecovariates'
dge_manual_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_manual_gvsd_morecovariates_patient+skin_layer+manual_gvsd_degresults.csv', sep = '/'), row.names = 1)

gsea_manual_gvsd_go <- run_gsea_analysis_deg('GO', dge_manual_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_manual_gvsd_kegg <- run_gsea_analysis_deg('KEGG', dge_manual_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_manual_gvsd_reactome <- run_gsea_analysis_deg('Reactome', dge_manual_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

## Leiden: Granuloma vs dermis ####
name_of_subset <- 'leiden_gvsd'
dge_leiden_gvsd <- read.csv(paste(deg_path, 'glmGamPoi_results', 'sce_leiden_gvsd_patient+leiden_gvsd_degresults.csv', sep = '/'), row.names = 1)

gsea_leiden_gvsd_go <- run_gsea_analysis_deg('GO', dge_leiden_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_leiden_gvsd_kegg <- run_gsea_analysis_deg('KEGG', dge_leiden_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)
gsea_leiden_gvsd_reactome <- run_gsea_analysis_deg('Reactome', dge_leiden_gvsd, filename = name_of_subset, number_of_top_pathways_up = 25, number_of_top_pathways_down = 5)

# Heatmaps ####
getwd()
gsea_ga_gvsd_go <- readRDS(paste0('GSEA/', 'ga_gvsd_gobp_gsea_results.RDS'))
gsea_sa_gvsd_go <- readRDS(paste0('GSEA/', 'sa_gvsd_gobp_gsea_results.RDS'))
gsea_nl_gvsd_go <- readRDS(paste0('GSEA/', 'nl_gvsd_gobp_gsea_results.RDS'))
gsea_manual_gvsd_go <- readRDS(paste0('GSEA/', 'manual_gvsd_morecovariates_gobp_gsea_results.RDS'))
gsea_leiden_gvsd_go <- readRDS(paste0('GSEA/', 'leiden_gvsd_gobp_gsea_results.RDS'))

## Granulomas vs dermis ####
granulomas_compare_go <- list(gsea_ga_gvsd_go, gsea_sa_gvsd_go, gsea_nl_gvsd_go, gsea_manual_gvsd_go, gsea_leiden_gvsd_go)
names(granulomas_compare_go) <- c('G.annulare', 'Sarcoidosis', 'Necrobiosis', 'Granuloma (manual)', 'Granuloma (Leiden)')
gsea_clusters_heatmap <- gsea_get_heatmap(granulomas_compare_go, file_name = 'all_granulomas_go', pathway_cutoff = 50)

# granulomas_compare_kegg <- list(gsea_ga_gvsd_kegg, gsea_sa_gvsd_kegg, gsea_nl_gvsd_kegg, gsea_manual_gvsd_kegg, gsea_leiden_gvsd_kegg)
# names(granulomas_compare_kegg) <- c('G.annulare', 'Sarcoidosis', 'Necrobiosis', 'Manual_granuloma', 'Leiden_granuloma')
# gsea_clusters_heatmap <- gsea_get_heatmap(granulomas_compare_kegg, file_name = 'all_granulomas_kegg', pathway_cutoff = 50)
# 
# granulomas_compare_reactome <- list(gsea_ga_gvsd_reactome, gsea_sa_gvsd_reactome, gsea_nl_gvsd_reactome, gsea_manual_gvsd_reactome, gsea_leiden_gvsd_reactome)
# names(granulomas_compare_reactome) <- c('G.annulare', 'Sarcoidosis', 'Necrobiosis', 'Manual_granuloma', 'Leiden_granuloma')
# gsea_clusters_heatmap <- gsea_get_heatmap(granulomas_compare_reactome, file_name = 'all_granulomas_reactome', pathway_cutoff = 50)
# 
# ## Disease-specific ####
# granulomas_compare_go <- list(gsea_ga_gvsd_go, gsea_sa_gvsd_go, gsea_nl_gvsd_go)
# names(granulomas_compare_go) <- c('G.annulare', 'Sarcoidosis', 'Necrobiosis')
# gsea_clusters_heatmap <- gsea_get_heatmap(granulomas_compare_go, file_name = '3_granulomas_go', pathway_cutoff = 50)
# 
# granulomas_compare_kegg <- list(gsea_ga_gvsd_kegg, gsea_sa_gvsd_kegg, gsea_nl_gvsd_kegg)
# names(granulomas_compare_kegg) <- c('G.annulare', 'Sarcoidosis', 'Necrobiosis')
# gsea_clusters_heatmap <- gsea_get_heatmap(granulomas_compare_kegg, file_name = '3_granulomas_kegg', pathway_cutoff = 50)
# 
# granulomas_compare_reactome <- list(gsea_ga_gvsd_reactome, gsea_sa_gvsd_reactome, gsea_nl_gvsd_reactome)
# names(granulomas_compare_reactome) <- c('G.annulare', 'Sarcoidosis', 'Necrobiosis')
# gsea_clusters_heatmap <- gsea_get_heatmap(granulomas_compare_reactome, file_name = '3_granulomas_reactome', pathway_cutoff = 50)
# 
