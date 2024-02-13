# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PEA tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Version: 1.0

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(123456)

# Set project library
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(mgsa) # for MGSA analysis
library(fgsea) # for GSEA analysis
library(clusterProfiler) # for PEA analysis
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

# Set relevant paths
getwd()
list.files()
in_path <- "Datasets/"
out_path <- "PEA/Results/"
bg_path <- "PEA/Background_genes/"

# Functions ------------------------------------------------
## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Read in data ===================================================
list.files(in_path)
df <- read.csv(paste0(in_path, 'severevshealthy_degresults.csv'), row.names = 1)

# ClusterProfiler ===================================================

## Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# # For clusterProfiler
# # Filter out the gmt files for KEGG, Reactome and GOBP
# genes_in_data <- df$gene_symbol
# gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE) # cluster_profiler needs a special parsing of the .gmt file, we cannot use pwl1 we created for MGSA analysis
# for (file in gmt_files){
#   #file <- gmt_files[1]
#   pwl2 <- read.gmt(file) # all genes from db (GOBP, Reactome, KEGG...)
#   pwl2 <- pwl2[pwl2$gene %in% genes_in_data,] # background gene sets (pathways) already filtered for those genes (15107) that are present in our data to avoid bias
# 
#   filename <- paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
#   saveRDS(pwl2, filename)
# }

## Prepare deg data -----------------------------------------------

# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
  
))

df <- df[df$diffexpressed != 'NO', ] # remove non-significant genes

# # Substitute names so they are annotated nicely in the heatmap later
# df$diffexpressed <- gsub('DOWN', 'Healthy', gsub('UP', 'Severe', df$diffexpressed))
# unique(df$diffexpressed)

# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
deg_results_list <- split(df, df$diffexpressed)

## Run ClusterProfiler -----------------------------------------------

# Settings
name_of_comparison <- 'severevshealthy'
background_genes <- 'reactome'
bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))
padj_cutoff <- 0.05
genecount_cutoff <- 5
filename <- paste0(out_path, 'clusterProfiler/', name_of_comparison, '_', background_genes) 
  
# Read in background genes (alternative)
if(background_genes == 'KEGG'){
  bg_genes <- readRDS(paste0(bg_path, 'kegg.RDS'))
} else if(background_genes == 'reactome'){
  bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))
} else if(background_genes == 'go.bp'){
  bg_genes <- readRDS(paste0(bg_path, 'go.bp.RDS'))
} else {
  stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
}

# Run clusterProfiler on each sub-dataframe
res <- lapply(names(deg_results_list),
                 function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
                                      TERM2GENE = bg_genes))
names(res) <- names(deg_results_list)

#Convert the go_all list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)

res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                                  diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(res_df)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df <- res_df[res_df$ID %in% target_pws, ]

print('Saving clusterprofiler results')
write.csv(res_df, paste0(filename, '_resclusterp.csv'), row.names = FALSE)

# ## Get heatmap -----------------------------------------------
# 
# # First convert the dataframe into the heatmap dataframe we need to build the heatmap
# heatmap_df <- tidyr::pivot_wider(res_df[,c('ID', 'minuslog10padj', 'diffexpressed')], values_from = minuslog10padj, names_from = diffexpressed, values_fill = 0)
# heatmap_df <- as.data.frame(heatmap_df)
# rownames(heatmap_df) <- heatmap_df$ID # to be able to set them as row names
# heatmap_df <- heatmap_df[!grepl('ID', colnames(heatmap_df))] # now remove the pathways column
# head(heatmap_df)
# #write.csv(heatmap_df, paste(path, 'clusterProfiler/heatmap_df.csv', sep = '/'), row.names = TRUE)
# 
# pathway_cutoff <- 30
# if(nrow(heatmap_df) > pathway_cutoff){
#   selected_pathways <- head(unique(res_df[order(-res_df$minuslog10padj),'ID']), pathway_cutoff) # order by minuslogpadj and get ones with larger minuslogpadj
#   heatmap_df <- subset(heatmap_df, rownames(heatmap_df) %in% selected_pathways) # subset to less pathways for visualisation purposes if there are too many
# }
# 
# ## Get heatmap 
# heatmap <- pheatmap(heatmap_df,
#                     show_colnames = T, show_rownames = T,
#                     cluster_rows = FALSE, cluster_cols = FALSE,
#                     clustering_distance_cols = 'euclidean',
#                     clustering_distance_rows = 'euclidean',
#                     clustering_method = 'ward.D',
#                     fontsize_col = 10,
#                     fontsize_row = 10,
# )
# 
# # Save it
# pdf(file = paste0(filename, '_clusterp_heatmap.pdf'), width = 15, height = 15)
# grid::grid.draw(heatmap$gtable)
# dev.off()

## Visualisations -----------------------------------------------
res_df <- read.csv(paste0(out_path, 'clusterProfiler/', 'severevshealthy_reactome_resclusterp.csv'))
bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))
# Convert clusterProfiler object to a new "enrichResult" object
# Select only upregulated genes in Severe
res_df <- res_df %>% filter(diffexpressed == 'UP') %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df) <- res_df$ID
# What the column names should be
colnames(res_df)

# For visualisation purposes, let's shorten the pathway names
res_df$Description <- gsub('(H|h)iv', 'HIV', 
                           gsub('pd 1', 'PD-1',
                                gsub('ecm', 'ECM', 
                                     gsub('(I|i)nterleukin', 'IL', 
                                          gsub('(R|r)na', 'RNA', 
                                               gsub('(D|d)na', 'DNA',
                                                    gsub(' i ', ' I ', 
                                                         gsub('(A|a)tp ', 'ATP ', 
                                                              gsub('(N|n)adh ', 'NADH ', 
                                                                   gsub('(N|n)ad ', 'NAD ',
                                                                        gsub('t cell', 'T cell',
                                                                             gsub('b cell', 'B cell',
                                                                                  gsub('built from .*', ' (...)',
                                                                                       gsub('mhc', 'MHC',
                                                                                            gsub('mhc class i', 'MHC I', 
                                                                                                 gsub('mhc class ii', 'MHC II', 
                                                                                                      stringr::str_to_sentence(
                                                                                                        gsub('_', ' ',  
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df$Description)))))))))))))))))))



# DE_genes is a vector of the names of your DE genes
# universe_vector is a vector of the names of all genes that are annotated with GO terms (your enrichment universe). 
# geneSets is a named list, where the names are enriched GO terms and the elements are DE genes annotated with that GO term.
enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_df,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = df$gene_symbol,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres)

# Barplot

# Most common method to visualize enriched terms. Shows enrichment scores (e.g. p values) and gene count or ratio as bar height and color 
p1 <- barplot(enrichres, showCategory = 20) 
p2 <- mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")

# Dotplot
# Dot plot is similar to bar plot with the capability to encode another score as dot size.
p3 <- dotplot(enrichres, showCategory = 15) + ggtitle("Severe vs Healthy")

# cnetplot
# Shows the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network.
# ## Convert gene ID to Symbol
# enrichres_eid <- setReadable(enrichres, 'org.Hs.eg.db', 'SYMBOL')
p4 <- cnetplot(enrichres)

# Heatplot
# The heatplot is similar to cnetplot, while displaying the relationships as a heatmap. 
# useful when there is a large number of significant terms
p5 <- heatplot(enrichres, showCategory = 5)

# Treeplot
# The treeplot() function performs hierarchical clustering of enriched terms. I
enrichres2 <- pairwise_termsim(enrichres) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
p6 <- treeplot(enrichres2)

# Enrichment map 
# organizes enriched pathways into a network with edges connecting overlapping gene sets. 
# In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module.
p7 <- emapplot(enrichres2)

# Upsetplot
# alternative to cnetplot. emphasizes the gene overlapping among different gene sets.
enrichres@result[["Description"]] <- gsub('Respiratory electron transport ATP synthesis by.*', 'Respiratory electron transport ATP synthesis by...', enrichres@result[["Description"]])
p8 <- upsetplot(enrichres)


plot_list <- list(p1, p2, p3, p4, p5, p6, p7, p8)
# Save all plots
new_path <- 'C:/Users/laura/Documents/Biostatsquid/Tutorials/PEA/clusterProfiler/'
for(i in 1:length(plot_list)){
  png(paste0(new_path, 'plot', i, '.png'), width = 480*1.5, height = 520 * 1.5)
  print(plot_list[i])
  
  dev.off()
}


# GSEA ===================================================

## Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# # For MGSA and GSEA
# # Filter out the gmt files for KEGG, Reactome and GOBP
# genes_in_data <- df$gene_symbol
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
for (file in gmt_files){
  # for debug
  #file <- gmt_files[1]

  # Read in gmt file
  gmt <- gmtPathways(file)
  hidden <- unique(unlist(gmt))
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  Ass <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(Ass)[2]){
    Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
  }

  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  Ass1 <- Ass[hidden1, colnames(Ass)[which(colSums(Ass[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  pwl1 <- matrix_to_list(Ass1) # for this we use the function we previously defined
  saveRDS(pwl1, file = paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RData', sep = ''))
}

#pwl1 <- readRDS(paste0(bg_path, 'kegg.RData')) # background gene sets (pathways) already filtered for those genes that are present in our data to avoid bias + gene sets with more than 5 genes annotated

## Prepare deg data -----------------------------------------------

# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
  
))

df <- df[df$diffexpressed != 'NO', ] # remove non-significant genes

# Substitute names
df$diffexpressed <- gsub('DOWN', 'Healthy', gsub('UP', 'Severe', df$diffexpressed))
unique(df$diffexpressed)

## Run GSEA -----------------------------------------------------------------

# Settings
background_genes <- 'go.bp'
bg_genes <- readRDS(paste0(bg_path, 'go.bp.RData'))
name_of_comparison <- 'severevshealthy'
number_of_top_pathways_up <- 25
number_of_top_pathways_down <- 5
filename <- paste0(out_path, 'GSEA/', name_of_comparison, '_', background_genes) 
  
# Read in background genes (alternative)
if(background_genes == 'kegg'){
  bg_genes <- readRDS(paste0(bg_path, 'kegg.RData'))
} else if(background_genes == 'reactome'){
  bg_genes <- readRDS(paste0(bg_path, 'reactome.RData'))
} else if(background_genes == 'go.bp'){
  bg_genes <- readRDS(paste0(bg_path, 'go.bp.RData'))
} else {
  stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
}

# Prepare input
# First, format the background genes
# bg_genes$term <- as.character(bg_genes$term) # fgsea doesn't like it when the entries are factors
# bg_genes$term <- gsub('GOBP_', '', gsub('KEGG_', '', gsub('REACTOME_', '', bg_genes$term)))

## FGSEA needs a named vector with the rankings and the genes as names
## We will use  signed p value as ranking: sign(df$log2fc)*(-log10(df$pval)))
# Note that you can also just use the log2FC - but it will not take into account 
# genes with a large FC but non-significant
rankings <- sign(df$log2fc)*(-log10(df$pval)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df$gene_symbol # genes as names
plot(rankings)
max(rankings)
# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])

rankings <- replace(rankings, rankings > max_ranking, max_ranking * 5)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 5)

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

# Easy peasy! Run fgsea with the pathways
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1)

saveRDS(GSEAres, file = paste0(filename, '_gsea_results.RDS'))
#data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))

# Check top pathways
topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
dev.off()

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres[order(pval)][padj < 0.01], bg_genes, rankings)
mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)
dev.off()



# MGSA ===================================================

## Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# # For MGSA and GSEA
# # Filter out the gmt files for KEGG, Reactome and GOBP
# genes_in_data <- df$gene_symbol
# gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE) 
# for (file in gmt_files){
#   # for debug
#   #file <- gmt_files[1]
#   
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
#   saveRDS(pwl1, file = paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RData', sep = ''))
# }

#pwl1 <- readRDS(paste0(bg_path, 'kegg.RData')) # background gene sets (pathways) already filtered for those genes that are present in our data to avoid bias + gene sets with more than 5 genes annotated

## Prepare deg data -----------------------------------------------

# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
  
))

df <- df[df$diffexpressed != 'NO', ] # remove non-significant genes

# Substitute names
df$diffexpressed <- gsub('DOWN', 'Healthy', gsub('UP', 'Severe', df$diffexpressed))
unique(df$diffexpressed)

## Run GSEA -----------------------------------------------------------------

# Settings
background_genes <- 'kegg'
bg_genes <- readRDS(paste0(bg_path, 'kegg.RData'))
name_of_comparison <- 'severevshealthy'
number_of_top_pathways_up <- 25
number_of_top_pathways_down <- 5
filename <- paste0(out_path, 'MGSA/', name_of_comparison, '_', background_genes) 

# Read in background genes (alternative)
if(background_genes == 'kegg'){
  bg_genes <- readRDS(paste0(bg_path, 'kegg.RData'))
} else if(background_genes == 'reactome'){
  bg_genes <- readRDS(paste0(bg_path, 'reactome.RData'))
} else if(background_genes == 'go.bp'){
  bg_genes <- readRDS(paste0(bg_path, 'go.bp.RData'))
} else {
  stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
}


##  Run mgsa for subset of genes in df
mgsa_fit_list <- mgsa(df$gene_symbol, bg_genes)
saveRDS(mgsa_fit_list, file = paste0(filename, '_mgsa_fit_list.RData'))
print('MGSA fit finished!')

## Get the enriched pathways 
res_list <- setsResults(mgsa_fit_list)
res_list$pathways <- gsub('^.*GOBP_|^.*KEGG_|^.*REACTOME_', '', rownames(res_list))
saveRDS(res_list, file = paste0(filename, '_mgsa_res_list.RData'))

target_pathways <- unique(res_list$pathways[res_list$estimate > 0.5]) # target pathways have an estimate higher than 0.5
res_all_df_filtered <- res_list[res_list$pathways %in% target_pathways, ]
write.csv(res_all_df_filtered, paste0(filename, '_mgsa_res_filtered.csv'))

