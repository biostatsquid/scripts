# ----------------------
# GSEA tutorial
# ----------------------

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(123456)

# Set project library
.libPaths('C:/Users/laura/Documents/Biostatsquid/Scripts/R4.2.3')
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)

# Set relevant paths
list.files()
in_path <- "Datasets/"
out_path <- "PEA/Results/"
bg_path <- "PEA/Background_genes/"

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  # gmt_file <- gmt_files[1]
  # genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# Analysis ====================================================

## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'severevshealthy_degresults.csv'), row.names = 1)

## 2. Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$gene_symbol
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
bg_genes <- prepare_gmt(gmt_files[2], my_genes, savefile = FALSE)

## 3. Prepare ranked list of genes -----------------------------------------------
## FGSEA needs a named vector with the rankings and the genes as names
## We will use  signed p value as ranking: sign(df$log2fc)*(-log10(df$pval)))
# Note that you can also just use the log2FC - but it will not take into account 
# genes with a large FC but non-significant
rankings <- sign(df$log2fc)*(-log10(df$pval)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df$gene_symbol # genes as names#
head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)
# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)
# ggplot(data.frame(gene_symbol = names(rankings), ranks = rankings), aes(gene_symbol, ranks)) + geom_point() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + geom_point() +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

## 5. Save the results -----------------------------------------------
name_of_comparison <- 'severevshealthy'
background_genes <- 'gobp'
filename <- paste0(out_path, 'GSEA/', name_of_comparison, '_', background_genes) 
# saveRDS(GSEAres, file = paste0(filename, '_gsea_results.RDS'))
# data.table::fwrite(GSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))


## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])
# Number of significant pathways at padj < 0.01
sum(GSEAres[, padj < 0.01])
sum(GSEAres[, pval < 0.05])

## Check top pathways
number_of_top_pathways_up <- 25
number_of_top_pathways_down <- 5

topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()

png(file = paste0(filename, '_gsea_top30pathways.png'), width = 1500, height = 800)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
dev.off()

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres[order(pval)][pval < 0.05], bg_genes, rankings)
mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
#pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)
#dev.off()

png(file = paste0(filename, '_gsea_mainpathways.png'), width = 1500, height = 800)
plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)
dev.off()

# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)

plotEnrichment(bg_genes[['REACTOME_FORMATION_OF_FIBRIN_CLOT_CLOTTING_CASCADE']],
               rankings) + 
  labs(title = 'Reactome pathway: Formation of fibrin clot / clotting cascade') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 






