# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ClusterProfiler tutorial
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

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(clusterProfiler) # for PEA analysis
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

# Set relevant paths
list.files()
in_path <- "C:/Users/laura/Documents/Biostatsquid/Projects/data/scRNAseq/"
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
                            pathway_name = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(res_df)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df <- res_df[res_df$ID %in% target_pws, ]

print('Saving clusterprofiler results')
write.csv(res_df, paste0(filename, '_resclusterp.csv'), row.names = FALSE)

##

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
