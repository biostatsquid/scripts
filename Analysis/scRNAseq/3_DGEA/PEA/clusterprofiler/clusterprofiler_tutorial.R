# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ClusterProfiler tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Version: 2.0

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
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(msigdbr)

# Set relevant paths
list.files()
in_path <- "C:/Users/laura/Documents/Biostatsquid/Projects/data/scRNAseq"
out_path <- "C:/Users/laura/Documents/Biostatsquid/Projects/PEA/clusterProfiler"
dir.create(out_path, showWarnings = F, recursive = T)

# Read in data ===================================================
list.files(in_path)
df <- read.csv(file.path(in_path, 'severevshealthy_degresults.csv'), row.names = 1)

# OVER-REPRESENTATION ANALYSIS (ORA) with ClusterProfiler ===================================================

## 1. Prepare pathways or gene sets -----------------------------------------------
# msigdbr R package provides Molecular Signatures Database (MSigDB) gene sets
# Read more here: https://igordot.github.io/msigdbr/articles/msigdbr-intro.html
# View available collections
# msigdbr_collections()
# Set gene set collection
name_of_gene_set <- 'H' # for hallmark gene sets
# Other options:
# 'C1' - positional gene sets
# 'C2' - curated gene sets (includes KEGG, Reactome, BioCarta, etc.)
# 'C3' - regulatory target gene sets  
# 'C4' - computational gene sets
# 'C5' - ontology gene sets (GO terms)
# 'C6' - oncogenic signature gene sets
# 'C7' - immunologic signature gene sets
# 'C8' - cell type signature gene sets
gene_sets_df <- msigdbr(species = 'Homo sapiens', category = name_of_gene_set)
#gene_sets_df <- msigdbr(species = 'Homo sapiens', category = 'C5', subcollection = 'GO:BP')
head(gene_sets_df)

print("Sample of gene sets:")
print(head(gene_sets_df))

# Convert to named list format required by clusterprofiler: the first column should have pathway names
# the second column should have gene symbols
gene_sets <- gene_sets_df %>%
  dplyr::select(gs_name, gene_symbol)

cat(paste("Loaded", length(unique(gene_sets$gs_name)), "gene sets\n"))

## 2. Prepare deg data -----------------------------------------------

# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
  
))
# Remove non-significant genes for ORA
df_sig <- df[df$diffexpressed != 'NO', ] # remove non-significant genes
cat(paste("Found", nrow(df_sig), "significantly differentially expressed genes\n"))
cat(paste("UP:", sum(df_sig$diffexpressed == 'UP'), "\n"))
cat(paste("DOWN:", sum(df_sig$diffexpressed == 'DOWN'), "\n"))

# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
deg_results_list <- split(df_sig, df_sig$diffexpressed)

## 4. Prepare background genes ----------------------------------------
# These should be ALL the genes you could have measured in your experiment,
# to avoid bias. I will use all the genes from the original dataset we used
# to compare severe vs healthy cells
data <- scRNAseq::BacherTCellData(filtered = TRUE, ensembl = FALSE, location = TRUE)
background_genes <- row.names(data)
background_genes[1:5]
rm(data) # to clear some space

## 5. Run ClusterProfiler -----------------------------------------------

# Settings
name_of_comparison <- 'severevshealthy'
padj_cutoff <- 0.05
genecount_cutoff <- 5
filename_prefix <- paste0(name_of_comparison, '_', name_of_gene_set)

# Run clusterProfiler on each sub-dataframe
res <- lapply(names(deg_results_list),
              function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
                                   TERM2GENE = gene_sets,
                                   universe = background_genes,
                                   # Let's use the cut-offs later to filter ORA results
                                   #pvalueCutoff = padj_cutoff,
                                   #minGSSize = genecount_cutoff,
                                   qvalueCutoff = 0.2))
names(res) <- names(deg_results_list)

# Check results
cat("Enrichment analysis results:\n")
for(i in names(res)) {
  cat(paste(i, ":", nrow(as.data.frame(res[[i]])), "enriched pathways\n"))
}

# Process enrichment results into dataframe
res_df_sig <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df_sig) <- names(res)
res_df_sig <- do.call(rbind, res_df_sig)
head(res_df_sig)
res_df_sig <- res_df_sig %>% mutate(minuslog10padj = -log10(p.adjust))
# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
res_df_sig <- res_df_sig %>% 
  dplyr::filter(p.adjust < padj_cutoff & Count > genecount_cutoff) # select only target pathways have p adjusted < 0.05 and at least 6 genes

print('Saving clusterprofiler results')
write.csv(res_df_sig, file.path(out_path, paste0(filename_prefix, '_resclusterp.csv')), row.names = FALSE)

# Also save the enrichresults object for plotting later
qs::qsave(res, file.path(out_path, paste0(filename_prefix, '_resclusterp_enrichres.qs')))
# You can also save it as RDS

##

# Visualisations -----------------------------------------------
# Visualize ORA results for UPregulated genes
results <- res$UP
# Most common method to visualize enriched terms. Shows enrichment scores (e.g. p values) and gene count or ratio as bar height and color 
p1 <- barplot(results, showCategory = 10, 
              title = paste("MSigDB", name_of_gene_set, "Enrichment - upregulated genes"))
p1
p2 <- mutate(results, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")
p2
# Dotplot
# Dot plot is similar to bar plot with the capability to encode another score as dot size.
p3 <- dotplot(results, showCategory = 15)
p3
# cnetplot
# Shows the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network.
p4 <- cnetplot(results, categorySize = "pvalue")
print(p4)
# Heatplot
# The heatplot is similar to cnetplot, while displaying the relationships as a heatmap. 
# useful when there is a large number of significant terms
p5 <- heatplot(results, showCategory = 12)
p5

# Treeplot
# The treeplot() function performs hierarchical clustering of enriched terms. I
enrichres2 <- pairwise_termsim(results) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
p6 <- treeplot(enrichres2)

# Enrichment map 
# organizes enriched pathways into a network with edges connecting overlapping gene sets. 
# In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module.
p7 <- emapplot(enrichres2)
p7
# Upsetplot
# alternative to cnetplot. emphasizes the gene overlapping among different gene sets.
p8 <- upsetplot(results)
p8

plot_list <- list(p1, p2, p3, p4, p5, p6, p7, p8)
# Save all plots
for(i in 1:length(plot_list)){
  png(file.path(out_path, paste0('clusterprofiler.plot', i, '.png')), width = 480*1.5, height = 520 * 1.5)
  print(plot_list[i])
  
  dev.off()
}

# Using other clusterProfiler functions ===============================
# GO Enrichment Analysis
go_bp <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",           # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)      # Convert IDs to symbols

go_mf <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",           # Molecular Function
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

go_cc <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",           # Cellular Component
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# KEGG Pathway Analysis
kegg_result <- enrichKEGG(gene = entrez_ids$ENTREZID,
                          organism = 'hsa',      # Human
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# Reactome Pathway Analysis
reactome_result <- enrichPathway(gene = entrez_ids$ENTREZID,
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.2,
                                 readable = TRUE)
