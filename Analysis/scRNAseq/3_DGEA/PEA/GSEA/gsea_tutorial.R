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

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)
library(msigdbr)

# Set relevant paths
project_path <- "C:/Users/laura/Documents/Biostatsquid/Projects"
in_path <- file.path(project_path, 'data/scRNAseq') # DGE results path
out_path <- file.path(project_path, "PEA/Results") # GSEA output path
dir.create(out_path, showWarnings = F, recursive = T)

# Analysis ====================================================

## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(file.path(in_path, 'severevshealthy_degresults.csv'), row.names = 1)
head(df)
# gene_symbol      pval      padj               log2fc
# 1 MIR1302-2HG 0.9998587 0.9999881 15.86284051374880022
# 2     FAM138A 0.9999881 0.9999881 -0.00000000000283586
# 3       OR4F5 0.9999881 0.9999881 -0.00000000000283586

## 2. Prepare background genes -----------------------------------------------
# msigdbr R package provides Molecular Signatures Database (MSigDB) gene sets
# Read more here: https://igordot.github.io/msigdbr/articles/msigdbr-intro.html
# msigdbr_collections()
gene_sets_df <- msigdbr(species = 'Homo sapiens', category = 'H')
#gene_sets_df <- msigdbr(species = 'Homo sapiens', category = 'C5', subcollection = 'GO:BP')
head(gene_sets_df)

# Convert to named list format required by fgsea
gene_sets <- gene_sets_df %>%
  split(x = .$gene_symbol, f = .$gs_name)

cat(paste("Loaded", length(gene_sets), "gene sets\n"))

## 3. Prepare ranked list of genes -----------------------------------------------
## FGSEA needs a named vector with the rankings and the genes as names
## We will use  signed p value as ranking: sign(df$log2fc)*(-log10(df$pval)))
# Note that you can also just use the log2FC only - but it will not take into account 
# genes with a large FC but non-significant
rankings <- sign(df$log2fc)*(-log10(df$pval)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df$gene_symbol # genes as names#
head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)
# Handle infinite values (very low p-values)
# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)

# Plot Top genes
top_genes <- head(rankings, 50)
ggplot(data.frame(gene_symbol = names(top_genes), ranks = top_genes), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
fgsea_res <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

# Add gene set descriptions
gene_set_info <- gene_sets_df %>%
  select(gs_name, gs_description) %>%
  distinct()

fgsea_res <- fgsea_res %>%
  left_join(gene_set_info, by = c("pathway" = "gs_name"))

head(fgsea_res)

## 5. Save the results -----------------------------------------------
filename_prefix <- 'severevshealthy.hallmark'
write.table(fgsea_res, file.path(out_path, paste0(filename_prefix, '_results.csv')), sep = ',', col.names = T, row.names = F)

## 6. Check results ------------------------------------------------------
# Select top pathways by absolute NES
top_pathways <- fgsea_res %>%
  arrange(padj) %>%
  head(30) %>% 
  as.data.frame()
p1 <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = -log10(padj))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
  scale_size_continuous(range = c(2, 8), name = "Gene Set Size") +
  theme_minimal() +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway") +
  theme(axis.text.y = element_text(size = 8))

# Top 6 enriched pathways (ordered by p-val)
head(fgsea_res[order(pval), ])
# Number of significant pathways at padj < 0.01
sum(fgsea_res[, padj < 0.01])
sum(fgsea_res[, pval < 0.05])

## Check top pathways
number_of_top_pathways_up <- 25
number_of_top_pathways_down <- 5

topPathwaysUp <- fgsea_res[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- fgsea_res[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = file.path(out_path, 'gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(gene_sets[topPathways], stats = rankings,  
              fgseaRes = fgsea_res, gseaParam = 0.5)
#dev.off()

png(file = file.path(out_path, paste0(filename_prefix, '_gsea_top30pathways.png')), width = 1500, height = 800)
plotGseaTable(gene_sets[topPathways], stats = rankings, fgseaRes = fgsea_res, gseaParam = 0.5)
dev.off()

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(fgsea_res[order(pval)][pval < 0.05], gene_sets, rankings)
mainPathways <- fgsea_res[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
#pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
plotGseaTable(gene_sets[mainPathways], rankings, fgsea_res, gseaParam = 0.5)
#dev.off()

png(file = file.path(out_path, paste0(filename_prefix, '_gsea_mainpathways.png')), width = 1500, height = 800)
plotGseaTable(gene_sets[mainPathways], rankings, fgsea_res, gseaParam = 0.5)
dev.off()

# plot the most significantly enriched pathway
plotEnrichment(gene_sets[[head(fgsea_res[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(fgsea_res[order(padj), ], 1)$pathway)

plotEnrichment(gene_sets[['REACTOME_FORMATION_OF_FIBRIN_CLOT_CLOTTING_CASCADE']],
               rankings) + 
  labs(title = 'Reactome pathway: Formation of fibrin clot / clotting cascade') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 


# Plotting several results ================================================





