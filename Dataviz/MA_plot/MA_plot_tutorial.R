# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MA plots in R tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Version: 1.0
# Description: MA plots with ggplot2

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Install libraries if needed
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# # For basic plotting
# install.packages(c("ggplot2", "dplyr", "gridExtra"))
# 
# # For bioinformatics applications
# BiocManager::install(c("limma", "edgeR", "DESeq2"))
# 
# # Additionally, we'll use the dataset "airway" for our analysis
# # https://bioconductor.org/packages/release/data/experiment/html/airway.html
# BiocManager::install('airway')


# Loading relevant libraries 
library(tidyverse)
library(gridExtra)
library(limma)
library(edgeR)
library(DESeq2)
library(airway)

# Set seed for reproducibility
set.seed(123)

# Get sample expression data =====================================
data("airway")
class(airway)
airway@colData %>% head()
table(airway$cell, airway$dex)

# Create DESeqDataSet
dds <- DESeqDataSet(airway, design = ~ dex)
# Run DESeq() to estimate size factors and normalize
dds <- DESeq(dds) 

# MA plot single sample =====================================
# Filter for N052611 cell line only
norm_counts <- counts(dds, normalized = TRUE)
norm_counts[1:5,1:5]

# Find the treated and untreated samples for N052611
sample_data <- colData(dds)

# For N052611, find treated and untreated samples (you can also figure it out by looking at the metadata directly)
n052611_samples <- rownames(sample_data)[sample_data$cell == "N052611"]
untreated_sample <- n052611_samples[sample_data[n052611_samples, "dex"] == "untrt"]
treated_sample <- n052611_samples[sample_data[n052611_samples, "dex"] == "trt"]

# Extract the counts as a dataframe
expr_df <- as.data.frame(norm_counts[, c(untreated_sample, treated_sample)])
head(expr_df)

# Let's change the column names to the dex column (treated or untreated) instead of the RunID to make it clearer
# Careful when replacing column names! Make sure you are matching the right name. 
expr_df <- expr_df %>% dplyr::rename(untrt = untreated_sample,
                                     trt = treated_sample)
head(expr_df)

# Calculate M and A values
M <- log1p(expr_df$trt) - log1p(expr_df$untrt)  # Log fold change
A <- (log1p(expr_df$untrt) + log1p(expr_df$trt)) / 2  # Average log expression

# Create data frame
ma_data <- data.frame(A = A,
                      M = M,
                      Gene = row.names(expr_df),
                      Differentially_expressed = abs(M) > 1)
head(ma_data)
table(ma_data$Differentially_expressed)
# Looks like we have 13285 different (abs(M) > 1) genes!

# Basic MA plot
ggplot(ma_data, aes(x = A, y = M)) +
  geom_point()

# Prettier MA plot
p1 <- ggplot(ma_data, aes(x = A, y = M)) +
  geom_point(aes(color = Differentially_expressed), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "violet")) +
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed", linewidth = 1) +
  labs(title = "MA Plot: treated vs untreated (N052611)",
       x = "A",
       y = "M",
       color = "High differential expression") +
  theme_bw() +
  theme(legend.position = "top")

# For example, let's label top 20 genes by fold change
genes_to_label <- ma_data %>%
  dplyr::filter(Differentially_expressed == TRUE) %>%
  dplyr::arrange(desc(abs(M))) %>%
  dplyr::slice(1:20)
print(genes_to_label)

# Now, use geom_text_repel to highlight these genes
p1 <- ggplot(ma_data, aes(x = A, y = M)) +
  geom_point(aes(color = Differentially_expressed), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "violet")) +
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed", linewidth = 1) +
  labs(title = "MA Plot: treated vs untreated (N052611)",
       x = "A",
       y = "M",
       color = "High differential expression") +
  ggrepel::geom_text_repel(data = genes_to_label,
                           aes(label = Gene),
                           size = 3, color = "black",
                           max.overlaps = Inf, box.padding = 0.4,
                           point.padding = 0.3, segment.color = "grey50"
  ) +
  theme_bw() +
  theme(legend.position = "top")

print(p1)

# MA plot multisample sample =====================================
# Get results comparing treated vs untreated across all samples
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
head(res_df)
# Remove NAs and infinite values
res_df <- res_df[complete.cases(res_df), ]

head(res_df)

# Add significance categories
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Significant (padj < 0.05, |FC| > 2)"
res_df$significance[res_df$padj < 0.05 & abs(res_df$log2FoldChange) <= 1] <- "Significant (padj < 0.05, |FC| ≤ 2)"
table(res_df$significance)

res_df <- res_df %>% mutate(significance_2 = case_when(
  res_df$padj < 0.05 & res_df$log2FoldChange > 1 ~ 'significantly upregulated (padj < 0.05, FC > 2)',
  res_df$padj < 0.05 & res_df$log2FoldChange < -1 ~ 'significantly upregulated (padj < 0.05, FC < -2)',
  .default = 'not significant'
))
table(res_df$significance_2)

p1 <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significance_2)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_x_log10() +
  # scale_color_manual(
  #   values = c("Not significant" = "gray60",
  #     "Significant (padj < 0.05, |FC| ≤ 2)" = "lightblue",
  #     "Significant (padj < 0.05, |FC| > 2)" = "darkblue")
  # ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = c(-1, 1), color = "blue", linetype = "dotted", linewidth = 0.5) +
  labs(title = "MA Plot: All Treated vs All Untreated Samples",
    x = "Mean of normalized counts (log10 scale)",
    y = "log2 fold change (Treated/Untreated)",
    color = "Differential Expression") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))


print(p1)


## Volcano plot
res_df <- res_df %>%
  mutate(significant = pvalue < 0.05)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values = c("TRUE" = "#D60093", "FALSE" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 50)) +
  theme_bw() +
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~p~value),
    title = "Volcano Plot (DESeq2 results)",
    color = "Significant (p < 0.05)"
  )

