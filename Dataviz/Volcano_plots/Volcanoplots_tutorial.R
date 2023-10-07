# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Volcano_plots tutorial
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: `Laura Twomey
# Date: October 2022
# Version: 1.0

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Set input path
path <- "C:/Users/laura/Documents/Biostatsquid/Scripts/Volcano_plots/Tutorial/"
setwd(path)

# Import DGE results
df <- read.csv(paste0(path, 'severevshealthy_degresults.csv'), row.names = 1)

# Create a basic volcano plot
ggplot(data = df, aes(x = log2fc, y = -log10(pval))) +
         geom_point()

# Add threshold lines
ggplot(data = df, aes(x = log2fc, y = -log10(pval))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() + 
  
# Theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

ggplot(data = df, aes(x = log2fc, y = -log10(pval))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 


# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df$diffexpressed[df$log2fc > 0.6 & df$pval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2fc < -0.6 & df$pval < 0.05] <- "DOWN"

# Explore a bit
head(df[order(df$padj) & df$diffexpressed == 'DOWN', ])

ggplot(data = df, aes(x = log2fc, y = -log10(pval), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)

# Edit axis labels and limits
ggplot(data = df, aes(x = log2fc, y = -log10(pval), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis
  
# Note. with coord_cartesian() even if we have genes with p-values or log2FC ourside our limits, they will still be plotted.

ggplot(data = df, aes(x = log2fc, y = -log10(pval), col = diffexpressed, label = gene_symbol)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Thf-like cells in severe COVID vs healthy patients') + # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 
  
# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
df$delabel <- ifelse(df$gene_symbol %in% head(df[order(df$padj), "gene_symbol"], 30), df$gene_symbol, NA)

ggplot(data = df, aes(x = log2fc, y = -log10(pval), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Thf-like cells in severe COVID vs healthy patients') + # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 

myvolcanoplot <- ggplot(data = df, aes(x = log2fc, y = -log10(pval), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Thf-like cells in severe COVID vs healthy patients') + # Plot title 
  geom_text_repel(max.overlaps = Inf) + # To show all labels 
  guides(col = guide_legend(override.aes = aes(label = ''))) # to remove the 'a' from the legend


# Open the file that will contain your plot (the name is up to you)
png(file = "myvolcanoplot.png") # you can change the size of the output file
# Execute the plot
myvolcanoplot
# Close the file that will contain the plot
dev.off()


# You can also save it as png, jpeg, tiff, bmp, svg, ps...
# For more on saving plots in R check https://biocorecrg.github.io/CRG_RIntroduction/with-the-console.html