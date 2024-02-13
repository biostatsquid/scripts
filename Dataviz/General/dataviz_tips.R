# ----------------------
# Dataviz
# ----------------------

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set project library
#.libPaths('C:/Users/laura/Documents/Biostatsquid/Scripts/R4.2.3')

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot

# Importing dataset ===================================================
data(iris)

# 1. Creating and using a custom theme ===================================================
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point()

ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() + 
  theme_classic() 

ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() + 
  theme_classic() + 
  theme(axis.text.x = element_text(colour = 'blue', size = 15),
        axis.line.y = element_blank(),
        legend.position = 'left')

biostatsquid_theme <- theme(plot.title = element_text(size = rel(2)),
                            panel.grid.major.y = element_line(colour = 'gray'),
                            panel.grid.minor.y = element_line(colour = 'gray'),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            plot.background = element_rect(fill = NULL, colour = 'white'),
                            panel.background = element_rect(fill = 'white'),
                            # Axis stuff
                            axis.line = element_line(colour = 'black', linewidth = 1),
                            axis.text = element_text(colour = "black", face = 'bold'),
                            axis.text.x = element_text(size = rel(1)),
                            axis.text.y = element_text(size = rel(1)),
                            axis.title = element_text(size = rel(1.2)),
                            axis.ticks = element_line(colour = 'black', linewidth = 1.2),
                            # Legend stuff
                            legend.position = "bottom",
                            legend.margin = margin(6, 6, 6, 6),
                            legend.title = element_text(face = 'bold'),
                            legend.background = element_blank(),
                            legend.box.background = element_rect(colour = "black"))

ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() + 
  biostatsquid_theme

theme_set(biostatsquid_theme)
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() +
  theme(axis.title.x = element_text(colour = 'orange'))

# 2. Combine plots ===================================================
library(cowplot)
p1 <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() 
p2 <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, fill = Species)) +
  geom_col() 

cowplot::plot_grid(p1, p2, ncol = 1, nrow = 2,
                   rel_heights = c(2, 3))
cowplot::plot_grid(p1, p2, ncol = 2, nrow = 1,
                   rel_widths = c(2, 3))

cowplot::plot_grid(p1, p2, ncol = 2, nrow = 1,
                   rel_widths = c(5, 3), labels = c('A.', 'B.'),
                   label_fontfamily = "Times", label_colour = "black")

library(patchwork)
plot_layout(p1 + p2 + gridExtra::tableGrob(iris[1:10, 1:2]), widths = c(2, 4, 2))

library(ggpubr)
ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = 'bottom')

# 3. Save plots ===================================================
png('C:/Users/laura/Documents/Biostatsquid/Plots/super_plot.png')
print(p1)
dev.off()

# Set the output directory where the plots will be saved to.
out_dir <-  "C:/Users/laura/Documents/Biostatsquid/Plots/"
dir.create(out_dir, showWarnings = FALSE)
# Let's be fancy and make the height and width depend on the number of columns and rows
# we want our plots distributed in
number_cols <- 2
number_rows <- 1
png(filename = paste0(out_dir, 'super_plot.png'), height = 750 * number_rows, width = 750 * number_cols)
print(p1 + p2)
dev.off()

number_cols <- 1
number_rows <- 1
combined_plot <- ggarrange(p1, p2, ncol = number_cols, nrow = number_rows, common.legend = TRUE, legend = 'bottom')
png(filename = paste0(out_dir, 'super_plot2.png'), height = 750 * number_rows, width = 750 * number_cols)
print(combined_plot)
dev.off()

ggsave(combined_plot, paste0(out_dir, 'super_plot.png'))










# 4. Italics and bold ===================================================

# 5. Colours ===================================================



https://epirhandbook.com/en/ggplot-tips.html