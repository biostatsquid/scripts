# ---------------------- #
# Stats_visualisations.R
# ---------------------- #
# Author: Laura Twomey
# Date of creation: 12.10.2023
# Version: 1.0

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(42)

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.

# Theme ====================================================================
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

# Load the data =============================================================
df <- data.frame('Cell' = paste(rep('Cell', 30), 1:30, sep = '_'),
                 'TP53_expr' = round(rnorm(30, mean = 40, sd = 2)))
#df <- df %>% mutate(TP53_expr_rel = TP53_expr/sum(TP53_expr))
head(df)
df %>% head()
ggplot(df, aes(x = TP53_expr)) +
  #geom_histogram(aes(y = after_stat(count / sum(count))), fill = '#0BD0D9', col = 'black', size = 0.9) +
  geom_density(size = 1.2, col = '#D60093') +
  ylab('Relative frequency') +
  xlab('TP53 expression levels') +
  scale_x_continuous(breaks = seq(38, 44, 1)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12))

pnorm(41, mean = 40, sd = 2)

density_values <- density(df$TP53_expr)
density_subset <- subset(data.frame(x = density_values$x, y = density_values$y), x <= 41)
ggplot(df, aes(x = TP53_expr)) +
  geom_density(size = 1.2, col = '#D60093') +
  ylab('Relative frequency') +
  xlab('TP53 expression levels') +
  scale_x_continuous(breaks = seq(35, 45, 1), limits = c(35, 45)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12))+
  geom_vline(xintercept = 41, linetype = "dashed", color = "black", size = 1.2) +  # Add a vertical line at x = 40
  geom_area(data = data.frame(x = density_subset$x, y = density_subset$y),
          aes(x = x, y = y), fill = "pink", alpha = 0.5)

probability_between_40_and_41 <- diff(pnorm(c(40, 41), mean = 40, sd = 2))
ggplot(df, aes(x = TP53_expr)) +
  geom_density(size = 1.2, col = '#D60093') +
  ylab('Relative frequency') +
  xlab('TP53 expression levels') +
  scale_x_continuous(breaks = seq(35, 45, 1), limits = c(35, 45)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12)) +
  geom_vline(xintercept = c(40, 41), linetype = "dashed", color = "black", size = 1.2) +
  geom_area(data = subset(density_subset, x > 40 & x <= 41), aes(x = x, y = y), fill = "lightblue", alpha = 0.5) 

probability_between_40_and_41 <- diff(pnorm(c(40, 41), mean = 40, sd = 2))
ggplot(df, aes(x = TP53_expr)) +
  geom_density(size = 1.2, col = '#D60093', fill = 'blue') +
  ylab('Relative frequency') +
  xlab('TP53 expression levels') +
  scale_x_continuous(breaks = seq(35, 45, 1), limits = c(35, 45)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12))
## Skewness

library(gendist)  # For generating skewed distributions

# Parameters for the distributions
size <- 1000
mean <- 0
std_dev <- 1

# Generate random data from a normal distribution
normal_data <- rnorm(size, mean, std_dev)

# Skewness parameters
right_skew <- 5
left_skew <- -5

# Create right skewed data
right_skewed_data <- rskew(size, mean, std_dev, right_skew)

# Create left skewed data
left_skewed_data <- rskew(size, mean, std_dev, left_skew)

# Combine data into a data frame
df <- data.frame(
  group = rep(c("Normal", "Right Skewed", "Left Skewed"), each = size),
  value = c(normal_data, right_skewed_data, left_skewed_data)
)

# Plot density plots
ggplot(df, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~group) +
  labs(title = "Density Plots of Normal and Skewed Distributions",
       x = "Value",
       y = "Density") +
  theme_minimal()


write.csv(df, "C:/Users/laura/Documents/Biostatsquid/Learn/Statistics_basics/density_plots/tp53_norm.csv", row.names = FALSE, quote = FALSE)


df1 <- data.frame('A' = round(rnorm(100, mean = 40, sd = 1)),
                 'B' = round(rnorm(100, mean = 40, sd = 10)),
                 'C' = round(rnorm(100, mean = 42, sd = 5)))
#df <- df %>% mutate(TP53_expr_rel = TP53_expr/sum(TP53_expr))
head(df)
df %>% head()
ggplot(df1, aes(x = A)) +
  #geom_histogram(aes(y = after_stat(count / sum(count))), fill = '#0BD0D9', col = 'black', size = 0.9) +
  geom_density(size = 1.2, col = '#D60093', fill = 'pink', alpha = 0.7) +
  ylab('Relative frequency') +
  xlab('Variable') +
  scale_x_continuous(breaks = seq(15, 75, 5), limits = c(15, 75)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12))

ggplot(df1, aes(x = B)) +
  #geom_histogram(aes(y = after_stat(count / sum(count))), fill = '#0BD0D9', col = 'black', size = 0.9) +
  geom_density(size = 1.2, col = '#D60093', fill = 'pink', alpha = 0.7) +
  ylab('Relative frequency') +
  xlab('Variable') +
  scale_x_continuous(breaks = seq(15, 75, 5), limits = c(15, 75)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12))


df2 <- as.data.frame(df1$A)
df2 <- rbind(df2, c(150))
df2 <- rbind(df2, c(153))

ggplot(df2, aes(x = df2$`df1$A`)) +
  #geom_histogram(aes(y = after_stat(count / sum(count))), fill = '#0BD0D9', col = 'black', size = 0.9) +
  geom_density(size = 1.2, col = '#D60093', fill = 'pink', alpha = 0.7) +
  ylab('Relative frequency') +
  xlab('Variable') +
  scale_x_continuous(breaks = seq(20, 160, 10), limits = c(20, 160)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1), limits = c(0, 0.4)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12))


df3 <- df1 %>% pivot_longer(values_to = 'Group', names_to = 'Sample', cols = c('A', 'B', 'C'))
ggplot(df3, aes(x = Group, group = Sample, fill = Sample)) +
  #geom_histogram(aes(y = after_stat(count / sum(count))), fill = '#0BD0D9', col = 'black', size = 0.9) +
  geom_density(size = 1.2, col = 'black', alpha = 0.7) +
  ylab('Relative frequency') +
  xlab('Gene expression') +
  scale_x_continuous(breaks = seq(20, 80, 10), limits = c(20, 80)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1), limits = c(0, 0.4)) +
  biostatsquid_theme +
  theme(axis.text = element_text(size = 12)) +
  facet_wrap(~ Sample)
