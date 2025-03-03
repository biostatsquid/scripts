# ======================== # 
# Violin_plots.R
# ======================== # 

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(reshape2) # To reshape data for point 2

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

# Importing dataset ===================================================
data(iris)
head(iris)
table(iris$Species)

# 0. Convert categorical variable to factor =================================
# This is especially important if it is numerical (e.g., doses of a drug, age... etc)
iris$Species <- as.factor(iris$Species)

# 1. Create violin plots ===================================================

# Basic plot
ggplot(data = iris, aes(x = Species, y = Sepal.Width)) +
  geom_violin(trim = FALSE) +
  coord_flip()

## Adding stats -----------------------------------------------------------

# Adding the mean or median as a point
# http://www.sthda.com/english/wiki/ggplot2-point-shapes
ggplot(data = iris, aes(x = Species, y = Sepal.Width)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, size = 2, colour = 'darkblue')
# Adding a boxplot
ggplot(data = iris, aes(x = Species, y = Sepal.Width)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, size = 2, colour = 'darkblue') +
  geom_boxplot(width = 0.1)
# Adding mean and standard deviation
# You can add it as a crossbar
ggplot(data = iris, aes(x = Species, y = Sepal.Width)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, size = 2, colour = 'darkblue') +
  stat_summary(fun.data = "mean_sdl", mult = 1.5, 
               geom = "crossbar", width = 0.2 )
# Or a poinrange
ggplot(data = iris, aes(x = Species, y = Sepal.Width)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, size = 2, colour = 'darkblue') +
  stat_summary(fun.data = "mean_sdl", mult = 1.5, 
               geom = "pointrange", width = 0.2 )

# Adding jitter
ggplot(data = iris, aes(x = Species, y = Sepal.Width)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, size = 3, colour = 'darkblue', fill = 'darkred') +
  geom_jitter(shape = 16, position = position_jitter(0.2))

## Adding colour -----------------------------------------------------------
# Adding colour
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8)

# Changing colours
# colour changes the contour of the violin plot
ggplot(data = iris, aes(x = Species, y = Sepal.Width, col = Species)) +
  geom_violin(trim = FALSE, linewidth = 0.8, fill = 'beige') +
  scale_color_brewer(palette = "Dark2")
# fill changes the inside of the violin plot
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, linewidth = 0.8) +
  scale_fill_brewer(palette = "Dark2")

# Custom colours
colour_dict <- list('setosa' = 'pink2', 'versicolor' = 'turquoise3', 'virginica' = 'orange3')
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict)

## Plot options -------------------------------------
# Axis titles
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "My amazing violin plot", x = "Species", y = "Sepal Width (cm)")

## Theme options -------------------------------------
# Everything is customisable!
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "My amazing violin plot", x = "Species", y = "Sepal Width (cm)") +
  theme_bw()

# Legend position
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "My amazing violin plot", x = "Species", y = "Sepal Width (cm)") +
  theme_bw() +
  theme(legend.position = 'bottom')

# Aspect ratio
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "My amazing violin plot", x = "Species", y = "Sepal Width (cm)") +
  theme_bw() +
  theme(legend.position = 'bottom',
        aspect.ratio = 0.5)

# Axis text
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "My amazing violin plot", x = "Species", y = "Sepal Width (cm)") +
  theme_bw() +
  theme(legend.position = 'bottom',
        aspect.ratio = 0.5,
        axis.text = element_text(size = 12, face = 'bold', colour = 'purple4'),
        axis.title = element_text(size = 16, face = 'bold', colour = 'darkblue'),
        plot.title = element_blank())

# Usings customisable themes
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(x = "Species", y = "Sepal Width (cm)") +
  biostatsquid_theme

## Multiple groups -----------------------------------------------------------
iris$Treatment <- rep(c('treated', 'untreated'), nrow(iris)/2)
head(iris)

colour_dict <- list('treated' = 'pink2', 'untreated' = 'turquoise3')
ggplot(data = iris, aes(x = Species, y = Sepal.Width, fill = Treatment)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  biostatsquid_theme 


# 2. Creating many violin plots at once: melt ===================================================
df_long <- melt(iris, id.vars = c('Species', 'Treatment'))
head(df_long)

# Use facet_wrap
colour_dict <- list('pink2', 'turquoise3', 'orange3')
ggplot(df_long, aes(x = Species, y = value, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "Violin Plots for All Variables", x = "", y = "") +
  facet_wrap(~variable) +
  biostatsquid_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold", colour = "green4"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1)) 

# Adding group by treatment
ggplot(df_long, aes(x = Species, y = value, fill = Treatment)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "Violin Plots for All Variables", x = "", y = "") +
  facet_wrap(~variable) +
  biostatsquid_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold", colour = "green4"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1)) 

# Use facet grid
ggplot(df_long, aes(x = Species, y = value, fill = Species)) +
  geom_violin(trim = FALSE, col = 'black', linewidth = 0.8) +
  scale_fill_manual(values = colour_dict) +
  labs(title = "Violin Plots for All Variables", x = "", y = "") +
  facet_grid(vars(Treatment), vars(variable)) +
  biostatsquid_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold", colour = "green4"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1)) 


# 3. Creating many violin plots at once: function ===================================================

# Basic plot
ggplot(data = iris, aes(x = Species, y = Sepal.Width, col = Species)) +
  geom_violin(trim = FALSE, linewidth = 0.8, fill = 'white') +
  scale_color_brewer(palette = "Dark2") +
  biostatsquid_theme

get_violin_plots <- function(data, x_variable, y_variable, col_variable){
  
  p1 <- ggplot(data = data, aes(x = !!sym(x_variable), y = !!sym(y_variable), col = !!sym(col_variable))) +
    geom_violin(trim = FALSE, linewidth = 0.8, fill = 'white') +
    scale_color_brewer(palette = "Dark2") +
    biostatsquid_theme
  
  return(p1)
}

get_violin_plots(data = iris, x_variable = 'Species', y_variable = 'Sepal.Width', col_variable = 'Species')

# Adding arguments to your function
get_violin_plots_v2 <- function(data, x_variable, y_variable, col_variable, title_x = NULL, title_y = NULL, custom_cols = brewer.pal(10, 'Set3')){
  
  p1 <- ggplot(data = data, aes(x = !!sym(x_variable), y = !!sym(y_variable), fill = !!sym(col_variable))) +
    geom_violin(trim = FALSE, linewidth = 0.8, colour = 'darkgrey') +
    scale_fill_manual(values = custom_cols) +
    labs(x = title_x, y = title_y) +
    biostatsquid_theme
  
  return(p1)
}

get_violin_plots_v2(data = iris, x_variable = 'Species', y_variable = 'Sepal.Width', col_variable = 'Species',
                    title_x = 'Species', title_y = 'Sepal.Width', 
                    custom_cols = list('pink', 'gold', 'brown'))

colour_dict <- list('pink2', 'turquoise3', 'orange3')
get_violin_plots_v2(data = iris, x_variable = 'Species', y_variable = 'Sepal.Width', col_variable = 'Species',
                    title_x = 'Species', title_y = 'Sepal.Width', 
                    custom_cols = colour_dict)                

get_violin_plots_v2(data = iris, x_variable = 'Species', y_variable = 'Petal.Length', col_variable = 'Species',
                    title_x = 'Species', title_y = 'Petal Length (cm)', 
                    custom_cols = colour_dict)   

# It works with any dataset!
data("diamonds")
head(diamonds)
# This won't work because you did not provide enough colours
get_violin_plots_v2(data = diamonds, x_variable = 'cut', y_variable = 'price', col_variable = 'color',
                    title_x = 'Diamond cut', title_y = 'Price (eur)', custom_cols = colour_dict)       

# This won't work because you did not provide enough colours
get_violin_plots_v2(data = diamonds, x_variable = 'cut', y_variable = 'price', col_variable = 'color',
                    title_x = 'Diamond cut', title_y = 'Price (eur)')  


# You can always customise it more!
get_violin_plots_v2(data = diamonds, x_variable = 'cut', y_variable = 'price', col_variable = 'color',
                    title_x = 'Diamond cut', title_y = 'Price (eur)')  +
  coord_flip() +
  scale_fill_manual('Colour of diamonds', values = brewer.pal(7, name = 'Set1')) +
  theme(axis.ticks = element_line(colour = 'orange'))
  




