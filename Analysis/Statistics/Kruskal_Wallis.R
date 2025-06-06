# ----------------------
#  Kruskal-Wallis & Dunn’s Test on Penguin Data
# ----------------------
# We’ll analyze penguin body mass across different species.

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(123456)

# Loading relevant libraries 
library(palmerpenguins)
library(FSA)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

# Load and clean data
data("penguins")
penguins_clean <- penguins %>%
  filter(!is.na(body_mass_g), !is.na(species))

# View summary
summary(penguins_clean$species)

# Kruskal-Wallis test
kruskal_test <- penguins_clean %>%
  kruskal_test(body_mass_g ~ species)
print(kruskal_test)
# If p < 0.05, proceed to Dunn’s test.

# Dunn’s Test (Post-hoc)
dunn_results <- penguins_clean %>%
  dunn_test(body_mass_g ~ species, p.adjust.method = "bonferroni")
# This shows pairwise comparisons between species (e.g., Adelie vs Gentoo).

# Get only significant comparisons
sig_comparisons <- dunn_results %>%
  filter(p.adj < 0.05) %>%
  select(group1, group2, p.adj) %>%
  mutate(p.adj = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*"
  ))

# Convert to list of pairs:
comparisons_list <- sig_comparisons %>%
  select(group1, group2) %>%
  as.list() %>%
  t()

# Check p-value
kw_p <- kruskal_test$p
# Create label based on significance
kw_label <- if (kw_p < 0.01) {
  "Kruskal-Wallis p < 0.01"
} else {
  paste0("Kruskal-Wallis p = ", signif(kw_p, 3))
}


ggplot(penguins_clean, aes(x = species, y = body_mass_g)) +
  geom_boxplot(aes(fill = species)) +
  scale_fill_manual(values = c('lightblue', 'lightpink', 'beige')) +
  theme_bw() +
  labs(title = "Penguin Body Mass by Species",
       y = "Body Mass (g)",
       x = "") +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold')
        ) +
  annotate("text", x = 1, y = 6700, label = kw_label, size = 5) +
  stat_pvalue_manual(sig_comparisons,
                     label = "p.adj",
                     y.position = seq(6200, 6600, length.out = nrow(sig_comparisons)),
                     tip.length = 0.01)

