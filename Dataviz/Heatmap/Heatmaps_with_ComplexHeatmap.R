# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ComplexHeatmap Tutorial: Comprehensive Guide with Patient Data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author: Laura Twomey
# Version: 1.0
# Description: ComplexHeatmap is a powerful R package for creating highly customizable 
# heatmaps with annotations. This tutorial demonstrates various features using simulated patient data.
# https://jokergoo.github.io/ComplexHeatmap-reference/book/

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 1. Create Sample Patient Data =====================
# Generate sample data for 50 patients
n_patients <- 50
n_genes <- 100

# Create gene expression matrix
gene_names <- paste0("Gene_", 1:n_genes)
patient_ids <- paste0("Patient_", sprintf("%02d", 1:n_patients))

# Simulate gene expression data (log2 fold changes)
expression_data <- matrix(
  rnorm(n_patients * n_genes, mean = 0, sd = 1.5),
  nrow = n_genes,
  ncol = n_patients,
  dimnames = list(gene_names, patient_ids)
)

# Add some structure to make heatmap more interesting
expression_data[1:20, 1:15] <- expression_data[1:20, 1:15] + 2    # Upregulated cluster
expression_data[21:40, 16:30] <- expression_data[21:40, 16:30] - 2 # Downregulated cluster
expression_data[41:60, 31:45] <- expression_data[41:60, 31:45] + 1.5 # Moderate upregulation

head(expression_data[1:5, 1:5])
dim(expression_data)

# 2 Create Patient Annotations =============================================
# Create patient annotation data
patient_annotations <- data.frame(
  Patient_ID = patient_ids,
  
  # Continuous variables
  Age = round(rnorm(n_patients, mean = 65, sd = 12)),
  BMI = round(rnorm(n_patients, mean = 25, sd = 4), 1),
  Tumor_Size = round(rnorm(n_patients, mean = 3.5, sd = 1.2), 1),
  
  # Discrete variables
  Gender = sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.6, 0.4)),
  Stage = sample(c("I", "II", "III", "IV"), n_patients, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2)),
  Treatment = sample(c("Surgery", "Chemo", "Radiation", "Combined"), n_patients, replace = TRUE),
  Response = sample(c("Complete", "Partial", "Stable", "Progressive"), n_patients, replace = TRUE),
  
  stringsAsFactors = FALSE
)

# Ensure age is within reasonable bounds
patient_annotations$Age <- pmax(18, pmin(90, patient_annotations$Age))
patient_annotations$BMI <- pmax(15, pmin(45, patient_annotations$BMI))
patient_annotations$Tumor_Size <- pmax(0.5, pmin(10, patient_annotations$Tumor_Size))

head(patient_annotations)

# 3. Basic Heatmap ====================================================

# Create basic heatmap
basic_heatmap <- Heatmap(
  expression_data,
  name = "Expression",
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  column_title = "Patient Gene Expression Heatmap"
)

# Draw the heatmap
draw(basic_heatmap)

# 4: Custom Color Schemes ===============================================

# Define custom color functions
col_fun1 <- colorRamp2(c(-3, 0, 3), c("darkgreen", "white", "orange3"))
col_fun2 <- colorRamp2(c(-3, 0, 3), c("navy", "lightgray", "darkred"))
col_fun3 <- RColorBrewer::brewer.pal(name = 'Spectral', n = 11)

# Heatmap with custom colors
custom_color_heatmap <- Heatmap(
  expression_data,
  name = "Expression",
  col = col_fun1,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "Custom Color Heatmap",
  heatmap_legend_param = list(
    title = "Log2 FC",
    legend_direction = "horizontal",
    legend_width = unit(6, "cm")
  )
)

draw(custom_color_heatmap, heatmap_legend_side = "bottom")

# 5. Adding Column Annotations =====================================
colour_list <- list(
  Age = colorRamp2(c(18, 90), c("lightblue", "darkblue")),
  BMI = colorRamp2(c(15, 45), c("lightgreen", "darkgreen")),
  Tumor_Size = colorRamp2(c(0.5, 10), c("lightyellow", "darkorange")),
  Gender = c("Male" = "lightblue", "Female" = "pink"),
  Stage = c("I" = "green", "II" = "yellow", "III" = "orange", "IV" = "red"),
  Treatment = c("Surgery" = "purple", "Chemo" = "blue", "Radiation" = "orange", "Combined" = "magenta4"),
  Response = c("Complete" = "darkgreen", "Partial" = "lightgreen", "Stable" = "yellow", "Progressive" = "red")
)

# Create column annotation object
col_annotation <- HeatmapAnnotation(
  # Continuous annotations
  Age = patient_annotations$Age,
  BMI = patient_annotations$BMI,
  Tumor_Size = patient_annotations$Tumor_Size,
  
  # Discrete annotations
  Gender = patient_annotations$Gender,
  Stage = patient_annotations$Stage,
  Treatment = patient_annotations$Treatment,
  Response = patient_annotations$Response,
  
  # Custom colors for discrete variables
  col = colour_list,
  
  # Annotation parameters
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Age = list(title = "Age (years)"),
    BMI = list(title = "BMI"),
    Tumor_Size = list(title = "Tumor Size (cm)")
  )
)

# Heatmap with annotations
annotated_heatmap <- Heatmap(
  expression_data,
  name = "Expression",
  col = col_fun1,
  top_annotation = col_annotation,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "Gene Expression with Patient Annotations"
)

draw(annotated_heatmap)

# Let's fix that a bit...
colour_list <- list(
  Age = colorRamp2(c(18, 90), c("white", "magenta4")),
  BMI = colorRamp2(c(15, 45), c("white", "magenta4")),
  Tumor_Size = colorRamp2(c(0.5, 10), c("white", "magenta4")),
  Gender = c("Male" = "aquamarine3", "Female" = "darkgoldenrod2"),
  Stage = c("I" = "white", "II" = "aquamarine3", "III" = "deepskyblue2", "IV" = "deepskyblue4"),
  Treatment = c("Surgery" = "white", "Chemo" = "lightblue", "Radiation" = "lightpink", "Combined" = "darkred"),
  Response = c("Complete" = "white", "Partial" = "lightblue", "Stable" = "lightpink", "Progressive" = "darkred")
)

# Create column annotation object for the top and one for the bottom
col_ann <- HeatmapAnnotation(
  Age = patient_annotations$Age,
  BMI = patient_annotations$BMI,
  Tumor_Size = patient_annotations$Tumor_Size,
  Gender = patient_annotations$Gender,

  # Custom colors for discrete variables
  col = colour_list,
  na_col = "gray",
  gp = gpar(col = "white"),
  
  # Annotation parameters
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Age = list(title = "Age (years)"),
    BMI = list(title = "BMI"),
    Tumor_Size = list(title = "Tumor Size (cm)")
  )
)

bottom_ann <- HeatmapAnnotation(
  
  Stage = patient_annotations$Stage,
  Treatment = patient_annotations$Treatment,
  Response = patient_annotations$Response,
  
  # Custom colors for discrete variables
  col = colour_list,
  gp = gpar(col = "white"),
  
  # Annotation parameters
  annotation_name_gp = gpar(fontsize = 10)
)

# Heatmap with annotations
annotated_heatmap <- Heatmap(
  expression_data,
  name = "Expression",
  
  # Clustering parameters
  clustering_distance_rows = "ward.D2", # check other clustering methods
  clustering_method_rows = "complete",
  clustering_distance_columns = "ward.D2",
  clustering_method_columns = "complete",
  
  col = col_fun2,
  
  top_annotation = col_ann,
  bottom_annotation = bottom_ann,
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  column_title = "Gene Expression with Patient Annotations"
)

draw(annotated_heatmap)

# What about row annotations?
# Create gene categories for row annotations
gene_categories <- data.frame(
  Gene = gene_names,
  Pathway = sample(c("Metabolism", "Cell Cycle", "Apoptosis", "DNA Repair"), 
                   n_genes, replace = TRUE),
  Function = sample(c("Oncogene", "Tumor Suppressor", "Metabolic", "Structural"), 
                    n_genes, replace = TRUE),
  stringsAsFactors = FALSE
)

# Row annotation
row_annotation <- rowAnnotation(
  Pathway = gene_categories$Pathway,
  Function = gene_categories$Function,
  
  col = list(
    Pathway = c("Metabolism" = "cyan2", "Cell Cycle" = "gold3", 
                "Apoptosis" = "magenta4", "DNA Repair" = "green4"),
    Function = c("Oncogene" = "gold3", "Tumor Suppressor" = "magenta4",
                 "Metabolic" = "cyan2", "Structural" = "magenta4")
  ),
  
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Pathway = list(title = "Biological Pathway"),
    Function = list(title = "Gene Function")
  )
)

# Heatmap with row and column annotations including regulation triangles
full_heatmap <- Heatmap(
  
  expression_data,
  name = "Expression",
  
  # Clustering parameters
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  
  col = col_fun2,
  
  top_annotation = col_ann,
  bottom_annotation = bottom_ann,
  left_annotation = row_annotation,
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45, column_names_side = "top",
  column_names_gp = gpar(fontsize = 15, face = 'bold', angle = 45),
  
  column_title = "Gene Expression with Patient Annotations",
  
  # Splitting
  column_split = patient_annotations$Stage,
  row_split = gene_categories$Pathway,
  
  row_title_rot = 0,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Gaps between splits
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm"),
  
  # Borders and spacing
  rect_gp = gpar(col = "white", lwd = 0.5),
  border = TRUE
)

draw(full_heatmap)

# 6. Bar Chart Annotations ============================================

# Generate cell composition data (percentages that sum to 100%)
cell_composition <- data.frame(
  Patient_ID = patient_ids,
  Tumor_cells = runif(n_patients, 20, 70),
  T_cells = runif(n_patients, 10, 40),
  Macrophages = runif(n_patients, 5, 25),
  B_cells = runif(n_patients, 2, 15),
  stringsAsFactors = FALSE
)
# Normalize to sum to 100%
row_sums <- rowSums(cell_composition[, 2:5])
cell_composition[, 2:5] <- cell_composition[, 2:5] / row_sums * 100

# Cell type colors
cell_colors <- c("Tumor_cells" = "#E74C3C", "T_cells" = "#3498DB", 
                 "Macrophages" = "#F39C12", "B_cells" = "#27AE60")

# Create summary statistics for bar chart annotations
treatment_counts <- table(patient_annotations$Treatment)
response_counts <- table(patient_annotations$Response)

# Cell composition barplot annotation
cell_barplot_ann <- HeatmapAnnotation(
  Cell_Composition = anno_barplot(
    cell_composition[, c("Tumor_cells", "T_cells", "Macrophages", "B_cells")],
    gp = gpar(fill = cell_colors, col = "white"),
    bar_width = 0.8,
    height = unit(3, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# Cell composition barplot annotation
cell_barplot_ann <- HeatmapAnnotation(
  Cell_Composition = anno_barplot(
    cell_composition[, c("Tumor_cells", "T_cells", "Macrophages", "B_cells")],
    gp = gpar(fill = cell_colors, col = "white"),
    bar_width = 0.8,
    height = unit(2, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

draw(cell_barplot_ann)
dev.off()

heatmap_with_cells <- Heatmap(
  expression_data,
  name = "Expression",
  
  # Clustering parameters
  clustering_distance_rows = "correlation",
  clustering_method_rows = "ward.D2",
  clustering_distance_columns = "correlation",
  clustering_method_columns = "ward.D2",
  
  col = col_fun2,
  
  top_annotation = c(col_ann, cell_barplot_ann),
  left_annotation = row_annotation,
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 10),
  
  column_title = "Gene Expression with Cell Composition Data",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Splitting
  column_split = patient_annotations$Stage,
  row_split = gene_categories$Pathway,
  
  row_title_rot = 0,
  
  # Gaps between splits
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm")
)

# Draw the heatmap with cell composition
draw(heatmap_with_cells)

# 7. Editing the legend and playing around with size ==========================
# Set width and height to make cells square
n_rows <- nrow(expression_data)
n_cols <- ncol(expression_data)

heatmap_with_cells <- Heatmap(
  expression_data,
  name = "Expression",
  
  # Clustering parameters
  clustering_distance_rows = "correlation",
  clustering_method_rows = "ward.D2",
  clustering_distance_columns = "correlation",
  clustering_method_columns = "ward.D2",
  
  col = col_fun2,
  
  top_annotation = c(col_ann, cell_barplot_ann),
  left_annotation = row_annotation,
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 10),
  
  column_title = "Gene Expression with Cell Composition Data",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Splitting
  column_split = patient_annotations$Stage,
  row_split = gene_categories$Pathway,
  
  row_title_rot = 0,
  
  # Gaps between splits
  column_gap = unit(0.7, "mm"),
  row_gap = unit(0.7, "mm"),
  
  # heatmap_width and heatmap_height control the width/height of the complete heatmap 
  # including all heatmap components (excluding the legends) while width and height 
  # only control the width/height of the heamtap body. All these four arguments can be set as absolute units.
  # width = unit(15, "cm"), height = unit(11, "cm")
  # heatmap_width = unit(12, "cm"), heatmap_height = unit(12, "cm")
  width = unit(n_cols/5, "cm"), height = unit(n_rows/5, "cm"),
  
  # heatmap_legend_param = list(
  #   title_gp = gpar(fontsize = 12, fontface = "bold"),
  #   labels_gp = gpar(fontsize = 10)
  # )
  
  # heatmap_legend_param = list(
  #   title = "Expression Level",
  #   title_gp = gpar(fontsize = 14, fontface = "bold", col = "navy"),
  #   labels_gp = gpar(fontsize = 10),
  #   legend_height = unit(6, "cm"),
  #   legend_width = unit(1.2, "cm"),
  #   color_bar = "continuous",
  #   border = "darkblue"
  # )
)

draw(heatmap_with_cells, heatmap_legend_side = "left") 
draw(heatmap_with_cells, heatmap_legend_side = "bottom")
# No heatmap legend, only annotation legends
draw(heatmap_with_cells, show_heatmap_legend = FALSE)


# 8. Advanced Annotation Types =============================================

# Text annotations with custom formatting
text_ann <- HeatmapAnnotation(
  Patient_Info = anno_text(
    paste0("Sample", gsub('Patient_', '', colnames(expression_data))),
    gp = gpar(fontsize = 8, col = "darkblue"),
    rot = 90,
    height = unit(1.5, "cm")
  ),
  
  # Points annotation for continuous data
  Age_Points = anno_points(
    patient_annotations$Age,
    pch = 16,
    size = unit(2, "mm"),
    gp = gpar(col = "darkred"),
    height = unit(2, "cm"),
    axis_param = list(at = c(20, 40, 60, 80), labels = c("20", "40", "60", "80"))
  ),
  
  # Line plot annotation
  BMI_Line = anno_lines(
    patient_annotations$BMI,
    gp = gpar(col = "blue", lwd = 2),
    height = unit(2, "cm")
  ),
  
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# Draw text annotation example
draw(text_ann)

# Let's add it to the heatmap, and order it by BMI... we need to turn off clustering by columns!

# First define the column order:
patient_annotations <- patient_annotations %>% dplyr::arrange(desc(BMI))
patient_order <- patient_annotations %>% pull(Patient_ID)
expression_data <- expression_data[,patient_order]

text_ann <- HeatmapAnnotation(
  Patient_Info = anno_text(
    paste0("Sample", gsub('Patient_', '', colnames(expression_data))),
    gp = gpar(fontsize = 8, col = "darkblue"),
    rot = 90,
    height = unit(1.5, "cm")
  ),
  
  # Points annotation for continuous data
  Age_Points = anno_points(
    patient_annotations$Age,
    pch = 16,
    size = unit(2, "mm"),
    gp = gpar(col = "darkred"),
    height = unit(2, "cm"),
    axis_param = list(at = c(20, 40, 60, 80), labels = c("20", "40", "60", "80"))
  ),
  
  # Line plot annotation
  BMI_Line = anno_lines(
    patient_annotations$BMI,
    gp = gpar(col = "blue", lwd = 2),
    height = unit(2, "cm")
  ),
  
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

heatmap_with_text_ann <- Heatmap(
  expression_data,
  name = "Expression",
  
  # Clustering parameters
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  cluster_columns = FALSE,
  
  col = col_fun2,
  
  bottom_annotation = text_ann,
  left_annotation = row_annotation,
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_rot = 90,
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 10),
  
  column_title = "Gene Expression with Cell Composition Data",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Splitting
  column_split = patient_annotations$Stage,
  row_split = gene_categories$Pathway,
  
  row_title_rot = 0,
  
  # Gaps between splits
  column_gap = unit(0.7, "mm"),
  row_gap = unit(0.7, "mm"),
  
  
  width = unit(n_cols/4, "cm"), height = unit(n_rows/5, "cm"),
  
)

# Try turninng on and off the column_split parameter
draw(heatmap_with_text_ann)


# 9. Combine heatmaps ================================================
combined_heatmap <- heatmap_with_cells + heatmap_with_text_ann
draw(combined_heatmap, column_title = "Multi-omics Patient Data")

# 10. Dendogram Customization ==========================================

# Custom dendrogram colors based on clusters
library(dendextend)

# Perform clustering
row_dend <- as.dendrogram(hclust(dist(expression_data)))
col_dend <- as.dendrogram(hclust(dist(t(expression_data))))

# Color dendrograms
row_dend <- color_branches(row_dend, k = 4)
col_dend <- color_branches(col_dend, k = 3)

dendrogram_heatmap <- Heatmap(
  expression_data,
  name = "Expression",
  col = col_fun2,
  cluster_rows = row_dend,
  cluster_columns = col_dend,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(3, "cm")
)

draw(dendrogram_heatmap)

# Saving it! ===================================================
# Save it
path_to_folder <- '/my/path/to/output/folder'
pdf(file.path(path_to_folder, "patient_heatmap.pdf"), width = 12, height = 15)
draw(bar_chart_heatmap)
dev.off()
