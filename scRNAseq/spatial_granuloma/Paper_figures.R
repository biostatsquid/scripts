# Name: Paper_figures.R
# Author: Laura Twomey
# Date of creation: 06 July 2022
# Figures for paper on spatial transcriptomics of skin granulomas

# Clean Environment ----------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects including hidden objects
gc() #free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F) #avoid truncated output in R console and scientific notation

# Set working directory ----------------------------------------------------
path <- "/Volumes/Drive/spatial_granuloma/output/Paper_figures"
setwd(path)

# Libraries  ------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(pheatmap)
library(anndata)
library(ggpubr)
library(png)

# Figure 3  ------------------------------------------------

## Data import ####
# The frequency table has genes x sample_SPECIMENs and 1 or 0 depending if the gene is present or not
freqtable_all <- readRDS(file = paste(path, 'Figure3_SpatialDE', 'freqtable_all.csv', sep = '/'))
# The annotations for the heatmap
col_annotations_df <- read.csv(paste(path, 'Figure3_SpatialDE', 'col_annotations_df.csv', sep = '/'), row.names = 1)
rownames(col_annotations_df) <- col_annotations_df$sample_SPECIMEN_pattern
col_annotations_df$KMeans5 <- factor(col_annotations_df$KMeans5)
# Select here the columns you want to plot as annotations in the heatmap
heatmap_ann <- col_annotations_df %>% dplyr::select(Disease, matches('signedFDR|KMeans5'), -matches('leiden_epidermis|leiden_granuloma'))
# Rename the columns so they appear nicely in the heatmap
heatmap_ann <- heatmap_ann %>% rename(`Core granuloma (Leiden)` = signedFDR.leiden_core_granuloma,
                                      `Border granuloma (Leiden)` = signedFDR.leiden_border_granuloma,
                                      #`Epidermis (Leiden)` = signedFDR.leiden_epidermis,
                                      #`Granuloma (Leiden)` = signedFDR.leiden_granuloma,
                                      `Granuloma (annotation)` = signedFDR.manual_granuloma, 
                                      `Epidermis (annotation)` = signedFDR.epidermis_interface, 
                                      `K-means clustering` = KMeans5)
# Set the order of the annotations in the heatmap
col_order <- c("Core granuloma (Leiden)", "Border granuloma (Leiden)", "Granuloma (annotation)", 
               "Epidermis (annotation)", "Disease", "K-means clustering")
heatmap_ann <- heatmap_ann[, col_order]

## Colours  ####

# Set patients for subsetting freqtable
ann_colors = list(
  patient = c("91253" = "#F46D43",#'red',#granuloma annulare
              "45703" = "#D53E4F",#'indianred',#granuloma annulare
              "50107" = "#9E0142",#'firebrick',#granuloma annulare
              "95096" = "orchid", # necrobiosis lipoidica
              "82301" = "mediumvioletred", #sarcoidosis-cutaneous
              "72859" = "blueviolet"), #sarcoidosis-suspected
  DISEASE = c("granuloma annulare" = "#9E0142",#'firebrick'
              "necrobiosis lipoidica" = "orchid", # necrobiosis lipoidica
              "sarcoidosis" = "blueviolet"), #sarcoidosis
  Disease = c("Granuloma annulare" = "#9E0142",#'firebrick'
              "Necrobiosis lipoidica" = "orchid", # necrobiosis lipoidica
              "Sarcoidosis" = "blueviolet"), #sarcoidosis
  leiden_core_granuloma = rev(brewer.pal(11, "RdBu")),
  leiden_border_granuloma = rev(brewer.pal(11, "RdBu")),
  leiden_epidermis = rev(brewer.pal(11, "RdBu")),
  leiden_granuloma = rev(brewer.pal(11, "RdBu")),
  manual_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_core_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_border_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_epidermis = rev(brewer.pal(11, "RdBu")),
  SignedFDR.leiden_granuloma = rev(brewer.pal(11, "RdBu")),
  SignedFDR.manual_granuloma = rev(brewer.pal(11, "RdBu")),
  `Core granuloma (Leiden)` = rev(brewer.pal(11, "RdBu")),
  `Border granuloma (Leiden)` = rev(brewer.pal(11, "RdBu")),
  `Epidermis (Leiden)` = rev(brewer.pal(11, "RdBu")),
  `Epidermis (annotation)` = rev(brewer.pal(11, "RdBu")),
  `Granuloma (annotation)` = rev(brewer.pal(11, "RdBu")),
  KMeans5 = c("1"="#88CFA4", "2"="#F88D52", "3"="skyblue",
              "4"="#9E0142", "5"="#5E4FA2"),
  `K-means clustering` = c("1"="#88CFA4", "2"="#F88D52", "3"="skyblue",
                           "4"="#9E0142", "5"="#5E4FA2")
  #sample_SPECIMEN_SAMPLE = sample(color, 20)
)

BuRd <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
GrBu <- colorRampPalette(colors = c(rev(brewer.pal(11, "PiYG"))[1:4],
                                    brewer.pal(11, "RdYlBu")[c(6, 8:11)]))

## Heatmap  ####
heatmap <- pheatmap(freqtable_all,
                    show_colnames = F, show_rownames = F,
                    clustering_method = 'ward.D',
                    annotation_col = heatmap_ann,
                    treeheight_row = 0,
                    treeheight_col = 10,
                    annotation_colors = ann_colors)

## UMAPS: K-Means clustering ####

# Getting patterns from a particular cluster:
to_python_format2(col_annotations_df$sample_SPECIMEN_pattern[col_annotations_df$KMeans5 == '3'])

kmeans <-  ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = KMeans5)) +
  geom_point() +
  scale_color_manual("K-means (5)",
                     values = c("1"="#88CFA4", "2"="#F88D52", "3"="dodgerblue",
                                "4"="#9E0142", "5"="#5E4FA2"),
                     labels = c("1", "Granuloma-2", "Interface", "Granuloma-1", "Epidermis")) +
  theme_bw() +
  ggtitle("K-Means clustering (K=5)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(), 
    legend.text = element_text(size=rel(1.06)),
    legend.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) 

disease <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = Disease)) +
  geom_point() +
  scale_color_manual("Disease",
                     values = ann_colors$Disease,
                     labels = c("Granuloma\nannulare", "Necrobiosis\nlipoidica", "Sarcoidosis")) +
  theme_bw() +
  ggtitle("Disease") +
  theme(#legend.position = "bottom",
    legend.key.height=unit(25, "pt"),
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.06)),
    legend.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)))

manualg <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.manual_granuloma)) +
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = BuRd(100)) +
  #scale_colour_gradientn(colours = BuRd(100), limits=c(-100, 100)) +
  ggtitle("Granuloma (annotation)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

epidermis <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.epidermis_interface)) +
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = GrBu(100)) +
  #scale_colour_gradientn(colours = GrBu(100), limits=c(-100, 100)) +
  ggtitle("Epidermis (annotation)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

core <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.leiden_core_granuloma)) + 
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = BuRd(100)) +
  #scale_colour_gradientn(colours = BuRd(100), limits=c(-100, 100)) +
  ggtitle("Granuloma core (Leiden cluster)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

border <- ggplot(col_annotations_df, aes(x = UMAP1, y = UMAP2, col = signedFDR.leiden_border_granuloma)) + 
  geom_point() +
  theme_bw() +
  scale_colour_gradientn(colours = BuRd(100)) +
  #scale_colour_gradientn(colours = BuRd(100), limits=c(-100, 100)) +
  ggtitle("Granuloma border (Leiden cluster)") +
  theme(#legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank())

### All UMAPs #####
all_UMAPs <- ggarrange(
  ggarrange(kmeans + theme(legend.position = "none"),
            as_ggplot(get_legend(kmeans)),
            nrow = 1, ncol = 2,
            widths = c(1, 0.5)),
  ggarrange(disease + theme(legend.position = "none"),
            as_ggplot(get_legend(disease)),
            nrow = 1, ncol = 2,
            widths = c(1, 0.5)),
  ggarrange(manualg + theme(legend.position = "none"),
            as_ggplot(get_legend(manualg)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ggarrange(epidermis + theme(legend.position = "none"),
            as_ggplot(get_legend(epidermis)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ggarrange(core + theme(legend.position = "none"),
            as_ggplot(get_legend(core)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ggarrange(border + theme(legend.position = "none"),
            as_ggplot(get_legend(border)), NULL,
            nrow = 1, ncol = 3,
            widths = c(1, 0.4,0.1)),
  ncol = 2, nrow = 3,
  labels = c("", "", "",
             "", "", "")) 

## Figures E and D merged ####
pdf('Figure3_SpatialDE/Heatmap_UMAPs_merged.pdf', width = 16, height = 12)
ggarrange(heatmap$gtable, all_UMAPs, 
          labels = c("D", "E"), font.label = list(size = 30, color = "black"),
          ncol = 2, nrow = 1, widths = c(1, 1), vjust = 0.05, hjust = 1.6) +
  theme(plot.margin = margin(0.9,0.1,0.1,1.5, "cm")) 

dev.off()


## Final figure (does not work yet, low quality) ####
figA <- readPNG(paste0('Figure3_SpatialDE/', 'lesional_zoom_figure.png'))
figB <- readPNG(paste0('Figure3_SpatialDE/', 'Top12_tight.png'))
figC <- readPNG(paste0('Figure3_SpatialDE/', 'patterns_sharedlegend_tight.png'))

im_A <- ggplot() + background_image(figA) # tried to insert the image as background.. there must be a better way
im_B <- ggplot() + background_image(figB) 
im_C <- ggplot() + background_image(figC) 

png('test.png')
ggarrange(im_B, im_C, 
          labels = c("B", "C"), vjust = 1, hjust = 1,
          ncol = 2, nrow = 1, widths = c(12, 10)) +
  theme(plot.margin = margin(0.1,0.1,2,0.1, "cm")) 

dev.off()


# Figure 2 paper --------------------------------------------------------------
## Data import #####
setwd("/Users/mendenlab/work/spatial_granuloma/scripts/summary")
df <- readRDS("../../results/current/final/Granuloma_QC_clustering_annotations_df.RDS")

## Colours #####
col.set <- list(spots = c("DERMIS" = '#E0EEE0', 
                          "INTERFACE" = 'deepskyblue',
                          "VESSEL" = 'darkgreen',
                          "EPIDERMIS" = 'blue',
                          "HAIR FOLLICLE" = "#543005",
                          "MUSCLE" = 'darkcyan',
                          "SEBACEOUS GLAND" = 'mistyrose',
                          "SWEAT GLAND" = 'yellow',
                          "GNL" = 'orchid',
                          "GSS" = 'blueviolet',
                          "GSC" = 'mediumvioletred',
                          "GA (1)" = '#F46D43', 
                          "GA (2)" = '#D53E4F', 
                          "GA (3)" = '#9E0142',
                          "UNDETERMINED" = 'black'),
                spots_general = c("DERMIS" = '#E0EEE0', 
                                  "INTERFACE" = 'deepskyblue',
                                  "VESSEL" = 'darkgreen',
                                  "EPIDERMIS" = 'blue',
                                  "HAIR FOLLICLE" = "#543005",
                                  "MUSCLE" = 'darkcyan',
                                  "SEBACEOUS GLAND" = 'mistyrose',
                                  "SWEAT GLAND" = 'yellow',
                                  "GA" = 'firebrick',  
                                  "GNL" = 'orchid',
                                  "GSS" = 'blueviolet',
                                  "GSC" = 'mediumvioletred',
                                  "UNDETERMINED" = 'black'),
                dermis = c("UNDETERMINED" = 'black',
                           "upper EPIDERMIS" = 'blue',
                           "middle EPIDERMIS" = 'dodgerblue',
                           "basal EPIDERMIS" = 'skyblue',
                           "DERdepth1" = '#006837',
                           "DERdepth2" = '#238443',
                           "DERdepth3" = '#41AB5D',
                           "DERdepth4" = '#78C679',
                           "DERdepth5" = '#ADDD8E',
                           "DERdepth6" = '#D9F0A3',
                           "DERdepth7" = '#F7FCB9'),
                lesions = c("LESIONAL" = 'tomato',
                            "NON LESIONAL" = 'darkolivegreen',
                            "UNDETERMINED" = 'lightgrey'),
                patients = c("91253" = "#F46D43",#'red',#granuloma annulare
                             "45703" = "#D53E4F",#'indianred',#granuloma annulare
                             "50107" = "#9E0142",#'firebrick',#granuloma annulare
                             "95096" = "orchid", # necrobiosis lipoidica
                             "82301" = "mediumvioletred", #sarcoidosis-cutaneous
                             "72859" = "blueviolet"), #sarcoidosis-suspected
                samples = c("P17851_1001" = '#F46D43',
                            "P17851_1002" = '#FDAE61',
                            "P17851_1003" = '#D53E4F',
                            "P17851_1004" = 'coral3',
                            "P18554_1001" = '#9E0142',
                            "P18554_1002" = 'firebrick3',
                            "P18554_1003" = 'orchid',
                            "P18554_1004" = 'magenta',
                            "P18554_1005" = 'blueviolet',
                            "P18554_1006" = 'mediumpurple',
                            "P18554_1007" = 'mediumvioletred',
                            "P18554_1008" = 'violet'),
                leiden_r0.5 = c("0" = "#E0EEE0",
                                "1" = "#41AB5D",
                                "2" = "#ADDD8E",
                                "3" = "yellow",
                                "4" = "orchid",
                                "5" = "dodgerblue",
                                "6" = "blueviolet",
                                "7" = "#D53E4F",
                                "8" = "#F46D43",
                                "9" = "mistyrose",
                                "10" = "blue",
                                "11" = "#9E0142",
                                "12" = "darkcyan"),
                leiden_r0.8 = c("0"  = 'yellow',
                                "1"  = '#ADDD8E',
                                "2"  = '#E0EEE0',
                                "3"  = '#D9F0A3',
                                "4"  = "#EAADE7",
                                "5"  = 'deepskyblue',
                                "6"  = "orchid",
                                "7"  = "#D53E4F",
                                "8"  = "darkgreen",
                                "9"  = "#F46D43",
                                "10" = "blueviolet",
                                "11" = 'skyblue',
                                "12" = 'blue',
                                "13" = "mediumvioletred",
                                "14" = 'mistyrose',
                                "15" = "#9E0142",
                                "16" = 'darkcyan'),
                leiden_r1.3 = c("0"  = 'darkolivegreen', #mainly contains non-lesional spots
                                "1"  = '#D9F0A3', # dermdepth6
                                "2"  = '#238443', # derdepth2
                                "3"  = 'firebrick', # granuloma
                                "4"  = '#78C679', # derdepth4
                                "5"  = '#78C679', # derdepth4
                                "6"  = "#41AB5D", #dermis depth 3
                                "7"  = "#006837", #derm depth1
                                "8"  = '#ADDD8E', # dermdepth5
                                "9"  = "#238443", #derm depth2
                                "10" = '#78C679', # derdepth4
                                "11" = 'blue', #epidermis
                                "12" = 'orchid', #purple for very low % granuloma --> this is mostly ga
                                "13" = '#F46D43', #orange for very low % granuloma --> this is mostly gnl
                                "14" = 'dodgerblue', #mostly interface, has mid epidermis colour
                                "15" = "deepskyblue", #mostly interface, has interface colour
                                "16" = '#cfafaf', #sebaceous gland. 
                                "17" = 'yellow', # sweat gland
                                "18" = 'darkcyan', # muscle
                                "19" = '#006837')) #derdepth1

## UMAPs using patient batch-corrected data coloured by spot_type ####
UMAP_BC_spottype <- ggplot(df,
                           aes(x = X_BC_umap1, y = X_BC_umap2, col = spot_type)) +
  geom_point(size = 0.5) +
  scale_color_manual("", values = col.set$spots) +
  labs(title = "Patient batch correction",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4)))

## UMAP using patient batch-corrected data coloured by leiden clusters with r=1.3 ######
UMAP_BC_leiden1.3 <- ggplot(df,
                            aes(x = X_BC_umap1, y = X_BC_umap2, col = leiden_r1.3_patient)) +
  geom_point(size = 0.5) +
  scale_color_manual("", values = col.set$leiden_r1.3) +
  labs(title = "Leiden clustering r = 1.3 (patient batch correction)",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4)))

## Spatial representation of manual annotations ####
manual_spatial_plots <- ggplot(data=df,
                               aes(x=spatial1, y=spatial2, col=spot_type)) +
  geom_point(size=0.5) +
  scale_color_manual('', values = col.set$spots) +
  scale_y_reverse() +
  labs(title = "Spatial representation of manual annotations") +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
        #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  facet_wrap(.~sample,
             scales = "fixed",
             nrow = 2)

# Leiden spatial 
## Spatial representation of Leiden clustering ----
leiden_spatial_plots <- ggplot(data=df,
                               aes(x=spatial1, y=spatial2, col=leiden_r1.3_patient)) +
  geom_point(size=0.5) +
  scale_color_manual('', values = col.set$leiden_r1.3) +
  scale_y_reverse() +
  labs(title = "Spatial representation of Leiden clusters") +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = NA, colour=NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #strip.background=element_blank(),#element_rect(fill=NA, colour="black"),
        #strip.text = element_blank(),#element_text(colour = "black", size = 12), #face="bold"
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  facet_wrap(.~sample,
             scales = "fixed",
             nrow = 2)

## Fractions and cluster size r1.3 ####
### Leiden legend
leiden_legend <- ggplot() +
  geom_point(aes(x=as.factor(seq(0:19)), y = rep(0,20), 
                 color = factor(as.character(seq(0, 19)),
                                ordered = TRUE,
                                levels = as.character(seq(0, 19))),
                 size = 10)) +
  scale_color_manual(values = col.set$leiden_r1.3) +
  labs(title = "Test") +
  theme_bw() + 
  coord_flip() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "white")) +
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color="white"),
        axis.title = element_blank(),
        plot.title = element_text(colour = 'white')) +
  scale_x_discrete(limits=rev)
#leiden_legend

### Cluster size Leiden 1.3
clustersize_leiden13 <- ggplot(data=df,
                               aes(x=factor(paste0(leiden_r1.3_patient, "     "),
                                            ordered = TRUE,
                                            levels = as.character(paste0(0:19, "     "))),
                                   fill = leiden_r1.3_patient)) +
  geom_bar(stat="count", position = "stack") +
  scale_fill_manual("", values = col.set$leiden_r1.3) +
  labs(title = "Cluster size") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(limits=rev) +
  coord_flip() 
#clustersize_leiden13

### Cluster fractions Leiden 1.3
fractions_leiden13 <- ggplot(data=df,
                             aes(x=leiden_r1.3_patient,
                                 fill=spot_type)) +
  geom_bar(stat="count", position = "fill") +
  scale_fill_manual("", values = col.set$spots) +
  labs(title = "Spot type fraction") +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) +
  scale_x_discrete(limits=rev, breaks = c(0, 0.1, 1)) + #TODO
  coord_flip() 

### Patient fractions Leiden 1.3
disease_leiden13 <- ggplot(data=df,
                             aes(x=leiden_r1.3_patient,
                                 fill=patient)) +
  geom_bar(stat="count", position = "fill") +
  scale_fill_manual("", values = col.set$patients) +
  labs(title = "Patient") +
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) +
  scale_x_discrete(limits=rev, breaks = c(0, 0.1, 1)) + #TODO
  coord_flip()

### Skin depth fractions Leiden 1.3
skin_leiden13 <- ggplot(data=df,
                           aes(x=leiden_r1.3_patient,
                               fill=skin_layer)) +
  geom_bar(stat="count", position = "fill") +
  scale_fill_manual("", values = col.set$dermis) +
  labs(title = "Skin depth") +
  theme_bw() +
  theme(axis.line.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) +
  scale_x_discrete(limits=rev, breaks = c(0, 0.1, 1)) + #TODO
  coord_flip()

#fractions_leiden13

### Combination
clustersize_fraction <- ggarrange(
  leiden_legend,
  clustersize_leiden13,
  fractions_leiden13,
  common.legend = T, legend = "none",
  widths = c(1, 3.5, 3),
  nrow = 1, ncol = 3)

annotate_figure(clustersize_fraction,
                top = text_grob("Leiden clusters (r = 1.3)", color = "black", size = 14),
                fig.lab = "C", fig.lab.face = "bold", fig.lab.size = 14
)


##Figure 2 all together #####

all_umaps_spatial <- ggarrange(
  ggarrange(
    ggarrange(
      UMAP_BC_spottype + theme(legend.position = "none") + labs(title = "\n Spot type"),
      manual_spatial_plots + theme(axis.text.x=element_blank()), 
      UMAP_BC_leiden1.3 + theme(legend.position = "none") + labs(title = "\n Leiden r = 1.3"),
      leiden_spatial_plots,
      nrow = 2, ncol = 2, widths = c(0.82,2), heights = c(1,1), labels = c("A", "B", 'C', 'D')),
    ggarrange(NULL, clustersize_fraction, NULL,
              nrow = 3, ncol = 1,
              heights = c(0.05, 1, 0.05), labels = c('E')),
    nrow = 1, ncol = 2,
    widths = c(6,2)),
  as_ggplot(get_legend(UMAP_BC_spottype + geom_point(size = 5) +
                         theme(legend.position = "bottom",
                               legend.direction = "horizontal") +
                         guides(colour=guide_legend(nrow=2, byrow=TRUE)))),
  nrow = 2, ncol = 1,
  heights = c(1, 0.05)
)


setwd("/Users/mendenlab/Documents/Master_Thesis_TwomeyL/Figures/")
pdf('UMAPs_spatial.pdf', width = 18, height = 10)
all_umaps_spatial
dev.off()

ggarrange(
  ggarrange(
  UMAP_BC_spottype + theme(legend.position = "none") + labs(title = "\n Spot type"),
  manual_spatial_plots + theme(axis.text.x=element_blank()), nrow = 1, ncol = 2, widths = c(0.82,2), heights = c(1)),
  as_ggplot(get_legend(UMAP_BC_spottype + geom_point(size = 5) +
                         theme(legend.position = "bottom",
                               legend.direction = "horizontal") +
                         guides(colour=guide_legend(nrow=2, byrow=TRUE)))), nrow = 2, ncol = 1,
  heights = c(1, 0.05)
)


### Combination of 4
pdf('clusters_leiden_13.pdf', width = 18, height = 10)
clustersize_fraction_disease <- ggarrange(
  leiden_legend,
  clustersize_leiden13,
  fractions_leiden13,
  disease_leiden13,
  skin_leiden13,
  common.legend = T, legend = "none",
  widths = c(0.5, 3, 2.5, 2.5, 2.5),
  nrow = 1, ncol = 5)

annotate_figure(clustersize_fraction_disease,
                top = text_grob("Leiden clusters (r = 1.3)", color = "black", size = 14),
                fig.lab = "C", fig.lab.face = "bold", fig.lab.size = 14
)
dev.off()


as_ggplot(get_legend(fractions_leiden13 +
                       theme(legend.position = "bottom",
                             legend.direction = "vertical") +
                       guides(colour=guide_legend(nrow=2, byrow=TRUE))))

as_ggplot(get_legend(disease_leiden13 +
                       theme(legend.position = "bottom",
                             legend.direction = "vertical") +
                       guides(colour=guide_legend(nrow=3, ncol = 2, byrow=FALSE))))

as_ggplot(get_legend(skin_leiden13 +
                       theme(legend.position = "bottom",
                             legend.direction = "vertical") +
                       guides(colour=guide_legend(nrow=2, byrow=TRUE))))

ggarrange(disease_leiden13, skin_leiden13, legend = 'right')


ggarrange(
  ggarrange(
    UMAP_BC_leiden1.3 + theme(legend.position = "none") + labs(title = "\n Leiden clustering"),
    leiden_spatial_plots + theme(axis.text.x=element_blank()), nrow = 1, ncol = 2, widths = c(0.82,2), heights = c(1)),
  as_ggplot(get_legend(UMAP_BC_leiden1.3 + geom_point(size = 5) +
                         theme(legend.position = "bottom",
                               legend.direction = "horizontal") +
                         guides(colour=guide_legend(nrow=2, byrow=TRUE)))), nrow = 2, ncol = 1,
  heights = c(1, 1)
)


# sAVE UMAPs and heatmap separately
setwd("/Users/mendenlab/Documents/Master_Thesis_TwomeyL/Figures/")
pdf('UMAPskmeans.pdf', width = 10, height = 10)
all_UMAPs
dev.off()

pdf('kmeans_heatmap.pdf', width = 6, height = 8)
ggarrange(heatmap$gtable)
dev.off()

pdf('Figure3_SpatialDE/Heatmap_UMAPs_merged.pdf', width = 16, height = 12)
ggarrange(heatmap$gtable, all_UMAPs, 
          labels = c("D", "E"), font.label = list(size = 30, color = "black"),
          ncol = 2, nrow = 1, widths = c(1, 1), vjust = 0.05, hjust = 1.6) +
  theme(plot.margin = margin(0.9,0.1,0.1,1.5, "cm")) 

dev.off()