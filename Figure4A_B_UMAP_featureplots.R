#Title : Script to generate UMAP and featureplots from Tabula muris senis dataset.
#Author : Onur Karakaslar and Neerja Katiyar

library("Seurat")
library("BiocFileCache")
library(ggplot2)
library(ggpubr)
library(paletteer)

palette <- RColorBrewer::brewer.pal(9, "Reds")
palette[1] <- "#adadad"

#Visualize data with Nebulosa.
tm_obj <- readRDS("../spleen_seurat_subset_rename_ident_Oct24.rds")
summary(tm_obj@meta.data)

pdf("UMAP_Tabula_muris_split_age.pdf")
DimPlot(tm_obj, reduction = "umap", split.by = "age")
dev.off()

pdf("spleen_Jun_Fos_Fosb.pdf", width = 9, height = 10)
FeaturePlot(tm_obj, features = c("Jun", "Fos", "Fosb"), cols = palette, pt.size = 1, order = TRUE,  split.by = "age")
dev.off()

