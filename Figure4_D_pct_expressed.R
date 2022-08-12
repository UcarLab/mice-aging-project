#Script to calculate heatmap of percent expressed cells.
#Author : Neerja Katiyar

library(ggpubr)
library(pheatmap)
library(tidyr)
library(Seurat)
library(RColorBrewer)

spleen.subset <- readRDS("../spleen_seurat_subset_rename_ident_Oct24.rds")
cell.types <- unique(Idents(spleen.subset))

spleen.subset$cluster_name <- Idents(spleen.subset)
cell.types <- unique(Idents(spleen.subset))

###################------------------------########################
df_combined <- data.frame()

for (i in cell.types) {
print(i)
data_celltype <- subset(spleen.subset@meta.data, cluster_name == i)
spleen.subset_celltype <- subset(spleen.subset, ident = i)

spleen_subset_celltype_young <- subset(spleen.subset_celltype, age == "3m")
spleen_subset_celltype_old <- subset(spleen.subset_celltype, age == "18m")
spleen_subset_celltype_young_gene <- subset(spleen_subset_celltype_young, Fosb > 0)
spleen_subset_celltype_old_gene <- subset(spleen_subset_celltype_old, Fosb > 0)

spleen_subset_celltype_young_gene_pct <- ncol(spleen_subset_celltype_young_gene)/ncol(spleen_subset_celltype_young)
spleen_subset_celltype_old_gene_pct <- ncol(spleen_subset_celltype_old_gene)/ncol(spleen_subset_celltype_old)

age_table <- data.frame(Young = spleen_subset_celltype_young_gene_pct, Aged = spleen_subset_celltype_old_gene_pct)
rownames(age_table) = i

df_combined <- rbind(df_combined, age_table)
}

print(df_combined)

df_mat <- as.matrix(df_combined)
pheatmap.val <- df_mat
breaksList = seq(0, max(pheatmap.val*100), by = .001)

gene = "Fosb"
print(gene)

pdf(paste0("Pheatmap_new_cell_type_vs_age_", gene, ".pdf"), width = 4)

pheatmap::pheatmap(pheatmap.val * 100, cluster_cols = FALSE, cluster_rows = FALSE,
        color =  colorRampPalette((brewer.pal(n = 9, name = "Blues")))(length(breaksList)),
        breaks = breaksList, border_color = "white", main = gene, display_numbers = T,
        fontsize_number = 9,
        number_color = "white", legend = F)
dev.off()

#######################------------------------------------#########################3

###################------------------------########################
df_combined <- data.frame()

for (i in cell.types) {

data_celltype <- subset(spleen.subset@meta.data, cluster_name == i)
spleen.subset_celltype <- subset(spleen.subset, ident = i)

spleen_subset_celltype_young <- subset(spleen.subset_celltype, age == "3m")
spleen_subset_celltype_old <- subset(spleen.subset_celltype, age == "18m")
spleen_subset_celltype_young_gene <- subset(spleen_subset_celltype_young, Fos > 0)
spleen_subset_celltype_old_gene <- subset(spleen_subset_celltype_old, Fos > 0)

spleen_subset_celltype_young_gene_pct <- ncol(spleen_subset_celltype_young_gene)/ncol(spleen_subset_celltype_young)
spleen_subset_celltype_old_gene_pct <- ncol(spleen_subset_celltype_old_gene)/ncol(spleen_subset_celltype_old)

age_table <- data.frame(Young = spleen_subset_celltype_young_gene_pct, Aged = spleen_subset_celltype_old_gene_pct)
rownames(age_table) = i

print(age_table)
df_combined <- rbind(df_combined, age_table)

}
print(df_combined)

df_mat <- as.matrix(df_combined)
pheatmap.val <- df_mat
breaksList = seq(0, max(pheatmap.val*100), by = .001)

gene = "Fos"
print(gene)

pdf(paste0("Pheatmap_new_cell_type_vs_age_", gene, ".pdf"), width = 4)

pheatmap::pheatmap(pheatmap.val * 100, cluster_cols = FALSE, cluster_rows = FALSE,
        color =  colorRampPalette((brewer.pal(n = 9, name = "Blues")))(length(breaksList)),
        breaks = breaksList, border_color = "white", main = gene, display_numbers = T,
        fontsize_number = 9,
        number_color = "white", legend = F)
dev.off()

#######################------------------------------------#########################3

###################------------------------########################
df_combined <- data.frame()

for (i in cell.types) {

data_celltype <- subset(spleen.subset@meta.data, cluster_name == i)
spleen.subset_celltype <- subset(spleen.subset, ident = i)

spleen_subset_celltype_young <- subset(spleen.subset_celltype, age == "3m")
spleen_subset_celltype_old <- subset(spleen.subset_celltype, age == "18m")
spleen_subset_celltype_young_gene <- subset(spleen_subset_celltype_young, Jun > 0)
spleen_subset_celltype_old_gene <- subset(spleen_subset_celltype_old, Jun > 0)

spleen_subset_celltype_young_gene_pct <- ncol(spleen_subset_celltype_young_gene)/ncol(spleen_subset_celltype_young)
spleen_subset_celltype_old_gene_pct <- ncol(spleen_subset_celltype_old_gene)/ncol(spleen_subset_celltype_old)

age_table <- data.frame(Young = spleen_subset_celltype_young_gene_pct, Aged = spleen_subset_celltype_old_gene_pct)
rownames(age_table) = i

print(age_table)
df_combined <- rbind(df_combined, age_table)

}

print(df_combined)

df_mat <- as.matrix(df_combined)
pheatmap.val <- df_mat
breaksList = seq(0, max(pheatmap.val*100), by = .001)

gene = "Jun"
print(gene)

pdf(paste0("Pheatmap_new_cell_type_vs_age_", gene, ".pdf"), width = 4)

pheatmap::pheatmap(pheatmap.val * 100, cluster_cols = FALSE, cluster_rows = FALSE,
        color =  colorRampPalette((brewer.pal(n = 9, name = "Blues")))(length(breaksList)),
        breaks = breaksList, border_color = "white", main = gene, display_numbers = T,
        fontsize_number = 9,
        number_color = "white", legend = F)
dev.off()

#######################------------------------------------#########################3

###################------------------------########################
df_combined <- data.frame()

for (i in cell.types) {

data_celltype <- subset(spleen.subset@meta.data, cluster_name == i)
spleen.subset_celltype <- subset(spleen.subset, ident = i)

spleen_subset_celltype_young <- subset(spleen.subset_celltype, age == "3m")
spleen_subset_celltype_old <- subset(spleen.subset_celltype, age == "18m")
spleen_subset_celltype_young_gene <- subset(spleen_subset_celltype_young, Fosl2 > 0)
spleen_subset_celltype_old_gene <- subset(spleen_subset_celltype_old, Fosl2 > 0)

spleen_subset_celltype_young_gene_pct <- ncol(spleen_subset_celltype_young_gene)/ncol(spleen_subset_celltype_young)
spleen_subset_celltype_old_gene_pct <- ncol(spleen_subset_celltype_old_gene)/ncol(spleen_subset_celltype_old)

age_table <- data.frame(Young = spleen_subset_celltype_young_gene_pct, Aged = spleen_subset_celltype_old_gene_pct)
rownames(age_table) = i

print(age_table)
df_combined <- rbind(df_combined, age_table)

}

print(df_combined)

df_mat <- as.matrix(df_combined)
pheatmap.val <- df_mat
breaksList = seq(0, max(pheatmap.val*100), by = .001)

gene = "Fosl2"
print(gene)

pdf(paste0("Pheatmap_new_cell_type_vs_age_", gene, ".pdf"), width = 4)

pheatmap::pheatmap(pheatmap.val * 100, cluster_cols = FALSE, cluster_rows = FALSE,
        color =  colorRampPalette((brewer.pal(n = 9, name = "Blues")))(length(breaksList)),
        breaks = breaksList, border_color = "white", main = gene, display_numbers = T,
        fontsize_number = 9,
        number_color = "white", legend = F)
dev.off()

#######################------------------------------------#########################3

