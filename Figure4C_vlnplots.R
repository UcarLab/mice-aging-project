
#Script to generate expression plots using seurat object.
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(gridBase)
library(gridExtra)
library(scater)

dir.create("Violin_plots")
args <- commandArgs(TRUE)

#args[1] = "/projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/RDS_files/seuratObject_rename_ident_May26_2021.rds"
#args[2] = "Pdk4"
mydata_orig = readRDS(args[1])
gene_name = args[2]
celltypes_list = args[3]

celltypes_subset <- c("Bcells", "Tcells", "Macrophages")

#mydata = readRDS("../seuratObject_rename_ident.rds")
mydata <- subset(mydata_orig, idents = celltypes_subset)

mydata$celltype.stim <- paste(Idents(mydata), mydata$age, sep = "_")
mydata$celltype <- Idents(mydata)
Idents(mydata) <- "celltype.stim"

#celltypes <- c("Bcells", "Tcells", "Macrophages", "NK_cells")

df_combined = data.frame()
plist_young <- list()
plist_old <- list()

for (i in celltypes_subset) {

print(i)
young_celltype = paste0(i, "_3m")
old_celltype = paste0(i, "_18m")

cell_young_old <- list(young_celltype, old_celltype)

mydata_celltype_Young <- subset(mydata, idents = young_celltype)
mydata_celltype_Old <- subset(mydata, idents = old_celltype)
mydata_new <- subset(mydata, subset = Fos > 0)
mydata_celltype_Young_no_zero <- subset(mydata_new, idents = young_celltype)
mydata_celltype_Old_no_zero <- subset(mydata_new, idents = old_celltype)

subset_celltype <- subset(mydata_new, idents = cell_young_old)

exprn_val <- subset(subset_celltype[["RNA"]]@data, rownames(subset_celltype[["RNA"]]@data) %in% gene_name)
exprn_val_df <- as.data.frame(exprn_val)
rownames(exprn_val_df) <- colnames(subset_celltype[["RNA"]]@data)
celltype_age <- subset_celltype@meta.data$celltype.stim
celltype <- subset_celltype@meta.data$celltype
#exprn_df_new <- as.data.frame(t(exprn_val_df))
df_final_AV <- cbind(exprn_val_df, celltype_age, celltype)
names(df_final_AV) <- c("exprn", "celltype_age", "celltype")

Celltype = i
var_18m <- paste0(Celltype, "_18m")
var_3m <- paste0(Celltype, "_3m")

df_final_AV_18m <- df_final_AV[df_final_AV$celltype_age == var_18m,]
df_final_AV_3m <- df_final_AV[df_final_AV$celltype_age == var_3m,]
test_sig <- t.test(df_final_AV_18m$exprn, df_final_AV_3m$exprn)
print(Celltype)
print(test_sig)
print(test_sig$p.value)

#######---------------------############
#df_table <- table(Idents(subset_LumAV))
m1 = ncol(mydata_celltype_Young_no_zero)
m2 = ncol(mydata_celltype_Old_no_zero)
m1_remain = ncol(mydata_celltype_Young) - m1
m2_remain = ncol(mydata_celltype_Old) - m2

df_3M_AV <- data.frame("Luminal AV cells" = c("Expressed", "Not expressed"), "Value" = c(m1,m1_remain))
df_18M_AV <- data.frame("Luminal AV cells" = c("Expressed", "Not expressed"), "Value" = c(m2,m2_remain))

##############################################################3
df_combined <- rbind(df_combined, df_final_AV)

bp1<- ggplot(df_3M_AV, aes(x="", y=Value, fill=Luminal.AV.cells)) + geom_bar(width = 0.2, stat = "identity") + theme_void()+theme(legend.position = "none")
pie1 <- bp1 + coord_polar("y", start=0) + scale_fill_grey()
bp2<- ggplot(df_18M_AV, aes(x="", y=Value, fill=Luminal.AV.cells)) + geom_bar(width = 0.2, stat = "identity") + theme_void()+ theme(legend.position = "none")
pie2 <- bp2 + coord_polar("y", start=0) + scale_fill_grey()

#print("pie1")
#print(class(pie1))
####################--------------------################

plist_young[[i]] <- pie1 
plist_old[[i]] <- pie2

}

library(gridExtra)
filename_piechart <- paste0("Piecharts_", gene_name, ".pdf", sep="")

pdf(filename_piechart, width=30)
grid.arrange(grobs=list(plist_young[[1]], plist_old[[1]], plist_young[[2]], plist_old[[2]], plist_young[[3]], plist_old[[3]], plist_young[[4]], plist_old[[4]]), nrow=1, top = "Percent of cells expressing gene")
dev.off()

#level_order <- c("Bcells_3m", "Bcells_18m", "Tcells_3m", "Tcells_18m", "Macrophages_3m", "Macrophages_18m", "NK_cells_3m", "NK_cells_18m")
write.table(df_combined, "DF_combined.txt", sep="\t", quote=FALSE)

filename_violinplot <- paste0("Violin_mice_celltypes_", gene_name, ".pdf", sep="")
pdf(filename_violinplot, width=166, height=40)
vlnplot_n <- df_combined %>%
	ggplot(aes(y=exprn, x = factor(celltype_age, level = level_order), fill = factor(celltype_age, level = level_order))) +
        geom_violin(position = position_dodge(0.8), alpha=0.5) +
        scale_fill_manual(values=c("lightskyblue1", "lightskyblue3", "lavenderblush1", "lavenderblush3", "indianred1", "indianred3", "mediumpurple1", "mediumpurple3")) + geom_jitter(shape=16, position=position_jitter(0.2),cex=5) +
        xlab("Sample")+ylab("Expression")+
        theme(plot.margin = unit(c(4,4,20,10), "cm"))+ggtitle(gene_name)+theme(plot.title = element_text(hjust = 0.5, size=80))+ theme(axis.text.x=element_text(size=60, angle=90), axis.text.y = element_text(size = 60), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60), legend.title=element_text(size=60), legend.text=element_text(size=60))+theme(legend.key.size = unit(2, "cm")) + theme(panel.background = element_blank(), axis.line = element_line(color = "black")) 

###################3--------------------################
ggdraw()+draw_plot(vlnplot_n, x=0, y=0, width=1, height=1)
dev.off()

#pdf("Plot_piecharts.pdf")
#do.call(grid.arrange, c(plist_young))
#dev.off()

##------------------------##################

