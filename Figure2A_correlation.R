
library("ggpubr")
library(VennDiagram)
library(gplots)
library(dplyr)

Celltypes <- c("MEMORY", "NAIVE", "PBL", "SPLEEN")

ens2gene <- function(genes){

  genome <- annotables::grcm38

  # ens to gene symbol mapping
  mapping <-
    base::subset(genome,
                 genome$ensgene %in% genes,
                 select = c('ensgene', 'symbol'))

    m <- match(genes, mapping$ensgene)

    sym.genes <- mapping$symbol[m]
    names(sym.genes) <- genes

    return(sym.genes)

}

ap1_genes = c("Jdp2", "Atf6", "Maf", "Batf3", "Mafb", "Mafa", "Mafg", "Atf5", "Atf3", "Batf2", "Atf2", "Atf7", "Atf6b", "Batf", "Junb", "Fos", "Jun", "Fosb", "Jund", "Atf4", "Maff", "Mafk")
#gzmk_genes = c("Gzmk", "Ccl5", "Ifng", "Pdcd1", "Tigit", "Lag3", "Tox")

for (i in Celltypes) {
	fileB6 <- paste0("../../", i, "_Age18vs3_B6_RNAseq.csv")
	fileNZO <- paste0("../../", i, "_Age18vs3_NZO_RNAseq.csv")
	
	B6_file_orig <- read.delim(fileB6, sep=",", quote = "\"", header=TRUE)
	#B6_file_orig$gene <- ens2gene(B6_file_orig$Gene.Ensembl)
	
	NZO_file_orig <- read.delim(fileNZO, sep=",", quote = "\"", header=TRUE)
	#NZO_file_orig$gene <- ens2gene(NZO_file_orig$Gene.Ensembl)
	
	print(dim(B6_file_orig))
	print(dim(NZO_file_orig))
	
	B6_NZO_file_orig <- merge(B6_file_orig, NZO_file_orig, by = "Gene.Ensembl")

        B6_NZO_file <- B6_NZO_file_orig[(((B6_NZO_file_orig$adj.P.Val.x<0.05) & (abs(B6_NZO_file_orig$logFC.x) > 1)) | ((B6_NZO_file_orig$adj.P.Val.y<0.05) & (abs(B6_NZO_file_orig$logFC.y)>1))),]
	print(head(B6_NZO_file))
	B6_NZO_file$gene <- ens2gene(B6_NZO_file$Gene.Ensembl)

	#B6_NZO_file$color_index <- if(B6_NZO_file$gene %in% ap1_genes == "TRUE", "red")
	#B6_NZO_file_subset <- subset(B6_NZO_file, color_index == "grey")
	#B6_NZO_file$color_index <- if(B6_NZO_file_subset$gene %in% gzmk_genes == "TRUE", "blue")
	B6_NZO_file$color_val <- "grey"
	
	B6_NZO_file <- B6_NZO_file%>%mutate(color_val = case_when(B6_NZO_file$gene %in% ap1_genes == "TRUE" ~ "red", TRUE~color_val))

	#ap1_gzmk_comb <- c(ap1_genes)
	
	#B6_NZO_file$color_val <- ifelse(B6_NZO_file$color_index=="TRUE", "red", "grey")
	
	print(i)
	print(dim(B6_NZO_file))
	print(head(B6_NZO_file))
	
	filename = paste0(i, "_B6_NZO_merged.txt")
	write.table(B6_NZO_file, filename, sep="\t", quote=FALSE)	

	cor_calc <- cor.test(B6_NZO_file$logFC.x, B6_NZO_file$logFC.y)
	print("Correlation test between B6 vs NZO for memory cells")
	print(cor_calc)
	plot_filename <- paste0("Plot_FC_", i, "_B6_NZO_RNASeq.pdf")
	print(plot_filename)
	pdf(plot_filename)
	print(ggscatter(B6_NZO_file, x = "logFC.x", y = "logFC.y", color = "color_val", label = "gene", repel = TRUE, label.select = ap1_genes, palette = c(grey = "grey", red = "red", blue = "blue"), add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", title = plot_filename, xlab = "FC (RNASeq) for B6", ylab = "FC (RNAseq) for NZO"))
	#print(ggscatter(B6_NZO_file, x = "logFC.x", y = "logFC.y", label = B6_NZO_file$gene, repel = TRUE, label.select = ap1_genes, color = "color_val", palette = c("red", "grey"), font.label = "black", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "FC (RNASeq) for B6", ylab = "FC (RNAseq) for NZO"))
	dev.off()
	#file_venn = paste0("Venn_diag_",i,"_B6_NZO.png")
	#venn.diagram(x=list(B6_file$Gene.Ensembl, NZO_file$Gene.Ensembl),
	#category.names = c("B6", "NZO"),
	#filename = file_venn,
	#output = TRUE)

	#x = list(B6_file$Gene.Ensembl, NZO_file$Gene.Ensembl)
	#v.table = venn(x)
	#intersections = attr(x = v.table, "intersections")
	#intersect_B6_only = intersections$`A`
	#intersect_NZO_only = intersections$`B`
	#intersect_B6_NZO = intersections$`A:B`
	#df <- data.frame(B6 = c(length(intersect_B6_NZO), length(intersect_B6_only)), nonB6 = c(length(intersect_NZO_only), 0))
	#rownames(df) <- c("NZO", "nonNZO")
	#print(df)
	#test_now <- chisq.test(df)
	#print(test_now)	
	#print(test_now$p.value)
}

