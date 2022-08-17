
#title - Gene set enrichment analysis using AP1 geneset.
#Author : Neerja Katiyar

library(fgsea)
library(data.table)
library(ggplot2)

# Modified function definition of plotEnrichment function

plotEnrichment_new <- function(pathway, stats, gseaParam = 1, ticksSize = 0.2) 
{
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    pathway_chk <- statsAdj[pathway]
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
        returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = "green", 
        size = 0.1) + geom_hline(yintercept = max(tops), colour = "red", 
        linetype = "dashed") + geom_hline(yintercept = min(bottoms), 
        colour = "red", linetype = "dashed") + geom_hline(yintercept = 0, 
        colour = "black") + geom_line(color = "green") + theme_bw() + 
        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
            y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + 
        theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
        labs(x = "rank", y = "enrichment score")
    return(g)
}

####################################################################
MAG_tissues <- read.table("../../MAG_RNAseq_analysis/MAG_tissues_all.txt", sep="\t", header=T)
MAG_tissues_sorted <- MAG_tissues[order(MAG_tissues$Score, decreasing=TRUE),]

####################################################################

#genelist_Taa <- c("Gzmk", "Ccl5", "Ifng", "Pdcd1", "Tigit", "Lag3", "Tox")
inflammatory_orig_file <- read.table("../Mice_pathways/nanostring_v2_mouse_inflammatory_markers.txt", sep="\t", header=T)
gene <- inflammatory_orig_file$GeneName
inflammatory_list <- as.list(as.data.frame(gene))
#names(inflammatory_list) <- "inflammatory_nanostring"

genelist_AP1 <- c("Jdp2", "Atf6", "Maf", "Batf3", "Mafb", "Mafa", "Mafg", "Atf5", "Atf3", "Batf2", "Atf2", "Atf7", "Atf6b", "Batf", "Junb", "Fos", "Jun", "Fosb", "Jund", "Atf4", "Maff", "Mafk")

genesets <- list(inflammatory_list$gene, genelist_AP1)
names(genesets) <- c("Inflammatory_list", "AP1_geneset")

genesets <- genesets[!duplicated(names(genesets))]

#MAG_tissues_sorted$gene <- ens2gene(file_mem_B6_sorted$Gene.Ensembl)
#file_mem_NZO_sorted$gene <- ens2gene(file_mem_NZO_sorted$Gene.Ensembl)

####################################################################

MAG_sorted_inp <- MAG_tissues_sorted$Score
names(MAG_sorted_inp) <- MAG_tissues_sorted$gene

print("Starting fgsea....")
####################################################################

fgsea_MAG <- fgsea(pathways = genesets, stats = MAG_sorted_inp)
print(length(MAG_sorted_inp))
print(head(MAG_sorted_inp))
fgsea_MAG_df <- as.data.frame(apply(fgsea_MAG,2,as.character))
fgsea_MAG_df_sorted <- fgsea_MAG_df[order(fgsea_MAG_df$padj),]
write.table(fgsea_MAG_df, "GSEA_fgsea_MAG_inflammatory_AP1_out.txt", sep="\t", quote=FALSE, row.names=FALSE)

####################################################################

pdf("GSEA_fgsea_MAG_inflammatory_enrichplot_top.pdf")
plotEnrichment(genesets[[head(fgsea_MAG_df_sorted, 1)$pathway]], MAG_sorted_inp) + labs(title=head(fgsea_MAG_df_sorted,1)$pathway)
dev.off()

pdf("GSEA_fgsea_MAG_AP1_enrichplot.pdf")
plotEnrichment(genesets[[tail(fgsea_MAG_df_sorted, 1)$pathway]], MAG_sorted_inp) + labs(title=tail(fgsea_MAG_df_sorted, 1)$pathway)
dev.off()

##################################################

