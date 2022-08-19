
library(ggplot2)
library(ggpubr)
library(rstatix)

#source("../../../../../../../code/color_values.R")
ens2gene <- function(genes){
  genome <- annotables::grcm38

  # ens to gene symbol mapping
  mapping <-
    base::subset(genome,
                 genome$ensgene %in% genes,
                 select = c('ensgene', 'symbol'))

    m <- match(genes, mapping$ensgene)

    ens.genes <- mapping$symbol[m]
    names(ens.genes) <- genes

    return(ens.genes)
}

data_AP1_all <- read.table("../Table_B6_NZO_FC_RNA_AP1_all_genes.txt", sep="\t", header=T, row.names=1)

data_AP1_family <- data_AP1_all[data_AP1_all$tissue == "PBL",]

data_AP1_family_Jun <- data_AP1_family[data_AP1_family$gene == "Jun",]
data_AP1_family_Junb <- data_AP1_family[data_AP1_family$gene == "Junb",]
data_AP1_family_Jund <- data_AP1_family[data_AP1_family$gene == "Jund",]
data_AP1_family_Fos <- data_AP1_family[data_AP1_family$gene == "Fos",]
data_AP1_family_Fosb <- data_AP1_family[data_AP1_family$gene == "Fosb",]
data_AP1_family_Fra1 <- data_AP1_family[data_AP1_family$gene == "Fra1",]
data_AP1_family_Fra2 <- data_AP1_family[data_AP1_family$gene == "Fra2",]
data_AP1_family_Atf <- data_AP1_family[data_AP1_family$gene == "Atf",]
data_AP1_family_Atf2 <- data_AP1_family[data_AP1_family$gene == "Atf2",]
data_AP1_family_Atf3 <- data_AP1_family[data_AP1_family$gene == "Atf3",]
data_AP1_family_Atf4 <- data_AP1_family[data_AP1_family$gene == "Atf4",]
data_AP1_family_Atf5 <- data_AP1_family[data_AP1_family$gene == "Atf5",]
data_AP1_family_Atf6 <- data_AP1_family[data_AP1_family$gene == "Atf6",]
data_AP1_family_Atf6b <- data_AP1_family[data_AP1_family$gene == "Atf6b",]
data_AP1_family_Atf7 <- data_AP1_family[data_AP1_family$gene == "Atf7",]
data_AP1_family_Batf <- data_AP1_family[data_AP1_family$gene == "Batf",]
data_AP1_family_Batf2 <- data_AP1_family[data_AP1_family$gene == "Batf2",]
data_AP1_family_Batf3 <- data_AP1_family[data_AP1_family$gene == "Batf3",]
data_AP1_family_Maf <- data_AP1_family[data_AP1_family$gene == "Maf",]
data_AP1_family_Mafa <- data_AP1_family[data_AP1_family$gene == "Mafa",]
data_AP1_family_Mafb <- data_AP1_family[data_AP1_family$gene == "Mafb",]
data_AP1_family_Maff <- data_AP1_family[data_AP1_family$gene == "Maff",]
data_AP1_family_Mafg <- data_AP1_family[data_AP1_family$gene == "Mafg",]
data_AP1_family_Mafk <- data_AP1_family[data_AP1_family$gene == "Mafk",]

######################################################################

data_AP1_family_Jun$AP1_subfamily <- c(rep("Jun"))
data_AP1_family_Junb$AP1_subfamily <- c(rep("Jun"))
data_AP1_family_Jund$AP1_subfamily <- c(rep("Jun"))
data_AP1_family_Fos$AP1_subfamily <- c(rep("Fos"))
data_AP1_family_Fosb$AP1_subfamily <- c(rep("Fos"))
#data_AP1_family_Fra1$AP1_subfamily <- c(rep("Fos"))
#data_AP1_family_Fra2$AP1_subfamily <- c(rep("Fos"))
#data_AP1_family_Atf$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Atf2$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Atf3$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Atf4$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Atf5$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Atf6$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Atf6b$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Atf7$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Batf$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Batf2$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Batf3$AP1_subfamily <- c(rep("Atf"))
data_AP1_family_Maf$AP1_subfamily <- c(rep("Maf"))
data_AP1_family_Mafa$AP1_subfamily <- c(rep("Maf"))
data_AP1_family_Mafb$AP1_subfamily <- c(rep("Maf"))
data_AP1_family_Maff$AP1_subfamily <- c(rep("Maf"))
data_AP1_family_Mafg$AP1_subfamily <- c(rep("Maf"))
data_AP1_family_Mafk$AP1_subfamily <- c(rep("Maf"))

df_AP1_family_all <- rbind(data_AP1_family_Jun, data_AP1_family_Junb, data_AP1_family_Jund, data_AP1_family_Fos, data_AP1_family_Fosb, data_AP1_family_Atf2, data_AP1_family_Atf3, data_AP1_family_Atf4, data_AP1_family_Atf5, data_AP1_family_Atf6, data_AP1_family_Atf6b, data_AP1_family_Atf7, data_AP1_family_Batf, data_AP1_family_Batf2, data_AP1_family_Batf3, data_AP1_family_Maf, data_AP1_family_Mafa, data_AP1_family_Mafb, data_AP1_family_Maff, data_AP1_family_Mafb, data_AP1_family_Mafg, data_AP1_family_Mafk)

write.table(df_AP1_family_all, "PBL_FC_AP1_all.txt", sep="\t", quote=FALSE)

#########################################################

df <- df_AP1_family_all

df_B6 <- df[df$strain == "B6",]
df_NZO <- df[df$strain == "NZO",]

data_Jun_B6 <- df_B6[df_B6$AP1_subfamily == "Jun",]
data_Jun_NZO <- df_NZO[df_NZO$AP1_subfamily == "Jun",]

data_Fos_B6 <- df_B6[df_B6$AP1_subfamily == "Fos",]
data_Fos_NZO <- df_NZO[df_NZO$AP1_subfamily == "Fos",]

data_Atf_B6 <- df_B6[df_B6$AP1_subfamily == "Atf",]
data_Atf_NZO <- df_NZO[df_NZO$AP1_subfamily == "Atf",]

data_Maf_B6 <- df_B6[df_B6$AP1_subfamily == "Maf",]
data_Maf_NZO <- df_NZO[df_NZO$AP1_subfamily == "Maf",]

Jun_B6_test <- t.test(data_Jun_B6$FC, mu=0, alternative = "greater")
Jun_NZO_test <- t.test(data_Jun_NZO$FC, mu=0, alternative = "greater")
Fos_B6_test <- t.test(data_Fos_B6$FC, mu=0, alternative = "greater")
Fos_NZO_test <- t.test(data_Fos_NZO$FC, mu=0, alternative = "greater")
Atf_B6_test <- t.test(data_Atf_B6$FC, mu=0, alternative = "greater")
Atf_NZO_test <- t.test(data_Atf_NZO$FC, mu=0, alternative = "greater")
Maf_B6_test <- t.test(data_Maf_B6$FC, mu=0, alternative = "greater")
Maf_NZO_test <- t.test(data_Maf_NZO$FC, mu=0, alternative = "greater")

################################################

print("Jun_B6_test...")
print(Jun_B6_test)
print("Jun_NZO_test...")
print(Jun_NZO_test)

print("Fos_B6_test...")
print(Fos_B6_test)
print("Fos_NZO_test...")
print(Fos_NZO_test)

print("Atf_B6_test...")
print(Atf_B6_test)
print("Atf_NZO_test...")
print(Atf_NZO_test)

print("Maf_B6_test...")
print(Maf_B6_test)
print("Maf_NZO_test...")
print(Maf_NZO_test)

#One sample t-test

Jun_B6_wilcox_test <- wilcox.test(data_Jun_B6$FC, mu=0, alternative = "greater")
Jun_NZO_wilcox_test <- wilcox.test(data_Jun_NZO$FC, mu=0, alternative = "greater")
Fos_B6_wilcox_test <- wilcox.test(data_Fos_B6$FC, mu=0, alternative = "greater")
Fos_NZO_wilcox_test <- wilcox.test(data_Fos_NZO$FC, mu=0, alternative = "greater")
Atf_B6_wilcox_test <- wilcox.test(data_Atf_B6$FC, mu=0, alternative = "greater")
Atf_NZO_wilcox_test <- wilcox.test(data_Atf_NZO$FC, mu=0, alternative = "greater")
Maf_B6_wilcox_test <- wilcox.test(data_Maf_B6$FC, mu=0, alternative = "greater")
Maf_NZO_wilcox_test <- wilcox.test(data_Maf_NZO$FC, mu=0, alternative = "greater")

print("Jun_B6_wilcox_test...")
print(Jun_B6_wilcox_test)
print("Jun_NZO_wilcox_test...")
print(Jun_NZO_wilcox_test)

print("Fos_B6_wilcox_test...")
print(Fos_B6_wilcox_test)
print("Fos_NZO_wilcox_test...")
print(Fos_NZO_wilcox_test)

print("Atf_B6_wilcox_test...")
print(Atf_B6_wilcox_test)
print("Atf_NZO_wilcox_test...")
print(Atf_NZO_wilcox_test)

print("Maf_B6_wilcox_test...")
print(Maf_B6_wilcox_test)
print("Maf_NZO_wilcox_test...")
print(Maf_NZO_wilcox_test)

p <- ggboxplot(df_AP1_family_all, x="AP1_subfamily", y="FC", color = "strain", palette = "jco", add="jitter", facet.by = "tissue", nrow=1, short.panel.labs = FALSE)  

pdf("Boxplot_PBL_FC_RNA_DE_not_filtered_all_24_AP1_genes.pdf", height=10, width=10)
p
dev.off()

