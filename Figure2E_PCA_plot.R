
#title - PCA plots from RNA seq data using top AP1 genes.
#Author - Onur, Neerja Katiyar

library(readr)
library(VennDiagram)
library(ggplot2)
library(grDevices)
library(edgeR)
library(dplyr)

source("../code/color_values.R")

gene2ens <- function(genes){
  
  genome <- annotables::grcm38
  
  # ens to gene symbol mapping
  mapping <-
    base::subset(genome,
                 genome$symbol %in% genes,
                 select = c('ensgene', 'symbol'))
  
    m <- match(genes, mapping$symbol)
    
    ens.genes <- mapping$ensgene[m]
    names(ens.genes) <- genes
    
    return(ens.genes)
}

ap1.genes <- c("Fos", "Fosb", "Jun", "Junb", "Maff")

length(ap1.genes)

ap1.genes.ens <- gene2ens(ap1.genes)

## Load Onur's data
count.matrix <- read.csv('../rna_count_matrix.csv', row.names = 1, check.names = F, stringsAsFactors = F)

# Order genes according to their standard deviation in decreasing order
count.matrix <- count.matrix [rev(order(apply(count.matrix, 1, sd))),]

# Remove duplicated genes
count.matrix <- count.matrix [!duplicated(rownames(count.matrix)),]


# Enforce all counts to be integers
count.matrix <- round(count.matrix, 0)

# remove low expressed genes      
count.matrix <- count.matrix [rowSums(cpm(count.matrix) >= 0.5) >= 2,]

# filter BM
BM.loc <- colnames(count.matrix) %>% sapply(function(x){grepl("BM", x, fixed = TRUE)})
count.matrix <- count.matrix[,!BM.loc]

# normalize with cpm
count.matrix.normalized <- cpm(count.matrix, log = T)

meta.data <- colnames(count.matrix) %>% strsplit("-", fixed = T) %>% do.call(rbind, .) %>% as.data.frame

colnames(meta.data) <- c("Strain", "Age", "Sex", "TCT", "SampleID")

pca.plot <- function(x, overlaid.info, sample.names = NULL, show.names = TRUE, color.vals = NULL){
  
  if(is.null(sample.names)){
    sample.names <- colnames(x)
  } else{
    if(length(sample.names) != ncol(x)){
      stop("The length of `sample.names` should be equal to number of samples.")
    }
  }
  
  # eliminate NaN values before-hand if there is any.
  pca <- stats::prcomp(t(stats::na.omit(x)), center = TRUE)

  d  <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)
  xl <- sprintf("PC 1: %.1f %%", d[1])
  yl <- sprintf("PC 2: %.1f %%", d[2])


  plot.df <- data.frame(PC1 = as.numeric(pca$x[,1]),
                   PC2 = as.numeric(pca$x[,2]),
                   overlaid.info = overlaid.info,
                   names = sample.names
                   )

  plot.pca <- ggplot2::ggplot(plot.df, ggplot2::aes(PC1, PC2, color = overlaid.info)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::labs(x=xl,y=yl) +
    ggplot2::theme_minimal() +
    ggplot2::labs(color = "Status") +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_light()

  if (typeof(overlaid.info) %in% c("character", "factor")){
    if (!is.null(color.vals)){
      plot.pca <- plot.pca +
      ggplot2::scale_color_manual(values = color.vals)
    } else{
      plot.pca <- plot.pca +
        ggplot2::scale_color_manual(values = RColorBrewer::brewer.pal(n = 9, name = "Set1"))
    }
  }

  if(show.names){
    plot.pca <- plot.pca + ggrepel::geom_text_repel(ggplot2::aes(label = names))
  }
  
  return(plot.pca)
}

color.vals <- c("memory-3mo" = "#45b575",
                "memory-18mo"= "#1a9850",
                  "naive-3mo"= "#99d594",
                  "naive-18mo" = "#80e378",
                  "PBL-3mo" = "#d73027",
                  "PBL-18mo" = "#a62019",
                  "spleen-3mo"="#6e93c2",
                  "spleen-12mo"="#4575b4",
                  "spleen-18mo"="#265591")
meta.data$TCT_Age<- paste0(meta.data$TCT, "-", meta.data$Age)
meta.data$TCT[meta.data$TCT == "memory"] <- "CD8 memory"
meta.data$TCT[meta.data$TCT == "naive"] <- "CD8 naive"

# eliminate NaN values before-hand if there is any.
pca <- stats::prcomp(t(stats::na.omit(count.matrix.normalized[(na.omit(ap1.genes.ens)),])), center = TRUE)

d  <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)
xl <- sprintf("PC 1: %.1f %%", d[1])
yl <- sprintf("PC 2: %.1f %%", d[2])

plot.df <- data.frame(PC1 = as.numeric(pca$x[,1]),
                      PC2 = as.numeric(pca$x[,2]),
                      TCT = meta.data$TCT,
                      Age = gsub("mo", "", meta.data$Age)
                      )

plot.pca.ap1.mice <- ggplot(plot.df, aes(PC1, PC2, color = Age, shape = TCT)) +
    geom_point(size = 4) +
    labs(x=xl,y=yl) +
    theme_minimal() +
    labs(shape = "Tissue/Cell Type") +
    #coord_fixed(ratio = 1) +
    theme_light(base_size = 16) + 
    scale_color_manual(values = color_values) + 
    scale_shape_manual(values=c(15:18))

ggsave(filename = "plot_PCA_AP1_mice.pdf", 
       plot = plot.pca.ap1.mice, 
       useDingbats = FALSE, width = 7, height = 5)

plot.pca.ap1.mice.wo12 <- ggplot(plot.df %>% filter(Age != 12), 
                            aes(PC1, PC2, color = Age, shape = TCT)) +
    geom_point(size = 4) +
    labs(x=xl,y=yl) +
    theme_minimal() +
    labs(shape = "Tissue/Cell Type") +
    #coord_fixed(ratio = 1) +
    theme_light(base_size = 16) + 
    scale_color_manual(values = color_values) + 
    scale_shape_manual(values=c(15:18))

pdf("plot_PCA_AP1_mice_wo12.pdf", width=7, height=5)
print(plot.pca.ap1.mice.wo12)
dev.off()
#ggsave(filename = "plot_PCA_AP1_mice_wo12.pdf", 
#       plot = plot.pca.ap1.mice.wo12, 
#       useDingbats = FALSE, width = 7, height = 20)
