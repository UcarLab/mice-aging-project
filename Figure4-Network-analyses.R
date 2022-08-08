library(ggplot2)
library(igraph)
setwd("/Users/karako/Dropbox (JAX)/MouseAging_clean/")


TCT <- "naive"
GSEA <- T
geneset2subset <- "wp"

if (GSEA){
  enrich.res <- read.csv(paste0("~/Dropbox (JAX)/MouseAging_clean/output/F3/Enrichment Files/GSEA/RNAseq/",TCT,"_gsea_summary.csv"))
  gene.loc <- "leadingEdge"
} else {
  enrich.res <- read.csv(paste0("~/Dropbox (JAX)/MouseAging_clean/output/F3/Enrichment Files/Hypergeometric/RNAseq/",TCT,"_er_summary.csv"))
  gene.loc <- "overlapping.genes"
}

gene.list  <- strsplit(enrich.res[,gene.loc],",", fixed =T)
ens.genes  <- do.call(c,gene.list)
ens.genes  <- unique(ens.genes)

genome <- annotables::grcm38

# entrez to gene symbol mapping
ens2gene <-
  base::subset(genome,
               genome$ensgene %in% ens.genes,
               select = c('ensgene', 'symbol'))

convert2symbol <- function(x){
  # Match to each annotation dataframe
  m <- match(x, ens2gene$ensgene)
  return(ens2gene$symbol[m])
}

df2graph <- function(graph.df, graph.meta, GSEA = F){

  g <- igraph::graph_from_data_frame(graph.df,
                                     directed = FALSE,
                                     vertices = NULL)
  sizes <- igraph::degree(g)
  size_label <- "# genes"

  # color of nodes

  if (GSEA){
    graph.meta$Status <- ifelse(graph.meta$NES > 0, "Up", "Down")
    V(g)$Status <- graph.meta$Status[match(V(g)$name, graph.meta$pathway)]
  } else {
    V(g)$Status <- graph.meta$Status[match(V(g)$name, graph.meta$module.name)]
  }
  V(g)$color <- ifelse(is.na(V(g)$Status),
                       "black",
                       ifelse(V(g)$Status == "Up",
                              "red", "blue"))


  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <- "gray"

  p <- ggraph::ggraph(g, layout = "stress")
  p <- p + ggraph::geom_edge_link(alpha = .8, colour = "darkgrey")
  p <- p + ggraph::geom_node_point(ggplot2::aes_(color = ~ I(color), size = ~size))
  p <- p + ggplot2::scale_size(range = c(4, 10),
                               breaks = round(seq(round(min(igraph::V(g)$size)),
                                                  round(max(igraph::V(g)$size)),
                                                  length.out = 4)),
                               name = size_label)

  p <- p + ggplot2::theme_void()
  p <- p + ggraph::geom_node_text(ggplot2::aes_(label = ifelse(V(g)$size > 0,
                                                               V(g)$name,"")),
                                  repel = T)

  return(p)
}

create_graphs <- function(enrich.res, geneset2subset = "vp2008", TCT = TCT, GSEA){

  if(GSEA){
    enrich.res.subset <- subset(enrich.res, Geneset == geneset2subset &
                                  padj < 0.05)
  } else {
    enrich.res.subset <- subset(enrich.res, geneset.name == geneset2subset)
  }


  # No p.val check needed for subsets since cut-offs are determined earlier
  enrich.res.subset.bl6 <- subset(enrich.res.subset, Contrast == "Age18vs3_B6")
  enrich.res.subset.nzo <- subset(enrich.res.subset, Contrast == "Age18vs3_NZO")

  prepare_subset_df <- function(df.subset, GSEA){
    if (GSEA){
      pathway.names <- df.subset$pathway
      overlapping.ens <- df.subset$leadingEdge
    } else {
      pathway.names <- df.subset$module.name
      overlapping.ens <- df.subset$overlapping.genes
    }


    overlapping.ens <- strsplit(overlapping.ens,",", fixed = T)
    overlapping.sym <- lapply(overlapping.ens, convert2symbol)
    names(overlapping.sym) <- pathway.names

    graph.mat <-
      do.call(rbind,
              lapply(pathway.names,
                     function(x) {
                       cbind(x, overlapping.sym[[x]])
                     }))

    graph.df <- as.data.frame(graph.mat)

    colnames(graph.df) <- c("Term", "Gene")

    return(graph.df)
  }

  if(GSEA){
    file.name.nzo <- paste0("./output/F4/Network_Analyses_GSEA_",
                            toupper(TCT), "_NZO_", geneset2subset,".pdf")
    file.name.bl6 <- paste0("./output/F4/Network_Analyses_GSEA_",
                            toupper(TCT), "_BL6_", geneset2subset,".pdf")
  } else {
    file.name.nzo <- paste0("./output/F4/Network_Analyses_HPEA_",
                            toupper(TCT), "_NZO_", geneset2subset,".pdf")
    file.name.bl6 <- paste0("./output/F4/Network_Analyses_HPEA_",
                            toupper(TCT), "_BL6_", geneset2subset,".pdf")
  }

  if(nrow(enrich.res.subset.nzo)>0){

    df.nzo <- prepare_subset_df(enrich.res.subset.nzo, GSEA)
    g.nzo <- df2graph(df.nzo, enrich.res.subset.nzo, GSEA = T)
    ggsave(plot = g.nzo,
           filename = file.name.nzo,
           width = 10, height = 6, useDingbats = F)
  }

  if(nrow(enrich.res.subset.bl6)>0){
    df.bl6 <- prepare_subset_df(enrich.res.subset.bl6, GSEA)
    g.bl6 <- df2graph(graph.df = df.bl6, graph.meta = enrich.res.subset.bl6,
                      GSEA= T)
    ggsave(plot = g.bl6,
           filename = file.name.bl6,
           width = 10, height = 6, useDingbats = F)
  }

}

create_graphs(enrich.res, geneset2subset = geneset2subset, TCT = TCT, GSEA = GSEA)


