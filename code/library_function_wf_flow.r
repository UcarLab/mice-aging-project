color_values <- c("PBL"="#d73027", ##d73027 # red
                  "pbl"="#d73027",
                  "red"="#d73027",
                  
                  
                  "spleen"="#4575b4", # blue
                  "spl"="#4575b4",
                  "blue"="#4575b4",
                  
                  "CD8.memory"="#1a9850", # dark green
                  "CD8 memory"="#1a9850",
                  "CD8+ memory" ="#1a9850",
                  "dgreen" ="#1a9850",
                  
                  "CD8.naive"="#99d594", # light green
                  "CD8 naive"="#99d594",
                  "CD8+ naive"="#99d594",
                  "CD8+\nnaive"="#99d594",
                  "lgreen"="#99d594",
                  
                  "B6" = "#ffa836", #"#F9A637", #"#ffa836", #orange #"#d73027",# red
                  "b6" = "#ffa836", #"#ffa836", #"#d73027", 
                  
                  
                  "NZO" = "#77B060", # "#900C3F", #"#86318c", #"#622466", #"#77B061", # dark red/purple-ish #"#4575b4", # dark blue
                  "nzo" = "#77B060", #"#900C3F", #"#86318c", #"#622466", #"#77B061", # "#4575b4",
                  
                  
                  "royalblue" = "#3878DB",
                  
                  
                  "3" = "#FFC300", # "#FCEAD6", #see through orange     #yel #FFC300  
                  "12" = "#EF843A", #"#FCBC87",  #medium orange     #orange #EF843A
                  "18" = "#900C3F", #"#F59537", #very orange      #dark red #900C3F
                  
                  "F" = "#F9C6D7", #pink
                  "female" = "#F9C6D7",
                  
                  "M" = "#67a3d9", # blue 
                  "male" = "#67a3d9",
                  
                  "age" = "#d73027", 
                  "sex" = "#4575b4", 
                  "gender" = "#4575b4",
                  "strain" = "#1a9850",
                  "Tissue" = "#91bfdb", # light blue
                  " unexplained" = "grey",
                  
                  "Age" = "#d73027", 
                  "Sex" = "#4575b4", 
                  "Gender" = "#4575b4",
                  "Strain" = "#1a9850",
                  "Tissue" = "#ffa836", 
                  "Unexplained" = "grey"
                  )

fomat_variable_names_cd_cells <- function(x){
  # x <- ggdf$variable
  if(!is.character(x)){
    x <- as.character(x)
  }
  
  # replace dots by space
  y <- sapply(x, function(x){
    new_x <- gsub("\\.", " ", x)
    return(new_x)
  })
  
  
  # add + after CD8 or CD4 or ILR7 or PD1
  z <- sapply(y, function(x){
    splitted <- strsplit(x, " ")[[1]]
    
    idx <- grep("CD4|CD8", splitted)
    
    if(length(idx)>0){
      splitted[idx] <- paste0(splitted[idx],"+")
    }
    
    if(any(grepl("IL7R|PD1", splitted))){
      splitted <- gsub("IL7R|PD1", "", splitted)
    }
    
    if(any(grepl("Naive", splitted))){
      splitted[grepl("Naive", splitted)] <- "naive"
    }
    return(paste0(splitted, collapse = " "))
  })
  
  return(z)
}



generate_cool_warm <- function(n){
  cool_warm <- colorspace::diverge_hcl(n, h = c(250, 10), c = 100, l = c(37, 88), 
                          power = c(0.7, 1.7))
  return(cool_warm)
}

coef.heatmap <- function(x, title_, path_to_File, plots_in_pdf, 
                         with_symbols = 16, main_p = NULL, main_c = NULL,
                         drop_third = FALSE){
  
  signs <- sign(x)
  
  txt2 = txt3 = txt4 = matrix(0,nr=nrow(x),nc=ncol(x))
  txt2[,1:2]=p.adjust(abs(x[,1:2]),"fdr")
  txt2[,3]=p.adjust(abs(x[,3]),"fdr")
  x <- txt2*signs
  
  MAX=max(abs(-sign(x)*log(abs(x))),na.rm=T)+1
  
  breaksList = seq(-MAX, MAX, by = 1)
  
  txt4 <- ifelse(txt2 < 0.05, formatC(txt2, format = "e", digits = 0), NA)
  
  if(!is.null(main_p)){
    # main_p <- coefficients_matrix$PBL$p
    # main_c <- coefficients_matrix$PBL$coef
    p_adj = c_after <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    
    p_adj[,1:2] <- p.adjust(main_p, "fdr")
    c_after[,1:2] <- main_c[,1:2]
    
    for(i in 1:nrow(p_adj)){
      for(j in 1:ncol(p_adj)){
        if(p_adj[i,j] < 0.05){
          c_after[i,j] <- round(c_after[i,j], digits = 2)
        }else{
          c_after[i,j] <- NA
        }
      }
    }
    if(drop_third){ # previously, there was no if(drop_third), but all c_after[,3] were set to NA
      c_after[,3] <- NA
    }
    #c_after[,3] <- NA
    txt4[,1:3] <- as.character(c_after)
  }
  
  rownames(x) <- fomat_variable_names_cd_cells(rownames(x))
  rownames(x)[nrow(x)] <- "CD4+ to CD8+ ratio"
  
  hm_flow <- -sign(x)*log(abs(x))
  
  rownames(txt3) = rownames(txt4) <- rownames(hm_flow)
  colnames(txt3) = colnames(txt4) <- colnames(hm_flow)
  
  sections <- factor(c(rep("Major cell population", 11),
                       rep("IL7R+ cell population", 11),
                       rep("PD1+ cell population", 11),
                       ""), 
                     levels = c("Major cell population",
                                "IL7R+ cell population",
                                "PD1+ cell population", 
                                ""))
  
  colvec_flow=colorRamp2(breaksList, generate_cool_warm(length(breaksList)))
  
  legend_flow <- Legend(col_fun = colvec_flow, 
                  title = "-log p * direction",
                  title_position = "topcenter",
                  at = pretty(trunc(breaksList), n = 3),
                  legend_height = unit(4, "cm"))
  
  if(drop_third == TRUE){
    hm_flow <- hm_flow[,c("B6", "NZO")]
    txt4 <- txt4[,c("B6", "NZO")]
  }
  
  if(plots_in_pdf){
    pdf(path_to_File, width = ifelse(drop_third, 5, 7), height = 12, useDingbats = F)
  }
  draw(Heatmap(hm_flow,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = sections,
          cluster_row_slices = FALSE,
          show_heatmap_legend = FALSE,
          col = colvec_flow,
          row_gap = unit(4, "mm"),
          rect_gp = gpar(col = "white", lwd = 2),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(!is.na(txt4[i, j])){
              if(with_symbols<0){
                grid.text(sprintf("%s", txt4[i, j]), x, y, 
                          gp = gpar(fontsize = 10))
              }else{
                grid.points(x, y, pch = with_symbols, size = unit(2, "mm"))
              }

            }
          }),
       padding = unit(c(3, 3, 3, 30), "mm"))
  draw(legend_flow, 
       just = c("left", "bottom"),
       x = unit(0.80, "npc"), 
       y = unit(0.825, "npc"))
  
  if(plots_in_pdf){
    dev.off()
  }
  
}

coef.heatmap_both <- function(x_both, title_ = "", path_to_File, plots_in_pdf, 
                              with_symbols = -1, drop_third = FALSE,
                              coef_both = NULL){
  
   # x_both <- store_hm_matrices
   # 
   # coef_both <- coefficients_matrix
  
  if(!is.null(coef_both)){
    x_both[["spleen"]] <- list(x_both[["spleen"]], coef_both[["spleen"]])
    x_both[["PBL"]] <- list(x_both[["PBL"]], coef_both[["PBL"]])
  }
  
  formatted_x <- lapply(x_both, function(x){
    #x <- x_both[[1]]
    if(!is.null(coef_both)){
      coefs <- x[[2]][["coef"]]
      pvals <- x[[2]][["p"]]
      x <- x[[1]]
    }
    
    signs <- sign(x)
    
    txt2 = txt3 = txt4 = matrix(0,nr=nrow(x),nc=ncol(x))
    txt2[,1:2]=p.adjust(abs(x[,1:2]),"fdr")
    txt2[,3]=p.adjust(abs(x[,3]),"fdr")
    x <- txt2*signs
    
    txt4 <- ifelse(txt2 < 0.05, formatC(txt2, format = "e", digits = 0), NA)
    
    if(!is.null(coef_both)){
      p_adj = c_after <- matrix(0, nrow = nrow(x), ncol = ncol(x))
      
      p_adj[,1:2] <- p.adjust(pvals, "fdr")
      c_after[,1:2] <- coefs[,1:2]
      
      for(i in 1:nrow(p_adj)){
        for(j in 1:ncol(p_adj)){
          if(p_adj[i,j] < 0.05){
            c_after[i,j] <- round(c_after[i,j], digits = 2)
          }else{
            c_after[i,j] <- NA
          }
        }
      }
      c_after[,3] <- NA
      txt4[,1:3] <- as.character(c_after)
    }
    
    rownames(x) <- fomat_variable_names_cd_cells(rownames(x))
    rownames(x)[nrow(x)] <- "CD4+ to CD8+ ratio"
    
    hm_flow <- -sign(x)*log(abs(x))
    
    rownames(txt3) = rownames(txt4) <- rownames(hm_flow)
    colnames(txt3) = colnames(txt4) <- colnames(hm_flow)
    
    
    if(drop_third){
      hm_flow <- hm_flow[,1:2]
      txt4 <- txt4[,1:2]
    }
    
    return(list(hm_flow, txt4))
  })
  
  names(formatted_x) <- names(x_both)
  
  cbound_x <- cbind(formatted_x[[1]][[1]], formatted_x[[2]][[1]])
  
  MAX=max(abs(cbound_x),na.rm=T)+1
  
  #max(abs(-sign(cbound_x)*log(abs(cbound_x))),na.rm=T)+1
  
  breaksList = seq(-MAX, MAX, by = 1)
  
  sections <- factor(c(rep("Major cell population", 11),
                       rep("IL7R+ cell population", 11),
                       rep("PD1+ cell population", 11),
                       ""), 
                     levels = c("Major cell population",
                                "IL7R+ cell population",
                                "PD1+ cell population", 
                                ""))
  
  
  colvec_flow=colorRamp2(breaksList, generate_cool_warm(length(breaksList)))
  
  legend_flow <- Legend(col_fun = colvec_flow, 
                        title = "-log p * direction",
                        title_position = "topcenter",
                        at = pretty(trunc(breaksList), n = 3),
                        legend_height = unit(4, "cm"))
  
  if(plots_in_pdf){
    pdf(path_to_File, width = ifelse(drop_third, 5, 7), height = 12, useDingbats = F)
  }
  draw(Heatmap(formatted_x[["spleen"]][[1]],
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_split = sections,
                cluster_row_slices = FALSE,
                show_heatmap_legend = FALSE,
                column_title = "Spleen",
                col = colvec_flow,
                row_gap = unit(4, "mm"),
                rect_gp = gpar(col = "white", lwd = 2),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(formatted_x[["spleen"]][[2]][i, j]))
                    if(with_symbols < 0){
                      grid.text(sprintf("%s", formatted_x[["spleen"]][[2]][i, j]), x, y,
                                gp = gpar(fontsize = 10))
                    }else{
                      grid.points(x, y, pch = with_symbols, size = unit(2, "mm"))
                    }

                }) +
    Heatmap(formatted_x[["PBL"]][[1]],
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_split = sections,
               cluster_row_slices = FALSE,
               show_heatmap_legend = FALSE,
               col = colvec_flow,
               column_title = "PBL",
               row_gap = unit(4, "mm"),
               rect_gp = gpar(col = "white", lwd = 2),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(formatted_x[["PBL"]][[2]][i, j]))
                   if(with_symbols < 0){
                     grid.text(sprintf("%s", formatted_x[["PBL"]][[2]][i, j]), x, y, 
                               gp = gpar(fontsize = 10))
                   }else{
                     grid.points(x, y, pch = with_symbols, size = unit(2, "mm"))
                   }
               }),
       padding = unit(c(3, 3, 3, 15), "mm"))
  draw(legend_flow, 
       just = c("left", "bottom"),
       x = unit(ifelse(drop_third,0.82, 0.88), "npc"), 
       y = unit(0.825, "npc"))
  
  if(plots_in_pdf){
    dev.off()
  }
}

library(dplyr)

coef.heatmap_both_horizontal <- function(x_both, title_ = "", path_to_File, plots_in_pdf, 
                              with_symbols = -1, drop_third = TRUE,
                              coef_both = NULL){
  
  # x_both <- store_hm_matrices
  # 
  # coef_both <- coefficients_matrix
  
  if(!is.null(coef_both)){
    x_both[["spleen"]] <- list(x_both[["spleen"]], coef_both[["spleen"]])
    x_both[["PBL"]] <- list(x_both[["PBL"]], coef_both[["PBL"]])
  }
  
  formatted_x <- lapply(x_both, function(x){
    #x <- x_both[[1]]
    if(!is.null(coef_both)){
      coefs <- x[[2]][["coef"]]
      pvals <- x[[2]][["p"]]
      x <- x[[1]]
    }
    
    signs <- sign(x)
    
    txt2 = txt3 = txt4 = matrix(0,nr=nrow(x),nc=ncol(x))
    txt2[,1:2]=p.adjust(abs(x[,1:2]),"fdr")
    txt2[,3]=p.adjust(abs(x[,3]),"fdr")
    x <- txt2*signs
    
    txt4 <- ifelse(txt2 < 0.05, formatC(txt2, format = "e", digits = 0), NA)
    
    if(!is.null(coef_both)){
      p_adj = c_after <- matrix(0, nrow = nrow(x), ncol = ncol(x))
      
      p_adj[,1:2] <- p.adjust(pvals, "fdr")
      c_after[,1:2] <- coefs[,1:2]
      
      for(i in 1:nrow(p_adj)){
        for(j in 1:ncol(p_adj)){
          if(p_adj[i,j] < 0.05){
            c_after[i,j] <- round(c_after[i,j], digits = 2)
          }else{
            c_after[i,j] <- NA
          }
        }
      }
      c_after[,3] <- NA
      txt4[,1:3] <- as.character(c_after)
    }
    
    rownames(x) <- fomat_variable_names_cd_cells(rownames(x))
    rownames(x)[nrow(x)] <- "CD4+ to CD8+ ratio"
    
    hm_flow <- -sign(x)*log(abs(x))
    
    rownames(txt3) = rownames(txt4) <- rownames(hm_flow)
    colnames(txt3) = colnames(txt4) <- colnames(hm_flow)
    
    
    if(drop_third){
      hm_flow <- hm_flow[,1:2]
      txt4 <- txt4[,1:2]
    }
    
    return(list(hm_flow, txt4))
  })
  
  names(formatted_x) <- names(x_both)
  
  cbound_x <- cbind(formatted_x[[1]][[1]], formatted_x[[2]][[1]])
  
  MAX=max(abs(cbound_x),na.rm=T)+1
  
  #max(abs(-sign(cbound_x)*log(abs(cbound_x))),na.rm=T)+1
  
  breaksList = seq(-MAX, MAX, by = 1)
  
  sections <- factor(c(rep("Major cell population", 11),
                       rep("IL7R+ cell population", 11),
                       rep("PD1+ cell population", 11),
                       ""), 
                     levels = c("Major cell population",
                                "IL7R+ cell population",
                                "PD1+ cell population", 
                                ""))
  
  
  colvec_flow=colorRamp2(breaksList, generate_cool_warm(length(breaksList)))
  
  legend_flow <- Legend(col_fun = colvec_flow, 
                        title = "-log p * direction",
                        title_position = "topcenter",
                        at = pretty(trunc(breaksList), n = 3),
                        legend_height = unit(4, "cm"))
  
  
  
  if(plots_in_pdf){
    pdf(path_to_File, width = ifelse(drop_third, 5, 7), height = 12, useDingbats = F)
  }
  draw(Heatmap(formatted_x[["spleen"]][[1]],
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               #row_split = sections,
               cluster_row_slices = FALSE,
               show_heatmap_legend = FALSE,
               column_title = "Spleen",
               col = colvec_flow,
               row_gap = unit(4, "mm"),
               rect_gp = gpar(col = "white", lwd = 2),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(formatted_x[["spleen"]][[2]][i, j]))
                   if(with_symbols < 0){
                     grid.text(sprintf("%s", formatted_x[["spleen"]][[2]][i, j]), x, y,
                               gp = gpar(fontsize = 10))
                   }else{
                     grid.points(x, y, pch = with_symbols, size = unit(2, "mm"))
                   }
                 
               }) +
         Heatmap(formatted_x[["PBL"]][[1]],
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 #row_split = sections,
                 cluster_row_slices = FALSE,
                 show_heatmap_legend = FALSE,
                 col = colvec_flow,
                 column_title = "PBL",
                 row_gap = unit(4, "mm"),
                 rect_gp = gpar(col = "white", lwd = 2),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(!is.na(formatted_x[["PBL"]][[2]][i, j]))
                     if(with_symbols < 0){
                       grid.text(sprintf("%s", formatted_x[["PBL"]][[2]][i, j]), x, y, 
                                 gp = gpar(fontsize = 10))
                     }else{
                       grid.points(x, y, pch = with_symbols, size = unit(2, "mm"))
                     }
                 }),
       padding = unit(c(3, 3, 3, 15), "mm"))
  draw(legend_flow, 
       just = c("left", "bottom"),
       x = unit(ifelse(drop_third,0.82, 0.88), "npc"), 
       y = unit(0.825, "npc"))
  
  if(plots_in_pdf){
    dev.off()
  }
}






# cool_warm <- function(n) {
#   colormap <- Rgnuplot:::GpdivergingColormap(seq(0,1,length.out=n),
#                                              rgb1 = colorspace::sRGB( 0.230, 0.299, 0.754),
#                                              rgb2 = colorspace::sRGB( 0.706, 0.016, 0.150),
#                                              outColorspace = "sRGB")
#   colormap[colormap>1] <- 1 # sometimes values are slightly larger than 1
#   colormap <- grDevices::rgb(colormap[,1], colormap[,2], colormap[,3])
#   colormap
# }



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if(numPlots==1){
    print(plots[[1]])

  }else{
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for(i in 1:numPlots){
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


PCA <- function(emat,color.var){
   emat=emat[rowSums(emat)>0,]

   pca=prcomp(t(emat),center=T,scale=T)

   d <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)

   xl <- sprintf("PC 1: %.1f %%", d[1])
   yl <- sprintf("PC 2: %.1f %%", d[2])


   dat=data.frame(PC1=as.numeric(pca$x[,1]),
                  PC2=as.numeric(pca$x[,2]),
                  Tissue=color.var[,"TISSUE"],
                  strain=color.var[,"STRAIN"],
                  age=color.var[,"AGE"],
                  gender=color.var[,"GENDER"])
 #p=vector("list",4)
   for(i in 3:ncol(dat)){
     p <- ggplot(dat, aes(PC1,PC2)) +
       geom_point(aes(color=dat[,i])) +
       labs(x=xl,y=yl,color=colnames(dat)[i])
     print(p)
   }
# #multiplot(p[[1]],p[[2]],p[[3]],p[[4]],cols=2) ## this doesn't work
}

regression <- function(dat){
  percent.vs.age=matrix(NA,nr=ncol(dat)-5,nc=4)
  for(i in 6:ncol(dat)){
    temp=summary(lm(dat[,i]~dat$Age))$coef
    if(nrow(temp)==2) percent.vs.age[i-5,]=temp[,c(1,4)]
  }

  percent.vs.age= percent.vs.age[,-ncol(percent.vs.age)+1]

  rownames(percent.vs.age)=colnames(dat)[6:ncol(dat)]
  colnames(percent.vs.age)=c("intercept","aging_coef","aging_p")

  return(percent.vs.age)
}

change.column.name <- function(x){
  ## y=x
  ## for(i in 1:length(x))
  ##   {
  ##     temp= strsplit(x[i],"[.]")[[1]]
  ##     temp[temp==
  ##     y[i]=paste(temp[-c(1:6,length(temp)-0:3)],collapse=" ")
  ##   }
  return(x)
}
pvalue.convert <- function(x){
  for(i in 1:ncol(x)){
    x[,i]=p.adjust(x[,i],"fdr")
  }
  y=x
  y[]=0
  y[x<0.05]=1
  return(y)
  
}

    
####### CD4T,8T subpopulation ggplot ##############



plotting <- function(dat0,datF,datM,datB6,datNZO){
  for(k in 6:ncol(dat0))
    { 
      rg=range(dat0[,k],na.rm=T)
      rgx=range(dat0$age,na.rm=T)

      pF <- ggplot(datF, aes(age, datF[,k] ))+geom_point(aes(color=type)) +ggtitle(paste(colnames(dat)[k],"Female"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
      pM <- ggplot(datM, aes(age, datM[,k] ))+geom_point(aes(color=type)) +ggtitle(paste(colnames(dat)[k],"Male"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

      pB6 <- ggplot(datB6, aes(age, datB6[,k] ))+geom_point(aes(color=Sex)) +ggtitle(paste(colnames(dat)[k],"B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
      pNZO <- ggplot(datNZO, aes(age, datNZO[,k] ))+geom_point(aes(color=Sex)) +ggtitle(paste(colnames(dat)[k],"NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

      multiplot(pF,pM,pB6,pNZO,cols=2)
    }
}

plotting2 <- function(dat0,datMB6,datMNZO,datFB6,datFNZO){
  for(k in 6:ncol(dat0)){
    rg=range(dat0[,k],na.rm=T)
    rgx=range(dat0$age,na.rm=T)

    pFB6 <- ggplot(datFB6, aes(age, datFB6[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"F B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
    pMNZO <- ggplot(datMNZO, aes(age, datMNZO[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"M NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

    pMB6 <- ggplot(datMB6, aes(age, datMB6[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"M B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
    pFNZO <- ggplot(datFNZO, aes(age, datFNZO[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"F NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

    
    multiplot(pFB6,pFNZO,pMB6,pMNZO,cols=2)
  }
}

plotting3 <- function(dat0,datB6,datNZO){
  for(k in 6:ncol(dat0)){
    rg=range(dat0[,k],na.rm=T)
    rgx=range(dat0$age,na.rm=T)

    pB6 <- ggplot(datB6, aes(age, datB6[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
    pNZO <- ggplot(datNZO, aes(age, datNZO[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
    
    multiplot(pB6,NA,pNZO,NA,cols=2)
  }
}

subplot <- function(ggdf){
  ggplot(ggdf, aes(x = new, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    facet_grid(Tissue~ type, scales = "free_x") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Age - sample ID") +
    ylab("percentage") #   +  scale_fill_manual(values = c("red","blue","green"))
}


# green: #1a9641
# blue: #2c7bb6
# mouse_colours <- c("B6" = "#2c7bb6",  # blue
#                    "NZO" = "#1a9641")  # green


subplot <- function(ggdf,title, mouse1 = TRUE){
  
  
  if(any(grepl("Strain", colnames(ggdf)))){
    column <- "Strain"
  }else{
    column <- "type"
  }
  
  ggdf <- ggdf %>% filter(new != 12) 
  
  
  gg1 <- ggplot(ggdf , aes(x = .data[[column]], y = value, 
                           fill = as.factor(new), group = as.factor(new))) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    geom_point(aes(shape = Sex), position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1), 
               na.rm = TRUE, color = "black") +  
    facet_grid(Tissue~ variable, scales = "free_x") + 
    stat_compare_means(aes(label =sprintf("%.3f", as.numeric(..p.format..))), 
                       label.y = 92, hide.ns = T, size = 3, method = "t.test") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(fill = "Age") +
    xlab("Strains") +
    ylab(title) +
    # geom_smooth(method = "lm", se = FALSE,size=1, na.rm = TRUE) +
    scale_y_continuous(breaks = c(0,50,100), limits = c(0,100))
    # scale_x_continuous(breaks = c(3,18)) # stat_smooth(method="lm",col="grey", se=FALSE)#   +  scale_fill_manual(values = c("red","blue","green"))

  
  if(!mouse1){
    gg1 <- gg1 +
      scale_fill_manual(values = color_values) +
      scale_color_manual(values = color_values)
      # stat_smooth(method="lm",col="grey", se=FALSE)#   
  }
  return(gg1)
}


subplot2 <- function(ggdf,marker){

  temp=as.matrix(data.frame(strsplit(as.character(ggdf$variable),marker)))
  temp=temp[nrow(temp),]
  ggdf$variable=temp
  
  ggplot(ggdf, aes(x = new, y = value, fill = variable)) +
    geom_point(aes(colour=variable)) +
    facet_grid(Tissue~ type, scales = "free_x") +
    theme_bw() + 
    scale_x_continuous(breaks = c(3,12,18)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Age - sample ID") +
    ylab(paste("percentage of",marker,"+")) #   +  scale_fill_manual(values = c("red","blue","green"))
}
  

subplot3 <- function(ggdf,marker){
  temp=as.matrix(data.frame(strsplit(as.character(ggdf$variable),marker)))
  temp=temp[nrow(temp),]
  ggdf$variable=temp

  ggplot(ggdf, aes(x = new, y = value,colour=variable)) +
    geom_jitter(width = 0.1, height = 0.1) +
    facet_grid(Tissue~ type, scales = "free_x") +
    xlab("Age") +
    ylab("") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values = c("red","blue","green","yellow"),
                       breaks=c("X..Viable","Granulocytes.Fre.Viable.Cells", 
                                "No.Granulocytes.Monocytes.Fre.Parent", 
                                "No.Granulocytes.Monocytes.Fre.Viable.Cells"),
                       labels=c("Viable.Cells.percentage",
                                "Granulocytes.\nFreq.of.Viable.Cells",
                                "No.Granulocytes.Monocytes.\nFreq.of.Parent", 
                                "No.Granulocytes.Monocytes.\nFreq.of.Viable.Cells")) #   +  scale_fill_manual(values = c("red","blue","green"))
}
  
subplot4 <- function(ggdf,marker,title){

  temp=as.matrix(data.frame(strsplit(as.character(ggdf$variable),paste(marker,".",sep=""))))
  temp=temp[nrow(temp),]
  ggdf$variable=temp
  
  ggplot(ggdf, aes(x = new, y = value, colour = type)) +
    geom_jitter(width = 0.1, height = 0.1) +
    facet_grid(Tissue~ variable, scales = "free_x") +
    theme_bw() +
    ggtitle(title) +
    geom_smooth(method = "lm", se = FALSE,size=0.3) +
    scale_x_continuous(breaks = c(3,12,18)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("age") +
    ylab(paste("percentage of",marker,"+"))#   +  scale_fill_manual(values = c("red","blue","green"))
}
