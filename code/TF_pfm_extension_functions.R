
split_jaspar_name_return_df_of_names <- function(concatenated_jaspar_ids){
  #concatenated_jaspar <- df_both_raw$`Transcription factor`
  standard_name <- lapply(concatenated_jaspar_ids, function(x){
    #x <- concatenated_jaspar[1]
    id_first <- substr(x, 1, 6)
    id_second <- substr(x, 7, 7)
    id_together <- paste0(id_first, ".",id_second)
    name_of <- toupper(substr(x, 8, nchar(x)))
    return(list(x, id_together, name_of))
  })
  
  df_of_names <- data.frame(matrix(unlist(standard_name), ncol = 3, 
                                   byrow = TRUE),
                            stringsAsFactors = F)
  colnames(df_of_names) <- c("Transcription factor", "Jaspar_matrix_ID", "name")
  return(df_of_names)
}


add_colnames_for_type <- function(df, addition, front = TRUE){
  if(front){
    colnames(df) <- paste0(addition, colnames(df))  
  }else{
    colnames(df) <- paste0(colnames(df), addition)
  }
  return(df)
}

library(JASPAR2018)

tf_name_to_jasp_id <- function(x){
  #x <- jasp_ID[1]
  first_eight <- strsplit(x, "")[[1]][1:8]
  if(!is.character(first_eight[8])){
    stop(paste("WARNING! JASPAR MATRIX ID IS POSSIBLY LONGER THAN 7:", x))
  }
  
  full_jasp_id <- paste0(paste0(first_eight[1:6], collapse = ""), ".", 
                         paste0(first_eight[7:7], collapse = ""))
  
  return(full_jasp_id)
}


db <- file.path(system.file("extdata", package="JASPAR2018"),
                "JASPAR2018.sqlite")

convert_ID_to_PM <- function(searchable_jasp_id, return_only_pfm = TRUE){
  PFM_matrix <- getMatrixByID(db, searchable_jasp_id)
  PWM_matrix <- toPWM(PFM_matrix, type = "prob")
  if(return_only_pfm){
    return(pfm = PFM_matrix)
  }else{
    return(list(pwm = PWM_matrix, pfm = PFM_matrix))  
  }
}

compute_similarity_return_clusters <- function(pm_matrices, 
                                               use_pfm = TRUE, 
                                               num_clusters){
  if(use_pfm){
    similarity_score_pfm <- matrix(0, nrow = length(pm_matrices), ncol = length(pm_matrices))
    
    for(i in 1:length(pm_matrices)){
      for(j in 1:length(pm_matrices)){
        #similarity_score_pwm[i,j] <- PWMSimilarity(PWM_matrices_from_IDs[[i]][[1]], PWM_matrices_from_IDs[[j]][[1]], method = "KL")
        similarity_score_pfm[i,j] <- PFMSimilarity(pm_matrices[[i]], pm_matrices[[j]])[2]
      }
    }
    
    colnames(similarity_score_pfm) <- rownames(similarity_score_pfm) <- names(pm_matrices)
    
    # if not rounded, the similarity matrix is not symmetrical
    similarity_score_pfm <- round(similarity_score_pfm, digits = 2)
    
    distance_object <- dist(similarity_score_pfm)
    
    hclust_object <- hclust(distance_object)
    
    cluster_ids <- stats::cutree(hclust_object, k = num_clusters)
    
    return(list(grouped_rows = cluster_ids, 
                as_cluster = hclust_object,
                distance_obj = distance_object,
                sim_score = similarity_score_pfm))
  }
  
  return(NULL)
  
}

concatenate_group_vector_names <- function(group_named_vector){
  distinct_groups <- length(unique(group_named_vector))
  names_per_group <- vector("list", length = distinct_groups)
  for(i in 1:distinct_groups){
    names_per_group[[i]] <- group_named_vector[which(group_named_vector == i)]
  }
  
  clean_unique_names_vector <- sapply(names_per_group, function(x){
    only_name <- sapply(names(x), function(y){
      return(substr(y, 8, nchar(y)))
    })
    
    unique_vector <- paste(sort(unique(toupper(only_name))), collapse = ",")
    return(unique_vector)
  })
  
  named_vector <- c(1:distinct_groups)
  names(named_vector) <- clean_unique_names_vector
  return(named_vector)
}

merge_heatmap_by_cluster_id <- function(some_heatmap, row_groups){
  
  #some_heatmap <- heatmaps_full_names
  #row_groups <- new_groups
  
  hm_with_group <- merge(as.data.frame(some_heatmap), 
                         as.data.frame(row_groups), 
                         by = 0)
  rownames(hm_with_group) <- hm_with_group[,"Row.names"]
  
  hm_with_group <- hm_with_group[,-1]
  
  # dplyr: groups the matrix by group ID and computes mean for each column
  averaged <- hm_with_group %>% group_by(row_groups) %>% summarise_all(mean)
  
  merged_row_names <- data.frame(row_groups = concatenate_group_vector_names(row_groups))
  merged_row_names$`Transcription factor` <- rownames(merged_row_names)
  
  heatmap_with_merged <- merge(merged_row_names, averaged, by = "row_groups", 
                               all = TRUE)
  
  #clean_heatmap <- subset(heatmap_with_merged, select = -c(row_groups))
  
  rownames(heatmap_with_merged) <- heatmap_with_merged[,"Transcription factor"]
  
  clean_heatmap <- subset(heatmap_with_merged, select = -c(`Transcription factor`, row_groups))
  
  matrix_hm <- as.matrix(clean_heatmap)
  
  return(matrix_hm)
}

# MAKE SURE THE COLNAMES OF THE HEATMAP YOU ENTER ARE UNIQUE!

cluster_heatmaps_on_TF_similarity <- function(heatmaps_full_names, 
                                              number_of_clusters = 7,
                                              return_dend = FALSE){
  # function takes in heatmap, returns merged rows per TF similarity clustering
  
  #heatmaps_full_names <- starting_heatmap
  
  rows <- rownames(heatmaps_full_names)
  # convert rownames to Jaspar Ids that can be searched in the database
  rows_japs_id <- sapply(rows, tf_name_to_jasp_id)
  
  # search the Jaspar database and return pfm matrices 
  pfm_matrices <- lapply(rows_japs_id, convert_ID_to_PM)
  
  # compute rel. % PFM similarity, do hierarchical clustering, return group IDs
  new_groups <- compute_similarity_return_clusters(pfm_matrices, 
                                                   num_clusters = number_of_clusters)
  
  # average heatmaps rows based on cluster membership
  merged_heatmap <- merge_heatmap_by_cluster_id(heatmaps_full_names, new_groups[["grouped_rows"]])
  if(!return_dend){
    return(merged_heatmap)
  }else{
    return(list(merged_hm = merged_heatmap, 
                groups = new_groups[["grouped_rows"]], 
                tree = new_groups[["as_cluster"]],
                distances = new_groups[[3]],
                similarities = new_groups[[4]]))
  }
  
}

introduce_new_line_to_rownames <- function(x, new_line_after = 40, splitSymbol = ","){
  elements <- strsplit(x, splitSymbol)[[1]]
  number_of_elements <- length(elements)
  with_new_line <- elements[1]
  if(number_of_elements > 1){
    current_line <- with_new_line
    # if more than n chars in current line --> introduce new line, reset counter
    for(i in 2:number_of_elements){
      if(nchar(current_line) > new_line_after){ #(i %% new_line_after) == 1
        if(splitSymbol != " "){
          with_new_line <- paste0(with_new_line, paste0(splitSymbol,"\n"), elements[i])
        }else{
          with_new_line <- paste0(with_new_line, "\n", elements[i])
        }
        
        current_line <- elements[i]
      }else{
        if(splitSymbol != " "){
          with_new_line <- paste0(with_new_line, paste0(splitSymbol," "), elements[i])
          current_line <- paste0(current_line, paste0(splitSymbol, " "), elements[i])
        }else{
          with_new_line <- paste0(with_new_line, " ", elements[i])
          current_line <- paste0(current_line, " ", elements[i])
        }

      }
    }
    
  }
  return(with_new_line)
}

library(grid)
library(gridExtra)
library(viridisLite)
library(viridis)
library(dendextend)

produce_dendrogram_plot <- function(some_dend,
                                    num_clusters, 
                                    path, 
                                    title_ = "Clustering of transcription factors\nbased on their PFM similarity", 
                                    use_colours = FALSE,
                                    subtitle_ = "",
                                    dark_background = FALSE,
                                    aspect_ratio = 0.6,
                                    y_pos = 0.96,
                                    box_border = 8,
                                    very_tall = 1,
                                    sub_offset = 0.035,
                                    png = FALSE){
  gg_pfm <- as.dendrogram(some_dend)
  
  num_c <- num_clusters
  
  colors_used <- inferno(num_c)
  
  labels(gg_pfm) <- sapply(labels(gg_pfm), function(x){
    splitted <- strsplit(x, "")[[1]]
    with_dot <- c(splitted[1:6], ".", splitted[7], " ", splitted[8:length(splitted)])
    as_string <- paste0(with_dot, collapse = "")
    return(toupper(as_string))
  })
  
  
  if(produce_plots_in_pdf){
    if(png){
      png(path, width = 10*(aspect_ratio), height = 10*very_tall, units = "in", res = 1200)
    }else{
      pdf(path, width = 10*(aspect_ratio), height = 10*very_tall)
    }
    
  }
  
  par(mar=c(1,1,3,12))
  
  if(is.character(dark_background)){
    par(bg = dark_background)
  }
  
  if(use_colours){
    gg_pfm %>% 
      set("labels_col", value = colors_used, k=num_c) %>%  
      plot(horiz = TRUE,
           axes = F)
  }else{
    gg_pfm %>% 
      plot(horiz = TRUE,
           axes = F)
  }
  
  gg_pfm %>%
    rect.dendrogram(k = num_c, 
                    horiz = TRUE, 
                    border = box_border, 
                    lty = 1, 
                    lwd = 2)
  
  grid.text(title_, gp = gpar(fontsize = 12, fontface = "bold"),
            x = 0.5, y = y_pos)
  
  grid.text(subtitle_, gp = gpar(fontsize = 10, fontface = "bold"),
            x = 0.5, y = (y_pos-sub_offset))
  
  if(produce_plots_in_pdf){
    dev.off()
  }
  
  return(1)
}

only_relevant_cols <- function(df, grep_term){
  #df <- df_of_FC
  #grep_term <- " O"
  if(grep_term == FALSE){
    return(df)
  }else{
    df_mod <- df[,c(TRUE, grepl(grep_term, colnames(df)[2:ncol(df)]))]
    return(df_mod)  
  }
}


keep_only_significant_FC_pv <- function(df_of_FC, df_of_pv, filter_by, threshold){
  # df_of_FC <- opening_subjects[["df_FC"]] #opening_FC
  # df_of_pv <- opening_subjects[["df_pv"]] #opening_pv
  # filter_by <- " M "
  # threshold <- 2
  
  # ASSUMING THE P-VALUES HAVE BEEN FDR ADJUSTED
  
  # only use columns relevant to the analysis, e.g. only old
  filtered_FC_cols <- only_relevant_cols(df_of_FC, filter_by)
  filtered_pv_cols <- only_relevant_cols(df_of_pv, filter_by)
  
  # only keep the rows in which there's at least one significant result
  if(ncol(filtered_pv_cols) > 2){
    only_significant_pv <- filtered_pv_cols[(rowSums(filtered_pv_cols[,2:ncol(filtered_pv_cols)] > threshold) > 0),]
    # subset the FC matrix to obtain the same rows in FC matrix
    only_significant_FC <- filtered_FC_cols %>% filter(`Transcription factor` %in% only_significant_pv$`Transcription factor`)
  }else{
    only_significant_pv <- filtered_pv_cols[(filtered_pv_cols[,2:ncol(filtered_pv_cols)] > threshold),]
    
    only_significant_FC <- filtered_FC_cols %>% filter(`Transcription factor` %in% only_significant_pv$`Transcription factor`)
  }
  
  return(list(df_FC = only_significant_FC, df_pv = only_significant_pv))
}


convert_to_df <- function(hm_matrix){
  hm_df <- cbind(`Transcription factor` = rownames(hm_matrix), 
                 as.data.frame(hm_matrix))
  return(hm_df)
}

convert_to_matrix <- function(hm_df, id_col_idx = 1){
  all_cols_but_id <- c(1:ncol(hm_df))
  all_cols_but_id <- all_cols_but_id[which(all_cols_but_id != id_col_idx)]
  new_matrix <- as.matrix(hm_df[,all_cols_but_id])
  rownames(new_matrix) <- hm_df[,id_col_idx]
  return(new_matrix)
}

testy_change <- function(x){
  return(1)
}

#cat("Hi, I'm new")

find_significant_per_gender <- function(df_of_fc, df_of_p_vals, gender, boundary = thr){
  # df_of_fc <- opening_subjects[["df_FC"]]
  # df_of_p_vals <- opening_subjects[["df_pv"]]
  # gender <- " M "
  # 
  
  filtered_cols <- keep_only_significant_FC_pv(df_of_fc, df_of_p_vals, gender,
                                               threshold = boundary)
  
  return(filtered_cols[["df_pv"]][,1])
}

convert_ln_to_log10 <- function(vector_of_nums){
  #pv <- df_pvals_raw[,2]
  as_p <- exp(-vector_of_nums)
  as_log10 <- -log10(as_p)
  return(as_log10)
}






