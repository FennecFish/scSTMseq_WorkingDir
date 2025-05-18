# read in any RDS with catch
safe_readRDS <- function(file_path) {
  tryCatch({
    # Attempt to read the RDS file
    data <- readRDS(file_path)
    return(data)
  }, error = function(e) {
    # Handle the error
    message(paste("Error reading RDS file:", file_path))
    message("Skipping to the next file.")
    return(NULL)  # Return NULL if an error occurs
  })
}

# select the top scSTM withhighest bound, and output its clusters
process_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    if(length(scSTMobj$bound) == 1){
      scSTMobj <- scSTMobj$runout
    }else{
      all_values <- unlist(scSTMobj$bound)
      max_value <- max(all_values, na.rm = T)
      if(length(which(all_values == max_value)) > 1){
        max_position_in_vector <- which(all_values == max_value)[1]
      } else{
        max_position_in_vector <- which(all_values == max_value)
      }
      scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
    }
  }
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}


select_top_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    if(length(scSTMobj$bound) == 1){
      scSTMobj <- scSTMobj$runout
    }else{
      all_values <- unlist(scSTMobj$bound)
      max_value <- max(all_values, na.rm = T)
      if(length(which(all_values == max_value)) > 1){
        max_position_in_vector <- which(all_values == max_value)[1]
      } else{
        max_position_in_vector <- which(all_values == max_value)
      }
      scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
    }
    }
    
  return(scSTMobj)
}

cluster_scSTM <- function(scSTMobj) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}


# The following code is to map topic back to original cell types
OneToOne_Mapping_Topics <- function(scSTMobj){
  sims <- scSTMobj$settings$sce
  max_topic_per_cell <- cluster_scSTM(scSTMobj) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Cell") %>%
    dplyr::rename(Topic = ".") 
  matched_data <- colData(sims) %>%
    as.data.frame %>%
    dplyr::right_join(max_topic_per_cell, by = "Cell") 
  topic_group_counts <- matched_data %>%
    group_by(Group, Topic) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Step 2: Calculate the total number of cells in each group
  group_totals <- matched_data %>%
    group_by(Group) %>%
    summarise(total = n(), .groups = 'drop')
  
  # Step 3: Calculate the likelihood of each topic belong to the given group
  likelihoods <- topic_group_counts %>%
    dplyr::left_join(group_totals, by = "Group") %>%
    dplyr::mutate(likelihood = count / total) %>%
    dplyr::select(Group, Topic, likelihood) %>%
    group_by(Group) %>%
    slice_max(likelihood, with_ties = FALSE) %>%
    ungroup()
  return(likelihoods)
}

plot_data_generate <- function(thresholds, pValueList, by_group = NULL){
  exist.change.list <- lapply(thresholds, function(thr) {
    res <- pValueList %>%
      mutate(wts.sig = wts < thr,
             mats.sig = mats < thr)
    return(res)
  })
  ProportionOfTrue <- data.frame()
  
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    if(is.null(by_group)){
      prop_truth <- exist.change.list[[i]] %>%
        as.data.frame() %>%
        summarize(
          wts_sig_proportion = mean(wts.sig),
          mats_sig_proportion = mean(mats.sig)
        ) %>%
        mutate(Threshold = threshold)
    }else{
      # if we want to group the power/ type I error by effect size
      prop_truth <- exist.change.list[[i]] %>%
        as.data.frame() %>%
        group_by(EffectSize) %>%
        summarize(
          wts_sig_proportion = mean(wts.sig),
          mats_sig_proportion = mean(mats.sig)
        ) %>%
        mutate(Threshold = threshold)
    }
    # Store the results in a dataframe
    ProportionOfTrue <- bind_rows(ProportionOfTrue, prop_truth)
  }
  return(ProportionOfTrue)
}