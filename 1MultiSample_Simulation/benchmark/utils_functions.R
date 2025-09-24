scSTMseq_plot_data <- function(thresholds, pValueList){
  exist.change.list <- lapply(thresholds, function(thr) {
    res <- pValueList %>%
      rownames_to_column("Model") %>%
      mutate(nSample = as.numeric(sub(".*nSample(\\d+)_.*", "\\1", Model)),
             nCellType = as.numeric(sub(".*nCellType(\\d+)_.*", "\\1", Model)),
             EffectSize = ifelse(grepl("Null", Model), 0, sub(".*HighVar([0-9.]+)$", "\\1", Model)),
             # EffectSize = as.numeric(sub(".*EffectSize(-?[0-9.]+)$", "\\1", Model)),
             Batch = unlist(strsplit(Model, "_"))[3],
             CancerCell = unlist(strsplit(Model, "_"))[4],
             Threshold = thr,
             # wts.sig = wts < thr,
             mats.sig = mats < thr)
    return(res)
  })
  
  ProportionOfTrue <- do.call(rbind, lapply(seq_along(exist.change.list), function(i){
    thr <- exist.change.list[[i]]$Threshold[1]
    prop_truth <- exist.change.list[[i]] %>%
      as.data.frame() %>%
      group_by(nSample, EffectSize, nCellType, Batch, CancerCell, Threshold) %>%
      summarize(
        # wts_sig_proportion = mean(wts.sig),
        mats_sig_proportion = mean(mats.sig),
        .groups = "drop")
  }))

  ProportionOfTrue <- ProportionOfTrue %>%
    mutate(Method = "scSTMseq.MATS") %>%
    pivot_wider(
      names_from = EffectSize,
      values_from = mats_sig_proportion,
      names_prefix = "EffectSize_"
    ) %>%
    dplyr::filter(nCellType != 15) %>%
    dplyr::rename(FDR = EffectSize_0) %>%
    pivot_longer(cols = starts_with("EffectSize_"),
                 names_to = "EffectSize",
                 values_to = "Power") %>%
    mutate(EffectSize = as.numeric(sub("EffectSize_", "", EffectSize)))
  return(ProportionOfTrue)
}

propeller_plot_dat <- function(dat, thresholds){
  thr.list <- do.call(rbind, lapply(thresholds, function(thr) {
    dat.fdr <- do.call(rbind, lapply(names(dat), function(name){
      x <- dat[[name]]
      asin.sig.change <- sapply(thr, function(t) x$fdr.asin < t)
      asin.sig.change <- colSums(asin.sig.change) != 0
      
      logit.sig.change <- sapply(thr, function(t) x$fdr.logit < t)
      logit.sig.change <- colSums(logit.sig.change) != 0
      
      temp <-  data.frame(
        nCellType = as.numeric(sub(".*nCellType(\\d+)_.*", "\\1", name)),
        nSample = as.numeric(sub(".*nSample(\\d+)_.*", "\\1", name)),
        Batch = unlist(strsplit(name, "_"))[3],
        CancerCell = unlist(strsplit(name, "_"))[4],
        EffectSize = ifelse(grepl("Null", name), 0, sub(".*HighVar([0-9.]+)$", "\\1", name)),
        Threshold = thr,
        propeller.asin = asin.sig.change,
        propeller.logit = logit.sig.change
        )
      return(temp)
    }))
  }))
  
  thr.list <- thr.list %>%
    pivot_longer(cols = starts_with("propeller"),
                 names_to = "Method",
                 values_to = "propeller") %>%
  group_by(nSample, EffectSize, nCellType, Batch, CancerCell, Threshold, Method) %>%
    summarize(propeller.res = mean(propeller),
              .groups = "drop") %>%
    pivot_wider(
      names_from = EffectSize,
      values_from = propeller.res,
      names_prefix = "EffectSize_") %>%
    dplyr::rename(FDR = EffectSize_0) %>%
    pivot_longer(cols = starts_with("EffectSize_"),
                 names_to = "EffectSize",
                 values_to = "Power") %>%
    mutate(EffectSize = as.numeric(sub("EffectSize_", "", EffectSize)))
  return(thr.list)
}

scSTM_clean <- function(scSTMseq, threshold){
  scSTMseq.pValue <- lapply(scSTMseq, function(x) {
    data.frame(mats = quantile(x$pValue.mats, probs = 0.50),
               wts = quantile(x$pValue.wts, probs = 0.50))
  })
  scSTMseq.pValue <- do.call(rbind, scSTMseq.pValue)
  rownames(scSTMseq.pValue) <- sub("\\.rds$", "", rownames(scSTMseq.pValue))
  scSTMseq.pValue$EffectSize <- as.numeric(sub(".*HighVar([0-9.]+)$", "\\1", rownames(scSTMseq.pValue)))
  scSTMseq.pValue$nCellType <- as.numeric(n_cell_types <- sub(".*nCellType(\\d+)_.*", "\\1", rownames(scSTMseq.pValue)))
  scSTMseq_plot <- plot_data_generate(threshold, pValueList = scSTMseq.pValue, by_group = "EffectSize")
  scSTMseq_plot <- scSTMseq_plot %>%
    dplyr::rename(WTS = wts_sig_proportion,
                  MATS = mats_sig_proportion) %>%
    pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "Power") %>%
    mutate(Methods = paste0("scSTMseq.", Methods)) %>%
    dplyr::select(EffectSize, nCellType, Threshold, Power, Methods)
  return(scSTMseq_plot)
}