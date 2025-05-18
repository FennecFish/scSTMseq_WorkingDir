setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
require(pals)
library(MASS)
library(tibble)
library(MANOVA.RM)

##### analysis ####'
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse"
# manova.dir <- "1000ManovaTheta_Pooled_noContent_Prevalence_Time"
manova.dir <- "ManovaTheta_Pooled_Content_Sample_Prevalence_Time"
paths <- Sys.glob(file.path(dir, "*nSample40*", manova.dir)) 
files <- unlist(lapply(paths, function(path) list.files(path, pattern = "Manova_pValue", full.names = TRUE)))

res <- vector(mode = "list")

for(file_name in files){
  set_level <- sub("^Manova_pValue_(.*)\\.rds$", "\\1", basename(file_name))
  design <- sub(".*/MultiResponse/([^/]+)/.*", "\\1", file_name)
  nCellType = unlist(strsplit(design, "_"))[2]
  ThetaManova <- readRDS(file_name)
  # rm.wts <- lapply(ThetaManova$manova.fit, function(x){
  #   dat <- x$resampling[1]
  #   return(dat)
  # })
  # 
  # rm.mats <- lapply(ThetaManova$manova.fit, function(x){
  #   dat <- x$resampling[2]
  #   return(dat)
  # })
  # 
  # dat.clean <- function(given_list){
  #   given_list <- do.call(rbind, given_list)
  #   rownames(given_list) <- paste0("replicate", 1:nrow(given_list))
  #   colnames(given_list) <- "pValue"
  #   return(given_list)
  # }
  ThetaManova <- ThetaManova %>% as.data.frame() %>% dplyr::filter(Type == "Interaction")
  # temp <- data.frame(pValue.Response = ThetaManova$Response,
  #                    pValue.Time = ThetaManova$Time,
  #                    pValue.Interaction = ThetaManova$`Response:Time`)
  # colnames(temp) <- c("pValue.wts", "pValue.mats")
  list.name <- paste0(nCellType, "_", set_level)
  res[[list.name]] <- ThetaManova
  cat(file_name, "\n")
}
name <- sub("^ManovaTheta_", "",basename(manova.dir))
saveRDS(res, file = paste0("res/1MultiSample_MultiResponse_Simulation/Manova_pValue_nSample20_", name,".rds"))

#################################################################################
################################## plot #########################################
#################################################################################
plot_data_generate <- function(thresholds, pValueList, by_group = NULL){
  exist.change.list <- lapply(thresholds, function(thr) {
    res <- pValueList %>%
      mutate(response.sig = Response < thr,
             time.sig = Time < thr,
             inter.sig = Response.Time < thr)
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
        group_by(EffectSize, nCellType) %>%
        summarize(
          response_sig_proportion = mean(response.sig),
          time_sig_proportion = mean(time.sig),
          interaction_sig_proportion = mean(inter.sig),
          .groups = "drop"
        ) %>%
        mutate(Threshold = threshold)
    }
    # Store the results in a dataframe
    ProportionOfTrue <- bind_rows(ProportionOfTrue, prop_truth)
  }
  return(ProportionOfTrue)
}

dat <- readRDS("/proj/milovelab/wu/scLDAseq/res/1MultiSample_MultiResponse_Simulation/Manova_pValue_nSample20_Pooled_Content_Sample_Prevalence_Time.rds")
pValue50 <- lapply(names(dat), function(name) {
  data <- dat[[name]]
  temp <- data.frame(
    nCellType = unlist(strsplit(name, "_"))[1],
    EffectSize = unlist(strsplit(name, "_"))[3],
    Response = quantile(data$Response, probs = 0.50),
    Time = quantile(data$Time, probs = 0.50),
    `Response:Time` = quantile(data$`Response:Time`, probs = 0.50)
  )
  return(temp)
})
pValue50 <- do.call(rbind, pValue50)
pValue50 <- pValue50 %>%
  mutate(
    nCellType = as.numeric(sub("nCellType", "", nCellType)),
    EffectSize = case_when(
      EffectSize == "NullModel" ~ 0,
      grepl("HighVar", EffectSize) ~ as.numeric(sub("HighVar", "", EffectSize))
    ))

threshold <- c(0.01, 0.05, 0.1)
topic_plot <- plot_data_generate(threshold, pValueList = pValue50, by_group = "EffectSize")
topic_plot <- topic_plot %>%
  mutate(Threshold = as.factor(Threshold)) %>%
  dplyr::rename("Time" = "time_sig_proportion",
                "Response" = "response_sig_proportion",
                "Time.Response.Inter" = "interaction_sig_proportion") %>%
  pivot_longer(cols = c("Time", "Response", "Time.Response.Inter"), names_to = "Covariate", values_to = "Power")

ggplot(topic_plot, aes(x = Covariate, y = Power, color = Threshold, group = Threshold)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize, 
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power/TypeI Error for Multi-Response Simulation",
       subtitle = "scSTMseq + Manova.RM",
       x = "Covariates",
       y = "Power/Type I Error") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )

################## MANOVA + Simulated Data ###########################
dat <- readRDS("/proj/milovelab/wu/scLDAseq/res/1MultiSample_MultiResponse_Simulation/Manova_sims_nSample20_Pooled_Content_Sample_Prevalence_Time.rds")

dat <- do.call(rbind, dat)
dat <- dat %>%
  as.data.frame() %>%
  mutate(
    nCellType = as.numeric(sub("nCellType", "", sapply(strsplit(rownames(dat), "_"), `[`, 1))),
    EffectSize = sapply(strsplit(rownames(dat), "_"), `[`, 3))
dat <- dat %>% 
  dplyr::rename("Response.Time" = "Response:Time") %>%
  mutate(EffectSize = case_when(
    EffectSize == "NullModel" ~ 0,
    grepl("HighVar", EffectSize) ~ as.numeric(sub("HighVar", "", EffectSize))
  ))

threshold <- c(0.01, 0.05, 0.1)
topic_plot <- plot_data_generate(threshold, pValueList = dat, by_group = "EffectSize")
topic_plot <- topic_plot %>%
  mutate(Threshold = as.factor(Threshold)) %>%
  dplyr::rename(Time = time_sig_proportion,
                Response = response_sig_proportion,
                Time.Response.Inter = interaction_sig_proportion) %>%
  pivot_longer(cols = c("Time", "Response", "Time.Response.Inter"), names_to = "Covariate", values_to = "Power")

ggplot(topic_plot, aes(x = Covariate, y = Power, color = Threshold, group = Threshold)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize, 
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "TypeIError/Power for Multi-Response Simulation",
       subtitle = "Manova with simulated data",
       x = "Covariates",
       y = "TypeIError/Power") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )
