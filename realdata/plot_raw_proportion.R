setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyverse)

collaposed.theta <- readRDS( "res/PD1/collaposed_theta.rds")
############### proportion change by E and NE ##############
average_by_group <- lapply(collaposed.theta, function(rep){
  avg_group <- rep %>%
    group_by(expansion, timepoint) %>%
    summarise(across(starts_with("V"), mean, .names = "mean_{.col}"), .groups = "drop")
})
average_by_group <- bind_rows(average_by_group)
average_by_group_long <- average_by_group %>%
  pivot_longer(cols = starts_with("mean"), names_to = "CellType", values_to = "Proportion")

ggplot(data = average_by_group_long, aes(x = CellType, y = Proportion, fill = timepoint)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~expansion)

# iRL
library(compositions)
ilr <- compositions::ilr(average_by_group[,grep("mean", colnames(average_by_group))])
average_by_group_ilr <- cbind(average_by_group[,!grepl("mean", colnames(average_by_group))], ilr)
average_by_group_ilr_long <- average_by_group_ilr %>%
  pivot_longer(cols = starts_with("V"), names_to = "CellType", values_to = "Proportion")

ggplot(data = average_by_group_ilr_long, aes(x = CellType, y = Proportion, fill = timepoint)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~expansion)
################### proportion change by sample ####################
average_by_sample <- lapply(collaposed.theta, function(rep){
  avg_group <- rep %>%
    group_by(expansion, timepoint, patient_id) %>%
    summarise(across(starts_with("V"), mean, .names = "mean_{.col}"), .groups = "drop")
})
average_by_sample <- bind_rows(average_by_sample)
average_by_sample_long <- average_by_sample %>%
  pivot_longer(cols = starts_with("mean"), names_to = "CellType", values_to = "Proportion")

ggplot(data = average_by_sample_long, aes(x = CellType, y = Proportion, fill = timepoint)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(expansion ~patient_id)

#################### relative difference #########################
collapsed.theta <- readRDS("res/PD1/collaposed_theta.rds")
### relative difference
names(collapsed.theta) <- paste0("Replicate", 1:length(collapsed.theta))
results <- lapply(names(collapsed.theta), function(name) {
  rep <- collapsed.theta[[name]]
  relative_difference <- rep %>%
    pivot_longer(cols = starts_with("V"), names_to = "Variable", values_to = "Value") %>%
    pivot_wider(names_from = timepoint, values_from = Value) %>%
    mutate(Relative_Difference = (On - Pre) / Pre) %>%
    dplyr::select(patient_id, expansion, Variable, Relative_Difference) %>%
    pivot_wider(names_from = Variable, values_from = Relative_Difference)

  colnames(relative_difference)[3:ncol(relative_difference)] <- paste0("Topic",1:(ncol(relative_difference)-2))

  # Summarize mean and SD for differences grouped by expansion
  # relative_difference %>%
  #   pivot_longer(cols = starts_with("Topic"), names_to = "Topics", values_to = "diff") %>%
  #   group_by(expansion, Topics) %>%
  #   summarise(mean_diff = mean(diff, na.rm = TRUE),.groups = "drop")
  relative_difference %>%
    pivot_longer(cols = starts_with("Topic"), names_to = "Topics", values_to = "diff") %>%
    mutate(Repliates = name)
}) %>% bind_rows()

results_mean <- results %>% 
  group_by(expansion, Topics, Repliates) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE), .groups = "drop")

mean_min <- function(x) mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)
mean_max <- function(x) mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
png("res/PD1/RelativeChange_Proportion_By_Replicates.png",
    width = 3000, height = 2000, res = 220)
ggplot(results_mean, aes(x = Topics, y = mean_diff, color = expansion)) +
  # Error bars for mean ± SD
  # scale_y_continuous(limits = c(-1, 1)) + 
  stat_summary(
    fun = "mean",
    fun.min = mean_min,
    fun.max = mean_max,
    geom = "errorbar",
    width = 0.5,
    position = position_dodge(width = 0.5),
    size = 1  # Make error bars thicker
  ) +
  # Points for individual values (jittered)
  geom_jitter(
    aes(group = expansion), 
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), 
    size = 1, 
    alpha = 0.5
  ) +
  # Points for the mean
  stat_summary(
    fun = "mean", 
    geom = "point", 
    size = 3,  # Larger mean points
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) + # Thicker reference line
  labs(
    title = "Inferred Topic/Cell Type Proportion Relative Difference Between Time points",
    x = "Inferred Topics",
    y = "Relative Difference Between Time",
    color = "Response Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Bigger, bold title
    axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 1),  # Bigger, bold x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Bigger, bold y-axis labels
    axis.title.x = element_text(size = 16, face = "bold"),  # Bigger, bold x-axis title
    axis.title.y = element_text(size = 16, face = "bold"),  # Bigger, bold y-axis title
    legend.text = element_text(size = 14),  # Bigger, bold legend text
    legend.title = element_text(size = 16),  # Bigger, bold legend title
    legend.position = "top"
  ) +
  scale_color_manual(
    values = c("E" = "skyblue", "NE" = "orange"), 
    labels = c("E" = "Expansion", "NE" = "Non-Expansion")
  )
dev.off()
#### time-response interaction #####
results <- lapply(collapsed.theta, function(rep) {
  relative_difference <- rep %>%
    pivot_longer(cols = starts_with("V"), names_to = "Variable", values_to = "Value") %>%
    pivot_wider(names_from = timepoint, values_from = Value) %>%
    mutate(Relative_Difference = (On - Pre) ) %>%
    dplyr::select(patient_id, expansion, Variable, Relative_Difference) %>%
    pivot_wider(names_from = Variable, values_from = Relative_Difference)
  
  colnames(relative_difference)[3:ncol(relative_difference)] <- paste0("Topic",1:(ncol(relative_difference)-2))
  
  # Summarize mean and SD for differences grouped by expansion
  relative_difference %>%
    group_by(expansion) %>%
    summarise(across(starts_with("Topic"), mean)) %>%
    pivot_longer(cols = starts_with("Topic"), names_to = "Topics", values_to = "diff") %>%
    pivot_wider(names_from = expansion, values_from = diff) %>%
    group_by(Topics) %>% 
    summarise(mean_diff = E - NE) 
}) %>% bind_rows()

ggplot(results, aes(x = Topics, y = mean_diff)) +
  
  # Error bars for mean ± SD
  stat_summary(
    fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x),
    geom = "errorbar", width = 1, position = position_dodge(width = 0.5)
  ) +
  # Points for the mean
  stat_summary(
    fun = "mean", geom = "point", size = 4, position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Reference line
  labs(
    title = "Proportion Differences by Topic and Expansion Group",
    x = "Topics",
    y = "Mean Proportion Difference (± SD)",
    color = "Expansion Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "top"
  ) +
  scale_color_manual(
    values = c("E" = "skyblue", "NE" = "orange"), 
    labels = c("E" = "Expansion", "NE" = "Non-expansion")
  )

################### ILR  relative diffference ##################

library(compositions)
results <- lapply(collapsed.theta, function(rep) {

  ilr <- compositions::ilr(rep[,grep("V", colnames(rep))])
  new_rep <- cbind(rep[,!grepl("V", colnames(rep))], ilr)
  
  relative_difference <- new_rep %>%
    pivot_longer(cols = starts_with("V"), names_to = "Variable", values_to = "Value") %>%
    pivot_wider(names_from = timepoint, values_from = Value) %>%
    mutate(Relative_Difference = (On - Pre)/On) %>%
    # mutate(Relative_Difference = (On - Pre)) %>%
    dplyr::select(patient_id, expansion, Variable, Relative_Difference) %>%
    pivot_wider(names_from = Variable, values_from = Relative_Difference)
  
  colnames(relative_difference)[3:ncol(relative_difference)] <- paste0("Topic",1:(ncol(relative_difference)-2))
  # Summarize mean and SD for differences grouped by expansion
  relative_difference %>%
    pivot_longer(cols = starts_with("Topic"), names_to = "Topics", values_to = "diff") %>%
    group_by(expansion, Topics) %>%
    summarise(mean_diff = mean(diff, na.rm = TRUE),.groups = "drop")
}) %>% bind_rows()

png("res/PD1/RelativeChange_ILR_Proportion_By_Replicates.png",
    width = 3000, height = 2000, res = 220)
ggplot(results, aes(x = Topics, y = mean_diff, color = expansion)) +
  scale_y_continuous(limits = c(-5, 5)) + 
  # Error bars for mean ± SD
  stat_summary(
    fun = "mean", fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x),
    geom = "errorbar", width = 0.5, position = position_dodge(width = 0.5)
  ) +
  # Points for individual values (outliers and all)
  geom_jitter(aes(group = expansion), position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), size = 1, alpha = 0.5) +
  # Points for the mean
  # stat_summary(
  #   fun = "mean", geom = "point", size = 4, position = position_dodge(width = 0.5)
  # ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Reference line
  labs(
    title = "ILR Transformed Inferred Topics Relative Differences Between Time points",
    x = "ILR Transformation of Topics",
    y = "Relative Differences Between Time points",
    color = "Response Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Bigger, bold title
    axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 1),  # Bigger, bold x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Bigger, bold y-axis labels
    axis.title.x = element_text(size = 16, face = "bold"),  # Bigger, bold x-axis title
    axis.title.y = element_text(size = 16, face = "bold"),  # Bigger, bold y-axis title
    legend.text = element_text(size = 14),  # Bigger, bold legend text
    legend.title = element_text(size = 16),  # Bigger, bold legend title
    legend.position = "top"
  ) +
  scale_color_manual(
    values = c("E" = "skyblue", "NE" = "orange"), 
    labels = c("E" = "Response", "NE" = "Non-Response")
  )
dev.off()
