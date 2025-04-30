setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mclust)

res <- readRDS("res/PD1/Manova_Res.rds")

rm.mats <- lapply(res, function(x){
  # x <- x$MANOVA.RM.res
  dat <- x$resampling[,2]
  return(dat)
})

dat.clean <- function(given_list){
  given_list <- do.call(rbind, given_list)
  rownames(given_list) <- paste0("replicate", 1:nrow(given_list))
  # colnames(given_list) <- "pValue"
  return(given_list)
}
dat <- dat.clean(rm.mats) %>% as.data.frame()

################################################################################
############################### p-value plot ###################################
################################################################################
summary_stats <- data.frame(
  median = apply(dat, 2, median),
  sd = apply(dat, 2, sd)) %>%  
  rownames_to_column("Variable") 
# summary_stats <- dat %>%
#   as.data.frame() %>%
#   dplyr::summarise(across(everything(), list(mean = mean, sd = sd))) %>%
#   pivot_longer(cols = everything(), 
#                names_to = c("Variable", "Statistic"), 
#                names_sep = "_")

# Separate the mean and SD into different columns
# summary_stats <- summary_stats %>%
#   pivot_wider(names_from = Statistic, values_from = value)
summary_stats$Variable <- c("Group", "Time", "Group:Time")
summary_stats$Variable <- factor(summary_stats$Variable, levels = c("Group", "Time", "Group:Time"))

dat_long <- dat %>% pivot_longer(cols = everything(), names_to = "Variable", values_to = "p_value")
dat_long$Variable <- factor(dat_long$Variable, levels = c("Response","Time","Response:Time"))

# Create the plot
png("res/PD1/Manova_pValue_Distribution.png", 
    width = 1500, height = 1000, res = 220)
ggplot(dat_long, aes(x = Variable, y = p_value)) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +  # Reference lines
  geom_boxplot(aes(group = Variable), fill = "lightgray", outlier.shape = NA, alpha = 0.5) +
  # Mark the median with a blue point
  stat_summary(fun = median, geom = "point", size = 3, color = "blue") +
  labs(
    title = "Distribution of pValues Across Replicates",
    x = "Covariates",
    y = "p-value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14)  # Increase facet label text size
  )
dev.off()

# previous plot
ggplot(summary_stats, aes(x = Variable, y = median, group = 1)) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  geom_point(color = "blue", size = 3) +       # Points for the mean values
  geom_errorbar(aes(ymin = median - sd, ymax = median + sd), width = 0.2, color = "black") +
  labs(title = "Distribution of pValues Across Replicates",
       x = "covariates",
       y = "p-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )


ggplot(summary_stats, aes(x = Variable, y = median, group = 1)) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  geom_point(color = "blue", size = 3) +       # Points for the mean values
  geom_errorbar(aes(ymin = median - sd, ymax = median + sd), width = 0.2, color = "black") +
  labs(title = "Distribution of pValues Across Replicates",
       x = "covariates",
       y = "p-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )
############################################################
############### proportion estimation ######################
############################################################
collapsed.theta <- readRDS("res/PD1/collaposed_theta.rds")

results <- lapply(collapsed.theta, function(df) {
  # Compute differences between On and Pre for each patient and expansion
  diff <- df %>%
    pivot_wider(names_from = timepoint, values_from = starts_with("V")) %>%
    mutate(across(matches("_On$"), ~ . - get(sub("_On$", "_Pre", cur_column())), 
                  .names = "diff_{.col}")) %>%
    dplyr::select(patient_id, expansion, starts_with("diff_V")) 
  
  colnames(diff)[3:ncol(diff)] <- paste0("Topic",1:(ncol(diff)-2))
  # Summarize mean and SD for differences grouped by expansion
  diff %>%
    pivot_longer(cols = starts_with("Topic"), names_to = "Topics", values_to = "diff") %>%
    group_by(expansion, Topics) %>%
    summarise(mean_diff = mean(diff, na.rm = TRUE),.groups = "drop")
}) %>% bind_rows()

res_plot <- results %>% 
  group_by(expansion, Topics) %>%
  summarise(mean = mean(mean_diff), sd = sd(mean_diff), .groups = "drop")

png("res/PD1/real_proportion_change_from_replicates.png",
    width = 3000, height = 2000, res = 220)
ggplot(res_plot, aes(x = Topics, y = mean, color = expansion)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5) +
  labs(
    title = "Inferred Topic/Cell Type Proportion Differences Between Timepoints By Group",
    x = "Inferred Topics",
    y = "Proportion Difference Across Time",
    color = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14),  # Increase facet label text size
    legend.text = element_text(size = 14),
    legend.position = "top") +
  scale_color_manual(
    values = c("E" = "skyblue", "NE" = "orange"), 
    labels = c("E" = "Expansion", "NE" = "Non-expansion"))
dev.off()

###############################################################################
############################# propeller #######################################
###############################################################################
prop.asin <- readRDS("res/PD1/propeller_BlockRep_asin.rds")
prop.logit <- readRDS("res/PD1/propeller_BlockRep_logit.rds")


prop.logit.long <- prop.logit %>%
  # Select only coefficient and FDR columns
  dplyr::select(starts_with("coefficients"), starts_with("fdr")) %>%
  # Reshape to long format for both coefficients and FDR
  pivot_longer(cols = starts_with("coefficients"),
               names_to = "Group", 
               values_to = "Coefficient") %>%
  mutate(FDR = case_when(
    Group == "coefficients.Response_Change" ~ fdr.Response_Change,
    Group == "coefficients.nonResponse_Change" ~ fdr.nonResponse_Change
  )) %>%
  # Add significance column
  mutate(Significance = ifelse(FDR < 0.05, "***", "")) %>%
  dplyr::select(- starts_with("fdr."))

# Create the plot
png("res/PD1/propeller.estimate.png",
    width = 1800, height = 1300, res = 220)
ggplot(prop.logit.long, aes(x = Coefficient, y = Cluster, color = Group)) +
  geom_point(size = 3) +
  geom_text(aes(label = Significance), vjust = -1.5, size = 5, color = "black") +
  scale_color_manual(values = c("coefficients.Response_Change" = "blue",
                                "coefficients.nonResponse_Change" = "red")) +
  labs(title = "Cluster Coefficients with Significance",
       x = "Coefficient",
       y = "Cluster",
       color = "Group") +
  theme_minimal() +
  theme(legend.position = "right")
dev.off()
