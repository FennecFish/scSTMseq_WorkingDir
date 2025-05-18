setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(SingleCellExperiment)
library(MASS)
library(MANOVA.RM)
##################################################################################
######################## direct plot simulated data ##############################
##################################################################################
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/nSample20_nCellType5_noBatch_StromalCell/1000sims"
files <- list.files(path = dir, pattern = "Null")
file_name <- files[1]
sims <- readRDS(paste0(dir, "/", file_name))

# collapse simulate cells by their group and sample
theta.collapsed <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Time, Sample, Group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Time, Sample) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  dplyr::select(-count) %>%
  pivot_wider(names_from = Group, values_from = proportion, values_fill = 0)

theta.collapsed %>% arrange(Sample)

##### Fit Manova.RM #######
response_vars <- grep("^Group", colnames(theta.collapsed), value = TRUE)
as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time"))
fit <- multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time")),
              data = theta.collapsed, subject = "Sample", within = "Time", iter = 1000)
summary(fit)
saveRDS(fit, "/proj/milovelab/wu/scLDAseq/res/1MultiSample_SingleResponse_Simulation/sims1000_Manova_pValue_Null_nSample20_nCellType5_noBatch_StromalCell.rds")

pValue <- readRDS("/proj/milovelab/wu/scLDAseq/res/1MultiSample_SingleResponse_Simulation/sims1000_Manova_pValue_Null_nSample20_nCellType5_noBatch_StromalCell.rds")
pValue <- as.numeric(do.call(rbind, pValue))
hist(pValue, main = "Distribution of MATS pValues under the Null with Simulated Data", breaks = 20)

# use qqplot
pValue <- sort(pValue)
# Generate the theoretical quantiles for a uniform(0,1) distribution
theoretical_quantiles <- qunif(ppoints(length(pValue)))

# Create the QQ-plot
qqplot(theoretical_quantiles, pValue, main = "QQ-Plot of P-values using direct simulations vs Uniform(0,1)",
       xlab = "Theoretical Quantiles (Uniform(0,1))", ylab = "Sample Quantiles (P-values)")
abline(0, 1, col = "red")  # Add a reference line

##################################################################################
################### distribution ad-hoc replicates ###############################
##################################################################################
pValue <- readRDS("/proj/milovelab/wu/scLDAseq/res/1MultiSample_SingleResponse_Simulation/1000ManovaTheta_Pooled_noContent_Prevalence_Time.rds")
pValue <- pValue[grep("NullModel",names(pValue))]
rep <- do.call(rbind, pValue)
rep_long <- pivot_longer(rep, cols = everything(), names_to = "methods", values_to = "pValue")
png("res/1MultiSample_SingleResponse_Simulation/1000sim_pValue_Histogram_allReplicates.png",
    height = 1200, width = 1800, res = 220)
ggplot(rep_long, aes(x = pValue)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), color = "black", aes(fill = methods)) +
  facet_wrap(~ methods) +
  theme_minimal() +
  labs(x = "p-value", y = "Frequency", title = "Histogram of p-values under the null across all replicates") +
  scale_fill_manual(values = c("pValue.wts" = "skyblue", "pValue.mats" = "coral"))
dev.off()

pValue.wts <- sort(rep$pValue.wts)
pValue.mats <- sort(rep$pValue.mats)
theoretical_quantiles <- qunif(ppoints(length(pValue.wts)))

png("res/1MultiSample_SingleResponse_Simulation/1000sim_pValue_qqplot_WTS_allReplicates.png",
    height = 1500, width = 1500, res = 220)
qqplot(theoretical_quantiles, pValue.wts, main = "QQ-Plot of P-values with WTS vs Uniform(0,1)",
       xlab = "Theoretical Quantiles (Uniform(0,1))", ylab = "Sample Quantiles (P-values)")
abline(0, 1, col = "red")  # Add a reference line
dev.off()

png("res/1MultiSample_SingleResponse_Simulation/1000sim_pValue_qqplot_MATS_allReplicates.png",
    height = 1500, width = 1500, res = 220)
qqplot(theoretical_quantiles, pValue.mats, main = "QQ-Plot of P-values with MATS vs Uniform(0,1)",
       xlab = "Theoretical Quantiles (Uniform(0,1))", ylab = "Sample Quantiles (P-values)")
abline(0, 1, col = "red")  # Add a reference line
dev.off()
##################################################################################
#################### distribution ad-hoc Medium #################################
##################################################################################
manova.pValue <- lapply(pValue, function(x) {
  data.frame(wts =quantile(x$pValue.wts, probs = 0.50),
             mats = quantile(x$pValue.mats, probs = 0.50))
})
manova.pValue <- do.call(rbind, manova.pValue)
manova.pValue.long <- pivot_longer(manova.pValue, cols = everything(), names_to = "methods", values_to = "pValue")

png("res/1MultiSample_SingleResponse_Simulation/1000sim_pValue_Histogram_Medium.png",
    height = 1200, width = 1800, res = 220)
ggplot(manova.pValue.long, aes(x = pValue)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), color = "black", aes(fill = methods)) +
  facet_wrap(~ methods) +
  theme_minimal() +
  labs(x = "p-value", y = "Frequency", title = "Histogram of P-values Under the Null Across Medium of Replicates") +
  scale_fill_manual(values = c("wts" = "skyblue", "mats" = "coral"))
dev.off()

pValue.wts.medium <- sort(manova.pValue$wts)
pValue.mats.medium <- sort(manova.pValue$mats)
theoretical_quantiles <- qunif(ppoints(length(pValue.wts.medium)))

png("res/1MultiSample_SingleResponse_Simulation/1000sim_pValue_qqplot_WTS_Medium.png",
    height = 1500, width = 1500, res = 220)
qqplot(theoretical_quantiles, pValue.wts.medium, main = "QQ-Plot of P-values with WTS vs Uniform(0,1)",
       xlab = "Uniform(0,1)", ylab = "P-Value of Medium Replicates")
abline(0, 1, col = "red")  # Add a reference line
dev.off()

png("res/1MultiSample_SingleResponse_Simulation/1000sim_pValue_qqplot_MATS_Medium.png",
    height = 1500, width = 1500, res = 220)
qqplot(theoretical_quantiles, pValue.mats.medium, main = "QQ-Plot of P-values with MATS vs Uniform(0,1)",
       xlab = "Uniform(0,1)", ylab = "P-Value of Medium Replicates")
abline(0, 1, col = "red")  # Add a reference line
dev.off()

check.na <- lapply(null.dat, function(x){
  y <- is.na(x[,1])
  return(y)
})
check.na <- do.call(rbind, check.na)
table(check.na)
null.dat.dist <- do.call(rbind, null.dat)
manova_long <- null.dat.dist %>%
  dplyr::select(starts_with("pValue")) %>%
  pivot_longer(cols = starts_with("pValue"),
               names_to = "Type",
               values_to = "pValue") %>%
  mutate(Type = sub("pValue\\.", "",Type))
# manova_long <- data.frame(
#   Value = c(null.dat.dist[,1], null.dat.dist[,2]),
#   Type = rep(c("WTS", "MATS"), each = nrow(null.dat.dist))
# )
png("res/1MultiSample_SingleResponse_Simulation/Histogram_pValue_AllReplicates_SingleResponse_noBatch_StromalCell_1000scSTM_Pooled_noContent_Prevalence_Time.png", width = 1500, height = 1500, res = 220)
ggplot(manova_long, aes(x = pValue, fill = Type)) +
  geom_histogram(bins = 20, color = "black", aes(fill = Type)) +
  labs(x = "p-value", y = "Frequency", title = "Histogram of p-values under the null across all replicates") +
  scale_fill_manual(values = c("wts" = "skyblue", "mats" = "coral")) + 
  facet_wrap(~ Type) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold the title
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Increase axis title size
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    legend.position = "none"  # Remove the legend as it's not necessary with the method labels on the x-axis
  )
dev.off()

manova.pValue <- lapply(null.dat, function(x) {
  data.frame(wts =quantile(x$pValue.wts, probs = 0.50),
             mats = quantile(x$pValue.mats, probs = 0.50))
})
manova.pValue <- do.call(rbind, manova.pValue)

manova_long <- manova.pValue %>%
  pivot_longer(cols = everything(),
               names_to = "Type",
               values_to = "pValue") %>%
  mutate(Type = sub("pValue\\.", "",Type))
# manova_long <- data.frame(
#   Value = c(null.dat.dist[,1], null.dat.dist[,2]),
#   Type = rep(c("WTS", "MATS"), each = nrow(null.dat.dist))
# )
png("res/1MultiSample_SingleResponse_Simulation/Histogram_pValue_50thQuantile_SingleResponse_noBatch_StromalCell_1000scSTM_Pooled_noContent_Prevalence_Time.png", width = 1500, height = 1500, res = 220)
ggplot(manova_long, aes(x = pValue, fill = Type)) +
  geom_histogram(bins = 20, color = "black", aes(fill = Type)) +
  labs(x = "p-value", y = "Frequency", title = "Histogram of p-values under the null using medium of replicates") +
  scale_fill_manual(values = c("wts" = "skyblue", "mats" = "coral")) + 
  facet_wrap(~ Type) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold the title
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Increase axis title size
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    legend.position = "none"  # Remove the legend as it's not necessary with the method labels on the x-axis
  )
dev.off()

