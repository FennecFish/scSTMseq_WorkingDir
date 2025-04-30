setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(stringr)
library(ggplot2)
library(dplyr)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_f_nc_null/", pattern = "scSTM*")
set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  files)
res <- data.frame()
for (level in set_level){
  alt <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test4/", 
                        "scSTM_", level, ".rds"))
  if(class(alt) == "selectModel"){
    all_values <- unlist(alt$bound)
    max_value <- max(all_values)
    max_position_in_vector <- which(all_values == max_value)
    alt <- alt$runout[[max_position_in_vector]]
  }
  null <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_f_nc_null/", 
                         "scSTM_", level, ".rds"))
  alt.l <- alt$convergence$bound[length(alt$convergence$bound)]
  null.l <- null$convergence$bound[length(null$convergence$bound)]
  test.stat <- -2*(null.l - alt.l)
  K <- alt$settings$dim$K
  p.value <- pchisq(test.stat, df= K - 1, lower.tail=FALSE)
  temp <- c(unlist(str_split(level, "_")), alt.l, null.l, K, round(test.stat,2), p.value)
  res <- rbind(res, temp)
}
colnames(res) <- c("seed", "control", "level", "nCellType", "alt.logLik", "null.logLik", "df", "statistics", "p.value")
write.csv(res, file = "res/LRT_single_V3.csv")

res <- read.csv("res/LRT_single_V3.csv")
power <- res %>%
  filter(control == "pos") %>%  # Subset where the null hypothesis is false
  summarise(
    total_cases = n(),  # Total cases where truth = 1
    successful_rejections = sum(p.value < 0.05),  # Cases where the null was correctly rejected
    power = successful_rejections / total_cases  # Proportion of successful rejections
  )

typeI <- res %>%
  filter(control == "neg") %>%  # Subset where the null hypothesis is false
  summarise(
    total_cases = n(),
    false_rejections = sum(p.value < 0.05),  # Cases where the null was correctly rejected
    alpha = false_rejections / total_cases  # Proportion of successful rejections
  )

# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_f_nc/", pattern = "scSTM*")
# set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  files)
# res <- data.frame()
# for (level in set_level){
#   res <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_f_nc/", 
#                         "scSTM_", level, ".rds"))
#   n_samples <- 1000  # Number of samples
#   samples <- MASS::mvrnorm(n_samples, mu = res$mu$gamma[2,], Sigma = diag(as.vector(res$mu$sn)))
#   # Calculate the credible intervals
#   alpha <- 0.05  # Significance level for 95% credible interval
#   credible_intervals <- apply(samples, 2, function(x) quantile(x, probs = c(alpha/2, 1 - alpha/2)))
#   
#   # Display the credible intervals
#   credible_intervals <- t(credible_intervals)  # Transpose for better readability
#   colnames(credible_intervals) <- c("Lower Bound", "Upper Bound")
#   credible_intervals
#   temp <- c(unlist(str_split(level, "_")), alt.l, null.l, K, round(test.stat,2), p.value)
#   res <- rbind(res, temp)
# }
# colnames(res) <- c("seed", "control", "level", "nCellType", "alt.logLik", "null.logLik", "df", "statistics", "p.value")
# write.csv(res, file = "res/LRT_single_V3.csv")



