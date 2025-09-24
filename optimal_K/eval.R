# setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(Matrix)
library(SingleCellExperiment)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(tibble)
library(stats)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/optimalK", pattern = "sims*")
sim_name <- sub("sims_(.*)\\.rds", "\\1", files)

dat <- data.frame()
# files[!file.exists(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/optimalK/combat_f_nc/optimalk_combat_f_nc_", sim_name, ".rds"))]
for (simname in sim_name){
  scSTM <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/optimalK/combat_f_nc/optimalk_combat_f_nc_", simname, ".rds")
  if(file.exists(scSTM)){
    cat(simname, " exists \n")
    res <- readRDS(scSTM)
    res <- res$results
    res$true_k <- strsplit(simname, "_")[[1]][1]
  } else(
    res = NA
  )
  dat <- rbind(dat,res)
  rm(scSTM)
}

write.csv(dat, file = "res/selectK_res.csv", row.names = F)

res <- read.csv(file = "res/selectK_res.csv")

res <- res %>% drop_na() 
mean_perplexity <- res %>%
  group_by(K, true_k) %>%
  summarise(mean_perplexity = mean(perplexity, na.rm = TRUE))

ggplot(mean_perplexity, aes(x = K, y = mean_perplexity)) +
  geom_line() +  
  geom_point() +  
  facet_wrap(~ true_k) +  
  geom_vline(xintercept = true_k, linetype = "dashed", color = "red") + 
  labs(title = "Perplexity by K with Facets for True K", x = "K", y = "Perplexity") +
  theme_minimal()

# run rpc
mean_rpc <- res %>%
  group_by(K, true_k) %>%
  summarise(mean_rpc = mean(rpc, na.rm = TRUE))
vline_data <- mean_rpc %>%
  distinct(true_k) 
true_label <- as_labeller(function(x) paste("True K =", x))

png("res/selectK_RPC.png", width = 2500, height = 1500, res = 300)
ggplot(mean_rpc, aes(x = K, y = mean_rpc)) +
  geom_line() +  
  geom_point() +  
  geom_vline(data = vline_data, aes(xintercept = true_k), linetype = "dashed", color = "red") +
  facet_wrap(~ true_k, scales = "free_x", labeller = true_label) +  
  scale_x_continuous(breaks = seq(3, 50, by = 5), limits = c(3, 50)) +
  labs(title = "Rate of Perplexity Change with increasing K", x = "K", y = "Rate of Perplexity Change") +
  theme_minimal()
dev.off()

#### check for residuals ####

mean_dat <- res %>%
  group_by(K, true_k) %>%
  summarise(mean_stats = mean(residual, na.rm = TRUE))
vline_data <- mean_dat %>%
  distinct(true_k) 
true_label <- as_labeller(function(x) paste("True K =", x))

png("res/selectK_RPC.png", width = 2500, height = 1500, res = 300)
ggplot(mean_dat, aes(x = K, y = mean_stats)) +
  geom_line() +  
  geom_point() +  
  geom_vline(data = vline_data, aes(xintercept = true_k), linetype = "dashed", color = "red") +
  facet_wrap(~ true_k, scales = "free_x", labeller = true_label) +  
  scale_x_continuous(breaks = seq(3, 50, by = 5), limits = c(3, 50)) +
  labs(title = "Residual with increasing K", x = "K", y = "Residual") +
  theme_minimal()
dev.off()

#### check for bound ####

mean_dat <- res %>%
  group_by(K, true_k) %>%
  summarise(mean_stats = mean(bound, na.rm = TRUE))
vline_data <- mean_dat %>%
  distinct(true_k) 
true_label <- as_labeller(function(x) paste("True K =", x))

png("res/selectK_RPC.png", width = 2500, height = 1500, res = 300)
ggplot(mean_dat, aes(x = K, y = mean_stats)) +
  geom_line() +  
  geom_point() +  
  geom_vline(data = vline_data, aes(xintercept = true_k), linetype = "dashed", color = "red") +
  facet_wrap(~ true_k, scales = "free_x", labeller = true_label) +  
  scale_x_continuous(breaks = seq(3, 50, by = 5), limits = c(3, 50)) +
  labs(title = "Bound with increasing K", x = "K", y = "Bound") +
  theme_minimal()
dev.off()



# test 
dat <- res %>% filter(true_k == 12)
dat <- dat[1:9,]

ggplot(dat, aes(x = K, y = rpc)) +
  geom_line() +  
  geom_point() +  
  geom_vline(data = dat, aes(xintercept = true_k), linetype = "dashed", color = "red") +
  # facet_wrap(~ true_k, scales = "free_x", labeller = true_label) +  
  scale_x_continuous(breaks = seq(3, 50, by = 5), limits = c(3, 50)) +
  labs(title = "Rate of Perplexity Change with increasing K", x = "K", y = "RPC") +
  theme_minimal()

ggplot(dat, aes(x = K, y = bound)) +
  geom_line() +  
  geom_point() +  
  geom_vline(data = dat, aes(xintercept = true_k), linetype = "dashed", color = "red") +
  # facet_wrap(~ true_k, scales = "free_x", labeller = true_label) +  
  scale_x_continuous(breaks = seq(3, 50, by = 5), limits = c(3, 50)) +
  labs(title = "Log Likelihood with increasing K", x = "K", y = "bound") +
  theme_minimal()

ggplot(dat, aes(x = K, y = residual)) +
  geom_line() +  
  geom_point() +  
  geom_vline(data = dat, aes(xintercept = true_k), linetype = "dashed", color = "red") +
  # facet_wrap(~ true_k, scales = "free_x", labeller = true_label) +  
  scale_x_continuous(breaks = seq(3, 50, by = 5), limits = c(3, 50)) +
  labs(title = "Residual with increasing K", x = "K", y = "Residual") +
  theme_minimal()
