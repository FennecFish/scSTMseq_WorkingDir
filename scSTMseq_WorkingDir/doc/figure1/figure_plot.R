setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(slam)
library(SingleCellExperiment)
library(Rcpp)
library(tidyverse)
library(mclust)
library(tibble)
library(stats)

res <- read.csv("res/res_fig1_adj_sims_V3.csv")

#change format from wide to long
res_long <- res %>% 
  gather(method, adjR, scSTM_combat_f_nc_cluster:sc3_cluster, factor_key=TRUE) %>%
  drop_na()

# Calculate the mean for each method per level
res_stats <- res_long %>%
  group_by(level, method) %>%
  summarize(mean_adjR = mean(adjR, na.rm = TRUE),
            sd = sd(adjR, na.rm = TRUE),  
            .groups = 'drop') 

# first demonstrate different version of scSTM
res_stats <- res_stats %>% #filter(grepl("^scSTM", method)) %>%
  mutate(line_type = ifelse(grepl("stm_cluster", method), "solid", "dashed")) %>%
  mutate(control = sapply(strsplit(level, "_"), `[`, 1),
         level = sapply(strsplit(level, "_"), `[`, 2))

# res_stats$control <- factor(res_stats$control)
png("res/figure1_scSTM_adjR.png", width = 2500, height = 1500, res = 300)
ggplot(res_stats, aes(x = level, y = mean_adjR, color = method, group = method, linetype = line_type)) +
  geom_line(linewidth = 1) +
  labs(title = "Mean Adjusted Rand Index by scSTM Methods and Level of Noise",
       x = "Level of Noise",
       y = "Mean adjRandIndex") +
  facet_grid(~control) +
 #  geom_errorbar(aes(ymin = mean_adjR - sd, ymax = mean_adjR + sd), width = 0.2) +
  theme_minimal() +
  # scale_color_brewer(palette = "Set5") +
  # scale_color_viridis_d() + 
  scale_color_brewer(palette = "Set2") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
  theme(axis.text.x = element_text(hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis text
        axis.title.x = element_text(size = 12, face = "bold"),  # Increase size of x-axis title
        axis.title.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis title
        plot.title = element_text(size = 12, face = "bold"),  # Increase size and bold the plot title
        strip.text = element_text(size = 12, face = "bold"))
dev.off()

# compare all methods, with top scSTM confg
all_res_stats <- res_stats %>% filter(!grepl("^scSTM_a", method)) %>%
  filter(!grepl("^scSTM_combat_f_nc", method)) %>%
  filter(!grepl("^scSTM_f_nc", method)) %>%
  mutate(line_type = ifelse(grepl("^scSTM", method), "dashed", "solid")) 
  
png("res/figure1_mean_adj.png", width = 2500, height = 1500, res = 300)
ggplot(all_res_stats, aes(x = level, y = mean_adjR, color = method, group = method, linetype = line_type)) +
  geom_line(linewidth = 0.8) +
  labs(title = "Mean Adjusted Rand Index by Method and Level of Noise",
       x = "Level of Noise",
       y = "Mean adjRandIndex") +
  theme_minimal() +
  # geom_errorbar(aes(ymin = mean_adjR - sd, ymax = mean_adjR + sd), width = 0.2) +
  # scale_color_brewer(palette = "Set5") +
  # scale_color_viridis_d() + 
  scale_color_brewer(palette = "Set2") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
  theme(axis.text.x = element_text(hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis text
        axis.title.x = element_text(size = 12, face = "bold"),  # Increase size of x-axis title
        axis.title.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis title
        plot.title = element_text(size = 12, face = "bold"),  # Increase size and bold the plot title
        strip.text = element_text(size = 12, face = "bold"))
dev.off()



filtered_data <- res_long %>%
  filter(level %in% paste0("L", 4:9) & 
           (method == "monocle3_cluster" | 
              method == "scSTM_combat_f_nc_cluster" | 
              method == "scSTM_f_c_cluster"))

png("res/figure1_boxplot.png", width = 3500, height = 3000, res = 300)
ggplot(filtered_data, aes(x=method, y=adjR)) +
  geom_boxplot() +
  ggtitle("Boxplot For Adjusted Rand Index For High Performing Methods") +
  facet_wrap(~level, ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
                axis.text.y = element_text(size = 15, face = "bold"),  # Increase size of y-axis text
        axis.title.x = element_text(size = 18, face = "bold"),  # Increase size of x-axis title
        axis.title.y = element_text(size = 18, face = "bold"),  # Increase size of y-axis title
        plot.title = element_text(size = 18, face = "bold"),  # Increase size and bold the plot title
        strip.text = element_text(size = 18, face = "bold"))
dev.off()
# 
# # methods <- unique(res_stats$method)
# res_stats <- res_stats %>% filter(method == "sc3_cluster" | 
#                                     method == "monocle3_cluster" | 
#                                     grepl("^scSTM", method))
# # plot
# png("res/figure1_mean_adj.png", width = 2500, height = 1500, res = 300)
# ggplot(res_stats, aes(x = level, y = mean_adjR, color = method, group = method, linetype = line_type)) +
#   geom_line(size = 1) +
#   labs(title = "Mean Adjusted Rand Index by Method and Level of Noise",
#        x = "Level of Noise",
#        y = "Mean adjRandIndex") +
#   theme_minimal() +
#   # scale_color_brewer(palette = "Set5") +
#   # scale_color_viridis_d() + 
#   scale_color_brewer(palette = "Set2") +
#   scale_linetype_manual(values = c("dashed", "solid", "dotted"), guide = "none")
# dev.off()
# 
# # boxplot
# sub_res_long <- res_long %>% 
#   filter(method != "scSTM_a_c_cluster") %>% 
#   filter(level %in% paste0("L",4:9))
# 
png("res/figure1_boxplot.png", width = 3500, height = 3000, res = 300)
ggplot(sub_res_long, aes(x=method, y=adjR)) +
  geom_boxplot() +
  ggtitle("Comparison of Clustering Accuracy among") +
  facet_wrap(~level, ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# 
# ## check to see the most DE genes
# scSTM_files <- list.files("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_f_c/")
# obj <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/scSTM_f_c/", scSTM_files[1]))
# obj$beta$logbeta
