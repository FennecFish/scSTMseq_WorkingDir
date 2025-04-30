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

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample20_nCellType5_noBatch_CancerCell/"
file_name <- list.files(paste0(dir, "scSTM_LinearMixed_noContent_Prevalence_Time"))[1]
set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", file_name)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

scSTM_name <- paste0(dir, "scSTM_LinearMixed_noContent_Prevalence_Time/scSTM_",
                     set_level,".rds")
if(file.exists(scSTM_name)){
  scSTMobj <- readRDS(scSTM_name)
  scSTMobj <- select_top_scSTM(scSTMobj)
}else{
  next
}
metadata <- colData(scSTMobj$settings$sce) %>% as.data.frame() %>% dplyr::select(Cell, Sample, Time, Response)
theta.old <- scSTMobj$theta
theta.old <- theta.old %>%
  as.data.frame() %>%
  rownames_to_column("Cell") %>%
  left_join(metadata, by = "Cell")
theta.collapsed <- theta.old %>%
  group_by(Sample, Time) %>%
  summarise(across(starts_with("V"), sum), .groups = 'drop') %>%
  ungroup()
theta.new <- theta.collapsed %>%
  rowwise() %>%
  mutate(across(starts_with("V"), ~ . / sum(c_across(starts_with("V"))))) %>%
  ungroup()
Y <- theta.new %>% 
  dplyr::select(starts_with("V"))
Y <- compositions::ilr(Y)
x <- theta.new %>%
  dplyr::select(-starts_with("V")) %>%
  dplyr::mutate(Time = as.factor(Time),
                Sample = as.factor(Sample))
x <- cbind(x, Y)
response_vars <- grep("^V", names(x), value = TRUE)

fit.temp <- multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time")),
                   data = x, subject = "Sample", within = "Time", iter = 1000)

n_patients <- 6
patients <- factor(rep(1:n_patients, each = 2))
time <- rep(c(0, 1), n_patients)
# Random intercept for each patient
random_intercept <- rnorm(n_patients, mean = 70, sd = 10) 
# assume weight increases by 2 unit, at time = 2
Group <- rep(1:3, each = 4)
weight <- random_intercept[as.numeric(patients)] + 2 * time*Group
TrtEffect <-  random_intercept[as.numeric(patients)] +  5* Group + rnorm(n_patients * 2, sd = 1)
data <- data.frame(Patient = patients, Time = time, Weight = weight, Group = Group, TrtEffect = TrtEffect)

# formula = as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time"))
formula = as.formula("cbind(Weight,TrtEffect) ~ Group + Time")
data = data
subject = "Patient"
# data = x
# subject = "Sample"
within = "Time"
iter = 1000
alpha = 0.05
resampling = "paramBS"
dec = 3
para = FALSE
seed = 123
############################### multRM #############################
input_list <- list(formula = formula, data = data,
                   iter = iter, alpha = alpha, resampling = resampling, seed = seed)
#--------------------------------------------------------------------------------#
# prepare data for the calculations
prepdt <- prepare.data(formula, data, subject, within)
# extract relevant info
p <- prepdt[["p"]]
EF <- prepdt[["EF"]]
nf <- prepdt[["nf"]]
split3 <- prepdt[["split3"]]
no.whole <- prepdt[["no.whole"]]
Yw2 <- prepdt[["data"]]
hypo_matrices <- prepdt[["hypo"]]
fac_names <- prepdt[["fac"]]
nind <- prepdt[["n"]]
fl <- prepdt[["fl"]]
lev_names <- prepdt[["lev_names"]]
no.subf <- prepdt[["no.subf"]]
#-------------------------------------------------------------------------#
WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
WTPS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 2)
MATS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 1)
quantiles <- matrix(NA, ncol = 2, nrow = length(hypo_matrices))
rownames(WTS_out) <- fac_names
rownames(WTPS_out) <- fac_names
rownames(MATS_out) <- fac_names
rownames(quantiles) <- fac_names
colnames(MATS_out) <- "Test statistic"
colnames(quantiles) <- c("WTS_resampling", "MATS_resampling")
# calculate results
for (i in 1:length(hypo_matrices)) {
  results <- multRM.Statistic(Yw2, nind, hypo_matrices[[i]], iter, alpha, resampling, 
                              para, CPU,
                              seed, p, t=prod(fl[within]))
  WTS_out[i, ] <- round(results$WTS, dec)
  WTPS_out[i, ] <- round(results$WTPS, dec)
  MATS_out[i] <- round(results$MATS, dec)
  quantiles[i, ] <- results$quantiles
}
mean_out <- matrix(round(results$Mean, dec), ncol = p, byrow = TRUE)
Var_out <- results$Cov
descriptive <- cbind(unique(lev_names), rep(nind, each =prod(fl[within])) , mean_out)
colnames(descriptive) <- c(EF, "n", split3)
rownames(descriptive) <- NULL
colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(MATS)"))
#WTPS_out[WTPS_out == 0] <- "<0.001"
colnames(MATS_out) <- "Test statistic"

# Output ------------------------------------------------------
output <- list()
output$input <- input_list
output$Descriptive <- descriptive
output$Covariance <- Var_out
output$Means <- mean_out
output$MATS <- MATS_out
output$WTS <- WTS_out
output$resampling <- WTPS_out
output$quantile <- quantiles
output$nf <- nf
output$H <- hypo_matrices
output$factors <- fac_names
output$p <- p
output$fl <- fl
output$BSMeans <- results$BSmeans
output$BSVar <- results$BSVar
output$levels <- lev_names
#output$nested <- nest
output$other <- list(no.subf = no.subf, no.whole = no.whole, p = p, within = within)
output$nested <- FALSE
output$modelcall <- multRM
output$modeltype <- "multRM"

# check for singular covariance matrix
test <- try(solve(output$Covariance), silent = TRUE)
if(!is.matrix(test)){
  warning("The covariance matrix is singular. The WTS provides no valid test statistic!")
}