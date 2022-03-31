#######################################
##
## Script name: NBMM-PL Simulation Analysis
##
## Purpose of script: Analyze simulated data using NBMM-PL
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-28
##
#######################################

## load up the packages we will need:  (comment as required)

library(DESeq2)
library(corrRNASeq)

## Source in function to calculate simulation metrics
source("run_methods/calculate_sim_metrics.R")


#######################################
## 01 - Read in data/contrasts and format
#######################################

############## Read in data/contrasts #################
sim_data <- readRDS("run_methods/sim_data.RDS")
contrasts <- readRDS("run_methods/contrasts.RDS")


############## Extract Counts/meta data #################
counts <- as.matrix(sim_data$counts)

meta_data = data.frame(
  patient = factor(sim_data$ids),
  time = factor(sim_data$time),
  group = factor(sim_data$groups)
)

de = sim_data$de

################################################
## 02 - Calculate offset (DESeq size factors)
################################################

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta_data,
                              design = ~ group * time)
dds <- DESeq(dds)
dds = estimateSizeFactors(dds)
meta_data$log_offset <- log(sizeFactors(dds))
meta_data$id = meta_data$patient

#######################################
## 03 - Fit models
#######################################

fit = corrSeq_fit(
  ~ group * time + offset(log_offset) + (1 | patient),
  expr_mat = counts,
  sample_data = meta_data,
  method = "nbmm_pl"
  ## Uncomment line to use parallelization
  #, parallel = T, cores = 8
)

#######################################
## 04 - Run Tests
#######################################


test_res = lapply(contrasts, function(contrast) {
  sum = corrSeq_summary(fit,
                        contrast = contrast,
                        df = "Satterthwaite",
                        sort_results = F)
  sum$summary_table$de = de
  sum
})


#######################################
## 05 - Calculate simulation metrics
#######################################

############## Calculate T1E for all contrast #################
T1E_vals = sapply(test_res, function(x) {
  tab = x$summary_table
  stat_calc(tab, pval_thresh = .05, stat = "T1E")
})

############## Calculate FDR for all contrast #################
FDR_vals = sapply(test_res, function(x) {
  tab = x$summary_table
  stat_calc(tab, pval_thresh = .05, stat = "FDR")
})

############## Calculate power for all contrast #################
power_vals = sapply(test_res, function(x) {
  tab = x$summary_table
  stat_calc(tab, pval_thresh = .05, stat = "power")
})

############## Calculate non-convergence (same for all contrasts)#########
non_conv = stat_calc(test_res$c1$summary_table, stat = "non_conv")
