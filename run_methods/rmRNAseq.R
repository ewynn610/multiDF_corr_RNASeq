#######################################
##
## Script name: rmRNAseq Simulation Analysis
##
## Purpose of script: Analyze simulated data using rmRNAseq
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-28
##
#######################################

## load up the packages we will need:  (comment as required)

library(rmRNAseq)

## Source in function to calculate simulation metrics
source("run_methods/calculate_sim_metrics.R")


##########################################
## 01 - Read in data/contrasts and format
##########################################

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

############## Format Contrasts for rmRNASeq #################
contrasts = lapply(contrasts, function(x) {
  x = t(rbind(x))
})

#######################################
## 02 - Fit models, run tests
#######################################

############## Fit models, run tests #################
dmat = model.matrix( ~ group * time, data = meta_data)
fit = TC_CAR1(
  counts = counts,
  design = dmat,
  Subject = as.factor(meta_data$patient),
  Time = as.numeric(meta_data$time),
  C.matrix = contrasts,
  Nboot = 100,
  ## Uncomment line to use parallelization
  #ncores = 8,
  print.progress = T
)


############## Make tables of results #################
test_res = lapply(paste0("c", 1:4), function(x) {
  p_val_adj <- p.adjust(fit$pqvalue$pv[, x], method = "BH")
  data.frame(
    p_val_raw = fit$pqvalue$pv[, x],
    p_val_adj = p_val_adj,
    de = de
  )
})
names(test_res) <- paste0("c", 1:4)


#######################################
## 03 - Calculate simulation metrics
#######################################

############## Calculate T1E for all contrast #################
T1E_vals = sapply(test_res, stat_calc, pval_thresh = .05, stat = "T1E")

############## Calculate FDR for all contrast #################
FDR_vals = sapply(test_res, stat_calc, pval_thresh = .05, stat = "FDR")

############## Calculate power for all contrast #################
power_vals = sapply(test_res, stat_calc, pval_thresh = .05, stat = "power")

############## Calculate non-convergence (same for all contrasts)#########
non_conv = stat_calc(test_res$c1, stat = "non_conv")
