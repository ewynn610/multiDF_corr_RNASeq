#######################################
##
## Script name: limma Simulation Analysis
##
## Purpose of script: Analyze simulated data using limma
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-28
##
#######################################

## load up the packages we will need:  (comment as required)

library(limma)

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

############## Format Contrasts for limma #################
contrasts = lapply(contrasts, function(x) {
  x = t(rbind(x))
})


#######################################
## 02 - Fit limma model with dupCorr
#######################################

design = model.matrix( ~ group * time, meta_data)

############## Run duplicateCorrelation #################
vobj_tmp = voom(counts, design)
dupcor = duplicateCorrelation(vobj_tmp, design, block = meta_data$patient)

vobj = voom(
  counts,
  design,
  plot = FALSE,
  block = meta_data$patient,
  correlation = dupcor$consensus
)

dupcor <-
  duplicateCorrelation(vobj, design, block = meta_data$patient)

############## Fit Model #################
fitDupCor <- lmFit(vobj,
                   design,
                   block = meta_data$patient,
                   correlation = dupcor$consensus)


#######################################
## 03 - Run Tests
#######################################

test_res <- lapply(contrasts, function(x) {
  cont <- contrasts.fit(fit = fitDupCor, contrasts = x)
  cont <- eBayes(cont)
  tab <- topTable(cont,
                  p.value = 1,
                  number = Inf,
                  sort.by = 'none')
  colnames(tab)[colnames(tab) == "P.Value"] = "p_val_raw"
  colnames(tab)[colnames(tab) == "adj.P.Val"] = "p_val_adj"
  tab$de = de
  tab
})

#######################################
## 04 - Calculate simulation metrics
#######################################

############## Calculate T1E for all contrast #################
T1E_vals = sapply(test_res, stat_calc, pval_thresh = .05, stat = "T1E")

############## Calculate FDR for all contrast #################
FDR_vals = sapply(test_res, stat_calc, pval_thresh = .05, stat = "FDR")

############## Calculate power for all contrast #################
power_vals = sapply(test_res, stat_calc, pval_thresh = .05, stat = "power")

############## Calculate non-convergence (same for all contrasts)#########
non_conv = stat_calc(test_res$c1, stat = "non_conv")
