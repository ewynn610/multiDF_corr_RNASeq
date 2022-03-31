#######################################
##
## Script name: edgeR* Simulation Analysis
##
## Purpose of script: Analyze simulated data using edgeR*
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-28
##
#######################################

## load up the packages we will need:  (comment as required)
library(edgeR)

## Source in function to calculate simulation metrics
source("run_methods/calculate_sim_metrics.R")


#######################################
## 01 - Read in data/contrasts and format
#######################################

############## Read in data #################
sim_data <- readRDS("run_methods/sim_data.RDS")

############## Extract Counts/meta data #################
counts <- as.matrix(sim_data$counts)

meta_data = data.frame(
  patient = factor(sim_data$ids),
  time = factor(sim_data$time),
  group = factor(sim_data$groups)
)

de = sim_data$de

############## Make Contrasts #################

## Add in 0's to contrasts. Number of patients-1
num_pat_design = length(unique(meta_data$patient)) - 1
pat_0s = rep(0, num_pat_design)


## Contrasts have to be made to incorporate patient as fixed effect
## Insert patient 0's after intercept, time variable positions
contrasts = list(
  # c1 is inestimable since group not included in the model
  # Difference between any timepoints in treatment group
  c2 = rbind(
    c(0, 1,  0, 0, pat_0s,  1,  0,  0),
    #Difference btwn. time0 and time1 in trt.
    c(0,  0,  1, 0, pat_0s,  0,  1,  0),
    #Difference btwn. time0 and time2 in trt.
    c(0,  0,  0, 1, pat_0s,  0,  0,  1)#, #Difference btwn. time0 and time3 in trt.
  ),
  # Are any of the interactions significant
  c3 = rbind(
    c(0, 0, 0,  0,  pat_0s,  1,  0, 0),
    c(0, 0, 0,  0,  pat_0s,  0,  1, 0),
    c(0, 0, 0,  0,  pat_0s,  0,  0, 1)
  ),
  # Any Significant Coefficients
  c4 = rbind(
    c(0,  0,  0,  0, pat_0s,  0,  0,  0),
    c(0,  1,  0,  0, pat_0s, 0,  0,  0),
    c(0,  0,  1,  0, pat_0s,  0,  0,  0),
    c(0,  0,  0,  1, pat_0s, 0, 0,  0),
    c(0,  0,  0,  0, pat_0s, 1,  0, 0),
    c(0,  0,  0,  0, pat_0s, 0,  1, 0),
    c(0,  0,  0,  0, pat_0s, 0,  0, 1)
  )
)

## Transpose contrasts so they work in edgeR
contrasts = lapply(contrasts, function(x) {
  x = t(rbind(x))
})


#######################################
## 02 - Run edgeR
#######################################

############## Fit edgeR models #################

design_mat <- model.matrix( ~ time + patient + group:time, data = meta_data)
## time0 interaction needs to be removed to make estimable
design_mat <- design_mat[, -grep("time0", colnames(design_mat))]

dge <- DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
dge <- edgeR::estimateDisp(dge, design = design_mat)

fit <- glmFit(dge, design_mat)

############## Run LRT test for each contrast #################
lrt_res = lapply(contrasts, function(x) {
  tab = glmLRT(fit, contrast = x)$table
  tab$p_val_adj = p.adjust(tab$PValue, method = "BH")
  tab$de = de
  colnames(tab)[colnames(tab) == "PValue"] = "p_val_raw"
  tab
  
})

#######################################
## 03 - Calculate simulation metrics
#######################################

############## Calculate T1E for all contrast #################
T1E_vals = sapply(lrt_res, stat_calc, pval_thresh = .05, stat = "T1E")

############## Calculate FDR for all contrast #################
FDR_vals = sapply(lrt_res, stat_calc, pval_thresh = .05, stat = "FDR")

############## Calculate power for all contrast #################
power_vals = sapply(lrt_res, stat_calc, pval_thresh = .05, stat = "power")

############## Calculate non-convergence (same for all contrasts)#########
non_conv = stat_calc(lrt_res$c2, stat = "non_conv")
