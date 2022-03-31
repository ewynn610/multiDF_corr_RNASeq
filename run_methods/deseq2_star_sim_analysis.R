#######################################
##
## Script name: DESeq2* Simulation Analysis
##
## Purpose of script: Analyze simulated data using DESeq2*
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-28
##
#######################################

## load up the packages we will need:  (comment as required)

library(DESeq2)

## Source in function to calculate simulation metrics
source("run_methods/calculate_sim_metrics.R")


#######################################
## 01 - Read in data and format
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


#######################################
## 02 - Fit DESeq2 Full model
#######################################

design_mat <-
  model.matrix( ~ time + group:time + patient, data = meta_data)

## Remove time0 term from design matrix so results will be estimable
design_mat <- design_mat[, -grep("time0", colnames(design_mat))]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta_data,
                              design = design_mat)
dds <- DESeq(dds)

#######################################
## 03 - Run Tests
#######################################

############## Make Reduced Model Matrices #################
reduced_design_mats = list()

## c1 is inestimable since group not included as fixed effect

## Difference between any timepoints in treatment group
design_mat_c2 = model.matrix( ~ group + group:time + patient, data = meta_data)
## Remove group 1 terms
reduced_design_mats$c2 = design_mat_c2[, -grep("group1", colnames(design_mat_c2))]

## Are any of the interactions significant
reduced_design_mats$c3 = model.matrix( ~ time + patient, data = meta_data)

## Any Significant Coefficients
reduced_design_mats$c4 = model.matrix( ~ patient, data = meta_data)

############## Run tests with reduced matrices #################
lrt_res = lapply(reduced_design_mats, function(x) {
  df = data.frame(results(
    DESeq(dds, test = 'LRT', reduced = x),
    independentFiltering = F,
    cooksCutoff = FALSE
  ))
  colnames(df)[colnames(df) == "pvalue"] = "p_val_raw"
  colnames(df)[colnames(df) == "padj"] = "p_val_adj"
  df$de = de
  df
})

#######################################
## 04 - Calculate simulation metrics
#######################################

############## Calculate T1E for all contrast #################
T1E_vals = sapply(lrt_res, stat_calc, pval_thresh = .05, stat = "T1E")

############## Calculate FDR for all contrast #################
FDR_vals = sapply(lrt_res, stat_calc, pval_thresh = .05, stat = "FDR")

############## Calculate power for all contrast #################
power_vals = sapply(lrt_res, stat_calc, pval_thresh = .05, stat = "power")

############## Calculate non-convergence (same for all contrasts)#########
non_conv = stat_calc(lrt_res$c2, stat = "non_conv")
