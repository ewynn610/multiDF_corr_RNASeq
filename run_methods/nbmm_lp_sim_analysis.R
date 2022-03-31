#######################################
##
## Script name: NBMM-LP Simulation Analysis
##
## Purpose of script: Analyze simulated data using NBMM-LP
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

## Save reversed group, time variables
## Use reversed vars because need trt group as reference for between time test
meta_data$group_rev <- rev(meta_data$group)
meta_data$time_rev<-rev(meta_data$time)

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

fit=corrSeq_fit(~group_rev*time_rev+offset(log_offset)+(1|patient),
                expr_mat=counts, 
                sample_data = meta_data, method="nbmm_lp"
                ## Uncomment line to use parallelization
                #, parallel = T, cores = 8
                )

#######################################
## 04 - Run Tests
#######################################

## For each test:
##    1. Fit reduced model
##    2. Run summary function
##    3. Remove reduced model object to save space
## I need this for contrast 4

test_res=list()

############## Difference Between groups at any time #################
fit_c1=corrSeq_fit(~time_rev+offset(log_offset)+(1|patient), expr_mat=counts,
                   sample_data = meta_data,method="nbmm_lp"
                   ## Uncomment line to use parallelization
                   #, parallel = T, cores = 8
                   )
test_res$c1=corrSeq_summary(fit, fit_c1, df=NA, sort_results = F)
rm(fit_c1)

########### Difference between any timepoints in treatment group ##############
mod_mat=data.frame(model.matrix(~group_rev*time_rev+log_offset, data=meta_data), 
                   patient=meta_data$patient)
fit_c2=corrSeq_fit(~group_rev1+group_rev1.time_rev1+group_rev1.time_rev2+
                     group_rev1.time_rev3+offset(log_offset)+(1|patient),
                   expr_mat=counts, sample_data = mod_mat, 
                   method="nbmm_lp"
                   ## Uncomment line to use parallelization
                   #, parallel = T, cores = 8
)
test_res$c2=corrSeq_summary(fit, fit_c2, df=NA, sort_results = F)
rm(fit_c2)

############## Difference Between groups at any time #################
fit_c3=corrSeq_fit(~group_rev+time_rev+offset(log_offset)+(1|patient), 
                   expr_mat=counts, sample_data = meta_data,
                   method="nbmm_lp"
                   ## Uncomment line to use parallelization
                   #, parallel = T, cores = 8
)
test_res$c3=corrSeq_summary(fit, fit_c3, df=NA, sort_results = F)
rm(fit_c3)

############## Difference Between groups at any time #################
fit_c4=corrSeq_fit(~offset(log_offset)+(1|patient), expr_mat=counts, 
                   sample_data = meta_data,method="nbmm_lp"
                   ## Uncomment line to use parallelization
                   #, parallel = T, cores = 8
                   )
test_res$c4=corrSeq_summary(fit, fit_c4, df=NA, sort_results = F)
rm(fit_c4)

############## Add DE to summary tables #################
test_res=lapply(test_res, function(x){
  x$summary_table$de=de
  x
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

############## Calculate non-convergence #########
non_conv = sapply(test_res, function(x) {
  tab = x$summary_table
  stat_calc(tab, stat = "non_conv")
})
