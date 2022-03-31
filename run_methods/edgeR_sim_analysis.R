#######################################
##
## Script name: edgeR Simulation Analysis
##
## Purpose of script: Analyze simulated data using edgeR (ignoring correlation)
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-28
##
#######################################

## load up the packages we will need:  (comment as required)

library(edgeR)

## Source in function to calculate simulation metrics
source("calculate_sim_metrics.R")


#######################################
## 01 - Read in data/contrasts and format
#######################################

############## Read in data/contrasts #################
sim_data <- readRDS("sim_data.RDS")
contrasts<-readRDS("contrasts.RDS")


############## Extract Counts/meta data #################
counts <- as.matrix(sim_data$counts)

meta_data=data.frame(patient=factor(sim_data$ids),
                     time=factor(sim_data$time),
                     group=factor(sim_data$groups))

de=sim_data$de

############## Format Contrasts for edgeR #################
contrasts=lapply(contrasts, function(x){
  x=t(rbind(x))
})


#######################################
## 02 - Run edgeR
#######################################

############## Fit edgeR models #################
design_mat <- model.matrix(~group*time, data=meta_data)

dge <- DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
dge <- edgeR::estimateDisp(dge, design=design_mat)

fit <- glmFit(dge, design_mat)

############## Run LRT test for each contrast #################
lrt_res=lapply(contrasts,function(x){
  tab=glmLRT(fit, contrast=x)$table
  tab$p_val_adj=p.adjust(tab$PValue, method="BH")
  tab$de=de
  colnames(tab)[colnames(tab)=="PValue"]="p_val_raw"
  tab
})

#######################################
## 03 - Calculate simulation metrics
#######################################

############## Calculate T1E for all contrast #################
T1E_vals=sapply(lrt_res, stat_calc, pval_thresh=.05, stat="T1E")

############## Calculate FDR for all contrast #################
FDR_vals=sapply(lrt_res, stat_calc, pval_thresh=.05, stat="FDR")

############## Calculate power for all contrast #################
power_vals=sapply(lrt_res, stat_calc, pval_thresh=.05, stat="power")

############## Calculate non-convergence (same for all contrasts)#########
non_conv=stat_calc(lrt_res$c1, stat="non_conv")
