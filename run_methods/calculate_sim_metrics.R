#######################################
##
## Script name: Simulation Metric Calculations
##
## Purpose of script: Calculate T1E, FDR, Power, and Convergence rate for simlations
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-30
##
#######################################


## load up the packages we will need:  (comment as required)



####################################################
## 01 - Define function to find T1e/FDR/power/conv
####################################################

############## Simulation Metric Function #################

## Function input:
## res_tab: table with column named raw p_val_raw (for T1E),
##          p_val_adj (for FDR, power), de (differential exprs.)
##          rows with NA p-vals represent non-converged models
## pval_thresh: threshold for calculating T1E, FDR, and power
## stat: simulation statistic to calculate. Options are "T1E", "FDR", "power",
##       and "non-conv" (proportion non-convergence)

stat_calc = function(res_tab, pval_thresh=.05, stat) {
  ## If stat is T1E, use raw p-values, otherwise use adjusted p-values
  pval_var = "p_val_adj"
  if (stat == "T1E")
    pval_var = "p_val_raw"
  else
    pval_var = "p_val_adj"
  
  ## Make two-by-two table of statistically sig. genes by DE genes
  two_by_two = table(res_tab[, pval_var] < pval_thresh, res_tab[, "de"])
  
  ## If table is not two by two:
  if (any(dim(two_by_two) != 2)) {
    ## If no rows/columns, stat=NA
    if (nrow(two_by_two) == 0) {
      ## If only non-sig genes, T1E, power, FDR are 0
      stat_val = NA
    } else if (rownames(two_by_two) == "FALSE") {
      stat_val = switch (
        stat,
        T1E = 0,
        power = 0,
        FDR = 0,
        non_conv = sum(is.na(res_tab[, pval_var]) / length(res_tab[, pval_var]))
      )
      ## Else calculate normally
    } else{
      stat_val = switch(
        stat,
        T1E = two_by_two["TRUE", "FALSE"] / sum(two_by_two[, "FALSE"]),
        FDR = two_by_two["TRUE", "FALSE"] / sum(two_by_two["TRUE", ]),
        power = two_by_two["TRUE", "TRUE"] / sum(two_by_two[, "TRUE"]),
        non_conv = sum(is.na(res_tab[, pval_var]) / length(res_tab[, pval_var]))
      )
    }
    ## Dimensions two-by-two, calculate normally
  } else
    stat_val = switch(
      stat,
      T1E = two_by_two["TRUE", "FALSE"] / sum(two_by_two[, "FALSE"]),
      FDR = two_by_two["TRUE", "FALSE"] / sum(two_by_two["TRUE", ]),
      power = two_by_two["TRUE", "TRUE"] / sum(two_by_two[, "TRUE"]),
      non_conv = sum(is.na(res_tab[, pval_var]) / length(res_tab[, pval_var]))
    )
  return(stat_val)
  
}
