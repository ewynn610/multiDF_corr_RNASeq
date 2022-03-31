#######################################
##
## Script name: Simulate data
##
## Purpose of script: Create Simulated Datasets
##
## Author: Elizabeth Wynn
##
## Date Created: 2022-03-31
##
#######################################


## load up the packages we will need:  (comment as required)

library(MASS) ## rnegbin function

## Read in simulation parameters
all_param_trip <- readRDS("simulate_data/all_param_trip.RDS")


#######################################
## 01 - Create Simulation Function
#######################################


sim_dat <- function(n_features = 500,      ## Number of genes to simulate
                    n_per_group = 5,       ## Number of subjects per group
                    n_times = 4,           ## Number of timepoints
                    n_groups = 2,          ## Number of groups
                    prop_DE = .2,          ## Proportion of DE
                    log_effect_size = 0.5, ## Log effect size for DE genes
                    ri_sd,                 ## Random intercept SD
                    beta0s,                ## Baseline means
                    dispersion            ##dispersions
){ 
  
  ## Generate ID's
  ids<-seq(1, n_groups*n_per_group, 1)
  
  ## Generate group labels
  groups<-NULL
  for(i in 0:(n_groups-1)){groups<-c(groups, rep(i, n_per_group))}
  
  ## Meta data
  temp<-data.frame(id=rep(ids,n_times), group=rep(groups,n_times),
                   time=rep(seq(0,n_times-1,1),each=n_groups*n_per_group))
  
  ## Assign DE genes and model parameters
  n_DE <- round(n_features * prop_DE)
  idx_DE <- sample(1:n_features, size = n_DE, replace = F)
  beta_groups <- 0 * numeric(n_features)
  beta_times <- 0 * numeric(n_features)
  beta_ints <- 0 * numeric(n_features)
  case<-sample(seq(1:2), n_DE, replace=T)
  beta_ints[idx_DE]<-ifelse(case==1, log_effect_size, ifelse(case==2, -log_effect_size, 0))
  beta_tot <- cbind(beta0s, beta_groups, beta_times, beta_ints)
  
  ## Make design matrix
  design_mat <-model.matrix(~temp$group+temp$time+temp$group*temp$time)
  
  ## Simulate counts
  counts <- do.call(rbind, lapply(1:n_features, function(i){
    
    ## Draw random intercepts
    r.ints<-rnorm(n_per_group*n_groups, mean=0, sd=ri_sd[i])
    ri<-rep(r.ints,n_times)
    
    ## Get means
    means_sub <- exp(design_mat %*% beta_tot[i, ] + ri)
    
    ## Draw counts
    counts_sub <- rnegbin(n_groups*n_times*n_per_group, mu = means_sub, theta = 1 / dispersion[i])
    return(counts_sub)
  }))
  
  ## Count sums
  cs_tmp <- colSums(counts / 1e6)
  
  ## Find samples with low counts
  idx=which(cs_tmp>40)
  
  ## Get other samples from individual who had samples with low counts
  idx_mate=which(temp$id %in% temp$id[idx])
  
  ## Which columns to resimulate
  idx=unique(sort(c(idx, idx_mate)))
  ## While there are still columns with less than 40 cpm, resimulate all samples from those individuals
  while(length(idx>0)){
    counts_idx <- do.call(rbind, lapply(1:n_features, function(i){
      r.ints<-rnorm(length(idx)/n_times, mean=0, sd=ri_sd[i]) #added in random intercepts
      ri<-rep(r.ints,n_times)
      means_sub <- exp(design_mat[idx,] %*% beta_tot[i, ] + ri)
      counts_sub <- rnegbin(length(idx), mu = means_sub, theta = 1 / dispersion[i])
      return(counts_sub)
    }))
    counts[,idx]<-counts_idx
    cs_tmp <- colSums(counts / 1e6)
    
    if(sum(complete.cases(counts)==0)) browser()
    
    
    ## Find samples with low counts
    idx=which(cs_tmp>40)
    
    ## Get other samples from individual who had samples with low counts
    idx_mate=which(temp$id %in% temp$id[idx])
    
    ## Which columns to resimulate
    idx=unique(sort(c(idx, idx_mate)))
    
  }
  
  ## Random intercept variance
  ri_var <- ri_sd^2
  
  ##
  de=rep(F, nrow(counts))
  de[idx_DE]=T
  
  ## Return simulation results
  return(list(counts = counts,
              groups = temp$group,
              time = temp$time,
              ids = temp$id,
              de=de, #row numbers of the differentially expressed genes
              param = cbind(beta_tot, dispersion, ri_var)))
}



#######################################
## 02 - #Simulate DE data 
#######################################


## Simulation settings
## Simulate 4 timepoints, 2 groups, and 3 samples per group
n_times=4
PD <- 0.2
n_gene <- 15000
n_sample_per_group <- 3
le_size <- 1/(n_times-1)
total_lib_size <- 25e6



## Get parameter values
idx_param <- sample(1:nrow(all_param_trip), size = n_gene, replace = T)
sim_param <- all_param_trip[idx_param, ]
sim_param$baseline_prop <- sim_param$mean_cpm / sum(sim_param$mean_cpm)
sim_baselines <- log(sim_param$baseline_prop * total_lib_size)

## Simulate data
suppressWarnings(
  sim_data <- sim_dat(n_features = n_gene,
                      prop_DE = PD,
                      log_effect_size = le_size,
                      n_per_group = n_sample_per_group,
                      n_times = n_times,
                      n_groups = 2,
                      beta0s = sim_baselines,
                      ri_sd = sqrt(sim_param$ri_var),
                      dispersion = sim_param$disp)
)


## Only keep genes with >1 cpm in >= n_sample_per Group
cpm_table <- apply(sim_data$counts, 2, function(x){
  return(x * 1e6 / sum(x))
})
summary(rowSums(1*(cpm_table > 1)) >= n_sample_per_group)
idx_keep <- which(rowSums(1*(cpm_table > 1)) >= n_sample_per_group)
sim_data$counts <- sim_data$counts[idx_keep, ]
sim_data$param <- sim_data$param[idx_keep, ]
sim_data$de <- sim_data$de[idx_keep]


