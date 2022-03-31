# multiDF_corr_RNASeq

This repository contains code for simulating longitudinal RNA-Seq data and performing multiple degree of freedom tests for differential expression using several different methods.
There are two folders:
* simulate_data contains code for simulating RNA-Seq data with two groups (control and treatment) and multiple timepoints per subject. The code includes options to choose the number of subjects, number of repeatedmeasures per subject and the amount and magnitude of differential expression. An .RDS object (all_param_trip.RDS) is also included with parameters estimated from real data that are used to simulated the data.
* run_methods contains scripts for running multi-DF tests on simulated data using several different methods as well as calculating type 1 error, false discover rate, power, and the proportion of converged models. An .RDS file with an example simulated data set is included to evaluate the methods. The following multi-df tests are implemented for each method:
    *   Are there differences in expression between the treatment and control at any of the time points?
    *   Is there a change in gene expression between any timepoints for the treatment group?
    *   Are there any significant interaction effects?
    *    Are there any significant model coefficients?

  
