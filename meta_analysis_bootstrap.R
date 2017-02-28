# This is a meta analysis boot strap using the structure of simulations-optimized as a template
# We wish to create intervals of confidence in the meta analysis through bootstrap resampling across
# a variable number of cohort studies

## Author: Paul Scherer and Tom Bishop
## Date: 28.02.2017

###############################################################################
###################### R Environment Settings #################################
###############################################################################
setwd('V:/Studies/InterConnect/Internal/Latent variable harmonisation/plots')
library(ggplot2)

###############################################################################
########################### DATA AND SETTINGS #################################
###############################################################################
seed <- 54  
set.seed(seed)

# Dataset properties (both validation and study)
set_beta <- 0.5
constant <- 20

stdDevs <- c(5,10,15)
final_output_list <- vector("list", length(stdDevs))


# each study will be based on a set of standard parameters, which will spawn a validation and study data set per cohort
# so might as well have a cohort object that contains:
# definitions for index gaussians
# study data:
## study data size 
## study data gaussians
# validation data:
## validation_data size

## Cohort1
indices1 <- data.frame(mean=c(30,40,50,60), std_dev = c(5,5,5,5), index_indicator = c(1,2,3,4))
validation_data1 <- data.frame(index = rep(x = indices$index_indicator, each = 25))
validation_data1$paee = unlist(unname(lapply(X = split(x = indices, f = as.factor(indices$index_indicator)), 
  FUN = function(x){
    output = rnorm(n = validation_index_size[i], mean = x$mean, sd = x$std_dev)
    return (output)
  })))
cohort1 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = validation_data1, study_index_size = 5000, study_data = study_data1)

## Cohort2
cohort2 <-

## Cohort3
cohort3 <- 

results_df <- data.frame()