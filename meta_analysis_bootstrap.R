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

# General Dataset properties (both validation and study)
set_beta <- 0.5
constant <- 20
stdDevs <- c(5,10,15)

# Validation Data specific properties

# Study Data Specific properties
study_index_size <- 5000

###############################################################################
########################### COHORT CREATION ###################################
###############################################################################
## Cohort1
indices1 <- data.frame(mean=c(30,40,50,60), std_dev = c(5,5,5,5), index_indicator = c(1,2,3,4))
# might need to make a for loop to go through the creation of multiple validation sets
validation_data1 <- data.frame(index = rep(x = indices$index_indicator, each = 25))
validation_data1$paee = unlist(unname(lapply(X = split(x = indices, f = as.factor(indices$index_indicator)), 
  FUN = function(x){
    output = rnorm(n = validation_index_size[i], mean = x$mean, sd = x$std_dev)
    return (output)
  })))
study_data1 = data.frame(index = rep(x = indices1$index_indicator,each= study_index_size))
# for each index, pick all the values (defined by study_index_size) at once
study_data1$paee = unlist(unname(lapply(X = split(x = indices1, f = as.factor(indices1$index_indicator)), 
                                      FUN = function(x){
                                          output = rnorm(n = study_index_size, mean = x$mean, sd = x$std_dev)
                                          return (output)
                                      })))
# generate our outcome variable plus some noise
study_data1$foo <-   rnorm(length(study_data$paee),(set_beta*study_data$paee) + constant,10)
cohort1 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = NULL, study_index_size = study_index_size, study_data = NULL)

## Cohort2
indices2 <- data.frame(mean=c(25,35,45,55,65), std_dev = c(5,5,5,5,5), index_indicator = c(1,2,3,4,5))
study_data2 = data.frame(index = rep(x = indices2$index_indicator,each= study_index_size))
# for each index, pick all the values (defined by study_index_size) at once
study_data2$paee = unlist(unname(lapply(X = split(x = indices2, f = as.factor(indices2$index_indicator)), 
                                      FUN = function(x){
                                          output = rnorm(n = study_index_size, mean = x$mean, sd = x$std_dev)
                                          return (output)
                                      })))
# generate our outcome variable plus some noise
study_data2$foo <-   rnorm(length(study_data$paee),(set_beta*study_data$paee) + constant,10)
cohort2 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = NULL, study_index_size = study_index_size, study_data = NULL)

## Cohort3
indices3 <- data.frame(mean=c(30,35,40,45,50,55), std_dev = c(5,5,5,5,5,5), index_indicator = c(1,2,3,4,5,6))
study_data3 = data.frame(index = rep(x = indices3$index_indicator,each= study_index_size))
# for each index, pick all the values (defined by study_index_size) at once
study_data3$paee = unlist(unname(lapply(X = split(x = indices3, f = as.factor(indices3$index_indicator)), 
                                      FUN = function(x){
                                          output = rnorm(n = study_index_size, mean = x$mean, sd = x$std_dev)
                                          return (output)
                                      })))
# generate our outcome variable plus some noise
study_data3$foo <-   rnorm(length(study_data$paee),(set_beta*study_data$paee) + constant,10)
cohort3 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = NULL, study_index_size = study_index_size, study_data = NULL)



