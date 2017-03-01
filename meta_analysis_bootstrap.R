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

# Resampling trial number
numtrials <- 1000


###############################################################################
########################### Functions #########################################
###############################################################################

validationCreate <- function(indices, stdDev, index_size){
  # creates a validation set with given index size and std accompanied by the indices dataframe from each cohort
  # Also comes with a paee
  validation_study = data.frame(index = rep(x = indices$index_indicator, each= index_size))
  validation_study$paee = unlist(unname(lapply(X = split(x = indices, f = as.factor(indices$index_indicator)), 
    FUN = function(x){
      output = rnorm(n = index_size, mean = x$mean, sd = stdDev)
      return (output)
  })))
  return validation_study
}

studyCreate <- function(indices, stdDev, index_size){
  # Creates a study dataset with indices presented, a standard deviation, and the size per index
  # Also comes with a paee and the foo condition
  study_data = data.frame(index = rep(x = indices$index_indicator,each= index_size))
  study_data$paee = unlist(unname(lapply(X = split(x = indices, f = as.factor(indices$index_indicator)), 
    FUN = function(x){
        output = rnorm(n = index_size, mean = x$mean, sd = stdDev)
        return (output)
    })))
  # generate our outcome variable plus some noise
  study_data$foo <-   rnorm(length(study_data$paee),(set_beta*study_data$paee) + constant, 10)
  return study_data
}

###############################################################################
########################### COHORT CREATION ###################################
###############################################################################
## Cohort1
indices1 <- data.frame(mean=c(30,40,50,60), std_dev = c(5,5,5,5), index_indicator = c(1,2,3,4))
study_data1 <- studyCreate(indices1, 5, study_index_size)
validation_data1 <- validationCreate(indices1, 5, 25)

## Cohort2
indices2 <- data.frame(mean=c(25,35,45,55,65), std_dev = c(5,5,5,5,5), index_indicator = c(1,2,3,4,5))
study_data2 <- studyCreate(indices2, 5, study_index_size)
validation_data2 <- validationCreate(indices2, 5, 25)

## Cohort3
indices3 <- data.frame(mean=c(30,35,40,45,50,55), std_dev = c(5,5,5,5,5,5), index_indicator = c(1,2,3,4,5,6))
study_data3 <- studyCreate(indices3, 5, study_index_size)
validation_data3 <- validationCreate(indices3, 5, 25)

## Cohort definition (should change to class if possible)
cohort1 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = validation_data1, study_index_size = study_index_size, study_data = study_data1)
cohort2 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = validation_data2, study_index_size = study_index_size, study_data = study_data2)
cohort3 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = validation_data3, study_index_size = study_index_size, study_data = study_data3)

###############################################################################
########################### Bootstrapping #####################################
###############################################################################
# We 'bootstrap' sampled distributions of the validation sets and perform regressions 
# with the means of those resampled distributions

