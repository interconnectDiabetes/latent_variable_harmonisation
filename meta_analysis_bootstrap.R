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

# Results dataframes and lists

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
  return (validation_study)
}

validationListCreate <- function(indices, stdDev_list, val_index_size_list){
  validationList <- vector("list", length(stdDev_list))
  for (s in 1:length(stdDev_list)){
    perStdDev <- vector("list", length(val_index_size_list))
    for (v in 1:length(val_index_size_list)){
      validationData <- validationCreate(indices, stdDev_list[s], val_index_size_list[v])
      perStdDev[v] = validation_data
    }
    validationList[s] = perStdDev
  }
  return(validationList)
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
  return (study_data)
}

studyListCreate <- function(indices, stdDev_list, index_size){
  studyList <- vector("list", length(stdDev_list))
  for (s in 1:length(stdDev_list)){
    studyData <- studyCreate(indices, stdDev_list[s], index_size)
    studyList[s] <- studyData
  }
  return (studylist)
}

###############################################################################
########################### COHORT CREATION ###################################
###############################################################################
## Cohort1
indices1 <- data.frame(mean=c(30,40,50,60), std_dev = stdDevs, index_indicator = c(1,2,3,4))
study_list1 <- vector("list", length(stdDevs))
validation_list1 <- vector("list", length(stdDevs))
study_data1 <- studyCreate(indices1, 5, study_index_size)
validation_data1 <- validationCreate(indices1, 5, 25)

## Cohort2
indices2 <- data.frame(mean=c(25,35,45,55,65), std_dev = stdDevs, index_indicator = c(1,2,3,4,5))

study_data2 <- studyCreate(indices2, 5, study_index_size)
validation_data2 <- validationCreate(indices2, 5, 25)

## Cohort3
indices3 <- data.frame(mean=c(30,35,40,45,50,55), std_dev = stdDevs, index_indicator = c(1,2,3,4,5,6))

study_data3 <- studyCreate(indices3, 5, study_index_size)
validation_data3 <- validationCreate(indices3, 5, 25)



## Cohort definition (should change to class if possible)
cohort1 <- list(indices = indices1, validation_index_size = c(25,50,100), validation_data = validation_list1, study_data = study_list1)
cohort2 <- list(indices = indices2, validation_index_size = c(25,50,100), validation_data = validation_list2, study_data = study_list2)
cohort3 <- list(indices = indices3, validation_index_size = c(25,50,100), validation_data = validation_list3, study_data = study_list3)

###############################################################################
########################### Bootstrapping #####################################
###############################################################################
# We 'bootstrap' sampled distributions of the validation sets and perform regressions 
# with the means of those resampled distributions
num_resample_samples <- 5

# per cohort
# per validation set of cohort

# Initialise some structures to store coeffcients
results_df <- data.frame()
betas <- vector("numeric")
std_errs <- vector("numeric")

# Resampling step (ie create the bootstrap) 
bootstrap_validation <- data.frame(index = rep(x = indices1$index_indicator, each= num_resample_samples))
bootstrap_validation$paee <- unlist(unname(lapply(X = split(x=validation_data1$paee, f= as.factor(validation_data1$index)), FUN = sample, size = num_resample_samples,replace=TRUE)))

# regression step (take the means of the bootstrap set and do regressions)
study_data1$paee_sample_ind_mean <- unlist(unname(lapply(X = split(x=bootstrap_validation$paee, f= as.factor(bootstrap_validation$index)),
                                                      FUN = function(paee_vals){
                                                        output = rep(x = mean(paee_vals), times = study_index_size)
                                                        return(output)
                                                      })))

reg_out_ind_mean <- lm(formula=foo~paee_sample_ind_mean, data=study_data1)
reg_coeff_ind_mean <- reg_out_ind_mean$coefficients["paee_sample_ind_mean"]
reg_std_ind_mean <- (summary(reg_out_ind_mean)$coefficients[,"Std. Error"])["paee_sample_ind_mean"]


























# Store values into dataframe
results <- as.data.frame(c(1:numtrials))
colnames(results) <- c("Trial")
results$valid_size <- rep(x=validation_index_size[i], times = numtrials)
results$reg_coeff_per_mean <- reg_coeff_ind_mean
results$reg_stdError_per_mean <- reg_std_ind_mean
results_df <- rbind(results_df, results)


# Summarizing the results dataframe
temp_output <- aggregate(results_df[,3:8], by=list(valid_size = results_df$valid_size), quantile, probs=c(0.05,0.5,0.95), names=TRUE)
final_output = data.frame(validation_size = temp_output[,1])

for (k in 2:ncol(temp_output)){
  temp = as.data.frame(temp_output[,k])
  colnames(temp) <- paste(colnames(temp_output)[k], colnames(temp), sep = "_")
  final_output = cbind(final_output, temp)
}