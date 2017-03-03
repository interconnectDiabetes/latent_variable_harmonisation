# BootStrap the Bootstrap
# Author: Paul Scherer and Tom Bishop
# Date: 03.02.2017


###############################################################################
###################### R Environment Settings #################################
###############################################################################
setwd('V:/Studies/InterConnect/Internal/Latent variable harmonisation/plots')
library(ggplot2)

###############################################################################
########################### DATA AND SETTINGS #################################
###############################################################################
seed <- 66
set.seed(seed)

# Base Data set properties
set_beta <- 0.5
constant <- 20

# Study Data properties
study_index_size <- 5000

# Validation Data properties
validation_index_sizes <- c(25,50,100)
validation_index_size <- 25 # temp

# Bootstrapping properties
num_trials <- 1000

###############################################################################
########################### Functions #########################################
###############################################################################
createStudyData <- function(coh_base, study_index_size){
	# given a cohort base, ie means, index, stdev it returns a study data dataframe
	# at this stage it can only accomodate the same index size for each index
	study_data = data.frame(index = rep(x = coh_base$indices, each=study_index_size))
	study_data$paee = unlist(unname(lapply(X = split(x = coh_base, f = as.factor(coh_base$indices)), 
		FUN = function(x){
			output = rnorm(n =study_index_size, mean = x$means, sd = x$std_dev)
			return (output)
		})))
	study_data$foo <- rnorm(length(study_data$paee), (set_beta*study_data$paee) + constant, 10)
	return(study_data)
}

createValidationData <- function(coh_base, validation_index_size) {
	# given a cohort base it generates the validation data
	validation_data = data.frame(index = rep(x = coh_base$indices, each = validation_index_size))
	validation_data$paee = unlist(unname(lapply(X = split(x = coh_base, f = as.factor(coh_base$indices)), 
		FUN=function(x){
			output = rnorm(n = validation_index_size, mean=x$mean, sd=x$std_dev)
		})))
	return (validation_data)
}

###############################################################################
########################### COHORT CREATION ###################################
###############################################################################
# Defining a base generator for one cohort (which spawns validation, study data)
coh_base = data.frame(means = c(30,40,50,60), indices = c(1,2,3,4), std_dev = 5) # change stddev and index size to variables later

# spawn the study data
study_data = createStudyData(coh_base, study_index_size)

# spawn the validation data
validation_data = createValidationData(coh_base, validation_index_size)

# spawn the bootstrap of the validation