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
bootstrap_trials <- 100
boostrap_index_size <- 10

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
# overall results dataframe
results_df <- data.frame()



# Defining a base generator for one cohort (which spawns validation, study data)
coh_base = data.frame(means = c(30,40,50,60), indices = c(1,2,3,4), std_dev = 5) # change stddev and index size to variables later

# spawn the study data
study_data = createStudyData(coh_base, study_index_size)

# spawn the validation data
validation_data = createValidationData(coh_base, validation_index_size)

# per cohort base result vector
betas <- vector("numeric")
std_errs <- vector("numeric")

# bootstrapping trials times
# spawn the bootstrap of the validation
# Resampling step (ie create the bootstrap)#
for (i in 1:bootstrap_trials){
	bootstrap_validation <- data.frame(index = rep(x = coh_base$indices, each= boostrap_index_size))
	bootstrap_validation$paee <- unlist(unname(lapply(X = split(x=validation_data$paee, f= as.factor(validation_data$index)), 
		FUN = sample, size = boostrap_index_size,replace=TRUE)))

	# Regression and storing of regression coefficients and stderrs
	study_data$paee_sample_ind_mean <- unlist(unname(lapply(X = split(x=bootstrap_validation$paee, f= as.factor(bootstrap_validation$index)),
	  FUN = function(paee_vals){
	    output = rep(x = mean(paee_vals), times = study_index_size)
	    return(output)
	  })))
	reg_out_ind_mean <- lm(formula=foo~paee_sample_ind_mean, data=study_data)
	reg_coeff_ind_mean <- reg_out_ind_mean$coefficients["paee_sample_ind_mean"]
	reg_std_ind_mean <- (summary(reg_out_ind_mean)$coefficients[,"Std. Error"])["paee_sample_ind_mean"]

	betas <- c(betas, reg_coeff_ind_mean)
	std_errs <- c(std_errs, reg_std_ind_mean)
}

results <- as.data.frame(c(1:bootstrap_trials))
colnames(results) <- c("BootStrap Trial")
results$valid_size <- rep(x=validation_index_size, times = bootstrap_trials)
results$reg_coeff_per_mean <- betas
results$reg_stdError_per_mean <- std_errs
results_df <- rbind(results_df, results)

