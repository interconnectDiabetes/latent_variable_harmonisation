###### Simulation Practice for Latent Variable Harmonisation ########

# The main objective of this program is to explore different data generation
# distributions and their relationship with regression calibration.

## Author: Paul Scherer
## Date: 27/09/2016

###############################################################################
########################### DATA AND SETTINGS #################################
###############################################################################
## Libraries
library("survival")
library(graphics)

## Seed
set.seed(66)

## Parameters
# The linear correlation coefficient between exposure and outcome.
set_beta <- 0.5

# Index properties (Assumption that they are Gaussian)
# Format : (mean, stdev)
index_mean1 <- 35.6
index_stdev1 <- 13.7
index_mean2 <- 43.7
index_stdev2 <- 15.2
index_mean3 <- 49.0
index_stdev3 <- 17.9
index_mean4 <- 56.2
index_stdev4 <- 18.4

mean_paee_difference <- abs(mean(c(index_mean1-index_mean2, index_mean2-index_mean3, index_mean3-index_mean4)))

# PAEE range
paee_range_min <- 20
paee_range_max <- 75

###############################################################################
############################# Functions #######################################
###############################################################################
data_generator <- function(exposure, beta=set_beta, constant=20, noise=2.5){
	# returns a datapoint given an exposure
	# relies on a set constant, gaussian noise component
	fbeta <- (beta*exposure) + constant
	data_point <- rnorm(1,fbeta,noise)
	return (data_point)
}

gaussian_index_sample <- function(x){
	# returns a datapoint sampled from the index distribution
	# This can be seen as an exposure to be used in the data_generator
	# :param: index = cambridge index
	if (x == 1){
		index_mean <- index_mean1
		index_stdev <- index_stdev1
	} else if (x == 2){
		index_mean <- index_mean2
		index_stdev <- index_stdev2
	} else if (x == 3){
		index_mean <- index_mean3
		index_stdev <- index_stdev3
	} else {
		index_mean <- index_mean4
		index_stdev <- index_stdev4
	}
	data_point <- rnorm(1, index_mean, index_stdev)
	return (data_point)
}


lambda_collector_cam_index <- function(data_set_outcome, data_set_exposure, data_set){
	rdr_regression_fit <- lm(formula=data_set_outcome~data_set_exposure, data=data_set)
	lambda <- rdr_regression_fit$coefficients["cam_index_means"]
	return (lambda)
}

calculate_index_from_paee <- function(paee){
	# take paee value and calculate probability of it belonging to index
	# random weighted assignment into that category
	prob1 <- pnorm(paee, mean=index_mean1, sd=index_stdev1)
	prob2 <- pnorm(paee, mean=index_mean2, sd=index_stdev2)
	prob3 <- pnorm(paee, mean=index_mean3, sd=index_stdev3)
	prob4 <- pnorm(paee, mean=index_mean4, sd=index_stdev4)

	probs_cam_index <- c(prob1, prob2, prob3, prob4)
	index_num <- sample(c(1:4), 1, replace=TRUE, prob=probs_cam_index )

	return (index_num)
}
####################################################################################################
################################ Creation of datasets ##############################################
####################################################################################################
# Data set of 1600 points for which we calculate the lambda, (validation)
validation_data <- as.data.frame(c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000)))
colnames(validation_data) <- c("cam_index")
validation_data$cam_index <- as.factor(validation_data$cam_index)
validation_data$cam_index <- validation_data$cam_index
validation_data$paee <- unlist(lapply(validation_data$cam_index, gaussian_index_sample))
validation_data$cam_index_means <- unlist(mapply(validation_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = index_mean1
    } else if (x == 2) {
      output = index_mean2
    } else if (x == 3) {
      output = index_mean3
    } else if (x == 4) {
      output = index_mean4
    } else {
      output = NA
  	}
  	return(output)
	}
))


# Data set of 20000 points for which we calculate the betas, (test)
study_data <- as.data.frame(c(rep(1,5000),rep(2,5000),rep(3,5000),rep(4,5000)))
colnames(study_data) <- c("cam_index")
study_data$cam_index <- as.factor(study_data$cam_index)
study_data$paee <- unlist(lapply(study_data$cam_index, gaussian_index_sample))
study_data$foo <- unlist(lapply(study_data$paee, data_generator, beta=set_beta))
study_data$cam_index_means <- unlist(mapply(study_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = index_mean1
    } else if (x == 2) {
      output = index_mean2
    } else if (x == 3) {
      output = index_mean3
    } else if (x == 4) {
      output = index_mean4
    } else {
      output = NA
  	}
  	return(output)
	}
))

################################### BASELINE with cam_index ###################################################
# calculation of correlation between foo and the cam_index as measurement
# coef = beta
foo_cam_study <- lm(formula=foo~cam_index, data=study_data)
coef_foo_cam_study <- foo_cam_study$coefficients[-1]
coef_foo_cam_study <- c(1,coef_foo_cam_study)
mean_coef <- mean(coef_foo_cam_study)
study_per_paee_estimate_cam_index <- mean_coef/mean_paee_difference # rough translation to per paee increase
## std error and variance
stdError_cc_study_data <- (summary(foo_cam_study)$coefficients[,"Std. Error"])[-1]
mean_stdError_coef_study <- mean(stdError_cc_study_data)
mean_stdError_coef_study <- mean_stdError_coef_study/mean_paee_difference
mean_variance_coef_study <- mean_stdError_coef_study^2
# confidence intervals
lci_beta_hat_star <- study_per_paee_estimate_cam_index - 1.96*mean_stdError_coef_study
uci_beta_hat_star <- study_per_paee_estimate_cam_index + 1.96*mean_stdError_coef_study


## coef = lambda
paee_cam_val <- lm(formula=paee~cam_index, data=validation_data)
coef_paee_cam_val <- paee_cam_val$coefficients[-1]
coef_paee_cam_val <- c(1,coef_paee_cam_val)
mean_coef <- mean(coef_paee_cam_val)
val_per_paee_coef <- mean_coef/mean_paee_difference
## std error and variance
stdError_cc_val_data <- (summary(paee_cam_val)$coefficients[,"Std. Error"])[-1]
mean_stdError_cc_val <- mean(stdError_cc_val_data)
mean_variance_cc_val <- mean_stdError_cc_val^2
# confidence intervals
lci_lambda_hat <- val_per_paee_coef - 1.96*mean_stdError_cc_val
uci_lambda_hat <- val_per_paee_coef + 1.96*mean_stdError_cc_val


## RC 
beta_hat_star <- study_per_paee_estimate_cam_index
lambda_hat <- val_per_paee_coef
beta_hat <- beta_hat_star/lambda_hat
var_beta <- mean_variance_coef_study/(lambda_hat^2) + (beta_hat_star/lambda_hat^2)^2*mean_variance_cc_val
# confidence intervals
lci_beta_hat <- beta_hat - 1.96*sqrt(var_beta)  
uci_beta_hat <- beta_hat + 1.96*sqrt(var_beta)   


# ___  ___ _____ _____ _   _ ___________   __  
# |  \/  ||  ___|_   _| | | |  _  |  _  \ /  | 
# | .  . || |__   | | | |_| | | | | | | | `| | 
# | |\/| ||  __|  | | |  _  | | | | | | |  | | 
# | |  | || |___  | | | | | \ \_/ / |/ /  _| |_
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \___/
###############################################################################################################
################################## WITH CAM_INDEX_MEANS #######################################################
## coef = beta
pc_study_data_means <- lm(formula=foo~cam_index_means, data=study_data)
cc_study_data_means <- pc_study_data_means$coefficients["cam_index_means"]
# stdError
stdError_cc_study_data_means <- (summary(pc_study_data_means)$coefficients[,"Std. Error"])[-1]
# confidence intervals
lci_beta_hat_star_means <- cc_study_data_means - 1.96*stdError_cc_study_data_means
uci_beta_hat_star_means <- cc_study_data_means + 1.96*stdError_cc_study_data_means


## coef = lambda
pc_val_data_means <- lm(formula=paee~cam_index_means, data=validation_data)
cc_val_data_means <- pc_val_data_means$coefficients["cam_index_means"]
# stdError
stdError_cc_val_data_means <- (summary(pc_val_data_means)$coefficients[,"Std. Error"])[-1]
# confidence intervals
lci_lambda_hat_means <- cc_val_data_means - 1.96*(stdError_cc_val_data_means)
uci_lambda_hat_means <- cc_val_data_means + 1.96*(stdError_cc_val_data_means)


# RC 
beta_hat_star_means <- cc_study_data_means
lambda_hat_means <- cc_val_data_means
beta_hat_means <- beta_hat_star_means/lambda_hat_means
var_beta_means <- stdError_cc_study_data_means/(lambda_hat_means^2) + (beta_hat_star_means/lambda_hat_means^2)^2*stdError_cc_val_data_means
# confidence intervals
uci_beta_hat_means <- beta_hat_means - 1.96*sqrt(var_beta_means)
lci_beta_hat_means <- beta_hat_means + 1.96*sqrt(var_beta_means)

# ___  ___ _____ _____ _   _ ___________   _____   ___  
# |  \/  ||  ___|_   _| | | |  _  |  _  \ / __  \ / _ \ 
# | .  . || |__   | | | |_| | | | | | | | `' / /'/ /_\ \
# | |\/| ||  __|  | | |  _  | | | | | | |   / /  |  _  |
# | |  | || |___  | | | | | \ \_/ / |/ /  ./ /___| | | |
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \_____/\_| |_/                                     
###############################################################################################################
################################## WITH SAMPLED PAEE VALS #####################################################
# Set lists of values for each PA categorization that will be used to draw random numbers for
# validation
cat1 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat1 <- cat1[!sapply(cat1,is.na)]

cat2 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat2 <- cat2[!sapply(cat2,is.na)]

cat3 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat3 <- cat3[!sapply(cat3,is.na)]

cat4 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
})
cat4 <- cat4[!sapply(cat4,is.na)]


# empty list initialization
beta_hat_star_sample_list <- vector('numeric')
lambda_hat_sample_list <- vector('numeric')
beta_hat_sample_list <- vector('numeric')
uci_beta_hat_star_sample_list <- vector('numeric')
lci_beta_hat_star_sample_list <- vector('numeric')
uci_lambda_hat_sample_list <- vector('numeric')
lci_lambda_hat_sample_list <- vector('numeric')
uci_beta_hat_sample_list <- vector('numeric')
lci_beta_hat_sample_list <- vector('numeric')

for (i in 1:10) {
	cat1_choice <- sample(cat1,1,replace=TRUE)
	cat2_choice <- sample(cat2,1,replace=TRUE)
	cat3_choice <- sample(cat3,1,replace=TRUE)
	cat4_choice <- sample(cat4,1,replace=TRUE)

	study_data$cam_index_means_sample <- unlist(mapply(study_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
	    if (is.na(x)) {
	      output = NA
	    } else if (x == 1){
	      output = cat1_choice
	    } else if (x == 2) {
	      output = cat2_choice
	    } else if (x == 3) {
	      output = cat3_choice
	    } else if (x == 4) {
	      output = cat4_choice
	    } else {
	      output = NA
	  	}
	  	return(output)
		}
	))

	validation_data$cam_index_means_sample <- unlist(mapply(validation_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
	    if (is.na(x)) {
	      output = NA
	    } else if (x == 1){
	      output = cat1_choice
	    } else if (x == 2) {
	      output = cat2_choice
	    } else if (x == 3) {
	      output = cat3_choice
	    } else if (x == 4) {
	      output = cat4_choice
	    } else {
	      output = NA
	  	}
	  	return(output)
		}
	))

	# calculation of beta hat star
	regression_study <- lm(formula=foo~cam_index_means_sample, data=study_data)
	beta_hat_star_sample <- regression_study$coefficients["cam_index_means_sample"]
	beta_hat_star_sample_list <- c(beta_hat_star_sample_list, beta_hat_star_sample)
	#stdError and confidence intervals
	stdError_beta_hat_star_sample <- (summary(regression_study)$coefficients[,"Std. Error"])[-1]
	uci_beta_hat_star_sample <- beta_hat_star_sample - 1.96*stdError_beta_hat_star_sample
	lci_beta_hat_star_sample <- beta_hat_star_sample + 1.96*stdError_beta_hat_star_sample
	uci_beta_hat_star_sample_list <- c(uci_beta_hat_star_sample_list, uci_beta_hat_star_sample)
	lci_beta_hat_star_sample_list <- c(lci_beta_hat_star_sample_list, lci_beta_hat_star_sample)

	# calculation of lambda hat
	regression_lambda <- lm(formula=paee~cam_index_means_sample, data=validation_data)
	lambda_hat_sample <- regression_lambda$coefficients["cam_index_means_sample"]
	lambda_hat_sample_list <- c(lambda_hat_sample_list, lambda_hat_sample)
	#stdError and confidence intervals
	stdError_lambda_hat_sample <- (summary(regression_lambda)$coefficients[,"Std. Error"])[-1]
	uci_lambda_hat_sample <- lambda_hat_sample + 1.96*(stdError_lambda_hat_sample)
	lci_lambda_hat_sample <- lambda_hat_sample + 1.96*(stdError_lambda_hat_sample)
	uci_lambda_hat_sample_list <- c(uci_lambda_hat_sample_list, uci_lambda_hat_sample)
	lci_lambda_hat_sample_list <- c(lci_lambda_hat_sample_list, lci_lambda_hat_sample)

	# calculation of calibrated beta
	beta_hat_sample <- beta_hat_star_sample/lambda_hat_sample
	beta_hat_sample_list <- c(beta_hat_sample_list,beta_hat_sample)
	#stdError and confidence intervals
	var_beta_sample <- stdError_beta_hat_star_sample/(lambda_hat_sample^2) + (beta_hat_star_sample/lambda_hat_sample^2)^2*stdError_lambda_hat_sample
	uci_beta_hat_sample <- beta_hat_sample - 1.96*sqrt(var_beta_sample)
	lci_beta_hat_sample <- beta_hat_sample + 1.96*sqrt(var_beta_sample)
	uci_beta_hat_sample_list <- c(uci_beta_hat_sample_list, uci_beta_hat_sample)
	lci_beta_hat_sample_list <- c(lci_beta_hat_sample_list, lci_beta_hat_sample)
}

mean_beta_hat_star_sample <- mean(beta_hat_star_sample_list)
mean_lambda_hat_sample <- mean(lambda_hat_sample_list)
mean_beta_hat_sample <- mean(beta_hat_sample_list)
mean_uci_beta_hat_star_sample <- mean(uci_beta_hat_star_sample_list)
mean_lci_beta_hat_star_sample <- mean(lci_beta_hat_star_sample_list)
mean_uci_lambda_hat_sample <- mean(uci_lambda_hat_sample_list)
mean_lci_lambda_hat_sample <- mean(lci_lambda_hat_sample_list)
mean_uci_beta_hat_sample <- mean(uci_beta_hat_sample_list)
mean_lci_beta_hat_sample <- mean(lci_beta_hat_sample_list)

# ___  ___ _____ _____ _   _ ___________   _____ ______ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \ / __  \| ___ \
# | .  . || |__   | | | |_| | | | | | | | `' / /'| |_/ /
# | |\/| ||  __|  | | |  _  | | | | | | |   / /  | ___ \
# | |  | || |___  | | | | | \ \_/ / |/ /  ./ /___| |_/ /
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \_____/\____/ 
                                                      

# ___  ___ _____ _____ _   _ ___________   _____ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \ |____ |
# | .  . || |__   | | | |_| | | | | | | |     / /
# | |\/| ||  __|  | | |  _  | | | | | | |     \ \
# | |  | || |___  | | | | | \ \_/ / |/ /  .___/ /
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \____/                                   
# Generating samples from fitted distributions based on maximum likelihood. Similar to the smapling method but
# we arent so reliant on the data distributions. 
# A distribution of paee values is fit for each index. Samples generated from this distribution are applied 
# to members of the study data
cat1 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat1 <- cat1[!sapply(cat1,is.na)]

cat2 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat2 <- cat2[!sapply(cat2,is.na)]

cat3 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat3 <- cat3[!sapply(cat3,is.na)]

cat4 <- mapply(validation_data$paee, validation_data$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
})
cat4 <- cat4[!sapply(cat4,is.na)]


# fitting of distributions via MLE
library(fitdistrplus)
# normal distributions
dist1 <- fitdist(cat1, "norm", method='mle')
dist2 <- fitdist(cat2, "norm", method='mle')
dist3 <- fitdist(cat3, "norm", method='mle')
dist4 <- fitdist(cat4, "norm", method='mle')

dist_list = list(dist1,dist2,dist3,dist4)

for (i in 1:4){
  plot(dist_list[[i]])
}


# take the sufficient statistics to recreate the distributions from which we can sample
# Index properties (Assumption that they are Gaussian)
# Format : (mean, stdev)
index_mean1 <- dist1[[1]][1]
index_stdev1 <- dist1[[1]][2]
index_mean2 <- dist2[[1]][1]
index_stdev2 <- dist2[[1]][2]
index_mean3 <- dist3[[1]][1]
index_stdev3 <- dist3[[1]][2]
index_mean4 <- dist4[[1]][1]
index_stdev4 <- dist4[[1]][2]
mean_paee_difference <- abs(mean(c(index_mean1-index_mean2, index_mean2-index_mean3, index_mean3-index_mean4)))

validation_data$cam_index_dfit <- unlist(mapply(validation_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else {
      output = gaussian_index_sample(x)
  	}
  	return(output)
	}
))

study_data$cam_index_dfit <- unlist(mapply(study_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else {
      output = gaussian_index_sample(x)
  	}
  	return(output)
	}
))

## coef = beta
pc_study_data_dfit <- lm(formula=foo~cam_index_dfit, data=study_data)
cc_study_data_dfit <- pc_study_data_dfit$coefficients["cam_index_dfit"]
# stdError
stdError_cc_study_data_dfit <- (summary(pc_study_data_dfit)$coefficients[,"Std. Error"])[-1]
# confidence intervals
lci_beta_hat_star_dfit <- cc_study_data_dfit - 1.96*stdError_cc_study_data_dfit
uci_beta_hat_star_dfit <- cc_study_data_dfit + 1.96*stdError_cc_study_data_dfit

## coef = lambda
pc_val_data_dfit <- lm(formula=paee~cam_index_dfit, data=validation_data)
cc_val_data_dfit <- pc_val_data_dfit$coefficients["cam_index_dfit"]
# stdError
stdError_cc_val_data_dfit <- (summary(pc_val_data_dfit)$coefficients[,"Std. Error"])[-1]
# confidence intervals
lci_lambda_hat_dfit <- cc_val_data_dfit - 1.96*(stdError_cc_val_data_dfit)
uci_lambda_hat_dfit <- cc_val_data_dfit + 1.96*(stdError_cc_val_data_dfit)

# RC 
beta_hat_star_dfit <- cc_study_data_dfit
lambda_hat_dfit <- cc_val_data_dfit
beta_hat_dfit <- beta_hat_star_dfit/lambda_hat_dfit
var_beta_dfit <- stdError_cc_study_data_dfit/(lambda_hat_dfit^2) + (beta_hat_star_dfit/lambda_hat_dfit^2)^2*stdError_cc_val_data_means
# confidence intervals
uci_beta_hat_dfit <- beta_hat_dfit - 1.96*sqrt(var_beta_dfit)
lci_beta_hat_dfit <- beta_hat_dfit + 1.96*sqrt(var_beta_dfit)

# ___  ___ _____ _____ _   _ ___________     ___ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \   /   |
# | .  . || |__   | | | |_| | | | | | | |  / /| |
# | |\/| ||  __|  | | |  _  | | | | | | | / /_| |
# | |  | || |___  | | | | | \ \_/ / |/ /  \___  |
# \_|  |_/\____/  \_/ \_| |_/\___/|___/       |_/
#                                                
# Density estimation for each index using kernel based estimation methods. Then sample directly using the kernels.
library(stats)

density1 <- density(cat1)
density2 <- density(cat2)
density3 <- density(cat3)
density4 <- density(cat4)

density_list = list(density1,density2,density3,density4)

for (i in 1:4){
  plot(density_list[[i]])
}