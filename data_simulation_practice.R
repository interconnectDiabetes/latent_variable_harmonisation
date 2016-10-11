###### Simulation Practice for Latent Variable Harmonisation ########

# The main objective of this program is to explore different data generation
# distributions and their relationship with regression calibrat+ion.


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
validation_data$foo <- unlist(lapply(validation_data$paee, data_generator, beta=set_beta))




# Data set of 20000 points for which we calculate the betas, (test)
test_data <- as.data.frame(c(rep(1,5000),rep(2,5000),rep(3,5000),rep(4,5000)))
colnames(test_data) <- c("cam_index")
test_data$cam_index <- as.factor(test_data$cam_index)
test_data$paee <- unlist(lapply(test_data$cam_index, gaussian_index_sample))
test_data$cam_index_means <- unlist(mapply(test_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
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
test_data$foo <- unlist(lapply(test_data$paee, data_generator, beta=set_beta))

####################################################################################################
################################ Attempt to recreate values ########################################
# now we work backworks to see if we can use the calculated foo value to find beta



## BASELINE with cam_index
# validation data
fc_val_data <- lm(formula=foo~cam_index, data=validation_data)
cc_val_data <- fc_val_data$coefficients[-1]
cc_val_data <- c(1,cc_val_data)
mean_cc <- mean(cc_val_data)
val_per_paee_cc <- mean_cc/mean_paee_difference

stdError_cc_val_data <- (summary(fc_val_data)$coefficients[,"Std. Error"])[-1]
mean_stdError_cc_val <- mean(stdError_cc_val_data)

# test data
fc_test_data <- lm(formula=foo~cam_index, data=test_data)
cc_test_data <- fc_test_data$coefficients[-1]
cc_test_data <- c(1,cc_test_data)
mean_cc <- mean(cc_test_data)
test_per_paee_cc <- mean_cc/mean_paee_difference

stdError_cc_test_data <- (summary(fc_test_data)$coefficients[,"Std. Error"])[-1]
mean_stdError_cc_test <- mean(stdError_cc_test_data)


## Replacing the cam_index with means
# validation
pc_val_data_means <- lm(formula=foo~cam_index_means, data=validation_data)
cc_val_data_means <- pc_val_data_means$coefficients["cam_index_means"]

stdError_cc_val_data_means <- (summary(pc_val_data_means)$coefficients[,"Std. Error"])[-1]
mean_stdError_cc_val_means <- mean(stdError_cc_val_data_means)

# test
pc_test_data_means <- lm(formula=foo~cam_index_means, data=test_data)
cc_test_data_means <- pc_test_data_means$coefficients["cam_index_means"]

stdError_cc_test_data_means <- (summary(pc_test_data_means)$coefficients[,"Std. Error"])[-1]
mean_stdError_cc_test_means <- mean(stdError_cc_test_data_means)


## Replacing the cam_index with number chosen for it at random and then monte carloing everything.

## Set lists of values for each PA categorization that will be used to draw random numbers for
#validation
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

test_per_paee_cc_list <- c(test_per_paee_cc)

for (i in 1:1000) {
	cat1_choice <- sample(cat1,1,replace=TRUE)
	cat2_choice <- sample(cat2,1,replace=TRUE)
	cat3_choice <- sample(cat3,1,replace=TRUE)
	cat4_choice <- sample(cat4,1,replace=TRUE)

	test_data$cam_index_means_sample <- unlist(mapply(test_data$cam_index, SIMPLIFY = FALSE, FUN=function(x){
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

	regression_fit <- lm(formula=foo~cam_index_means_sample, data=test_data)
	cc_test_data_sample <- regression_fit$coefficients["cam_index_means_sample"]

	lambda_men_var <- (summary(regression_fit)$coefficients["cam_index_means_sample","Std. Error"])^2
	test_per_paee_cc_list <- c(test_per_paee_cc_list, cc_test_data_sample)
}
mean_test_per_paee_cc <- mean(test_per_paee_cc_list)
