# Jan 30 meeting preparation
# Simulated data run 1000 runs, each individual recieved a raodonm value drawin 
# from a distribution defined for their PA category.

# 1. Generate data for a “large” validation study (e.g. N=1000, 250 individuals in each pa_index category).   
# The PAEE for individuals with pa_index=1 are values sampled from a Normal distribution, mean 30, SD 2.5, 
# i.e. N(30,2.5); for pa_index=2, sample from N(40,2.5); for pa_index=3, sample from N(50,2.5); for pa_index=4, 
# sample from N(60,2.5).  This way, the sampled distributions should be almost entirely non-overlapping.

# 2. For each individual, randomly select a value from the appropriate artificial validation study distribution
#  and estimate the PAEE/T2D association using the real InterAct data. Repeat this 1000 (or some other large 
#  number) times.

# 3. Repeat (1) and (2) for a “small” validation study (e.g. N=100, 25 in each pa_index category), sampling 
# values from the same 4 underlying distributions.

## Author: Paul Scherer
## Date: 13.01.2017

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