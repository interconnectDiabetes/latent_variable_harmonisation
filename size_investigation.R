# Jan 30 meeting preparation
# Simulated data run 1000 runs, each individual recieved a raodonm value drawin 
# from a distribution defined for their PA category.

# 1000 runs, each individual received a random value drawn from a distribution defined for their PA category. So Stephen is proposing:

# 1. Generate data for a "large"	 validation study (e.g. N=1000, 250 individuals in each pa_index category).   
# The PAEE for individuals with pa_index=1 are values sampled from a Normal distribution, mean 30, SD 2.5, 
# i.e. N(30,2.5); for pa_index=2, sample from N(40,2.5); for pa_index=3, sample from N(50,2.5); for pa_index=4, 
# sample from N(60,2.5).  This way, the sampled distributions should be almost entirely non-overlapping.

# 2. For each individual, randomly select a value from the appropriate artificial validation study distribution
#  and estimate the PAEE/T2D association using the real InterAct data. Repeat this 1000 (or some other large 
#  number) times.

# 3. Repeat (1) and (2) for a "small" validation study (e.g. N=100, 25 in each pa_index category), sampling 
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
index_mean1 <- 30
index_stdev1 <- 2.5
index_mean2 <- 40
index_stdev2 <- 2.5
index_mean3 <- 50
index_stdev3 <- 2.5
index_mean4 <- 60
index_stdev4 <- 2.5

# PAEE range
paee_range_min <- 20
paee_range_max <- 70

###############################################################################
############################# Functions #######################################
###############################################################################
data_generator <- function(exposure, beta=set_beta, constant=20, noise=2.5){
  # Gives the y in y = mx+b given the other parameters 
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

# 'Large' Validation Dataset 
large_val <- as.data.frame(c(rep(1,250),rep(2,250),rep(3,250),rep(4,250)))
colnames(large_val) <- c("cam_index")
large_val$cam_index <- as.factor(large_val$cam_index)
large_val$paee <- unlist(lapply(large_val$cam_index, gaussian_index_sample))
large_val$cam_index_means <- unlist(mapply(large_val$cam_index, SIMPLIFY = FALSE, FUN=function(x){
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

# 'Small' Validation Dataset 
small_val <- as.data.frame(c(rep(1,25),rep(2,25),rep(3,25),rep(4,25)))
colnames(small_val) <- c("cam_index")
small_val$cam_index <- as.factor(small_val$cam_index)
small_val$paee <- unlist(lapply(small_val$cam_index, gaussian_index_sample))
small_val$cam_index_means <- unlist(mapply(small_val$cam_index, SIMPLIFY = FALSE, FUN=function(x){
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


# Study data set of 20000 points for which we calculate the betas, (test)
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

# Create lists to draw random samples from for predictions this is for the small validation set
s_cat1 <- mapply(small_val$paee, small_val$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
s_cat1 <- s_cat1[!sapply(s_cat1,is.na)]

s_cat2 <- mapply(small_val$paee, small_val$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
s_cat2 <- s_cat2[!sapply(s_cat2,is.na)]

s_cat3 <- mapply(small_val$paee, small_val$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
s_cat3 <- s_cat3[!sapply(s_cat3,is.na)]

s_cat4 <- mapply(small_val$paee, small_val$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
})
s_cat4 <- s_cat4[!sapply(s_cat4,is.na)]


# similarly for the large validation set
cat1 <- mapply(large_val$paee, large_val$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat1 <- cat1[!sapply(cat1,is.na)]

cat2 <- mapply(large_val$paee, large_val$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat2 <- cat2[!sapply(cat2,is.na)]

cat3 <- mapply(large_val$paee, large_val$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat3 <- cat3[!sapply(cat3,is.na)]

cat4 <- mapply(large_val$paee, large_val$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
})
cat4 <- cat4[!sapply(cat4,is.na)]


