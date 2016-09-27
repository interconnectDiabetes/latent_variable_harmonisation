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

## Seed
set.seed(66)

## Parameters
# The linear correlation coefficient between exposure and outcome.
set_beta <- 0.5 

# Index properties (Assumption that they are Gaussian)
# Format : (mean, stdev)
index1 <- c(10,5)
index2 <- c(20,5)
index3 <- c(30,5)
index4 <- c(40,5)

###############################################################################
############################# Functions #######################################
###############################################################################
data_generator <- function(exposure, set_beta){
	# returns a datapoint given an exposure with some gaussian noise
	return (data_point)
}

gaussian_index_sample <- function(index){
	# returns a datapoint sampled from the index distribution
	# :param: index = cambridge index
	index_mean <- index[1]
	index_stdev <- index[2]
	data_point <- rnorm(1, index_mean, index_stdev)

	return (data_point)
}




