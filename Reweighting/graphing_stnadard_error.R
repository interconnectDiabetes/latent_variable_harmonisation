# Investigation into semi monotonic behaviour of standard error
## Author: Paul Scherer
##         Tom Bishop
## Institute: University of Cambridge MRC Epidemiology Unit 
## Date: 6.09.2017

# I have had my results for a long time; but I do not yet know how I am to arrive at them - Ma boi Gauss

# In this file we want to investigate whether we can replicate the semi monotonic behaviour exhibited by 
# standard error dependent on standard measurement error

library(ggplot2)
library(reshape2)
library(parallel)
library(stats)
library(metafor)

# Raw Data Set Properties
trueBeta = 0.5
constant = 20

study_size = 20000
raw_mean = 40
raw_stdev = 15

# Create Dataframe with x, create y from x, then add normal random noise to x
raw_data = data.frame(x = rnorm(n = study_size, mean = raw_mean, sd = raw_stdev))
raw_data$y = (trueBeta*raw_data$x) + constant

# Add measurement error to x and then calculate the lm between x with error and y to get the standard error
error_upperbound = 1000
std_error = vector("numeric", length = error_upperbound)
estimates = vector("numeric", length = error_upperbound)

for (measurement_error in 1:error_upperbound){
	raw_copy = raw_data
	x_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
	raw_copy$measured_x = raw_copy$x + x_error
	fmla = as.formula(y ~ measured_x)
	linear_model = lm(formula = fmla, data = raw_copy)
	x_estimate = linear_model$coefficients["measured_x"]
	x_stdErr = summary(linear_model)$coefficients["measured_x", "Std. Error"]
	std_error[measurement_error] = x_stdErr
	estimates[measurement_error] = x_estimate
}

plot (x = (1:error_upperbound), y = std_error, xlab = "Measurement Error", ylab = "Standard Error")
