# Investigation into semi monotonic behaviour of standard error in Continuous Variables
## Author: Paul Scherer
##         Tom Bishop
## Institute: University of Cambridge MRC Epidemiology Unit 
## Date: 6.09.2017

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

###############################################################################
########################### Functions #########################################
###############################################################################

createStudyData <- function(raw_data, measurement_error, number_of_indices){
    x_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
    raw_data$measured_x = raw_data$x + x_error
    return(raw_data)
}

createValidationData <- function(val_size, measurement_error, number_of_indices) {
    validation_data = data.frame(x =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
    x_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
    validation_data$measured_x = validation_data$x + x_error
    return (validation_data)
}

#######################################################################################
########################### Motivation of the Problem #################################
#######################################################################################
## Graphing Standard Error as a function of Measurement Error
error_upperbound = 100
std_error = vector("numeric", length = error_upperbound)
estimates = vector("numeric", length = error_upperbound)

for (measurement_error in 1:error_upperbound){
	x_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
	raw_data$measured_x = raw_data$x + x_error
	fmla = as.formula(y ~ measured_x)
	linear_model = lm(formula = fmla, data = raw_data)
	x_estimate = linear_model$coefficients["measured_x"]
	x_stdErr = summary(linear_model)$coefficients["measured_x", "Std. Error"]
	std_error[measurement_error] = x_stdErr
	estimates[measurement_error] = x_estimate
}
plot (x = (1:error_upperbound), y = std_error, xlab = "Measurement Error", ylab = "Standard Error")
plot (x = (1:error_upperbound), y = estimates, xlab = "Measurement Error", ylab = "estimates")

plot (x = estimates, y = std_error, xlab = "Estimates", ylab = "Standard Error")

#######################################################################################
######################## Exhibiting the Effects on Studies ############################
#######################################################################################
# General Parameters
numLevels = 4
validation_size = 400

# Study A
measureError_A = 5
studyData_A = createStudyData(raw_data = raw_data, measurement_error = measureError_A, number_of_indices = numLevels)
validation_data_A = createValidationData(val_size = validation_size, measurement_error = measureError_A, number_of_indices = numLevels )

# Study B
measureError_B = 10
studyData_B = createStudyData(raw_data = raw_data, measurement_error = measureError_B, number_of_indices = numLevels)
validation_data_B = createValidationData(val_size = validation_size, measurement_error = measureError_B, number_of_indices = numLevels )

# Study C
measureError_C = 60
studyData_C = createStudyData(raw_data = raw_data, measurement_error = measureError_C, number_of_indices = numLevels)
validation_data_C = createValidationData(val_size = validation_size, measurement_error = measureError_C, number_of_indices = numLevels )

## Modeling x ~ measured_x in validation
fmla_harmonisation = as.formula(x~measured_x)
# Study A
lambda_lm_A <- lm(formula=fmla_harmonisation, data=validation_data_A)
lambda_A = lambda_lm_A$coefficients["measured_x"]
stdError_lambda_A = summary(lambda_lm_A)$coefficients["measured_x","Std. Error"]

# Study B
lambda_lm_B <- lm(formula=fmla_harmonisation, data=validation_data_B)
lambda_B = lambda_lm_B$coefficients["measured_x"]
stdError_lambda_B = summary(lambda_lm_B)$coefficients["measured_x","Std. Error"]

# Study C
lambda_lm_C <- lm(formula=fmla_harmonisation, data=validation_data_C)
lambda_C = lambda_lm_C$coefficients["measured_x"]
stdError_lambda_C = summary(lambda_lm_C)$coefficients["measured_x","Std. Error"]


## Transformation and Regression using Index Means
# Study A
new = data.frame(measured_x = studyData_A$measured_x)
studyData_A$harmonised_x = predict(lambda_lm_A, newdata = new)
lm_A <- lm(formula=y~harmonised_x, data=studyData_A)
estimate_A = lm_A$coefficients["harmonised_x"]
stdError_A = summary(lm_A)$coefficients["harmonised_x","Std. Error"]

# Study B
new = data.frame(measured_x = studyData_B$measured_x)
studyData_B$harmonised_x = predict(lambda_lm_B, newdata = new)
lm_B <- lm(formula=y~harmonised_x, data=studyData_B)
estimate_B = lm_B$coefficients["harmonised_x"]
stdError_B = summary(lm_B)$coefficients["harmonised_x","Std. Error"]

# Study C
new = data.frame(measured_x = studyData_C$measured_x)
studyData_C$harmonised_x = predict(lambda_lm_C, newdata = new)
lm_C <- lm(formula=y~harmonised_x, data=studyData_C)
estimate_C = lm_C$coefficients["harmonised_x"]
stdError_C = summary(lm_C)$coefficients["harmonised_x","Std. Error"]


## Random Effects Model Forest Plot Before Reweighting
estimates = c(estimate_A, estimate_B, estimate_C)
stand_errs = c(stdError_A, stdError_B, stdError_C)
labels = c("A","B","C")
res <- rma(yi = estimates, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~harmonised_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0.5, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1)
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
