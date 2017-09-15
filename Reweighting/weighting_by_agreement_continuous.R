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

# Set save point for the plots
setwd("U:/workspace/paee_harmonisation/Reweighting/plots")

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

createStudyData <- function(raw_data, measurement_error){
    x_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
    raw_data$measured_x = raw_data$x + x_error
    return(raw_data)
}

createValidationData <- function(val_size, measurement_error) {
    validation_data = data.frame(x =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
    x_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
    validation_data$measured_x = validation_data$x + x_error
    return (validation_data)
}

#######################################################################################
########################### Motivation of the Problem #################################
#######################################################################################
## Graphing Standard Error as a function of Measurement Error
error_upperbound = 50
std_error = vector("numeric", length = error_upperbound)
estimates = vector("numeric", length = error_upperbound)

set.seed(1)
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
plot (x = (1:error_upperbound), y = std_error, xlab = "Measurement Error", ylab = "Standard Error", main = "Standard Error based on Measurement Error \n Before Changes")
dev.copy(png,'continous_plot1.png')
dev.off()
plot (x = (1:error_upperbound), y = estimates, xlab = "Measurement Error", ylab = "estimates", main = "Estimates (Accuracy based on Measurement Error) \n Before Changes")
dev.copy(png,'continous_plot2.png')
dev.off()
plot (x = estimates, y = std_error, xlab = "Estimates", ylab = "Standard Error", main = "Non Monotonic Relationship Between \n Accuracy and Standard Error \n Before Changes")
dev.copy(png,'continous_plot3.png')
dev.off()


## Random Effects Model Forest Plot After Regression Calibration
estimates_forRMA = estimates
stand_errs = std_error
labels = 1:length(estimates)
res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~harmonised_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Before Any Changes")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'continous_plot_big_no_changes.png')
dev.off()

#######################################################################################
######################## Before Any Changes to the Process ############################
#######################################################################################
# General Parameters
validation_size = 400

# Study A
measureError_A = 5
studyData_A = createStudyData(raw_data = raw_data, measurement_error = measureError_A)
validation_data_A = createValidationData(val_size = validation_size, measurement_error = measureError_A )

# Study B
measureError_B = 10
studyData_B = createStudyData(raw_data = raw_data, measurement_error = measureError_B)
validation_data_B = createValidationData(val_size = validation_size, measurement_error = measureError_B )

# Study C
measureError_C = 200
studyData_C = createStudyData(raw_data = raw_data, measurement_error = measureError_C)
validation_data_C = createValidationData(val_size = validation_size, measurement_error = measureError_C )


lm_A <- lm(formula=y~measured_x, data=studyData_A)
estimate_A = lm_A$coefficients["measured_x"]
stdError_A = summary(lm_A)$coefficients["measured_x","Std. Error"]

lm_B <- lm(formula=y~measured_x, data=studyData_B)
estimate_B = lm_B$coefficients["measured_x"]
stdError_B = summary(lm_B)$coefficients["measured_x","Std. Error"]

lm_C <- lm(formula=y~measured_x, data=studyData_C)
estimate_C = lm_C$coefficients["measured_x"]
stdError_C = summary(lm_C)$coefficients["measured_x","Std. Error"]


## Random Effects Model Forest Plot Before Regression Calibration
estimates = c(estimate_A, estimate_B, estimate_C)
stand_errs = c(stdError_A, stdError_B, stdError_C)
labels = c("A","B","C")
res <- rma(yi = estimates, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~measured_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Before Any Changes to Process")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'continous_plot4.png')
dev.off()


#######################################################################################
######################## Using Regression Calibration Classic #########################
#######################################################################################
## Graphing Standard Error as a function of Measurement Error
std_error = vector("numeric", length = error_upperbound)
estimates = vector("numeric", length = error_upperbound)

set.seed(1)
for (measurement_error in 1:error_upperbound){
    studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measurement_error)
    validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measurement_error)

    # calculate the lambda
    fmla_harmonisation = as.formula(x~measured_x)
    lambda_lm_graph <- lm(formula=fmla_harmonisation, data=validation_data_graph)
    lambda_graph = lambda_lm_graph$coefficients["measured_x"]
    stdError_lambda_graph = summary(lambda_lm_graph)$coefficients["measured_x","Std. Error"]

    # calculate the estimates using the normal study data
    lm_graph <- lm(formula=y~measured_x, data=studyData_graph)
    estimate_graph = lm_graph$coefficients["measured_x"]
    stdError_graph = summary(lm_graph)$coefficients["measured_x","Std. Error"]

    # Propagate Errors in Standard Error
    var_lambda = (sqrt(validation_size) * summary(lambda_lm_graph)$coefficients["measured_x","Std. Error"])^2
    var_Beta = (sqrt(validation_size) * summary(lm_graph)$coefficients["measured_x","Std. Error"])^2
    #var_Beta = summary(lm_A)$sigma**2
    # var_Beta = (summary(lm_A)$coefficients["measured_x","Std. Error"])^2
    lambda_pure = unlist(unname(lambda_lm_graph$coefficients["measured_x"]))
    beta_lambda_div_sq = (unname(unlist(lm_graph$coefficients["measured_x"]/(lambda_lm_graph$coefficients["measured_x"])^2)))^2
    delta_variance = (var_Beta / (lambda_pure)^2) + beta_lambda_div_sq * var_lambda
    delta_stdError_B = sqrt(delta_variance)/sqrt(validation_size) 
    delta_stdError_B = sqrt(delta_variance)

    std_error[measurement_error] = delta_stdError_B
    estimates[measurement_error] = estimate_graph/lambda_graph
}
plot (x = (1:error_upperbound), y = std_error, xlab = "Measurement Error", ylab = "Standard Error", main = "Standard Error based on Measurement Error \n Regression Calibration")
dev.copy(png,'continous_plot5.png')
dev.off()
plot (x = (1:error_upperbound), y = estimates, xlab = "Measurement Error", ylab = "estimates", main = "Estimates (Accuracy based on Measurement Error) \n Regression Calibration")
dev.copy(png,'continous_plot6.png')
dev.off()
plot (x = estimates, y = std_error, xlab = "Estimates", ylab = "Standard Error", main = "Non Monotonic Relationship Between \n Accuracy and Standard Error \n Regression Calibration")
dev.copy(png,'continous_plot7.png')
dev.off()

## Random Effects Model Forest Plot After Regression Calibration
estimates_forRMA = estimates
stand_errs = std_error
labels = 1:length(estimates)
res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~harmonised_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "After Regression Calibration")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'continous_plot_big_regression_calibration.png')
dev.off()




# Using real regression calibration.
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


lm_A <- lm(formula=y~measured_x, data=studyData_A)
estimate_A = lm_A$coefficients["measured_x"]
stdError_A = summary(lm_A)$coefficients["measured_x","Std. Error"]

lm_B <- lm(formula=y~measured_x, data=studyData_B)
estimate_B = lm_B$coefficients["measured_x"]
stdError_B = summary(lm_B)$coefficients["measured_x","Std. Error"]

lm_C <- lm(formula=y~measured_x, data=studyData_C)
estimate_C = lm_C$coefficients["measured_x"]
stdError_C = summary(lm_C)$coefficients["measured_x","Std. Error"]

corrected_estimate_A = estimate_A / lambda_A
corrected_estimate_B = estimate_B / lambda_B
corrected_estimate_C = estimate_C / lambda_C

# calculate the corrected standard error
var_lambda = (summary(lambda_lm_A)$coefficients["measured_x","Std. Error"])^2
var_Beta = (summary(lm_A)$coefficients["measured_x","Std. Error"])^2
#var_Beta = summary(lm_A)$sigma**2
lambda_pure = unlist(unname(lambda_lm_A$coefficients["measured_x"]))
beta_lambda_div_sq = (unname(unlist(lm_A$coefficients["measured_x"]/(lambda_lm_A$coefficients["measured_x"])^2)))^2
delta_variance = (var_Beta / (lambda_pure)^2) + (beta_lambda_div_sq * var_lambda)
delta_stdError_A = sqrt(delta_variance)/sqrt(validation_size) 
delta_stdError_A = sqrt(delta_variance)

var_lambda = (summary(lambda_lm_B)$coefficients["measured_x","Std. Error"])^2
var_Beta = (sqrt(validation_size) * summary(lm_B)$coefficients["measured_x","Std. Error"])^2
#var_Beta = summary(lm_A)$sigma**2
lambda_pure = unlist(unname(lambda_lm_B$coefficients["measured_x"]))
beta_lambda_div_sq = (unname(unlist(lm_B$coefficients["measured_x"]/(lambda_lm_B$coefficients["measured_x"])^2)))^2
delta_variance = (var_Beta / (lambda_pure)^2) + beta_lambda_div_sq * var_lambda
delta_stdError_B = sqrt(delta_variance)/sqrt(validation_size) 
delta_stdError_B = sqrt(delta_variance)


var_lambda = (summary(lambda_lm_C)$coefficients["measured_x","Std. Error"])^2
var_Beta = (sqrt(validation_size) * summary(lm_C)$coefficients["measured_x","Std. Error"])^2
#var_Beta = summary(lm_A)$sigma**2
lambda_pure = unlist(unname(lambda_lm_C$coefficients["measured_x"]))
beta_lambda_div_sq = (unname(unlist(lm_C$coefficients["measured_x"]/(lambda_lm_C$coefficients["measured_x"])^2)))^2
delta_variance = (var_Beta / (lambda_pure)^2) + beta_lambda_div_sq * var_lambda
delta_stdError_C = sqrt(delta_variance)/sqrt(validation_size) 
delta_stdError_C = sqrt(delta_variance)

## Random Effects Model Forest Plot After Regression Calibration
estimates = c(corrected_estimate_A, corrected_estimate_B, corrected_estimate_C)
stand_errs = c(delta_stdError_A, delta_stdError_B, delta_stdError_C)
labels = c("A","B","C")
res <- rma(yi = estimates, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~harmonised_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "After Regression Calibration")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'continous_plot8.png')
dev.off()


#######################################################################################
######################## USING THE MODEL OF ERROR *PREDICT* ###########################
#######################################################################################
## PREDICT (ERROR IN VARIABLES MODELS)

## Graphing Standard Error as a function of Measurement Error
std_error = vector("numeric", length = error_upperbound)
estimates = vector("numeric", length = error_upperbound)

set.seed(1)

for (measurement_error in 1:error_upperbound){
    studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measurement_error)
    validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measurement_error)

    # calculate the Error Model
    error_model_graph <- lm(formula=fmla_harmonisation, data=validation_data_graph)
    lambda_graph = error_model_graph$coefficients["measured_x"]
    stdError_lambda_graph = summary(error_model_graph)$coefficients["measured_x","Std. Error"]

    # calculate the estimates using the normal study data
    new = data.frame(measured_x = studyData_graph$measured_x)
    studyData_graph$harmonised_x = predict(error_model_graph, newdata = new)
    lm_graph <- lm(formula=y~harmonised_x, data=studyData_graph)
    estimate_graph = lm_graph$coefficients["harmonised_x"]
    stdError_graph = summary(lm_graph)$coefficients["harmonised_x","Std. Error"]


    std_error[measurement_error] = stdError_graph
    estimates[measurement_error] = estimate_graph
}
plot (x = (1:error_upperbound), y = std_error, xlab = "Measurement Error", ylab = "Standard Error", main = "Standard Error based on Measurement Error \n Error Model Pred")
dev.copy(png,'continous_plot9.png')
dev.off()
plot (x = (1:error_upperbound), y = estimates, xlab = "Measurement Error", ylab = "estimates", main = "Estimates (Accuracy based on Measurement Error) \n Error Model Pred")
dev.copy(png,'continous_plot10.png')
dev.off()
plot (x = estimates, y = std_error, xlab = "Estimates", ylab = "Standard Error", main = "Monotonic Relationship Between \n Accuracy and Standard Error \n Error Model Pred")
dev.copy(png,'continous_plot11.png')
dev.off()

## Random Effects Model Forest Plot After Regression Calibration
estimates_forRMA = estimates
stand_errs = std_error
labels = 1:length(estimates)
res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~harmonised_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Error Model Predict Big Regression")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'continous_plot_big_error_model.png')
dev.off()


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


## Transformation and Regression
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


## Random Effects Model Forest Plot After Regression Modeling
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
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Using Predict aka 'Error Model'")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'continous_plot12.png')
dev.off()


# #######################################################################################
# ######################## USING THE MODEL OF ERROR *PREDICT* ###########################
# ################## DONE MANUALLY WITH MODEL PARAMETER ESTIMATES #######################
# #######################################################################################
# ## PREDICT (ERROR IN VARIABLES MODELS)

# ## Graphing Standard Error as a function of Measurement Error
# std_error = vector("numeric", length = error_upperbound)
# estimates = vector("numeric", length = error_upperbound)

# set.seed(1)

# for (measurement_error in 1:error_upperbound){
#     studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measurement_error)
#     validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measurement_error)

#     # calculate the Error Model
#     error_model_graph <- lm(formula=fmla_harmonisation, data=validation_data_graph)
#     lambda_graph = error_model_graph$coefficients["measured_x"]
#     intercept_lambda_graph = error_model_graph$coefficients["(Intercept)"]

#     # calculate the estimates using the normal study data
#     studyData_graph$lambdadx = studyData_graph$measured_x*lambda_graph + intercept_lambda_graph
#     lm_graph <- lm(formula=y~lambdadx, data=studyData_graph)
#     estimate_graph = lm_graph$coefficients["lambdadx"]
#     stdError_graph = summary(lm_graph)$coefficients["lambdadx","Std. Error"]

#     estimates[measurement_error] = estimate_graph
#     std_error[measurement_error] = stdError_graph
# }
# plot (x = (1:error_upperbound), y = std_error, xlab = "Measurement Error", ylab = "Standard Error", main = "Standard Error based on Measurement Error \n Error Model Pred")
# dev.copy(png,'continous_plot13.png')
# dev.off()

# plot (x = (1:error_upperbound), y = estimates, xlab = "Measurement Error", ylab = "estimates", main = "Estimates (Accuracy based on Measurement Error) \n Error Model Pred")
# dev.copy(png,'continous_plot14.png')
# dev.off()

# plot (x = estimates, y = std_error, xlab = "Estimates", ylab = "Standard Error", main = "Monotonic Relationship Between \n Accuracy and Standard Error \n Error Model Pred")
# dev.copy(png,'continous_plot15.png')
# dev.off()


# ## By manually predicting values using the model created
# ## Modeling x ~ measured_x in validation
# fmla_harmonisation = as.formula(x~measured_x)
# # Study A
# lambda_lm_A <- lm(formula=fmla_harmonisation, data=validation_data_A)
# lambda_A = lambda_lm_A$coefficients["measured_x"]
# intercept_lambda_A = lambda_lm_A$coefficients["(Intercept)"]
# stdError_lambda_A = summary(lambda_lm_A)$coefficients["measured_x","Std. Error"]

# # Study B
# lambda_lm_B <- lm(formula=fmla_harmonisation, data=validation_data_B)
# lambda_B = lambda_lm_B$coefficients["measured_x"]
# intercept_lambda_B = lambda_lm_B$coefficients["(Intercept)"]
# stdError_lambda_B = summary(lambda_lm_B)$coefficients["measured_x","Std. Error"]

# # Study C
# lambda_lm_C <- lm(formula=fmla_harmonisation, data=validation_data_C)
# lambda_C = lambda_lm_C$coefficients["measured_x"]
# intercept_lambda_C = lambda_lm_C$coefficients["(Intercept)"]
# stdError_lambda_C = summary(lambda_lm_C)$coefficients["measured_x","Std. Error"]

# studyData_A$lambdadx = studyData_A$measured_x*lambda_A + intercept_lambda_A
# lm_A <- lm(formula=y~lambdadx, data=studyData_A)
# estimate_A = lm_A$coefficients["lambdadx"]
# stdError_A = summary(lm_A)$coefficients["lambdadx","Std. Error"]

# studyData_B$lambdadx = studyData_B$measured_x*lambda_B + intercept_lambda_B
# lm_B <- lm(formula=y~lambdadx, data=studyData_B)
# estimate_B = lm_B$coefficients["lambdadx"]
# stdError_B = summary(lm_B)$coefficients["lambdadx","Std. Error"]

# studyData_C$lambdadx = studyData_C$measured_x*lambda_C + intercept_lambda_C
# lm_C <- lm(formula=y~lambdadx, data=studyData_C)
# estimate_C = lm_C$coefficients["lambdadx"]
# stdError_C = summary(lm_C)$coefficients["lambdadx","Std. Error"]


# ## Random Effects Model Forest Plot After Regression Calibration in all members of measured x
# estimates = c(estimate_A, estimate_B, estimate_C)
# stand_errs = c(stdError_A, stdError_B, stdError_C)
# labels = c("A","B","C")
# res <- rma(yi = estimates, sei = stand_errs, method='DL', slab = labels)
# weights_res <- weights.rma.uni(res)

# # Forest Plot
# res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
# fmla = as.formula(y~harmonised_x)
# forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
#     .(sprintf("%.3f", round(res$QEp,3))),')')),
# xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0.5, p = ',
#     .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Using 'Error Model'")
# usr <- par("usr")
# text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
# text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
# abline(v = 0.5, col = "lightgray")
# dev.copy(png,'continous_plot16.png')
# dev.off()
