# Investigation into Data Agreement as Weighting System
## Author: Paul Scherer
##         Tom Bishop
## Institute: University of Cambridge MRC Epidemiology Unit 
## Date: 4.09.2017

# Set up the Data
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
set.seed(666)
study_size = 20000
raw_mean = 40
raw_stdev = 15

# Simple creation of the raw data and the related outcome without error
raw_data = data.frame(exp_true = rnorm(n = study_size, mean = raw_mean, sd = raw_stdev))
raw_data$outcome = (trueBeta*raw_data$exp_true) + constant


###############################################################################
########################### Functions #########################################
###############################################################################
createStudyData <- function(raw_data, measurement_error, number_of_indices){
  exposure_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
  raw_data$exposure_measured = raw_data$exp_true + exposure_error
  raw_data$index = cut_number(x = raw_data$exposure_measured,n = number_of_indices, labels =FALSE)
  raw_data$index = as.factor(raw_data$index)
  raw_data = raw_data[with(raw_data, order(index)),]
  return(raw_data)
}

createValidationData <- function(val_size, measurement_error, number_of_indices) {
  validation_data = data.frame(exp_true =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
  exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
  validation_data$exposure_measured = validation_data$exp_true + exposure_error
  validation_data$index = cut_number(x = validation_data$exposure_measured,n = number_of_indices, labels =FALSE)
  validation_data$index = as.factor(validation_data$index)
  return (validation_data)
}

createMeansList <- function(data, number_of_indices) {
  meansList = vector(mode="list", length = number_of_indices)
  for (i in 1:number_of_indices) {
    meansList[i] = mean(unname(unlist((split(x=data$exp_true, f= as.factor(data$index)))[i])))
  }
  return(meansList)
}




#######################################################################################
######################## Establish a baseline #####################
#######################################################################################

numbaselines = 10
numLevels = 4
validation_size = 400
upperbound = 30
lowerbound = 10
set.seed(666)

per_error_weights = list()
per_error_est = list()
per_error_se = list()

for (i in 1:numbaselines){
  
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    # Find Error Model
    mapping_formula = as.formula(exp_true ~ exposure_measured)
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    new = data.frame(exposure_measured = studyData$exposure_measured)
    studyData$exposure_corrected = predict(mapping_model, newdata = new)
    
    # Use Error Model to Scale Index to Gold
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    stdError = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    #assume lambda is 1
    lambda = summary(mapping_model)$coefficients["exposure_measured","Estimate"]
    lambda_stdError = summary(mapping_model)$coefficients["exposure_measured","Std. Error"]
    
    # recalibrate the standard error using delta function
    
    variance_beta = stdError^2
    beta_lambda_div_sq = (unname(unlist(estimate/(lambda)^2)))^2
    var_lambda = (lambda_stdError)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_stdError_graph = sqrt(delta_variance)
    
    estimates[measurement_error_counter-lowerbound+1] = estimate
    std_errors[measurement_error_counter-lowerbound+1] = delta_stdError_graph
  }
  
  
  ## Plot meta analysis for many studies
  
  estimates_forRMA = estimates[seq(1, length(estimates), 1)]
  stand_errs = std_errors[seq(1, length(std_errors), 1)]
  
  labels = c(lowerbound:upperbound)
  #labels = 1:length(estimates_forRMA)
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  
}

#Example meta analysis

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~harmonised_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                              .(sprintf("%.3f", round(res$QEp,3))),')')),
       xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                         .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Example baseline")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'baseline 10 30 example.png')
dev.off()

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)

mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='Baseline - average over 10 studies')
dev.copy(png,'baseline average weights.png')
dev.off()
#mean_est = colMeans(per_error_est_final)
#plot(labels, mean_est)
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='Baseline - average over 10 studies')
dev.copy(png,'baseline average se.png')
dev.off()
