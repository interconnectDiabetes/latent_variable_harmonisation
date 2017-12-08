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
  validation_data = validation_data[order(validation_data$index),]
  return (validation_data)
}

createValidationDataRep <- function(val_size, measurement_error, number_of_indices) {
  validation_data = data.frame(exp_true =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
  exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
  exposure_error2 = rnorm(n = val_size, mean = 0, sd = measurement_error)
  validation_data$exposure_measured = validation_data$exp_true + exposure_error
  validation_data$exposure_measured2 = validation_data$exp_true + exposure_error2
  validation_data$index = cut_number(x = validation_data$exposure_measured,n = number_of_indices, labels =FALSE)
  validation_data$index = as.factor(validation_data$index)
  validation_data$index2 = cut_number(x = validation_data$exposure_measured2,n = number_of_indices, labels =FALSE)
  validation_data$index2 = as.factor(validation_data$index2)
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
######################## Establish a baseline - no RC #####################
#######################################################################################


numLevels = 4
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 25
#set.seed(666)

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'baseline no RC multi.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

for (i in 1:numbaselines){
  
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    
    # Raw model, no correction
    corrected_model <- lm(formula=outcome~exposure_measured, data=studyData)
    estimate = corrected_model$coefficients["exposure_measured"]
    stdError = summary(corrected_model)$coefficients["exposure_measured","Std. Error"]
    

    
    estimates[measurement_error_counter-lowerbound+1] = estimate
    std_errors[measurement_error_counter-lowerbound+1] = stdError
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
  
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "baseline")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()
par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)

mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='Baseline - average over 25 studies')
dev.copy(png,'baseline no RC average weights.png')
dev.off()
mean_est = colMeans(per_error_est_final)
plot(labels, mean_est, xlab='Amount of error in instrument', ylab='Mean est given to study', main='Baseline - average over 25 studies')
dev.copy(png,'baseline no RC average est.png')
dev.off()
mean_se = colMeans(per_error_se_final)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='Baseline - average over 25 studies')
dev.copy(png,'baseline no RC average se.png')
dev.off()


#######################################################################################
######################## Establish a baseline #####################
#######################################################################################


numLevels = 4
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 25
#set.seed(666)

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'baseline multi.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

for (i in 1:numbaselines){
  
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    
    # Find Error Model
    mapping_formula = as.formula(exp_true ~ exposure_measured)
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    new = data.frame(exposure_measured = studyData$exposure_measured)
    studyData$exposure_corrected = predict(mapping_model, newdata = new)
    
    #get lambda
    lambda = summary(mapping_model)$coefficients["exposure_measured","Estimate"]
    lambda_stdError = summary(mapping_model)$coefficients["exposure_measured","Std. Error"]
    
    # Use Error Model to correct and estimate beta
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate_corrected = corrected_model$coefficients["exposure_corrected"]
    stdError_corrected = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    #raw model
    raw_model <- lm(formula=outcome~exposure_measured, data=studyData)
    estimate_uncorrected = raw_model$coefficients["exposure_measured"]
    stdError_uncorrected = summary(raw_model)$coefficients["exposure_measured","Std. Error"]
    
    #get lambda
    lambda = summary(mapping_model)$coefficients["exposure_measured","Estimate"]
    lambda_stdError = summary(mapping_model)$coefficients["exposure_measured","Std. Error"]
    
    # recalibrate the standard error using delta function
    
    variance_beta = stdError_uncorrected^2
    beta_lambda_div_sq = (unname(unlist(estimate_uncorrected/(lambda)^2)))^2
    var_lambda = (lambda_stdError)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_stdError_graph = sqrt(delta_variance)
    
    
    estimates[measurement_error_counter-lowerbound+1] = estimate_corrected
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
  
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "baseline")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()
par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)

mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='Baseline - average over 25 studies')
dev.copy(png,'baseline average weights.png')
dev.off()
mean_est = colMeans(per_error_est_final)
plot(labels, mean_est, xlab='Amount of error in instrument', ylab='Mean est given to study', main='Baseline - average over 25 studies')
dev.copy(png,'baseline average est.png')
dev.off()
mean_se = colMeans(per_error_se_final)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='Baseline - average over 25 studies')
dev.copy(png,'baseline average se.png')
dev.off()

#######################################################################################
######################## Sub in means only - no RC #####################
#######################################################################################


numLevels = 4
validation_size = 400
lowerbound = 10
upperbound = 30
set.seed(666)
numbaselines = 25
## Plotting Bit
## Graphing Standard Error as a function of Measurement Error
per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'sub in means only.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))


for (i in 1:numbaselines){
  
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    # Find Error Model
    mapping_formula = as.formula(exp_true ~ index + 0)
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    new = data.frame(index = studyData$index)
    studyData$index_mapped = predict(mapping_model, newdata = new)
    
    # Use Error Model to Scale Index to Gold
    corrected_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate = corrected_model$coefficients["index_mapped"]
    stdError = summary(corrected_model)$coefficients["index_mapped","Std. Error"]
    
    estimates[measurement_error_counter-lowerbound+1] = estimate
    std_errors[measurement_error_counter-lowerbound+1] = stdError
  }
  
  
  ## Plot meta analysis for many studies
  estimates_forRMA = estimates[seq(1, length(estimates), 1)]
  stand_errs = std_errors[seq(1, length(std_errors), 1)]
  
  labels = c(lowerbound:upperbound)
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "sub means only")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
  
}


dev.off()
par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)


mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='Map to means no RC - average over 25 studies')
dev.copy(png,'sub in means no RC average weights.png')
dev.off()
mean_est = colMeans(per_error_est_final)
plot(x=labels, y=mean_est,  xlab='Amount of error in instrument', ylab='Mean estimate for study', main='Map to means no RC - average over 25 studies')
dev.copy(png,'sub in means no RC average est.png')
dev.off()
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='Map to means no RC - average over 25 studies')
dev.copy(png,'sub in means no RC average se.png')
dev.off()


#######################################################################################
######################## Map to means with RC #####################
#######################################################################################

# Concern is that there is no estimate of error in mapping, only at the calibration stage

numLevels = 4
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 25

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'sub in means with RC.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))


for (i in 1:numbaselines){
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    # Find mapping Model
    mapping_formula = as.formula(exp_true ~ index + 0)
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    new_study = data.frame(index = studyData$index)
    new_validation = data.frame(index = validation_data$index)
    validation_data$index_mapped = predict(mapping_model, newdata = new_validation)
    studyData$index_mapped = predict(mapping_model, newdata = new_study)
    
    # Use mapped values in model - estimate beta
    corrected_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate = corrected_model$coefficients["index_mapped"]
    std_error = summary(corrected_model)$coefficients[2,2]
    
    #RC model
    
    RC_model = lm(formula=exp_true~index_mapped, data=validation_data)
    lambda = RC_model$coefficients["index_mapped"]
    lambda_std_error = summary(RC_model)$coefficients[2,2]
    
    # recalibrate the standard error using delta function
    
    variance_beta = std_error^2
    beta_lambda_div_sq = (unname(unlist(estimate/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_stdError_graph = sqrt(delta_variance)
    
    estimates[measurement_error_counter-lowerbound+1] = estimate/lambda
    std_errors[measurement_error_counter-lowerbound+1] = delta_stdError_graph
  }
  
  
  ## Plot meta analysis for many studies
  estimates_forRMA = estimates[seq(1, length(estimates), 1)]
  stand_errs = std_errors[seq(1, length(std_errors), 1)]
  
  labels = c(lowerbound:upperbound)
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Example map to means with RC")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()
par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)

mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='Map to means with RC - average over 25 studies')
dev.copy(png,'map to means with RC average weights.png')
dev.off()
#mean_est = colMeans(per_error_est_final)
#plot(labels, mean_est)
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='Map to means with RC - average over 25 studies')
dev.copy(png,'map to means with RC average se.png')
dev.off()

#######################################################################################
######################## Bootstrapped direct mapping  #####################
#######################################################################################
# Obviously the same as the means.

numLevels = 4
validation_size = 400
validation_index_size = validation_size/numLevels
upperbound = 30
num_boots = 100
lowerbound = 10
numbaselines=25
set.seed(666)

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'bootstrap means multi.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

for (i in 1:numbaselines){
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    boot_est = vector()
    
    for (j in 1:num_boots){
      
      bootstrap_validation <- data.frame(index = as.factor(rep(x = c(1:numLevels), each= validation_index_size)))
      bootstrap_validation$exp_true <- unlist(unname(lapply(X = split(x=validation_data$exp_true, f= as.factor(validation_data$index)), 
                                                            FUN = sample, size = validation_index_size, replace=TRUE)))
 
      # Find mapping Model
      mapping_formula = as.formula(exp_true ~ index + 0)
      mapping_model <- lm(formula=mapping_formula, data=bootstrap_validation)
      new_study = data.frame(index = studyData$index)
      new_validation = data.frame(index = bootstrap_validation$index)
      bootstrap_validation$index_mapped = predict(mapping_model, newdata = new_validation)
      studyData$index_mapped = predict(mapping_model, newdata = new_study)
      
      # Use mapped values in model
      corrected_model <- lm(formula=outcome~index_mapped, data=studyData)
      estimate = corrected_model$coefficients["index_mapped"]
      
      boot_est = c(boot_est, estimate)
    }
    
    estimates[measurement_error_counter-lowerbound+1] = mean(boot_est)
    std_errors[measurement_error_counter-lowerbound+1] = var(boot_est)^0.5
    print(measurement_error_counter)
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
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Pure bootstrap")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
  
}  

dev.off()
par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)

mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='Bootstrap means - average over 25 studies')
dev.copy(png,'bs mean average weights.png')
dev.off()
mean_est = colMeans(per_error_est_final)
plot(labels, mean_est, xlab='Amount of error in instrument', ylab='Mean estimate for to study', main='Bootstrap means - average over 25 studies')
dev.copy(png,'bs mean average se.png')
dev.off()
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='Bootstrap means - average over 25 studies')
dev.copy(png,'bs mean average se.png')
dev.off()


#######################################################################################
######################## Direct mapping with RC #####################
#######################################################################################


numLevels = 4
validation_size = 400
lowerbound = 10
upperbound = 30
set.seed(666)
numbaselines = 25
## Plotting Bit
## Graphing Standard Error as a function of Measurement Error
per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'SEM multi.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

for (i in 1:numbaselines){
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    # Find Error Model
    mapping_formula = as.formula(exp_true ~ index + 0)
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    new = data.frame(index = studyData$index)
    studyData$index_mapped = predict(mapping_model, newdata = new)
    
    # Use Error Model to Scale Index to Gold
    corrected_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate = corrected_model$coefficients["index_mapped"]
    stdError = summary(corrected_model)$coefficients["index_mapped","Std. Error"]
    
    #assume lambda is 1
    lambda = 1
    lambda_stdError = summary(mapping_model)$coefficients[1,"Std. Error"]
    
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
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "direct mapping with RC")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}

#dev.copy(png,'direct index mapping multi.png')
dev.off()

par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)


mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='SEM with RC - average over 25 studies')
dev.copy(png,'SEM with RC average weights.png')
dev.off()
#mean_est = colMeans(per_error_est_final)
#plot(labels, mean_est)
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='SEM with RC - average over 25 studies')
dev.copy(png,'SEM with RC average se.png')
dev.off()


#######################################################################################
######################## Map to means with RC and repeated measures #####################
#######################################################################################

# Concern is that there is no estimate of error in mapping, only at the calibration stage

numLevels = 4
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 25

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'double measure.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))


for (i in 1:numbaselines){
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationDataRep(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    #measurement error
    
    meas_form1 = as.formula(exp_true ~ index + 0)
    meas_model1 <- lm(formula=meas_form1, data=validation_data)
    
    meas_form2 = as.formula(exp_true ~ index2 + 0)
    meas_model2 <- lm(formula=meas_form2, data=validation_data)

    validation_data$index_mapped = predict(meas_model1, newdata = validation_data)
    validation_data$index_mapped2 = predict(meas_model2, newdata = validation_data)
    studyData$index_mapped = predict(meas_model1, newdata = studyData)
    
    meas_err = lm(formula=index_mapped2~index_mapped, data=validation_data)
    
    reg_calib = lm(formula = exp_true ~ index_mapped, data=validation_data)
    
    alpha = summary(meas_err)$coefficients[2,1]
    alpha_var = summary(meas_err)$coefficients[2,2]^2
    lambda = 1
    lambda_var = summary(reg_calib)$coefficients[2,2]^2
    
    lambda_var_corr = lambda_var/alpha^2 + (lambda/alpha^2)*alpha_var
    lambda_std_error = lambda_var_corr^0.5
    
    # Use mapped values in model - estimate beta
    corrected_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate = corrected_model$coefficients["index_mapped"]
    std_error = summary(corrected_model)$coefficients[2,2]
    
    # recalibrate the standard error using delta function
    
    variance_beta = std_error^2
    beta_lambda_div_sq = (unname(unlist(estimate/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_stdError_graph = sqrt(delta_variance)
    
    estimates[measurement_error_counter-lowerbound+1] = estimate/lambda
    std_errors[measurement_error_counter-lowerbound+1] = delta_stdError_graph
  }
  
  
  ## Plot meta analysis for many studies
  estimates_forRMA = estimates[seq(1, length(estimates), 1)]
  stand_errs = std_errors[seq(1, length(std_errors), 1)]
  
  labels = c(lowerbound:upperbound)
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "double measure")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()
par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)

mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='double measure - average over 25 studies')
dev.copy(png,'map to means with RC average weights.png')
dev.off()
#mean_est = colMeans(per_error_est_final)
#plot(labels, mean_est)
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='double measure - average over 25 studies')
dev.copy(png,'map to means with RC average se.png')
dev.off()


#######################################################################################
######################## Bootstrapped individuals  #####################
#######################################################################################
# Obviously the same as the means.

numLevels = 4
validation_size = 400
validation_index_size = validation_size/numLevels
study_index_size = study_size/numLevels
upperbound = 30
num_boots = 1000
lowerbound = 10
numbaselines=25
set.seed(666)

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'bootstrap ind multi.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

for (i in 1:numbaselines){
  std_errors = vector("numeric", length = upperbound - lowerbound)
  estimates = vector("numeric", length = upperbound - lowerbound)
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    boot_est = vector()
    
    for (j in 1:num_boots){
      
      bootstrap_validation = validation_data[,c('exp_true', 'index')]
      bootstrap_validation$exp_measured <- unlist(unname(lapply(X = split(x=validation_data$exp_true, f= validation_data$index), 
                                                                FUN = sample, size = validation_index_size, replace=TRUE)))
      studyData$exp_measured <- unlist(unname(lapply(X = split(x=validation_data$exp_true, f= as.factor(validation_data$index)), 
                                                     FUN = sample, size = study_index_size, replace=TRUE)))
      
      # Find calibration Model
      calib_formula = as.formula(exp_true ~ exp_measured)
      calib_model <- lm(formula=calib_formula, data=bootstrap_validation)
      studyData$exp_calib = predict(calib_model, newdata = studyData)
      
      # Use mapped values in model
      corrected_model <- lm(formula=outcome~exp_calib, data=studyData)
      estimate = corrected_model$coefficients["exp_calib"]
      
      boot_est = c(boot_est, estimate)
    }
    
    estimates[measurement_error_counter-lowerbound+1] = mean(boot_est)
    std_errors[measurement_error_counter-lowerbound+1] = var(boot_est)^0.5
    print(measurement_error_counter)
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
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
  fmla = as.formula(y~harmonised_x)
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Per person bootstrap")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
  
}  

dev.off()
par(old.par)

per_error_weights_final  = do.call(rbind, per_error_weights)
per_error_est_final  = do.call(rbind, per_error_est)
per_error_se_final  = do.call(rbind, per_error_se)

mean_weights = colMeans(per_error_weights_final)
plot(x=labels, y=mean_weights, xlab='Amount of error in instrument', ylab='Mean weight given to study', main='Bootstrap inds - average over 25 studies')
dev.copy(png,'bs ind average weights.png')
dev.off()
mean_est = colMeans(per_error_est_final)
plot(labels, mean_est, xlab='Amount of error in instrument', ylab='Mean estimate for to study', main='Bootstrap inds - average over 25 studies')
dev.copy(png,'bs ind average se.png')
dev.off()
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels, y=mean_se, xlab='Amount of error in instrument', ylab='Mean se for study', main='Bootstrap inds - average over 25 studies')
dev.copy(png,'bs ind average se.png')
dev.off()
mean_se = colMeans(per_error_se_final)
plot(labels, mean_se)
plot(x=labels[1:19], y=mean_se[1:19], xlab='Amount of error in instrument', ylab='Mean se for study', main='Bootstrap inds - average over 25 studies')
dev.copy(png,'bs ind average se without 30.png')
dev.off()


