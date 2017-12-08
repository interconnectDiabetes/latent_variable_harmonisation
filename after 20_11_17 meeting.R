
# Set up the Data
library(ggplot2)
library(reshape2)

library(stats)
library(metafor)

# Set save point for the plots
setwd("U:/workspace/paee_harmonisation/Reweighting/plots")

# Raw Study Data Set Properties
trueBeta = 0.5
constant = 20
set.seed(666)
study_size = 20000
raw_mean = 40
raw_stdev = 15
outcome_error = 0

# Simple creation of the true exposure data and the related outcome with error - indices etc created later
raw_data = data.frame(exp_true = rnorm(n = study_size, mean = raw_mean, sd = raw_stdev))
raw_data$outcome = (trueBeta*raw_data$exp_true) + constant + rnorm(n = study_size, mean = 0, sd = outcome_error)


###############################################################################
########################### Functions #########################################
###############################################################################

#Use the raw data to generate the indices for a study
createStudyData <- function(raw_data, measurement_error, cut_points){
  exposure_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
  raw_data$exposure_measured = raw_data$exp_true + exposure_error
  raw_data$index = cut(x=raw_data$exposure_measured, breaks = cut_points, labels = FALSE)
  raw_data$index = as.factor(raw_data$index)
  raw_data = raw_data[with(raw_data, order(index)),]
  return(raw_data)
}

# Create a validation data set
createValidationData <- function(val_size, measurement_error, cut_points) {
  validation_data = data.frame(exp_true =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
  exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
  validation_data$exposure_measured = validation_data$exp_true + exposure_error
  validation_data$index = cut(x=validation_data$exposure_measured, breaks = cut_points, labels = FALSE)
  validation_data$index = as.factor(validation_data$index)
  validation_data = validation_data[order(validation_data$index),]
  return (validation_data)
}

#create a validation data set with a repeat measure
createValidationDataRep <- function(val_size, measurement_error, cut_points) {
  validation_data = data.frame(exp_true =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
  exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
  exposure_error2 = rnorm(n = val_size, mean = 0, sd = measurement_error)
  validation_data$exposure_measured = validation_data$exp_true + exposure_error
  validation_data$exposure_measured2 = validation_data$exp_true + exposure_error2
  validation_data$index = cut(x=validation_data$exposure_measured, breaks = cut_points, labels = FALSE)
  validation_data$index = as.factor(validation_data$index)
  validation_data$index2 = cut(x=validation_data$exposure_measured2, breaks = cut_points, labels = FALSE)
  validation_data$index2 = as.factor(validation_data$index2)
  return (validation_data)
}


# Create a list of means per index level
createMeansList <- function(data, number_of_indices) {
  meansList = vector(mode="list", length = number_of_indices)
  for (i in 1:number_of_indices) {
    meansList[i] = mean(unname(unlist((split(x=data$exp_true, f= as.factor(data$index)))[i])))
  }
  return(meansList)
}

######################################




#######################################################################################
####### Continuous and REPEAT categorical (index mean subst) with RC ##################
#######################################################################################
#
# Code to run continuous and categorical (sub in index mean) measured exposures
# RC is used in both scenarios, with a SE correction using the delta method
# Both cont and cat versions are plotted on the same forest plot
# 
# The code is run in a loop, with the error in the exposure increasing from lowerbound to upperbound in increments of 1
# At the end, the results are meta-analysed with a random effects model (to represent how they might be used)
# This is saved as a forest plot
#
# The REPEAT means that the validation study has 2 measures (ie index values) per person. The first measures is used to map from indices
# to means. The second measure also undergoes mapping and is used to perform regression calibration between the index means
# 


# set up parameters
#numLevels = 2
#cut_points = c(-Inf,40,Inf)
numLevels = 4
cut_points = c(-Inf,30,40,50,Inf)
#numLevels = 8
#cut_points = c(-Inf,10,20,30,40,50,60,70,Inf)
#numLevels = 16
#cut_points = c(-Inf,0,10,16,22,26,30,35,40,45,50,54,58,64,70,78,Inf)
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('dualRC_repeat_valid_cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  # structures for storage
  std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  estimates_cat = vector("numeric", length = upperbound - lowerbound)
  lambda_cat = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cat = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cat = vector("numeric", length = upperbound - lowerbound)
  std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  estimates_cont = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cont = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cont = vector("numeric", length = upperbound - lowerbound)
  
  means = data.frame()
  sds = data.frame()
  sizes = data.frame()
  
  for (measurement_error_counter in lowerbound:upperbound) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    # generate 2 validation studies of half the intended size to simulate splitting the validation study in half
    validation_data = createValidationDataRep(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points )
    
    
    
    #Store information about indices
    #means = rbind(means,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=mean)[,2]))
    #sds = rbind(sds,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=sd)[,2]))
    #sizes = rbind(sizes,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=length)[,2]))
    
    ######### CONTINUOUS
    
    #naive model - outcome with measured exposure
    naive_model <- lm(formula=outcome~exposure_measured, data=studyData)
    estimate_naive = naive_model$coefficients["exposure_measured"]
    std_error_naive = summary(naive_model)$coefficients["exposure_measured","Std. Error"]
    
    # Find RC Model (relationship between measured exposure and true exposure)
    # Then use this to substitute in the corrected exposure in the main study and validation study
    RC_formula = as.formula(exp_true ~ exposure_measured)
    RC_model <- lm(formula=RC_formula, data=validation_data)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)
    
    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting SE
    
    lambda = summary(RC_model)$coefficients["exposure_measured","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["exposure_measured","Std. Error"]
    
    # recalibrate the standard error using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    estimates_cont[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cont[measurement_error_counter-lowerbound+1] = delta_std_error
    
    # Store lambda SE and beta SE in case we need them (not used)
    lambda_SE_cont[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cont[measurement_error_counter-lowerbound+1] = std_error_naive
    
    
    ##### CATEGORICAL
    
    # Additional step for categorical - find relationship between index and true exposure
    # Then use this to substitute in the index means in the main study and validation study
    # Do this for both measurements
    mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    validation_data$index_mapped = predict(mapping_model, newdata = validation_data)
    studyData$index_mapped = predict(mapping_model, newdata = studyData)


    mapping_formula = as.formula(exp_true ~ index2 + 0) # 0 forces no intercept
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    validation_data$index_mapped2 = predict(mapping_model, newdata = validation_data)


    #naive model - outcome with mapped exposure
    #although not that naive since we did the mapping...

    naive_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate_naive = naive_model$coefficients["index_mapped"]
    std_error_naive = summary(naive_model)$coefficients["index_mapped","Std. Error"]

    # Find RC Model (relationship between two sets of index means to estimate error)
    # Then use this to substitute in the corrected exposure in the main study

    RC_formula = as.formula(index_mapped2 ~ index_mapped)
    RC_model <- lm(formula=RC_formula, data=validation_data)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)

    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    # Note that this is the same as the naive model since lambda is 1 no correction occured

    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]

    # Extract the coefficients from the RC model for correcting estimate SE

    lambda = summary(RC_model)$coefficients["index_mapped","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["index_mapped","Std. Error"]

    # refine the standard error of the main model estimate using delta function

    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)

    # Store estimate and corrected std error

    estimates_cat[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cat[measurement_error_counter-lowerbound+1] = delta_std_error

    # Store lambda SE and beta SE in case we need them (not used)

    lambda_cat[measurement_error_counter-lowerbound+1] = lambda
    lambda_SE_cat[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cat[measurement_error_counter-lowerbound+1] = std_error

    
    
    #################################################################################
    
    
    
    # ## new first step
    # # Additional step for categorical - find relationship between index and true exposure
    # # Then use this to substitute in the index means in the main study and validation study
    # 
    # mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
    # mapping_model <- lm(formula=mapping_formula, data=validation_data)
    # validation_data$index_mapped = predict(mapping_model, newdata = validation_data)
    # studyData$index_mapped = predict(mapping_model, newdata = studyData)
    # 
    # 
    # mapping_formula2 = as.formula(exp_true ~ index2 + 0) # 0 forces no intercept
    # mapping_model2 <- lm(formula=mapping_formula2, data=validation_data)
    # validation_data$index_mapped2 = predict(mapping_model2, newdata = validation_data)
    # 
    # 
    # #naive model - outcome with mapped exposure
    # #although not that naive since we did the mapping...
    # 
    # naive_model <- lm(formula=outcome~index_mapped, data=studyData)
    # estimate_naive = naive_model$coefficients["index_mapped"]
    # std_error_naive = summary(naive_model)$coefficients["index_mapped","Std. Error"]
    # 
    # #First stage RC
    # 
    # RC_formula_1 = as.formula(exp_true ~ index_mapped)
    # RC_model_1 <- lm(formula=RC_formula_1, data=validation_data)
    # #already done
    # #studyData$exposure_corrected = predict(RC_model_1, newdata = studyData)
    # 
    # RC_formula_2 = as.formula(index_mapped2 ~ index_mapped)
    # RC_model_2 <- lm(formula=RC_formula_2, data=validation_data)
    # studyData$exposure_corrected = predict(RC_model_2, newdata = studyData)
    # 
    # # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # # store the estimate and its standard error
    # # Note that this is the same as the naive model since lambda is 1 no correction occured
    # 
    # corrected_model_1 <- lm(formula=outcome~index_mapped, data=studyData)
    # estimate_1 = corrected_model_1$coefficients["index_mapped"]
    # std_error_1 = summary(corrected_model_1)$coefficients["index_mapped","Std. Error"]
    # 
    # corrected_model_2 <- lm(formula=outcome~exposure_corrected, data=studyData)
    # estimate_2 = corrected_model_2$coefficients["exposure_corrected"]
    # std_error_2 = summary(corrected_model_2)$coefficients["exposure_corrected","Std. Error"]
    # 
    # # Extract the coefficients from the RC model for correcting estimate SE
    # 
    # lambda_1 = summary(RC_model_1)$coefficients["index_mapped","Estimate"]
    # lambda_std_error_1 = summary(RC_model_1)$coefficients["index_mapped","Std. Error"]
    # 
    # lambda_2 = summary(RC_model_2)$coefficients["index_mapped","Estimate"]
    # lambda_std_error_2 = summary(RC_model_2)$coefficients["index_mapped","Std. Error"]
    # 
    # # refine the standard error of the main model estimate using delta function
    # 
    # variance_beta = std_error_naive^2
    # beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda_1)^2)))^2
    # var_lambda = (lambda_std_error_1)^2
    # delta_variance = (variance_beta / (lambda_1)^2) + (beta_lambda_div_sq * var_lambda)
    # delta_std_error_1 = sqrt(delta_variance)
    # 
    # variance_beta = delta_std_error_1^2
    # beta_lambda_div_sq = (unname(unlist(estimate_1/(lambda_2)^2)))^2
    # var_lambda = (lambda_std_error_2)^2
    # delta_variance = (variance_beta / (lambda_2)^2) + (beta_lambda_div_sq * var_lambda)
    # delta_std_error_2 = sqrt(delta_variance)
    # 
    # # Store estimate and corrected std error
    # 
    # estimates_cat[measurement_error_counter-lowerbound+1] = estimate_2
    # std_errors_cat[measurement_error_counter-lowerbound+1] = delta_std_error_2
    # 
    # # Store lambda SE and beta SE in case we need them (not used)
    # 
    # lambda_SE_cat[measurement_error_counter-lowerbound+1] = lambda_std_error
    # beta_SE_cat[measurement_error_counter-lowerbound+1] = std_error
  }
  
  
  ## Plot meta analysis for many studies, both cont and cat
  estimates_forRMA = c(estimates_cat[seq(1, length(estimates_cat), 1)],estimates_cont[seq(1, length(estimates_cont), 1)])
  stand_errs = c(std_errors_cat[seq(1, length(std_errors_cat), 1)], std_errors_cont[seq(1, length(std_errors_cont), 1)])
  lambda_SE = c(lambda_SE_cat[seq(1, length(lambda_SE_cat), 1)], lambda_SE_cont[seq(1, length(lambda_SE_cont), 1)])
  
  labels = c(paste0(c(lowerbound:upperbound), '_cat'), paste0(c(lowerbound:upperbound), '_cont'))
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%) ",round(stand_errs,digits=3))
  
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association H'[0]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = paste0('cat and cont with modelled RC outcome error ',outcome_error, ' levels ', numLevels, ' validation ', validation_size ))
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], "Study error, weight, SE", adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()




#######################################################################################
####### Continuous and split measurement categorical (index mean subst) with RC ########
#######################################################################################
#
# Code to run continuous and categorical (sub in index mean) measured exposures
# RC is used in both scenarios, with a SE correction using the delta method
# Both cont and cat versions are plotted on the same forest plot
# 
# The code is run in a loop, with the error in the exposure increasing from lowerbound to upperbound in increments of 1
# At the end, the results are meta-analysed with a random effects model (to represent how they might be used)
# This is saved as a forest plot
#
# The SPLIT means that the validation study is split half. One half is used to map from indices
# to means. The other half is used to perform regression calibration between the index means
# and the true measured exposures


# set up parameters
#numLevels = 2
#cut_points = c(-Inf,40,Inf)
#numLevels = 4
#cut_points = c(-Inf,30,40,50,Inf)
#numLevels = 8
#cut_points = c(-Inf,10,20,30,40,50,60,70,Inf)
#numLevels = 16
#cut_points = c(-Inf,0,10,16,22,26,30,35,40,45,50,54,58,64,70,78,Inf)
numLevels = 32
cut_points = c(-Inf,0,5,7,10,12,16,19,22,23,26,28,30,32,35,37,40,43,45,47,50,52,54,56,58,61,64,67,70,72,74,78,Inf)
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('split_valid_cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  # structures for storage
  std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  estimates_cat = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cat = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cat = vector("numeric", length = upperbound - lowerbound)
  std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  estimates_cont = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cont = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cont = vector("numeric", length = upperbound - lowerbound)
  
  means = data.frame()
  sds = data.frame()
  sizes = data.frame()
  
  for (measurement_error_counter in lowerbound:upperbound) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    
    # split the validation study in half
    
    # validation_data$split = unlist(lapply(X = split(x = validation_data, validation_data$index),FUN = function(X){
    #   y = cut_number(x = as.numeric(X$index), n = 2, labels=FALSE)
    #   return(y)
    # }))
    
    validation_data$split = unlist(lapply(X = split(x = validation_data, validation_data$index),FUN = function(X){
      y1 = rep(x=1,times=ceiling(length(X$index)/2))
      y2 = rep(x=2,times=length(X$index)-ceiling(length(X$index)/2))
      y = c(y1,y2)
      return(y)
    }))
    
    
    validation_data_A = validation_data[validation_data$split == 1,]
    validation_data_B = validation_data[validation_data$split == 2,]
    
    
    #Store information about indices
    #means = rbind(means,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=mean)[,2]))
    #sds = rbind(sds,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=sd)[,2]))
    #sizes = rbind(sizes,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=length)[,2]))
    
    ######### CONTINUOUS
    
    #naive model - outcome with measured exposure
    naive_model <- lm(formula=outcome~exposure_measured, data=studyData)
    estimate_naive = naive_model$coefficients["exposure_measured"]
    std_error_naive = summary(naive_model)$coefficients["exposure_measured","Std. Error"]
    
    # Find RC Model (relationship between measured exposure and true exposure)
    # Then use this to substitute in the corrected exposure in the main study and validation study
    RC_formula = as.formula(exp_true ~ exposure_measured)
    RC_model <- lm(formula=RC_formula, data=validation_data)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)
    
    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting SE
    
    lambda = summary(RC_model)$coefficients["exposure_measured","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["exposure_measured","Std. Error"]
    
    # recalibrate the standard error using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    estimates_cont[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cont[measurement_error_counter-lowerbound+1] = delta_std_error
    
    # Store lambda SE and beta SE in case we need them (not used)
    lambda_SE_cont[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cont[measurement_error_counter-lowerbound+1] = std_error_naive
    
    
    ##### CATEGORICAL
    
    # Additional step for categorical - find relationship between index and true exposure
    # Then use this to substitute in the index means in the main study and validation study
    mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
    mapping_model <- lm(formula=mapping_formula, data=validation_data_A)
    validation_data_B$index_mapped = predict(mapping_model, newdata = validation_data_B)
    studyData$index_mapped = predict(mapping_model, newdata = studyData)
    
    #naive model - outcome with mapped exposure
    #although not that naive since we did the mapping...
    
    naive_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate_naive = naive_model$coefficients["index_mapped"]
    std_error_naive = summary(naive_model)$coefficients["index_mapped","Std. Error"]
    
    # Find RC Model (relationship between measured exposure and true exposure)
    # Then use this to substitute in the corrected exposure in the main study and validation study
    # However we have already done the mapping and exposure_corrected is just the same as
    # index_mapped. This is because lambda is ~1
    
    RC_formula = as.formula(exp_true ~ index_mapped)
    RC_model <- lm(formula=RC_formula, data=validation_data_B)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)
    
    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    # Note that this is the same as the naive model since lambda is 1 no correction occured
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting estimate SE
    
    lambda = summary(RC_model)$coefficients["index_mapped","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["index_mapped","Std. Error"]
    
    # refine the standard error of the main model estimate using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    
    estimates_cat[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cat[measurement_error_counter-lowerbound+1] = delta_std_error
    
    # Store lambda SE and beta SE in case we need them (not used)
    
    lambda_SE_cat[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cat[measurement_error_counter-lowerbound+1] = std_error
    
  }
  
  
  ## Plot meta analysis for many studies, both cont and cat
  estimates_forRMA = c(estimates_cat[seq(1, length(estimates_cat), 1)],estimates_cont[seq(1, length(estimates_cont), 1)])
  stand_errs = c(std_errors_cat[seq(1, length(std_errors_cat), 1)], std_errors_cont[seq(1, length(std_errors_cont), 1)])
  lambda_SE = c(lambda_SE_cat[seq(1, length(lambda_SE_cat), 1)], lambda_SE_cont[seq(1, length(lambda_SE_cont), 1)])
  
  labels = c(paste0(c(lowerbound:upperbound), '_cat'), paste0(c(lowerbound:upperbound), '_cont'))
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%) ",round(stand_errs,digits=3))
  
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association H'[0]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = paste0('cat and cont with modelled RC outcome error ',outcome_error, ' levels ', numLevels, ' validation ', validation_size ))
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], "Study error, weight, SE", adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()



#######################################################################################
######## Baseline Continuous and categorical (index mean subst) with RC ###############
#######################################################################################
#
# Code to run continuous and categorical (sub in index mean) measured exposures
# RC is used in both scenarios, with a SE correction using the delta method
# Both cont and cat versions are plotted on the same forest plot
# 
# The code is run in a loop, with the error in the exposure increasing from lowerbound to upperbound in increments of 1
# At the end, the results are meta-analysed with a random effects model (to represent how they might be used)
# This is saved as a forest plot
#


# set up parameters
#numLevels = 2
#cut_points = c(-Inf,40,Inf)
#numLevels = 4
#cut_points = c(-Inf,30,40,50,Inf)
#numLevels = 8
#cut_points = c(-Inf,10,20,30,40,50,60,70,Inf)
#numLevels = 16
#cut_points = c(-Inf,0,10,16,22,26,30,35,40,45,50,54,58,64,70,78,Inf)
numLevels = 32
cut_points = c(-Inf,0,5,7,10,12,16,19,22,23,26,28,30,32,35,37,40,43,45,47,50,52,54,56,58,61,64,67,70,72,74,78,Inf)
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('baseline_cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  # structures for storage
  std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  estimates_cat = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cat = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cat = vector("numeric", length = upperbound - lowerbound)
  std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  estimates_cont = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cont = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cont = vector("numeric", length = upperbound - lowerbound)
  
  means = data.frame()
  sds = data.frame()
  sizes = data.frame()
  
  for (measurement_error_counter in lowerbound:upperbound) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    #Store information about indices
    means = rbind(means,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=mean)[,2]))
    sds = rbind(sds,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=sd)[,2]))
    sizes = rbind(sizes,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=length)[,2]))
    
    ######### CONTINUOUS
    
    #naive model - outcome with measured exposure
    naive_model <- lm(formula=outcome~exposure_measured, data=studyData)
    estimate_naive = naive_model$coefficients["exposure_measured"]
    std_error_naive = summary(naive_model)$coefficients["exposure_measured","Std. Error"]
    
    # Find RC Model (relationship between measured exposure and true exposure)
    # Then use this to substitute in the corrected exposure in the main study and validation study
    RC_formula = as.formula(exp_true ~ exposure_measured)
    RC_model <- lm(formula=RC_formula, data=validation_data)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)
    
    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting SE
    
    lambda = summary(RC_model)$coefficients["exposure_measured","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["exposure_measured","Std. Error"]
    
    # recalibrate the standard error using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    estimates_cont[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cont[measurement_error_counter-lowerbound+1] = delta_std_error
    
    # Store lambda SE and beta SE in case we need them (not used)
    lambda_SE_cont[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cont[measurement_error_counter-lowerbound+1] = std_error_naive
    
    
    ##### CATEGORICAL
    
    # Additional step for categorical - find relationship between index and true exposure
    # Then use this to substitute in the index means in the main study and validation study
    mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    validation_data$index_mapped = predict(mapping_model, newdata = validation_data)
    studyData$index_mapped = predict(mapping_model, newdata = studyData)
    
    #naive model - outcome with mapped exposure
    #although not that naive since we did the mapping...
    
    naive_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate_naive = naive_model$coefficients["index_mapped"]
    std_error_naive = summary(naive_model)$coefficients["index_mapped","Std. Error"]
    
    # Find RC Model (relationship between measured exposure and true exposure)
    # Then use this to substitute in the corrected exposure in the main study and validation study
    # However we have already done the mapping and exposure_corrected is just the same as
    # index_mapped. This is because lambda is ~1
    
    RC_formula = as.formula(exp_true ~ index_mapped)
    RC_model <- lm(formula=RC_formula, data=validation_data)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)
    
    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    # Note that this is the same as the naive model since lambda is 1 no correction occured
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting estimate SE
    
    lambda = summary(RC_model)$coefficients["index_mapped","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["index_mapped","Std. Error"]
    
    # refine the standard error of the main model estimate using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    
    estimates_cat[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cat[measurement_error_counter-lowerbound+1] = delta_std_error
    
    # Store lambda SE and beta SE in case we need them (not used)
    
    lambda_SE_cat[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cat[measurement_error_counter-lowerbound+1] = std_error
    
  }
  
  
  ## Plot meta analysis for many studies, both cont and cat
  estimates_forRMA = c(estimates_cat[seq(1, length(estimates_cat), 1)],estimates_cont[seq(1, length(estimates_cont), 1)])
  stand_errs = c(std_errors_cat[seq(1, length(std_errors_cat), 1)], std_errors_cont[seq(1, length(std_errors_cont), 1)])
  lambda_SE = c(lambda_SE_cat[seq(1, length(lambda_SE_cat), 1)], lambda_SE_cont[seq(1, length(lambda_SE_cont), 1)])
  
  labels = c(paste0(c(lowerbound:upperbound), '_cat'), paste0(c(lowerbound:upperbound), '_cont'))
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%) ",round(stand_errs,digits=3))
  
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association H'[0]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = paste0('cat and cont with modelled RC outcome error ',outcome_error, ' levels ', numLevels, ' validation ', validation_size ))
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], "Study error, weight, SE", adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()



#######################################################################################
####### Continuous and REPEAT categorical (index mean subst) with RC ##################
#######################################################################################
#
# Code to run continuous and categorical (sub in index mean) measured exposures
# RC is used in both scenarios, with a SE correction using the delta method
# Both cont and cat versions are plotted on the same forest plot
# 
# The code is run in a loop, with the error in the exposure increasing from lowerbound to upperbound in increments of 1
# At the end, the results are meta-analysed with a random effects model (to represent how they might be used)
# This is saved as a forest plot
#
# The REPEAT means that the validation study has 2 measures (ie index values) per person. The first measures is used to map from indices
# to means. The second measure also undergoes mapping and is used to perform regression calibration between the index means
# 


# set up parameters
#numLevels = 2
#cut_points = c(-Inf,40,Inf)
#numLevels = 4
#cut_points = c(-Inf,30,40,50,Inf)
#numLevels = 8
#cut_points = c(-Inf,10,20,30,40,50,60,70,Inf)
numLevels = 16
cut_points = c(-Inf,0,10,16,22,26,30,35,40,45,50,54,58,64,70,78,Inf)
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('alt_repeat_valid_cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  # structures for storage
  std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  estimates_cat = vector("numeric", length = upperbound - lowerbound)
  lambda_cat = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cat = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cat = vector("numeric", length = upperbound - lowerbound)
  std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  estimates_cont = vector("numeric", length = upperbound - lowerbound)
  lambda_SE_cont = vector("numeric", length = upperbound - lowerbound)
  beta_SE_cont = vector("numeric", length = upperbound - lowerbound)
  
  means = data.frame()
  sds = data.frame()
  sizes = data.frame()
  
  for (measurement_error_counter in lowerbound:upperbound) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    # generate 2 validation studies of half the intended size to simulate splitting the validation study in half
    validation_data = createValidationDataRep(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points )
    
    validation_data_A = validation_data[,c('exp_true', 'index')]
    validation_data_B = validation_data[,c('exp_true', 'index2')]
    colnames(validation_data_B) = c('exp_true', 'index')
    
    #Store information about indices
    #means = rbind(means,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=mean)[,2]))
    #sds = rbind(sds,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=sd)[,2]))
    #sizes = rbind(sizes,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=length)[,2]))
    
    ######### CONTINUOUS
    
    #naive model - outcome with measured exposure
    naive_model <- lm(formula=outcome~exposure_measured, data=studyData)
    estimate_naive = naive_model$coefficients["exposure_measured"]
    std_error_naive = summary(naive_model)$coefficients["exposure_measured","Std. Error"]
    
    # Find RC Model (relationship between measured exposure and true exposure)
    # Then use this to substitute in the corrected exposure in the main study and validation study
    RC_formula = as.formula(exp_true ~ exposure_measured)
    RC_model <- lm(formula=RC_formula, data=validation_data)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)
    
    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting SE
    
    lambda = summary(RC_model)$coefficients["exposure_measured","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["exposure_measured","Std. Error"]
    
    # recalibrate the standard error using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    estimates_cont[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cont[measurement_error_counter-lowerbound+1] = delta_std_error
    
    # Store lambda SE and beta SE in case we need them (not used)
    lambda_SE_cont[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cont[measurement_error_counter-lowerbound+1] = std_error_naive
    
    
    ##### CATEGORICAL
    
    # Additional step for categorical - find relationship between index and true exposure
    # Then use this to substitute in the index means in the main study and validation study
    mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
    mapping_model <- lm(formula=mapping_formula, data=validation_data_A)
    validation_data_B$index_mapped = predict(mapping_model, newdata = validation_data_B)
    studyData$index_mapped = predict(mapping_model, newdata = studyData)
    
    #naive model - outcome with mapped exposure
    #although not that naive since we did the mapping...
    
    naive_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate_naive = naive_model$coefficients["index_mapped"]
    std_error_naive = summary(naive_model)$coefficients["index_mapped","Std. Error"]
    
    # Find RC Model (relationship between measured exposure and true exposure)
    # Then use this to substitute in the corrected exposure in the main study and validation study
    # However we have already done the mapping and exposure_corrected is just the same as
    # index_mapped. This is because lambda is ~1
    
    RC_formula = as.formula(exp_true ~ index_mapped)
    RC_model <- lm(formula=RC_formula, data=validation_data_B)
    studyData$exposure_corrected = predict(RC_model, newdata = studyData)
    
    # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    # Note that this is the same as the naive model since lambda is 1 no correction occured
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting estimate SE
    
    lambda = summary(RC_model)$coefficients["index_mapped","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["index_mapped","Std. Error"]
    
    # refine the standard error of the main model estimate using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    
    estimates_cat[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cat[measurement_error_counter-lowerbound+1] = delta_std_error
    
    # Store lambda SE and beta SE in case we need them (not used)
    
    lambda_SE_cat[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cat[measurement_error_counter-lowerbound+1] = std_error
    
  
  }
  
  
  ## Plot meta analysis for many studies, both cont and cat
  estimates_forRMA = c(estimates_cat[seq(1, length(estimates_cat), 1)],estimates_cont[seq(1, length(estimates_cont), 1)])
  stand_errs = c(std_errors_cat[seq(1, length(std_errors_cat), 1)], std_errors_cont[seq(1, length(std_errors_cont), 1)])
  lambda_SE = c(lambda_SE_cat[seq(1, length(lambda_SE_cat), 1)], lambda_SE_cont[seq(1, length(lambda_SE_cont), 1)])
  
  labels = c(paste0(c(lowerbound:upperbound), '_cat'), paste0(c(lowerbound:upperbound), '_cont'))
  res <- rma(yi = estimates_forRMA, sei = stand_errs, method='DL', slab = labels)
  weights_res <- weights.rma.uni(res)
  
  per_error_weights[[i]] = weights_res
  per_error_est[[i]] = estimates_forRMA
  per_error_se[[i]] = stand_errs
  
  # Forest Plot
  res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%) ",round(stand_errs,digits=3))
  
  forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                .(sprintf("%.3f", round(res$QEp,3))),')')),
         xlab=bquote(paste('Test of Association H'[0]*': true beta association = 0, p = ',
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = paste0('cat and cont with modelled RC outcome error ',outcome_error, ' levels ', numLevels, ' validation ', validation_size ))
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], "Study error, weight, SE", adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
}


dev.off()


