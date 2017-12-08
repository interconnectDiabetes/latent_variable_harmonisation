
# Set up the Data
library(ggplot2)
library(reshape2)
library(parallel)
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
createStudyData <- function(raw_data, measurement_error, number_of_indices){
  exposure_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
  raw_data$exposure_measured = raw_data$exp_true + exposure_error
  raw_data$index = cut_number(x = raw_data$exposure_measured,n = number_of_indices, labels =FALSE)
  raw_data$index = as.factor(raw_data$index)
  raw_data = raw_data[with(raw_data, order(index)),]
  return(raw_data)
}

#Use the raw data to generate the indices for a study for fixed levels
createStudyDataFix <- function(raw_data, measurement_error, cut_points){
  exposure_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
  raw_data$exposure_measured = raw_data$exp_true + exposure_error
  raw_data$index = cut(x=raw_data$exposure_measured, breaks = cut_points, labels = FALSE)
  raw_data$index = as.factor(raw_data$index)
  raw_data = raw_data[with(raw_data, order(index)),]
  return(raw_data)
}

# Create a validation data set
createValidationData <- function(val_size, measurement_error, number_of_indices) {
  validation_data = data.frame(exp_true =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
  exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
  validation_data$exposure_measured = validation_data$exp_true + exposure_error
  validation_data$index = cut_number(x = validation_data$exposure_measured,n = number_of_indices, labels =FALSE)
  validation_data$index = as.factor(validation_data$index)
  validation_data = validation_data[order(validation_data$index),]
  return (validation_data)
}

# Create a validation data set for fixed levels
createValidationDataFix <- function(val_size, measurement_error, cut_points) {
  validation_data = data.frame(exp_true =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
  exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
  validation_data$exposure_measured = validation_data$exp_true + exposure_error
  validation_data$index = cut(x=validation_data$exposure_measured, breaks = cut_points, labels = FALSE)
  validation_data$index = as.factor(validation_data$index)
  validation_data = validation_data[order(validation_data$index),]
  return (validation_data)
}

# Function not normally used - for generating a validation study with repeat measures
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

# Create a list of means per index level
createMeansList <- function(data, number_of_indices) {
  meansList = vector(mode="list", length = number_of_indices)
  for (i in 1:number_of_indices) {
    meansList[i] = mean(unname(unlist((split(x=data$exp_true, f= as.factor(data$index)))[i])))
  }
  return(meansList)
}




#######################################################################################
############# Continuous and categorical (index mean subst) with RC #####################
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
numLevels = 16
validation_size = 800
upperbound = 30
lowerbound = 10
numbaselines = 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
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
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
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
######################## Bootstrapped continuous and categorical  #####################
#######################################################################################
# In this section we generate bootstrap values for the estimated association and standard error
# This helps us to validate whether the estimate and SE from using the delta method above
# are correct (this was Ian White's adivce).
#
# In the continuous example, we randomly draw from the validation study and main study, 
# and performing the classical RC process 1000 times 
# (i.e. calculate lambda, correct error-prone measurements, fit corrected model to estimate 
# association). The Delta correction is not applied. Our bootstrapped estimate is the mean
# of the estimates, and our bootstrapped SE is the SD of the estimates.
#
# In the categorical example, we randomly draw from the validation study and main study. 
# In the sampled validation study we calculate the index means and substitute into the sampled
#main study. Then we fit the corrected model. This is done 1000 times. The Delta correction is
# not applied. Our bootstrapped estimate is the mean of the estimates, and our bootstrapped SE
# is the SD of the estimates.


# set up parameters

numLevels = 4
validation_size = 400
validation_index_size = validation_size/numLevels
study_index_size = study_size/numLevels
upperbound = 30
num_boots = 100
lowerbound = 10
numbaselines= 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

# set up plot
svg(filename = paste0('bootstrapped_',num_boots,'_cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  
  std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  estimates_cat = vector("numeric", length = upperbound - lowerbound)

  std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  estimates_cont = vector("numeric", length = upperbound - lowerbound)
  
  for (measurement_error_counter in lowerbound:upperbound) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study - these remain constant for 
    # each bootstrap
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    boot_est_cat = vector()
    boot_est_cont = vector()
    
    for (j in 1:num_boots){
      
      
      # bootstrap sample from the validation study, note that we do the sample on a
      # per index basis. For example, if we have a validation study of 400 people, 4 levels and 100
      # people in each level. We draw 100 samples with replacement from index level 1, 100 from
      # level 2 and so on, as opposed to just drawing 400 random samples with replacement from the
      # validation study, which might yield 95 in level 1, 93 in level 2, 109 in level 3 and
      # 99 in level 4.
      # The reason for doing this is that it guarantees that a level is not empty when dealing
      # with smaller validation studies / lots of levels
      
      bootstrap_validation <- data.frame(index = as.factor(rep(x = c(1:numLevels), each= validation_index_size)))
      bootstrap_validation[,c('exp_true','exposure_measured')] <- unname(do.call('rbind',lapply(X = split(x=validation_data[,c('exp_true','exposure_measured')], f= as.factor(validation_data$index)), 
                                                            FUN = 
                                                            function(x){
                                                              output = x[sample(1:validation_index_size, replace=TRUE),]
                                                              return (output)
                                                            }
                                                            
                                                            )))
      
      # we also bootstrap the main study (as advised by Ian)
      
      bootstrap_study <- data.frame(index = as.factor(rep(x = c(1:numLevels), each= study_index_size)))
      bootstrap_study[, c("exp_true", "outcome", "exposure_measured", "index")] =
        unname(do.call('rbind', lapply(
          X = split(x = studyData[, c("exp_true", "outcome", "exposure_measured", "index")], f = as.factor(studyData$index)),
          FUN =
            function(x) {
              output = x[sample(1:study_index_size, replace = TRUE), ]
              return (output)
            }
          
        )))
      
      ##### CONTINUOUS
      
      # Find RC Model (relationship between measured exposure and true exposure)
      # Then use this to substitute in the corrected exposure in the main study and validation study
      RC_formula = as.formula(exp_true ~ exposure_measured)
      RC_model <- lm(formula=RC_formula, data=bootstrap_validation)
      bootstrap_study$exposure_corrected = predict(RC_model, newdata = bootstrap_study)
      
      # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
      # store the estimate and its standard error
      corrected_model <- lm(formula=outcome~exposure_corrected, data=bootstrap_study)
      estimate = corrected_model$coefficients["exposure_corrected"]
      
      # store estimate
      boot_est_cont = c(boot_est_cont, estimate)
      
      
      ##### CATEGORICAL

      # Additional step for categorical - find relationship between index and true exposure
      # Then use this to substitute in the index means in the main study and validation study
      mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
      mapping_model <- lm(formula=mapping_formula, data=bootstrap_validation)
      bootstrap_validation$index_mapped = predict(mapping_model, newdata = bootstrap_validation)
      bootstrap_study$index_mapped = predict(mapping_model, newdata = bootstrap_study)
      
      # Find RC Model (relationship between measured exposure and true exposure)
      # Then use this to substitute in the corrected exposure in the main study and validation study
      # However for categorical case exposure_corrected should just be the same as index_mapped
      # This is because lambda is 1
      
      RC_formula = as.formula(exp_true ~ index_mapped)
      RC_model <- lm(formula=RC_formula, data=bootstrap_validation)
      bootstrap_study$exposure_corrected = predict(RC_model, newdata = bootstrap_study)
      
      # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
      # store the estimate and its standard error
      # Note that this is the same as the naive model since lambda is 1 no correction occured
      
      corrected_model <- lm(formula=outcome~exposure_corrected, data=bootstrap_study)
      estimate = corrected_model$coefficients["exposure_corrected"]

      
      boot_est_cat = c(boot_est_cat, estimate)
      
    }
    
    # generate the bootstrap values for association estimate and standard error
    
    estimates_cont[measurement_error_counter-lowerbound+1] = mean(boot_est_cont)
    std_errors_cont[measurement_error_counter-lowerbound+1] = var(boot_est_cont)^0.5
    
    estimates_cat[measurement_error_counter-lowerbound+1] = mean(boot_est_cat)
    std_errors_cat[measurement_error_counter-lowerbound+1] = var(boot_est_cat)^0.5
    print(measurement_error_counter)
  }
  
  
  ## Plot meta analysis for many studies, both cont and cat
  estimates_forRMA = c(estimates_cat[seq(1, length(estimates_cat), 1)],estimates_cont[seq(1, length(estimates_cont), 1)])
  stand_errs = c(std_errors_cat[seq(1, length(std_errors_cat), 1)], std_errors_cont[seq(1, length(std_errors_cont), 1)])

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
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = paste0('Bootstrapped ', num_boots, ' cat and cont with modelled RC outcome error ',outcome_error, ' levels ', numLevels, ' validation ', validation_size ))
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], "Study error, weight, SE", adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
  
}  

dev.off()








#######################################################################################
####### Fixed level Continuous and categorical (index mean subst) with RC #############
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
numLevels = 4
cut_points = c(-Inf,30,40,50,Inf)
#numLevels = 8
#cut_points = c(-Inf,10,20,30,40,50,60,70,Inf)
#numLevels = 16
#cut_points = c(-Inf,0,10,16,22,26,30,35,40,45,50,54,58,64,70,78,Inf)
#numLevels = 32
#cut_points = c(-Inf,0,5,7,10,12,16,19,22,23,26,28,30,32,35,37,40,43,45,47,50,52,54,56,58,61,64,67,70,72,74,78,Inf)
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
svg(filename = paste0('fixed_cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
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
    studyData = createStudyDataFix(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    validation_data = createValidationDataFix(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points)
    
  
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
############# Fixed level Bootstrapped continuous and categorical  ####################
#######################################################################################
# In this section we generate bootstrap values for the estimated association and standard error
# This helps us to validate whether the estimate and SE from using the delta method above
# are correct (this was Ian White's adivce).
#
# In the continuous example, we randomly draw from the validation study and main study, 
# and performing the classical RC process 1000 times 
# (i.e. calculate lambda, correct error-prone measurements, fit corrected model to estimate 
# association). The Delta correction is not applied. Our bootstrapped estimate is the mean
# of the estimates, and our bootstrapped SE is the SD of the estimates.
#
# In the categorical example, we randomly draw from the validation study and main study. 
# In the sampled validation study we calculate the index means and substitute into the sampled
#main study. Then we fit the corrected model. This is done 1000 times. The Delta correction is
# not applied. Our bootstrapped estimate is the mean of the estimates, and our bootstrapped SE
# is the SD of the estimates.


# set up parameters
#numLevels = 2
#cut_points = c(-Inf,40,Inf)
numLevels = 4
cut_points = c(-Inf,30,40,50,Inf)
#numLevels = 8
#cut_points = c(-Inf,10,20,30,40,50,60,70,Inf)
#numLevels = 16
#cut_points = c(-Inf,0,10,16,22,26,30,35,40,45,50,54,58,64,70,78,Inf)
#numLevels = 32
#cut_points = c(-Inf,0,5,7,10,12,16,19,22,23,26,28,30,32,35,37,40,43,45,47,50,52,54,56,58,61,64,67,70,72,74,78,Inf)
validation_size = 400
upperbound = 30
lowerbound = 10
numbaselines = 1
set.seed(12345)
num_boots = 1000

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

# set up plot
svg(filename = paste0('fixed_bootstrapped_',num_boots,'_cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  
  std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  estimates_cat = vector("numeric", length = upperbound - lowerbound)
  
  std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  estimates_cont = vector("numeric", length = upperbound - lowerbound)
  
  for (measurement_error_counter in lowerbound:upperbound) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study - these remain constant for 
    # each bootstrap
    studyData = createStudyDataFix(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    validation_data = createValidationDataFix(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    boot_est_cat = vector()
    boot_est_cont = vector()
    
    for (j in 1:num_boots){
      
      
      # bootstrap sample from the validation study
      bootstrap_validation <- validation_data[sample(1:validation_size, replace = TRUE), ]
      
      # we also bootstrap the main study (as advised by Ian)
      
      bootstrap_study <- studyData[sample(1:study_size, replace = TRUE), ]
      
      ##### CONTINUOUS
      
      # Find RC Model (relationship between measured exposure and true exposure)
      # Then use this to substitute in the corrected exposure in the main study and validation study
      RC_formula = as.formula(exp_true ~ exposure_measured)
      RC_model <- lm(formula=RC_formula, data=bootstrap_validation)
      bootstrap_study$exposure_corrected = predict(RC_model, newdata = bootstrap_study)
      
      # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
      # store the estimate and its standard error
      corrected_model <- lm(formula=outcome~exposure_corrected, data=bootstrap_study)
      estimate = corrected_model$coefficients["exposure_corrected"]
      
      # store estimate
      boot_est_cont = c(boot_est_cont, estimate)
      
      
      ##### CATEGORICAL
      
      # Additional step for categorical - find relationship between index and true exposure
      # Then use this to substitute in the index means in the main study and validation study
      mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
      mapping_model <- lm(formula=mapping_formula, data=bootstrap_validation)
      bootstrap_validation$index_mapped = predict(mapping_model, newdata = bootstrap_validation)
      bootstrap_study$index_mapped = predict(mapping_model, newdata = bootstrap_study)
      
      # Find RC Model (relationship between measured exposure and true exposure)
      # Then use this to substitute in the corrected exposure in the main study and validation study
      # However for categorical case exposure_corrected should just be the same as index_mapped
      # This is because lambda is 1
      
      RC_formula = as.formula(exp_true ~ index_mapped)
      RC_model <- lm(formula=RC_formula, data=bootstrap_validation)
      bootstrap_study$exposure_corrected = predict(RC_model, newdata = bootstrap_study)
      
      # Use mapped values to generate a corrected model - estimate beta using the corrected exposure and outcome
      # store the estimate and its standard error
      # Note that this is the same as the naive model since lambda is 1 no correction occured
      
      corrected_model <- lm(formula=outcome~exposure_corrected, data=bootstrap_study)
      estimate = corrected_model$coefficients["exposure_corrected"]
      
      
      boot_est_cat = c(boot_est_cat, estimate)
      
    }
    
    # generate the bootstrap values for association estimate and standard error
    
    estimates_cont[measurement_error_counter-lowerbound+1] = mean(boot_est_cont)
    std_errors_cont[measurement_error_counter-lowerbound+1] = var(boot_est_cont)^0.5
    
    estimates_cat[measurement_error_counter-lowerbound+1] = mean(boot_est_cat)
    std_errors_cat[measurement_error_counter-lowerbound+1] = var(boot_est_cat)^0.5
    print(measurement_error_counter)
  }
  
  
  ## Plot meta analysis for many studies, both cont and cat
  estimates_forRMA = c(estimates_cat[seq(1, length(estimates_cat), 1)],estimates_cont[seq(1, length(estimates_cont), 1)])
  stand_errs = c(std_errors_cat[seq(1, length(std_errors_cat), 1)], std_errors_cont[seq(1, length(std_errors_cont), 1)])
  
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
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = paste0('Bootstrapped ', num_boots, ' cat and cont with modelled RC outcome error ',outcome_error, ' levels ', numLevels, ' validation ', validation_size ))
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], "Study error, weight, SE", adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
  
}  

dev.off()

