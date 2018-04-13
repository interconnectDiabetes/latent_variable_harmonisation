
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
####### Baselines         ########
#######################################################################################
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
error_vect = c(10,15,20,25,30)
#error_vect = c(10:30)
numbaselines = 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('baselines_for_pub_', numLevels,'_level_',validation_size,'_valid.svg'), width = 6, height = 8)
old.par <- par(mfrow=c(2, 1))
par(mar=c(4.1,4.1,2.1,2.1))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  # structures for storage

  std_errors_cat = numeric()
  estimates_cat = numeric()

  std_errors_cont = numeric()
  estimates_cont = numeric()
  
  ######### CONTINUOUS
  
  #for (measurement_error_counter in lowerbound:upperbound) {
  for (measurement_error_counter in error_vect) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    #naive model - outcome with measured exposure
    naive_model <- lm(formula=outcome~exposure_measured, data=studyData)
    estimate_naive = naive_model$coefficients["exposure_measured"]
    std_error_naive = summary(naive_model)$coefficients["exposure_measured","Std. Error"]
    
    # Store estimate and corrected std error
    estimates_cont = c(estimates_cont, estimate_naive)
    std_errors_cont = c(std_errors_cont, std_error_naive)

  }
  ## Plot meta analysis for cont 
  
  labels = paste0( error_vect)
  
  
  forest(x = estimates_cont, sei = std_errors_cont,
         xlab='Estimate of association', cex=1, cex.lab=1, cex.axis=1, 
         refline = 0.5, slab = labels, psize= c(1,1,1,1,1),
         alim=c(0.05,0.55), xlim = c(-0.1,0.85))
  usr <- par("usr")
  text(usr[2], usr[4], "Estimate [95% CI]", adj = c(1, 6),cex=1)
  text(usr[1], usr[4], "Instrument error size", adj = c( 0, 6 ),cex=1)
  title(main = paste0('Continuous instr. study size ', study_size, ' levels ', numLevels, ' validation ', validation_size ), line = -0.5, adj = 1)
  #abline(v = 0.5, col = "lightgray")
  
  
  
  ##### CATEGORICAL
  
  for (measurement_error_counter in error_vect) {  
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points)
    

    
    # Do the mapping using data from validation sample A, and apply to A, B and study
    
    mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
    mapping_model <- lm(formula=mapping_formula, data=validation_data)
    validation_data$index_mapped = predict(mapping_model, newdata = validation_data)
    studyData$index_mapped = predict(mapping_model, newdata = studyData)
    
    #get mapped exposure to outcome relationship (naive model)
    
    mapped_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate_mapped = mapped_model$coefficients["index_mapped"]
    std_error_mapped = summary(mapped_model)$coefficients["index_mapped","Std. Error"]
    
    
    # Store estimate and corrected std error
    
    estimates_cat = c(estimates_cat, estimate_mapped)
    std_errors_cat = c(std_errors_cat, std_error_mapped)
    
  }
  
  ## Plot  
  
  labels = paste0( error_vect)


  forest(x = estimates_cat, sei = std_errors_cat,
         xlab='Estimate of association', cex=1, cex.lab=1, cex.axis=1, 
         refline = 0.5, slab = labels, psize= c(1,1,1,1,1),
         alim=c(0.05,0.55), xlim = c(-0.1,0.85))
  usr <- par("usr")
  text(usr[2], usr[4], "Estimate [95% CI]", adj = c(1, 6),cex=1)
  text(usr[1], usr[4], "Instrument error size", adj = c( 0, 6 ),cex=1)
  title(main = paste0('Categorical instr. study size ', study_size, ' levels ', numLevels, ' validation ', validation_size ), line = -0.5, adj = 1)
  
}


dev.off()

#######################################################################################
####### new notation Continuous and split measurement categorical         ########
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
#error_vect = c(10,15,20,25,30)
error_vect = c(10,15,20,25,30)
numbaselines = 5
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('split_cat_and_cont_for_pub_', numLevels,'_level_',validation_size,'_valid.svg'), width = 6, height = 8)
old.par <- par(mfrow=c(2, 1))
par(mar=c(4.1,4.1,2.1,2.1))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1

estimates_all = data.frame()
std_error_all = data.frame()

for (i in 1:numbaselines){

  
  std_errors_cat = numeric()
  estimates_cat = numeric()
  lambda_SE_cat = numeric()
  beta_SE_cat = numeric()
  std_errors_cont = numeric()
  estimates_cont = numeric()
  lambda_SE_cont = numeric()
  beta_SE_cont = numeric()
  
  lambda_cat = numeric()
  
  means = data.frame()
  sds = data.frame()
  sizes = data.frame()
  
  
  ######### CONTINUOUS
  
  #for (measurement_error_counter in lowerbound:upperbound) {
  for (measurement_error_counter in error_vect) {
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    #Store information about indices
    #means = rbind(means,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=mean)[,2]))
    #sds = rbind(sds,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=sd)[,2]))
    #sizes = rbind(sizes,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=length)[,2]))
    
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
    estimates_cont = c(estimates_cont, estimate)
    std_errors_cont = c(std_errors_cont, delta_std_error)
    
    # Store lambda SE and beta SE in case we need them (not used)
    lambda_SE_cont = c(lambda_SE_cont, lambda_std_error)
    beta_SE_cont = c(beta_SE_cont, std_error_naive)
    
  }
  estimates_all = rbind(estimates_all, estimates_cont)
  std_error_all = rbind(std_error_all, std_errors_cont)
}

    
  ## Plot  
  
  labels = paste0( error_vect)
  
  
  forest(x = estimates_cont, sei = std_errors_cont,
         xlab='Estimate of association', cex=1, cex.lab=1, cex.axis=1, 
         refline = 0.5, slab = labels, psize= c(1,1,1,1,1),
         alim=c(0.3,0.7), xlim = c(0.2,0.9))
  usr <- par("usr")
  text(usr[2], usr[4], "Estimate [95% CI]", adj = c(1, 6),cex=1)
  text(usr[1], usr[4], "Instrument error size", adj = c( 0, 6 ),cex=1)
  title(main = paste0('Continuous instr. study size ', study_size, ' levels ', numLevels, ' validation ', validation_size ), line = -0.5, adj = 0.5)
  
estimates_all = data.frame()
std_error_all = data.frame()
  
for (i in 1:numbaselines){
    
    
    std_errors_cat = numeric()
    estimates_cat = numeric()
    lambda_SE_cat = numeric()
    beta_SE_cat = numeric()
    std_errors_cont = numeric()
    estimates_cont = numeric()
    lambda_SE_cont = numeric()
    beta_SE_cont = numeric()
    
    lambda_cat = numeric()
    
    means = data.frame()
    sds = data.frame()
    sizes = data.frame()  
  
  ##### CATEGORICAL
  
  for (measurement_error_counter in error_vect) {  
    
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
    
    
    # Do the mapping using data from validation sample A, and apply to A, B and study
    
    mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
    mapping_model <- lm(formula=mapping_formula, data=validation_data_A)
    validation_data_A$index_mapped = predict(mapping_model, newdata = validation_data_A)
    validation_data_B$index_mapped = predict(mapping_model, newdata = validation_data_B)
    studyData$index_mapped = predict(mapping_model, newdata = studyData)
    
    #get mapped exposure to outcome relationship (naive model)
    
    mapped_model <- lm(formula=outcome~index_mapped, data=studyData)
    estimate_mapped = mapped_model$coefficients["index_mapped"]
    std_error_mapped = summary(mapped_model)$coefficients["index_mapped","Std. Error"]
    
    #error in mapping (no calibration, "lambda" is 1) - use validation A
    
    mapping_error_formula = as.formula(exp_true ~ index_mapped)
    mapping_error_model <- lm(formula=mapping_error_formula, data=validation_data_A)
    estimate_mapping = mapping_error_model$coefficients["index_mapped"] # should be 1
    std_error_mapping = summary(mapping_error_model)$coefficients["index_mapped","Std. Error"]
    
    #calibration & associated error - use validation B and apply calibration to study
    
    calib_formula = as.formula(exp_true ~ index_mapped)
    calib_model <- lm(formula=calib_formula, data=validation_data_B)
    studyData$exposure_corrected = predict(calib_model, newdata = studyData)
    estimate_calib = calib_model$coefficients["index_mapped"]
    std_error_calib = summary(calib_model)$coefficients["index_mapped","Std. Error"]
    
    # Use mapped and calibrated values to generate a corrected model - estimate beta using the corrected exposure and outcome
    # store the estimate and its standard error
    
    corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
    estimate_corrected = corrected_model$coefficients["exposure_corrected"]
    std_error_corrected = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # refine the standard error of the main model estimate using delta function
    
    variance_beta = std_error_mapped^2
    beta_lambda_div_sq = (unname(unlist(estimate_mapped/(estimate_calib)^2)))^2
    var_lambda = (std_error_calib)^2
    delta_variance = (variance_beta / (estimate_calib)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error_corrected = sqrt(delta_variance)
    
    
    # Store estimate and corrected std error
    
    estimates_cat = c(estimates_cat, estimate_corrected)
    std_errors_cat = c(std_errors_cat, delta_std_error_corrected)
    
    # Store lambda SE and beta SE in case we need them (not used)
    
    #lambda_SE_cat[measurement_error_counter-lowerbound+1] = lambda_std_error
    #beta_SE_cat[measurement_error_counter-lowerbound+1] = std_error
    
    lambda_cat = c(lambda_cat, estimate_calib)
    
  }
  estimates_all = data.frame()
  std_error_all = data.frame()
}
  
  
  ## Plot  
  
  labels = paste0( error_vect)
  
  
  forest(x = estimates_cat, sei = std_errors_cat,
         xlab='Estimate of association', cex=1, cex.lab=1, cex.axis=1, 
         refline = 0.5, slab = labels, psize= c(1,1,1,1,1),
         alim=c(0.3,0.7), xlim = c(0.2,0.9))
  usr <- par("usr")
  text(usr[2], usr[4], "Estimate [95% CI]", adj = c(1, 6),cex=1)
  text(usr[1], usr[4], "Instrument error size", adj = c( 0, 6 ),cex=1)
  title(main = paste0('Categorical instr. study size ', study_size, ' levels ', numLevels, ' validation ', validation_size ), line = -0.5, adj = 0.5)
  



dev.off()


#######################################################################################
####### new notation REPEAT measurement categorical         ########
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
#numLevels = 32
#cut_points = c(-Inf,0,5,7,10,12,16,19,22,23,26,28,30,32,35,37,40,43,45,47,50,52,54,56,58,61,64,67,70,72,74,78,Inf)
validation_size = 200
upperbound = 30
lowerbound = 10
#error_vect = c(10,15,20,25,30)
error_vect = c(10,15,20,25,30)
numbaselines = 1
set.seed(12345)

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('repeat_cat_for_pub_', numLevels,'_level_',validation_size,'_valid.svg'), width = 6, height = 4)

par(mar=c(4.1,4.1,2.1,2.1))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  # structures for storage
  # std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  # estimates_cat = vector("numeric", length = upperbound - lowerbound)
  # lambda_SE_cat = vector("numeric", length = upperbound - lowerbound)
  # beta_SE_cat = vector("numeric", length = upperbound - lowerbound)
  # std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  # estimates_cont = vector("numeric", length = upperbound - lowerbound)
  # lambda_SE_cont = vector("numeric", length = upperbound - lowerbound)
  # beta_SE_cont = vector("numeric", length = upperbound - lowerbound)
  
  std_errors_cat = numeric()
  estimates_cat = numeric()
  lambda_SE_cat = numeric()
  beta_SE_cat = numeric()
  
  lambda_cat = numeric()
  
  means = data.frame()
  sds = data.frame()
  sizes = data.frame()
  
  
  ##### CATEGORICAL
  
  for (measurement_error_counter in error_vect) {  
    
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
    estimate_corrected = corrected_model$coefficients["exposure_corrected"]
    std_error = summary(corrected_model)$coefficients["exposure_corrected","Std. Error"]
    
    # Extract the coefficients from the RC model for correcting estimate SE
    
    lambda = summary(RC_model)$coefficients["index_mapped","Estimate"]
    lambda_std_error = summary(RC_model)$coefficients["index_mapped","Std. Error"]
    
    # refine the standard error of the main model estimate using delta function
    
    variance_beta = std_error_naive^2
    beta_lambda_div_sq = (unname(unlist(estimate_naive/(lambda)^2)))^2
    var_lambda = (lambda_std_error)^2
    delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
    delta_std_error_corrected = sqrt(delta_variance)
    
    # Store estimate and corrected std error
    
    estimates_cat = c(estimates_cat, estimate_corrected)
    std_errors_cat = c(std_errors_cat, delta_std_error_corrected)
    
    
  }
  
  ## Plot  
  
  labels = paste0(error_vect)
  
  
  forest(x = estimates_cat, sei = std_errors_cat,
         xlab='Estimate of association', cex=1, cex.lab=1, cex.axis=1, 
         refline = 0.5, slab = labels, psize= c(1,1,1,1,1),
         alim=c(0.3,0.7), xlim = c(0.2,0.9))
  usr <- par("usr")
  text(usr[2], usr[4], "Estimate [95% CI]", adj = c(1, 6),cex=1)
  text(usr[1], usr[4], "Instrument error size", adj = c( 0, 6 ),cex=1)
  title(main = paste0('Categorical instr. study size ', study_size, ' levels ', numLevels, ' validation ', validation_size ), line = -0.5, adj = 0.5)
  
}


dev.off()


#######################################################################################
####### BOOTSTRAP REPEAT measurement categorical         ########
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
#numLevels = 32
#cut_points = c(-Inf,0,5,7,10,12,16,19,22,23,26,28,30,32,35,37,40,43,45,47,50,52,54,56,58,61,64,67,70,72,74,78,Inf)
validation_size = 400
upperbound = 30
lowerbound = 10
#error_vect = c(10,15,20,25,30)
error_vect = c(10,15,20,25,30)
numbaselines = 1
set.seed(12345)
num_boots = 2000

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('repeat_boot_cat_for_pub_', numLevels,'_level_',validation_size,'_valid.svg'), width = 6, height = 4)

par(mar=c(4.1,4.1,2.1,2.1))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){

  std_errors_cat = numeric()
  estimates_cat = numeric()

  
  ##### CATEGORICAL
  
  for (measurement_error_counter in error_vect) {  
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    # generate 2 validation studies of half the intended size to simulate splitting the validation study in half
    validation_data = createValidationDataRep(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points )
    
    boot_est_cat = vector()
    
    for (j in 1:num_boots){
      
      
      bootstrap_validation <- validation_data[sample(1:validation_size, replace = TRUE), ]
      
      # we also bootstrap the main study (as advised by Ian)
      
      bootstrap_study <- studyData[sample(1:study_size, replace = TRUE), ]
      
      
      validation_data_A = bootstrap_validation[,c('exp_true', 'index')]
      validation_data_B = bootstrap_validation[,c('exp_true', 'index2')]
      colnames(validation_data_B) = c('exp_true', 'index')
      
      
      mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
      mapping_model <- lm(formula=mapping_formula, data=validation_data_A)
      validation_data_A$index_mapped = predict(mapping_model, newdata = validation_data_A)
      validation_data_B$index_mapped = predict(mapping_model, newdata = validation_data_B)
      bootstrap_study$index_mapped = predict(mapping_model, newdata = bootstrap_study)
      
      #calibration & associated error - use validation B and apply calibration to study
      
      calib_formula = as.formula(exp_true ~ index_mapped)
      calib_model <- lm(formula=calib_formula, data=validation_data_B)
      bootstrap_study$exposure_corrected = predict(calib_model, newdata = bootstrap_study)
      
      # Use mapped and calibrated values to generate a corrected model - estimate beta using the corrected exposure and outcome
      # store the estimate and its standard error
      
      corrected_model <- lm(formula=outcome~exposure_corrected, data=bootstrap_study)
      estimate_corrected = corrected_model$coefficients["exposure_corrected"]
      
      boot_est_cat = c(boot_est_cat, estimate_corrected)
    }
    
    estimates_cat = c(estimates_cat, mean(boot_est_cat))
    std_errors_cat = c(std_errors_cat, var(boot_est_cat)^0.5)
    print(measurement_error_counter)
    
  }
  
  ## Plot  
  
  labels = paste0(error_vect)
  
  
  forest(x = estimates_cat, sei = std_errors_cat,
         xlab='Estimate of association', cex=1, cex.lab=1, cex.axis=1, 
         refline = 0.5, slab = labels, psize= c(1,1,1,1,1),
         alim=c(0.3,0.7), xlim = c(0.2,0.9))
  usr <- par("usr")
  text(usr[2], usr[4], "Estimate [95% CI]", adj = c(1, 6),cex=1)
  text(usr[1], usr[4], "Instrument error size", adj = c( 0, 6 ),cex=1)
  title(main = paste0('Categorical instr. study size ', study_size, ' levels ', numLevels, ' validation ', validation_size ), line = -0.5, adj = 0.5)
  
}


dev.off()


#######################################################################################
####### BOOTSTRAP SPLIT measurement categorical         ########
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
#error_vect = c(10,15,20,25,30)
error_vect = c(10,15,20,25,30)
numbaselines = 1
set.seed(12345)
num_boots = 1000

#structures to store results
per_error_weights = list()
per_error_est = list()
per_error_se = list()

#set up plot
svg(filename = paste0('split_boot_cat_for_pub_', numLevels,'_level_',validation_size,'_valid.svg'), width = 6, height = 4)

par(mar=c(4.1,4.1,2.1,2.1))

# numbaselines parameter allows multiple runs of the whole loop - in case there are concerns that set of results
# are just due to the initial seed chosen. Most of the time just set it to 1
for (i in 1:numbaselines){
  
  std_errors_cat = numeric()
  estimates_cat = numeric()
  
  
  ##### CATEGORICAL
  
  for (measurement_error_counter in error_vect) {  
    
    # Seed is set for each measurement error to ensure consistency when generating
    # main study and validation study
    
    set.seed(i*measurement_error_counter)
    
    # generate the main study and validation study
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, cut_points = cut_points)
    
    # generate 2 validation studies of half the intended size to simulate splitting the validation study in half
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, cut_points = cut_points )
    
    validation_data$split = unlist(lapply(X = split(x = validation_data, validation_data$index),FUN = function(X){
      y1 = rep(x=1,times=ceiling(length(X$index)/2))
      y2 = rep(x=2,times=length(X$index)-ceiling(length(X$index)/2))
      y = c(y1,y2)
      return(y)
    }))
    
    boot_est_cat = vector()
    
    for (j in 1:num_boots){
      
      
      bootstrap_validation <- validation_data[sample(1:validation_size, replace = TRUE), ]
      
      # we also bootstrap the main study (as advised by Ian)
      
      bootstrap_study <- studyData[sample(1:study_size, replace = TRUE), ]
      
      
      validation_data_A = bootstrap_validation[bootstrap_validation$split == 1,]
      validation_data_B = bootstrap_validation[bootstrap_validation$split == 2,]
      
      
      
      mapping_formula = as.formula(exp_true ~ index + 0) # 0 forces no intercept
      mapping_model <- lm(formula=mapping_formula, data=validation_data_A)
      validation_data_A$index_mapped = predict(mapping_model, newdata = validation_data_A)
      validation_data_B$index_mapped = predict(mapping_model, newdata = validation_data_B)
      bootstrap_study$index_mapped = predict(mapping_model, newdata = bootstrap_study)
      
      #calibration & associated error - use validation B and apply calibration to study
      
      calib_formula = as.formula(exp_true ~ index_mapped)
      calib_model <- lm(formula=calib_formula, data=validation_data_B)
      bootstrap_study$exposure_corrected = predict(calib_model, newdata = bootstrap_study)
      
      # Use mapped and calibrated values to generate a corrected model - estimate beta using the corrected exposure and outcome
      # store the estimate and its standard error
      
      corrected_model <- lm(formula=outcome~exposure_corrected, data=bootstrap_study)
      estimate_corrected = corrected_model$coefficients["exposure_corrected"]
      
      boot_est_cat = c(boot_est_cat, estimate_corrected)
    }
    
    estimates_cat = c(estimates_cat, mean(boot_est_cat))
    std_errors_cat = c(std_errors_cat, var(boot_est_cat)^0.5)
    print(measurement_error_counter)
    
  }
  
  ## Plot  
  
  labels = paste0(error_vect)
  
  
  forest(x = estimates_cat, sei = std_errors_cat,
         xlab='Estimate of association', cex=1, cex.lab=1, cex.axis=1, 
         refline = 0.5, slab = labels, psize= c(1,1,1,1,1),
         alim=c(0.3,0.7), xlim = c(0.2,0.9))
  usr <- par("usr")
  text(usr[2], usr[4], "Estimate [95% CI]", adj = c(1, 6),cex=1)
  text(usr[1], usr[4], "Instrument error size", adj = c( 0, 6 ),cex=1)
  title(main = paste0('Categorical instr. study size ', study_size, ' levels ', numLevels, ' validation ', validation_size ), line = -0.5, adj = 0.5)
  
}


dev.off()

######
# test splitting

numLevels = 4
cut_points = c(-Inf,30,40,50,Inf)
#numLevels = 8
#cut_points = c(-Inf,10,20,30,40,50,60,70,Inf)
#numLevels = 16
#cut_points = c(-Inf,0,10,16,22,26,30,35,40,45,50,54,58,64,70,78,Inf)
#numLevels = 32
#cut_points = c(-Inf,0,5,7,10,12,16,19,22,23,26,28,30,32,35,37,40,43,45,47,50,52,54,56,58,61,64,67,70,72,74,78,Inf)
validation_size = 800


measurement_error_counter = 15

A_frame = data.frame()
B_frame = data.frame()

for (i in 1:10) {  
  
  # Seed is set for each measurement error to ensure consistency when generating
  # main study and validation study
  
  set.seed(i)

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
  
  tempA = aggregate(validation_data_A$exp_true, by=list(validation_data_A$index), mean, na.rm = TRUE)
  tempB = aggregate(validation_data_B$exp_true, by=list(validation_data_B$index), mean, na.rm = TRUE)
  
  A_frame = rbind(A_frame, tempA[,2])
  B_frame = rbind(B_frame, tempB[,2])
}