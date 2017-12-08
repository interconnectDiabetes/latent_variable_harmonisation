
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
outcome_error = 0

# Simple creation of the raw data and the related outcome with error
raw_data = data.frame(exp_true = rnorm(n = study_size, mean = raw_mean, sd = raw_stdev))
raw_data$outcome = (trueBeta*raw_data$exp_true) + constant + rnorm(n = study_size, mean = 0, sd = outcome_error)


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
######################## Combined continuous and mean subst with RC #####################
#######################################################################################


numLevels = 16
validation_size = 1600
upperbound = 30
lowerbound = 10
numbaselines = 1
set.seed(12345)

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = paste0('cat_and_cont_with_RC_outcome_error_',outcome_error,'_', numLevels,'_level_',validation_size,'_valid.svg'), width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))


for (i in 1:numbaselines){
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
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    #Get information about indices
    means = rbind(means,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=mean)[,2]))
    sds = rbind(sds,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=sd)[,2]))
    sizes = rbind(sizes,t(aggregate(x = validation_data$exp_true, by= list(validation_data$index), FUN=length)[,2]))
    
    ##### CATEGORICAL
    
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
    
    estimates_cat[measurement_error_counter-lowerbound+1] = estimate
    std_errors_cat[measurement_error_counter-lowerbound+1] = delta_stdError_graph
    lambda_SE_cat[measurement_error_counter-lowerbound+1] = lambda_std_error
    beta_SE_cat[measurement_error_counter-lowerbound+1] = std_error
    
    ######### CONTINUOUS
    
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
    
    
    estimates_cont[measurement_error_counter-lowerbound+1] = estimate_corrected
    std_errors_cont[measurement_error_counter-lowerbound+1] = delta_stdError_graph
    lambda_SE_cont[measurement_error_counter-lowerbound+1] = lambda_stdError
    beta_SE_cont[measurement_error_counter-lowerbound+1] = stdError_uncorrected
    
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
######################## Bootstrapped continuous and cat  #####################
#######################################################################################
# Obviously the same as the means.

numLevels = 32
validation_size = 400
validation_index_size = validation_size/numLevels
study_index_size = study_size/numLevels
upperbound = 30
num_boots = 1000
lowerbound = 10
numbaselines= 1
set.seed(12345)

per_error_weights = list()
per_error_est = list()
per_error_se = list()

svg(filename = 'bootstrap cont & cont error 30 32 levels.svg', width = 10, height = 15)
old.par <- par(mfrow=c(numbaselines^0.5, numbaselines^0.5))

for (i in 1:numbaselines){
  
  std_errors_cat = vector("numeric", length = upperbound - lowerbound)
  estimates_cat = vector("numeric", length = upperbound - lowerbound)

  std_errors_cont = vector("numeric", length = upperbound - lowerbound)
  estimates_cont = vector("numeric", length = upperbound - lowerbound)
  
  for (measurement_error_counter in lowerbound:upperbound) {
    set.seed(i*measurement_error_counter)
    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    
    boot_est_cat = vector()
    boot_est_cont = vector()
    
    for (j in 1:num_boots){
      
      bootstrap_validation <- data.frame(index = as.factor(rep(x = c(1:numLevels), each= validation_index_size)))
      bootstrap_validation[,c('exp_true','exposure_measured')] <- unname(do.call('rbind',lapply(X = split(x=validation_data[,c('exp_true','exposure_measured')], f= as.factor(validation_data$index)), 
                                                            FUN = 
                                                            function(x){
                                                              output = x[sample(1:validation_index_size, replace=TRUE),]
                                                              return (output)
                                                            }
                                                            
                                                            )))
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
      
      # Find Error Model
      mapping_formula = as.formula(exp_true ~ exposure_measured)
      mapping_model <- lm(formula=mapping_formula, data=bootstrap_validation)
      new = data.frame(exposure_measured = bootstrap_study$exposure_measured)
      bootstrap_study$exposure_corrected = predict(mapping_model, newdata = new)
      
      # Use Error Model to Scale Index to Gold
      corrected_model <- lm(formula=outcome~exposure_corrected, data=bootstrap_study)
      estimate = corrected_model$coefficients["exposure_corrected"]
      
      boot_est_cont = c(boot_est_cont, estimate)
      
      
      ##### CATEGORICAL

      
      # Find mapping Model
      mapping_formula = as.formula(exp_true ~ index + 0)
      mapping_model <- lm(formula=mapping_formula, data=bootstrap_validation)
      new_study = data.frame(index = bootstrap_study$index)
      new_validation = data.frame(index = bootstrap_validation$index)
      bootstrap_validation$index_mapped = predict(mapping_model, newdata = new_validation)
      bootstrap_study$index_mapped = predict(mapping_model, newdata = new_study)
      
      # Use mapped values in model
      corrected_model <- lm(formula=outcome~index_mapped, data=bootstrap_study)
      estimate = corrected_model$coefficients["index_mapped"]
      
      boot_est_cat = c(boot_est_cat, estimate)
      
    }
    
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
                           .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Cont & Cat exposures bootstrapped error 30")
  usr <- par("usr")
  text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
  text(usr[1], usr[4], "Study error, weight, SE", adj = c( 0, 1 ),cex=1)
  abline(v = 0.5, col = "lightgray")
  
  
}  

dev.off()

#######################################################################################
######################## comparing SEs from empirical and mdoel  #####################
#######################################################################################
# Obviously the same as the means.

numLevels = 4
validation_size = 400
validation_index_size = validation_size/numLevels
measurement_error = 30
num_boots = 1000
set.seed(30)

    studyData = createStudyData(raw_data = raw_data, measurement_error = measurement_error, number_of_indices = numLevels)
    validation_data = createValidationData(val_size = validation_size, measurement_error = measurement_error, number_of_indices = numLevels )
    
    boot_est = vector()
    set.seed(30)
    for (j in 1:num_boots){
      
      bootstrap_validation <- data.frame(index = as.factor(rep(x = c(1:numLevels), each= validation_index_size)))
      bootstrap_validation[,c('exp_true','exposure_measured')] <- unname(do.call('rbind',lapply(X = split(x=validation_data[,c('exp_true','exposure_measured')], f= as.factor(validation_data$index)), 
                                                                                                FUN = 
                                                                                                  function(x){
                                                                                                    output = x[sample(1:validation_index_size, replace=TRUE),]
                                                                                                    return (output)
                                                                                                  }
                                                                                                
      )))

      mapping_formula = as.formula(exp_true ~ exposure_measured)
      mapping_model <- lm(formula=mapping_formula, data=bootstrap_validation)
      new = data.frame(exposure_measured = studyData$exposure_measured)
      studyData$exposure_corrected = predict(mapping_model, newdata = new)

      corrected_model <- lm(formula=outcome~exposure_corrected, data=studyData)
      estimate = corrected_model$coefficients["exposure_corrected"]
      
      boot_est = c(boot_est, estimate)
    }
    
    estimate_cont_boot = mean(boot_est)
    std_error_cont_boot = var(boot_est)^0.5

    set.seed(30)
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
estimate_cat_boot = mean(boot_est)
std_error_cat_boot = var(boot_est)^0.5

  ##### CATEGORICAL
  
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
  
  #RC model since the mapping model is just a mapping
  
  RC_model = lm(formula=exp_true~index_mapped, data=validation_data)
  lambda = RC_model$coefficients["index_mapped"]
  lambda_std_error = summary(RC_model)$coefficients[2,2]
  
  # recalibrate the standard error using delta function
  
  variance_beta = std_error^2
  beta_lambda_div_sq = (unname(unlist(estimate/(lambda)^2)))^2
  var_lambda = (lambda_std_error)^2
  delta_variance = (variance_beta / (lambda)^2) + (beta_lambda_div_sq * var_lambda)
  delta_stdError_graph = sqrt(delta_variance)
  
  estimate_cat_model = estimate/lambda
  std_error_cat_model = delta_stdError_graph
  lambda_SE_cat_model = lambda_std_error
  beta_SE_cat_model = std_error
  lambda_cat_model = lambda
  
  ######### CONTINUOUS
  
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
  
  estimate_cont_model = estimate_corrected
  std_error_cont_model = delta_stdError_graph
  lambda_SE_cont_model = lambda_stdError
  beta_SE_cont_model = stdError_uncorrected
  lambda_cont_model = lambda
  
  #######################################################################################
  ######################## Bootstrapped direct mapping  #####################
  #######################################################################################
  # Obviously the same as the means.
  
  numLevels = 4
  validation_size = 400
  validation_index_size = validation_size/numLevels
  study_index_size = study_size/numLevels
  upperbound = 30
  num_boots = 1000
  lowerbound = 10
  numbaselines = 1
  set.seed(666)
  
  per_error_weights = list()
  per_error_est = list()
  per_error_se = list()
  
  svg(filename = 'bootstrap cat 1000 multi.svg', width = 4*numbaselines^0.5, height = 4*numbaselines^0.5)
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
      print(paste0(i, " ", measurement_error_counter))
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
           xlab=bquote(paste('Test of Association H'[0]*': true beta association = 0, p = ',
                             .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Pure bootstrap")
    usr <- par("usr")
    text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
    text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
    abline(v = 0.5, col = "lightgray")
    
    
  }  
  
  dev.off()
