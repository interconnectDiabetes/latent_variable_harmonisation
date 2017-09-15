# Set up the Data
library(ggplot2)
library(reshape2)
library(parallel)
library(stats)
library(metafor)
library(Amelia)

# Raw Data Set Properties
trueBeta = 0.5
constant = 20

study_size = 20000
raw_mean = 40
raw_stdev = 15

# Simple creation of the raw data and the related outcome without error
raw_data = data.frame(gold = rnorm(n = study_size, mean = raw_mean, sd = raw_stdev))
raw_data$outcome = (trueBeta*raw_data$gold) + constant
###############################################################################
########################### Functions #########################################
###############################################################################
createStudyData <- function(raw_data, measurement_error, number_of_indices){
  exposure_error = rnorm(n = study_size, mean = 0, sd = measurement_error)
  raw_data$exposure_measured = raw_data$gold + exposure_error
  raw_data$index = cut_number(x = raw_data$exposure_measured,n = number_of_indices, labels =FALSE)
  raw_data$index = as.ordered(raw_data$index)
  raw_data = raw_data[with(raw_data, order(index)),]
  return(raw_data)
}

createValidationData <- function(val_size, measurement_error, number_of_indices) {
  validation_data = data.frame(gold =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
  exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
  validation_data$exposure_measured = validation_data$gold + exposure_error
  validation_data$index = cut_number(x = validation_data$exposure_measured,n = number_of_indices, labels =FALSE)
  validation_data$index = as.ordered(validation_data$index)
  return (validation_data)
}

createMeansList <- function(data, number_of_indices) {
  meansList = vector(mode="list", length = number_of_indices)
  for (i in 1:number_of_indices) {
    meansList[i] = mean(unname(unlist((split(x=data$gold, f= as.factor(data$index)))[i])))
  }
  return(meansList)
}



#######################################################################################
######################## Crazy Amelia ############################
#######################################################################################
## Graphing Standard Error as a function of Measurement Error
error_upperbound = 200
std_error = vector("numeric", length = error_upperbound)
estimates = vector("numeric", length = error_upperbound)
imputations = 5
numLevels = 4
validation_size = 400

for (measurement_error_counter in 1:error_upperbound){
  studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
  validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
  study_data_cp = studyData_graph[,c(2,4)]
  study_data_cp$exp = NA
  
  validation_data_cp = validation_data_graph[,c(1,3)]
  names(validation_data_cp) = c("exp", "index")
    validation_data_cp$outcome = NA

  
  combined = rbind(study_data_cp,validation_data_cp)
  
  a.out <- amelia(combined, m=imputations, noms = 'index')

  imp_est = vector()
  imp_se = vector()
  for (i in 1:imputations){
    
    results = lm(formula = outcome ~ exp, data = a.out$imputations[[i]])
    output = summary(results)
    imp_est = c(imp_est, output$coefficients[2,1])
    imp_se = c(imp_se, output$coefficients[2,2])
  }

  rubin_est = sum(imp_est)/imputations
  W = sum(imp_se)/imputations
  B = sum((imp_est-rubin_est)^2)/(imputations-1)
  rubin_var = W + (1+1/imputations)*B
  rubin_se = rubin_var^0.5
  
  
  std_error[measurement_error_counter] = rubin_se
  estimates[measurement_error_counter] = rubin_est
}
plot (x = (1:error_upperbound), y = std_error, xlab = "Measurement Error", ylab = "Standard Error", main = "Standard Error based on Measurement Error \n Regression Calibration")
dev.copy(png,'myplot5.png')
dev.off()
plot (x = (1:error_upperbound), y = estimates, xlab = "Measurement Error", ylab = "estimates", main = "Estimates (Accuracy based on Measurement Error) \n Regression Calibration")
dev.copy(png,'myplot6.png')
dev.off()
plot (x = estimates, y = std_error, xlab = "Estimates", ylab = "Standard Error", main = "Non Monotonic Relationship Between \n Accuracy and Standard Error \n Regression Calibration")
dev.copy(png,'myplot7.png')
dev.off()


## Random Effects Model Forest Plot After Regression Modeling
estimates_REMA = estimates[seq(1, length(estimates), 100)]
#estimates_REMA = estimates[1:10]
stand_errs = std_error[seq(1, length(std_error), 100)]
#stand_errs = std_error[1:10]
labels = as.character(c(1:length(estimates_REMA)))
res <- rma(yi = estimates_REMA, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~harmonised_x)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                              .(sprintf("%.3f", round(res$QEp,3))),')')),
       xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0.5, p = ',
                         .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Using Predict aka 'Error Model'")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")
dev.copy(png,'myplot12.png')
dev.off()



