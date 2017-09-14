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
    raw_data$index = as.factor(raw_data$index)
    raw_data = raw_data[with(raw_data, order(index)),]
    return(raw_data)
}

createValidationData <- function(val_size, measurement_error, number_of_indices) {
    validation_data = data.frame(gold =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
    exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
    validation_data$exposure_measured = validation_data$gold + exposure_error
    validation_data$index = cut_number(x = validation_data$exposure_measured,n = number_of_indices, labels =FALSE)
    validation_data$index = as.factor(validation_data$index)
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
######################## Before Any Changes to the Process (MEANS) ####################
#######################################################################################
# General Parameters
numLevels = 4
validation_size = 400
upperbound = 50
## Graphing Standard Error as a function of Measurement Error
std_error = vector("numeric", length = upperbound)
estimates_graph = vector("numeric", length = upperbound)
for (measurement_error_counter in 1:upperbound) {
    studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    validation_means_graph = createMeansList(validation_data_graph, numLevels)

    studyData_graph$ind_mean <- unlist(lapply(X=studyData_graph$index, FUN=function(index_val){
        output =  validation_means_graph[index_val]
    }))
    lm_graph <- lm(formula=outcome~ind_mean, data=studyData_graph)
    estimate_graph = lm_graph$coefficients["ind_mean"]
    stdError_graph = summary(lm_graph)$coefficients["ind_mean","Std. Error"]

    estimates_graph[measurement_error_counter] = estimate_graph
    std_error[measurement_error_counter] = stdError_graph

}
plot(x = (1:upperbound), y = std_error, xlab = "measurement_error", ylab = "stderr", main = "Using Index Means")
plot(x = (1:upperbound), y = estimates_graph, xlab = "measurement_error", ylab = "estimates", main = "Using Index Means")
plot(x = estimates_graph, y = std_error, xlab = "estimates", ylab = "stderror", main = "Using Index Means")


# Study A
measureError_A = 5
studyData_A = createStudyData(raw_data = raw_data, measurement_error = measureError_A, number_of_indices = numLevels)
validation_data_A = createValidationData(val_size = validation_size, measurement_error = measureError_A, number_of_indices = numLevels )
validation_means_A = createMeansList(validation_data_A, numLevels)

# Study B
measureError_B = 10
studyData_B = createStudyData(raw_data = raw_data, measurement_error = measureError_B, number_of_indices = numLevels)
validation_data_B = createValidationData(val_size = validation_size, measurement_error = measureError_B, number_of_indices = numLevels )
validation_means_B = createMeansList(validation_data_B, numLevels)

# Study C
measureError_C = 15
studyData_C = createStudyData(raw_data = raw_data, measurement_error = measureError_C, number_of_indices = numLevels)
validation_data_C = createValidationData(val_size = validation_size, measurement_error = measureError_C, number_of_indices = numLevels )
validation_means_C = createMeansList(validation_data_C, numLevels)

## Transformation and Regression using Index Means
# Study A
studyData_A$ind_mean <- unlist(lapply(X=studyData_A$index, FUN=function(index_val){
    output =  validation_means_A[index_val]
}))
lm_A <- lm(formula=outcome~ind_mean, data=studyData_A)
estimate_A = lm_A$coefficients["ind_mean"]
stdError_A = summary(lm_A)$coefficients["ind_mean","Std. Error"]

# Study B
studyData_B$ind_mean <- unlist(lapply(X=studyData_B$index, FUN=function(index_val){
    output =  validation_means_B[index_val]
}))
lm_B <- lm(formula=outcome~ind_mean, data=studyData_B)
estimate_B = lm_B$coefficients["ind_mean"]
stdError_B = summary(lm_B)$coefficients["ind_mean","Std. Error"]

# Study C
studyData_C$ind_mean <- unlist(lapply(X=studyData_C$index, FUN=function(index_val){
    output =  validation_means_C[index_val]
}))
lm_C <- lm(formula=outcome~ind_mean, data=studyData_C)
estimate_C = lm_C$coefficients["ind_mean"]
stdError_C = summary(lm_C)$coefficients["ind_mean","Std. Error"]

## Random Effects Model Forest Plot Before Reweighting
estimates = c(estimate_A, estimate_B, estimate_C)
stand_errs = c(stdError_A, stdError_B, stdError_C)
labels = c("A","B","C")
res <- rma(yi = estimates, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(outcome~ind_mean)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0.5, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1)
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1, main = "Using Index Means")
abline(v = 0.5, col = "lightgray")






#######################################################################################
######################## Using Error Model Esque Scale it to Gold #####################
#######################################################################################
# Obviously the same as the means.

numLevels = 4
validation_size = 400

## Plotting Bit
## Graphing Standard Error as a function of Measurement Error
std_error = vector("numeric", length = upperbound)
estimates_graph = vector("numeric", length = upperbound)
for (measurement_error_counter in 1:upperbound) {
    studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )

    # Find Error Model
    index_transformation = as.formula(gold ~ index)
    lambda_lm_graph <- lm(formula=index_transformation, data=validation_data_graph)
    new = data.frame(index = studyData_graph$index)
    studyData_graph$index_scaled = predict(lambda_lm_graph, newdata = new)

    # Use Error Model to Scale Index to Gold
    lm_graph <- lm(formula=outcome~index_scaled, data=studyData_graph)
    estimate_graph = lm_graph$coefficients["index_scaled"]
    stdError_graph = summary(lm_graph)$coefficients["index_scaled","Std. Error"]

    estimates_graph[measurement_error_counter] = estimate_graph
    std_error[measurement_error_counter] = stdError_graph
}
plot(x = (1:upperbound), y = std_error, xlab = "measurement_error", ylab = "stderr", main = "Using Error Model Pred")
plot(x = (1:upperbound), y = estimates_graph, xlab = "measurement_error", ylab = "estimates", main = "Using Error Model Pred")
plot(x = estimates_graph, y = std_error, xlab = "estimates", ylab = "stderror", main = "Using Error Model Pred")



# # Study A
# measureError_A = 5
# studyData_A = createStudyData(raw_data = raw_data, measurement_error = measureError_A, number_of_indices = numLevels)
# validation_data_A = createValidationData(val_size = validation_size, measurement_error = measureError_A, number_of_indices = numLevels )
# validation_means_A = createMeansList(validation_data_A, numLevels)

# # Study B
# measureError_B = 10
# studyData_B = createStudyData(raw_data = raw_data, measurement_error = measureError_B, number_of_indices = numLevels)
# validation_data_B = createValidationData(val_size = validation_size, measurement_error = measureError_B, number_of_indices = numLevels )
# validation_means_B = createMeansList(validation_data_B, numLevels)

# # Study C
# measureError_C = 60
# studyData_C = createStudyData(raw_data = raw_data, measurement_error = measureError_C, number_of_indices = numLevels)
# validation_data_C = createValidationData(val_size = validation_size, measurement_error = measureError_C, number_of_indices = numLevels )
# validation_means_C = createMeansList(validation_data_C, numLevels)


## Transformation of Indexes to Gold Scale
index_transformation = as.formula(gold ~ index)
# Study A
lambda_lm_A <- lm(formula=index_transformation, data=validation_data_A)
new = data.frame(index = studyData_A$index)
studyData_A$index_scaled = predict(lambda_lm_A, newdata = new)

studyData_A$ind_mean <- unlist(lapply(X=studyData_A$index, FUN=function(index_val){
    output =  validation_means_A[index_val]
}))

# Study B
lambda_lm_B <- lm(formula=index_transformation, data=validation_data_B)
new = data.frame(index = studyData_B$index)
studyData_B$index_scaled = predict(lambda_lm_B, newdata = new)

studyData_B$ind_mean <- unlist(lapply(X=studyData_B$index, FUN=function(index_val){
    output =  validation_means_B[index_val]
}))

# Study C
lambda_lm_C <- lm(formula=index_transformation, data=validation_data_C)
new = data.frame(index = studyData_C$index)
studyData_C$index_scaled = predict(lambda_lm_C, newdata = new)

studyData_C$ind_mean <- unlist(lapply(X=studyData_C$index, FUN=function(index_val){
    output =  validation_means_C[index_val]
}))


# Study A
lm_A <- lm(formula=outcome~index_scaled, data=studyData_A)
estimate_A = lm_A$coefficients["index_scaled"]
stdError_A = summary(lm_A)$coefficients["index_scaled","Std. Error"]

# Study B
lm_B <- lm(formula=outcome~index_scaled, data=studyData_B)
estimate_B = lm_B$coefficients["index_scaled"]
stdError_B = summary(lm_B)$coefficients["index_scaled","Std. Error"]

# Study C
lm_C <- lm(formula=outcome~index_scaled, data=studyData_C)
estimate_C = lm_C$coefficients["index_scaled"]
stdError_C = summary(lm_C)$coefficients["index_scaled","Std. Error"]

## Random Effects Model Forest Plot Before Reweighting
estimates = c(estimate_A, estimate_B, estimate_C)
stand_errs = c(stdError_A, stdError_B, stdError_C)
labels = c("A","B","C")
res <- rma(yi = estimates, sei = stand_errs, method='DL', slab = labels)
weights_res <- weights.rma.uni(res)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(outcome~index_scaled)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0.5, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1)
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1, main = "Using Error Model of Gold~Index")
abline(v = 0.5, col = "lightgray")



# plot(x = studyData_A$outcome, y = studyData_A$ind_mean, xlab = "outcome", ylab = "indexGolded", main = "Comparing Means and Pred A", col = "red")
# points(x = studyData_A$outcome, y = studyData_A$index_scaled, col = "green")

# plot(x = studyData_B$outcome, y = studyData_B$ind_mean, xlab = "outcome", ylab = "indexGolded", main = "Comparing Means and Pred B", col = "red")
# points(x = studyData_B$outcome, y = studyData_B$index_scaled, col = "green")

# plot(x = studyData_C$outcome, y = studyData_C$ind_mean, xlab = "outcome", ylab = "indexGolded", main = "Comparing Means and Pred C", col = "red")
# points(x = studyData_C$outcome, y = studyData_C$index_scaled, col = "green")

#######################################################################################
######################## Using Index Means Esque Scale it to Gold #####################
############################## Also Regression Calibration ############################
#######################################################################################
# Obviously the same as the means.

numLevels = 4
validation_size = 400

## Plotting Bit
## Graphing Standard Error as a function of Measurement Error
std_error = vector("numeric", length = upperbound)
estimates_graph = vector("numeric", length = upperbound)
for (measurement_error_counter in 1:upperbound) {
    studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measurement_error_counter, number_of_indices = numLevels)
    validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measurement_error_counter, number_of_indices = numLevels )
    validation_means_graph = createMeansList(validation_data_graph, numLevels)

    studyData_graph$ind_mean <- unlist(lapply(X=studyData_graph$index, FUN=function(index_val){
        output =  validation_means_graph[index_val]
    }))

    validation_data_graph$ind_mean <- unlist(lapply(X=validation_data_graph$index, FUN=function(index_val){
        output =  validation_means_graph[index_val]
    }))

    # calculate lambda
    lambda_lm_graph = lm(formula = gold~ind_mean, data = validation_data_graph)
    lambda_graph = lambda_lm_graph$coefficients["ind_mean"]
    lambda_stdError_graph = summary(lambda_lm_graph)$coefficients["ind_mean","Std. Error"]

    # calculate beta hat
    lm_graph <- lm(formula=outcome~ind_mean, data=studyData_graph)
    estimate_graph = lm_graph$coefficients["ind_mean"]
    stdError_graph = summary(lm_graph)$coefficients["ind_mean","Std. Error"]

    # recalibrate the standard error using delta function
    variance_beta = (sqrt(validation_size) * summary(lm_graph)$coefficients["ind_mean","Std. Error"])^2
    #variance_beta = summary(lm_graph)$sigma**2
    # variance_beta = (summary(lm_graph)$coefficients["ind_mean","Std. Error"])^2
    lambda_pure = unlist(unname(lambda_lm_graph$coefficients["ind_mean"]))
    beta_lambda_div_sq = (unname(unlist(lm_graph$coefficients["ind_mean"]/(lambda_lm_graph$coefficients["ind_mean"])^2)))^2
    var_lambda = (sqrt(validation_size) * summary(lambda_lm_graph)$coefficients["ind_mean","Std. Error"])^2
    delta_variance = (variance_beta / (lambda_pure)^2) + (beta_lambda_div_sq * var_lambda)
    delta_stdError_graph = sqrt(delta_variance)/sqrt(validation_size) 
    delta_stdError_graph = sqrt(delta_variance)

    estimates_graph[measurement_error_counter] = estimate_graph/lambda_graph
    std_error[measurement_error_counter] = delta_stdError_graph
}

plot(x = (1:upperbound), y = std_error, xlab = "measurement_error", ylab = "stderr", main = "Using Index Means Regression Calibration")
plot(x = (1:upperbound), y = estimates_graph, xlab = "measurement_error", ylab = "estimates", main = "Using Index Means Regression Calibration")
plot(x = estimates_graph, y = std_error, xlab = "estimates", ylab = "stderror", main = "Using Index Means Regression Calibration")


## Transformation and Regression using Index Means
# Study A
studyData_A$ind_mean <- unlist(lapply(X=studyData_A$index, FUN=function(index_val){
    output =  validation_means_A[index_val]
}))
lm_A <- lm(formula=outcome~ind_mean, data=studyData_A)
estimate_A = lm_A$coefficients["ind_mean"]
stdError_A = summary(lm_A)$coefficients["ind_mean","Std. Error"]

# Study B
studyData_B$ind_mean <- unlist(lapply(X=studyData_B$index, FUN=function(index_val){
    output =  validation_means_B[index_val]
}))
lm_B <- lm(formula=outcome~ind_mean, data=studyData_B)
estimate_B = lm_B$coefficients["ind_mean"]
stdError_B = summary(lm_B)$coefficients["ind_mean","Std. Error"]

# Study C
studyData_C$ind_mean <- unlist(lapply(X=studyData_C$index, FUN=function(index_val){
    output =  validation_means_C[index_val]
}))
lm_C <- lm(formula=outcome~ind_mean, data=studyData_C)
estimate_C = lm_C$coefficients["ind_mean"]
stdError_C = summary(lm_C)$coefficients["ind_mean","Std. Error"]


## Calculate the lambdas
# Study A
# all have study indmean already
lambda_lm_A = lm(formula = gold~ind_mean, data = validation_data_A)
lambda_A = lambda_lm_A$coefficients["ind_mean"]

# calculate the corrected standard error
var_Beta = (sqrt(validation_size) * summary(lm_A)$coefficients["ind_mean","Std. Error"])^2
#var_Beta = summary(lm_A)$sigma**2
# var_Beta = (summary(lm_A)$coefficients["ind_mean","Std. Error"])^2
lambda_pure = unlist(unname(lambda_lm_A$coefficients["ind_mean"]))
beta_lambda_div_sq = (unname(unlist(lm_A$coefficients["ind_mean"]/(lambda_lm_A$coefficients["ind_mean"])^2)))^2
var_lambda = (sqrt(validation_size) * summary(lambda_lm_A)$coefficients["ind_mean","Std. Error"])^2
delta_variance = (var_Beta / (lambda_pure)^2) + (beta_lambda_div_sq * var_lambda)
delta_stdError_A = sqrt(delta_variance)/sqrt(validation_size) 
delta_stdError_A = sqrt(delta_variance)


# Study B
lambda_lm_B = lm(formula = gold~ind_mean, data = validation_data_B)
lambda_B = lambda_lm_B$coefficients["ind_mean"]

var_lambda = (sqrt(validation_size) * summary(lambda_lm_B)$coefficients["ind_mean","Std. Error"])^2
var_Beta = (sqrt(validation_size) * summary(lm_B)$coefficients["ind_mean","Std. Error"])^2
#var_Beta = summary(lm_A)$sigma**2
# var_Beta = (summary(lm_A)$coefficients["ind_mean","Std. Error"])^2
lambda_pure = unlist(unname(lambda_lm_B$coefficients["ind_mean"]))
beta_lambda_div_sq = (unname(unlist(lm_B$coefficients["ind_mean"]/(lambda_lm_B$coefficients["ind_mean"])^2)))^2
delta_variance = (var_Beta / (lambda_pure)^2) + beta_lambda_div_sq * var_lambda
delta_stdError_B = sqrt(delta_variance)/sqrt(validation_size) 
delta_stdError_B = sqrt(delta_variance)


# Study C
lambda_lm_C = lm(formula = gold~ind_mean, data = validation_data_C)
lambda_C = lambda_lm_C$coefficients["ind_mean"]

var_lambda = (sqrt(validation_size) * summary(lambda_lm_C)$coefficients["ind_mean","Std. Error"])^2
var_Beta = (sqrt(validation_size) * summary(lm_C)$coefficients["ind_mean","Std. Error"])^2
#var_Beta = summary(lm_A)$sigma**2
# var_Beta = (summary(lm_A)$coefficients["ind_mean","Std. Error"])^2
lambda_pure = unlist(unname(lambda_lm_C$coefficients["ind_mean"]))
beta_lambda_div_sq = (unname(unlist(lm_C$coefficients["ind_mean"]/(lambda_lm_C$coefficients["ind_mean"])^2)))^2
delta_variance = (var_Beta / (lambda_pure)^2) + beta_lambda_div_sq * var_lambda
delta_stdError_C = sqrt(delta_variance)/sqrt(validation_size) 
delta_stdError_C = sqrt(delta_variance)

corrected_estimate_A = estimate_A / lambda_A
corrected_estimate_B = estimate_B / lambda_B
corrected_estimate_C = estimate_C / lambda_C


## Random Effects Model Forest Plot After Regression Calibration
estimates = c(corrected_estimate_A, corrected_estimate_B, corrected_estimate_C)
stand_errs = c(delta_stdError_A, delta_stdError_B, delta_stdError_C)
labels = c("A","B","C")
res <- rma(yi = estimates, sei = stand_errs, method='DL', slab = labels)

# Forest Plot
res$slab <- paste(res$slab, " (", round(weights.rma.uni(res),digits=1), "%)")
fmla = as.formula(y~ind_mean)
forest(res, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
    .(sprintf("%.3f", round(res$QEp,3))),')')),
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0.5, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "After Regression Calibration")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")