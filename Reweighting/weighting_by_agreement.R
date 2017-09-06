# Investigation into Data Agreement as Weighting System
## Author: Paul Scherer
##         Tom Bishop
## Institute: University of Cambridge MRC Epidemiology Unit 
## Date: 4.09.2017

# Title: Battling Regression Dilution with Spear(man's) Correlation.

# # Setting up experiment.
# 1. Define 'true' relationship to generate the data.
# 2. Generate base and validation data. (Problem Point 1 for implementation and assumptions). With different measurement errors. More random measurement error in the exposure causes regression dilution (more wrong estimates), and surer standard errors destroying all that is good with meta analyses. 
#   - Doing the meta analysis with the heavily diluted studies will give them the higher weighting in a REM, for all things controlled.

# # Method
# 1. Calculate agreement of the ordinal base variable and gold standard in the validation study. Using Pearson, Spearman, Whitney Mann or Kendall-Tau
# 2. Use the agreement as a rating of how good the study is. Divide number of study participants as part of the whole and normalize this value to rest of group.
# 3. Use normalized value of 

# # Results
# There is no way to judge except by eye or some automation

# # Conclusion
# Meta Analyses in particular individual level (repooled analyses) and federated harmonisation exercises should be done in this new way to battle regression dilution via bad measurement.





# Set up the Data

## We pull raw gold standard data from a normal distribution representing the amount of Physical Activity a person does
## A study attempts to record this with a device that has a set amount of random measurement error, and then
## puts this measurement into one of n indices of exposure, this is done for sample of raw data.

## the gold standard data has a set true linear relationship with some continous "outcome". We know that the backtransformation
## of index values to the level of the gold standard variable introduces regression dilution due to the random measurement error
## introduced at the measurement stage. We investigate, all things kept constant, whether changing the number of n
## indices introduces more regression dilution to the correclation coefficient between the backtransformed activity level.

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
    # raw_data$index = kmeans(x = raw_data$exposure_measured, centers = number_of_indices, iter.max = 100)
    raw_data = raw_data[with(raw_data, order(index)),]
    return(raw_data)
}

createValidationData <- function(val_size, measurement_error, number_of_indices) {
    validation_data = data.frame(gold =  rnorm(n = val_size, mean = raw_mean, sd = raw_stdev))
    exposure_error = rnorm(n = val_size, mean = 0, sd = measurement_error)
    validation_data$exposure_measured = validation_data$gold + exposure_error
    validation_data$index = cut_number(x = validation_data$exposure_measured,n = number_of_indices, labels =FALSE)
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
########################### Motivation of the Problem #################################
#######################################################################################
## Setting up the Studies
# General Parameters
numLevels = 4
validation_size = 400


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

## Graphing Standard Error as a function of Measurement Error
upperbound = 80
stdErrorNuggets = vector("numeric", length = upperbound)
estimatesNuggets = vector("numeric", length = upperbound)
for (measureError_Nuggets in 1:upperbound) {
    studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measureError_Nuggets, number_of_indices = numLevels)
    validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measureError_Nuggets, number_of_indices = numLevels )
    validation_means_graph = createMeansList(validation_data_graph, numLevels)

    studyData_graph$ind_mean <- unlist(lapply(X=studyData_graph$index, FUN=function(index_val){
        output =  validation_means_graph[index_val]
    }))
    lm_graph <- lm(formula=outcome~ind_mean, data=studyData_graph)
    estimate_graph = lm_graph$coefficients["ind_mean"]
    stdError_graph = summary(lm_graph)$coefficients["ind_mean","Std. Error"]
    stdErrorNuggets[measureError_Nuggets] = stdError_graph
    estimatesNuggets[measureError_Nuggets] = estimate_graph
}
plot(x = (1:upperbound), y = stdErrorNuggets, xlab = "measurement_error", ylab = "stderr")
plot(x = (1:upperbound), y = estimatesNuggets, xlab = "measurement_error", ylab = "estimates")

plot(x = (1:upperbound), y = stdErrorNuggets, xlab = "measurement_error", ylab = "stderr", col = "red")
points(x = (1:upperbound), y = estimatesNuggets, xlab = "measurement_error", ylab = "estimates", col = "green")


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
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)


## Calculating Tau coefficient for weighting

## Reweighting and Weight Renormalizing.

## Random Effects Model Forest Plot After Reweighting
