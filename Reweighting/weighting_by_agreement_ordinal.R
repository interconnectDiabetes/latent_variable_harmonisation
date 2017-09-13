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

set.seed(123)

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
    raw_data$index = as.ordered(raw_data$index)
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
# General Parameters
numLevels = 4
validation_size = 400

## Graphing Standard Error as a function of Measurement Error when using ordinal independent variable index
upperbound = 100
stdErrorLinear = vector("numeric", length = upperbound)
estimatesLinear = vector("numeric", length = upperbound)
stdErrorQuadratic = vector("numeric", length = upperbound)
estimatesQuadratic = vector("numeric", length = upperbound)
stdErrorCubic = vector("numeric", length = upperbound)
estimatesCubic = vector("numeric", length = upperbound)
for (measureError_status in 1:upperbound) {
    studyData_graph = createStudyData(raw_data = raw_data, measurement_error = measureError_status, number_of_indices = numLevels)
    validation_data_graph = createValidationData(val_size = validation_size, measurement_error = measureError_status, number_of_indices = numLevels )

    lm_graph_linear <- lm(formula=outcome~index, data=studyData_graph)
    estimate_graph_linear = lm_graph_linear$coefficients["index.L"]
    stdError_graph_linear = summary(lm_graph_linear)$coefficients["index.L","Std. Error"]
    stdErrorLinear[measureError_status] = stdError_graph_linear
    estimatesLinear[measureError_status] = estimate_graph_linear

    lm_graph_quadratic <- lm(formula=outcome~index, data=studyData_graph)
    estimate_graph_quadratic = lm_graph_quadratic$coefficients["index.Q"]
    stdError_graph_quadratic = summary(lm_graph_quadratic)$coefficients["index.Q","Std. Error"]
    stdErrorQuadratic[measureError_status] = stdError_graph_quadratic
    estimatesQuadratic[measureError_status] = estimate_graph_quadratic

    lm_graph_cubic <- lm(formula=outcome~index, data=studyData_graph)
    estimate_graph_cubic = lm_graph_cubic$coefficients["index.C"]
    stdError_graph_cubic = summary(lm_graph_cubic)$coefficients["index.C","Std. Error"]
    stdErrorCubic[measureError_status] = stdError_graph_cubic
    estimatesCubic[measureError_status] = estimate_graph_cubic

}

# Measurement vs. Std.Error
plot(x = (1:upperbound), y = stdErrorLinear, xlab = "Measurement Error", ylab = "Standard Error", main = "Assuming Linear Relationship")
# plot(x = (1:upperbound), y = stdErrorQuadratic, xlab = "Measurement Error", ylab = "Standard Error", col = "green", main = "Assuming Quadratic Relationship")
plot(x = (1:upperbound), y = stdErrorCubic, xlab = "Measurement Error", ylab = "Standard Error", col = "red", main = "Assuming Cubic Relationship")

# Measurement vs. Estimates
plot(x = (1:upperbound), y = estimatesLinear, xlab = "Measurement Error", ylab = "Estimates", main = "Assuming Linear Relationship")
# plot(x = (1:upperbound), y = estimatesQuadratic, xlab = "Measurement Error", ylab = "Estimates", col = "green", main = "Assuming Quadratic Relationship")
plot(x = (1:upperbound), y = estimatesCubic, xlab = "Measurement Error", ylab = "Estimates", col = "red", main = "Assuming Cubic Relationship")

# Estimates vs. Std.Error
plot(x = estimatesLinear, y = stdErrorLinear, xlab = "Estimates", ylab = "Standard Error", main = "Non Monotonic Relationship Between \n Accuracy and Standard Error")
# plot(x = estimatesQuadratic, y = stdErrorQuadratic, xlab = "Estimates", ylab = "Standard Error", col = "green", main = "Assuming Quadratic Relationship")
plot(x = estimatesCubic, y = stdErrorCubic, xlab = "Estimates", ylab = "Standard Error", col = "red", main = "Assuming Cubic Relationship")


## Setting up the Studies
# General Parameters
numLevels = 4
validation_size = 400

# Study A
measureError_A = 5
studyData_A = createStudyData(raw_data = raw_data, measurement_error = measureError_A, number_of_indices = numLevels)

# Study B
measureError_B = 10
studyData_B = createStudyData(raw_data = raw_data, measurement_error = measureError_B, number_of_indices = numLevels)

# Study C
measureError_C = 60
studyData_C = createStudyData(raw_data = raw_data, measurement_error = measureError_C, number_of_indices = numLevels)

## Transformation and Regression using Index Means
# Study A
lm_A <- lm(formula=outcome~index, data=studyData_A)
estimate_A = lm_A$coefficients["index.L"]
stdError_A = summary(lm_A)$coefficients["index.L","Std. Error"]

# Study B
lm_B <- lm(formula=outcome~index, data=studyData_B)
estimate_B = lm_B$coefficients["index.L"]
stdError_B = summary(lm_B)$coefficients["index.L","Std. Error"]

# Study C
lm_C <- lm(formula=outcome~index, data=studyData_C)
estimate_C = lm_C$coefficients["index.L"]
stdError_C = summary(lm_C)$coefficients["index.L","Std. Error"]


## Random Effects Model Forest Plot Before Reweighting
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
xlab=bquote(paste('Test of Association'[0.5]*': true beta association = 0.5, p = ',
    .(sprintf("%.3f", round(res$pval,3))))), cex=1, cex.lab=0.75, cex.axis=1, main = "Before Any Changes to Process")
usr <- par("usr")
text(usr[2], usr[4], "Beta [95% CI]", adj = c(1, 4),cex=1)
text(usr[1], usr[4], paste0(gsub(paste0("Study Data","\\$"),"", deparse(fmla)),collapse="\n"), adj = c( 0, 1 ),cex=1)
abline(v = 0.5, col = "lightgray")


##########################################################################################
##########################################################################################
## Same as Basic Process except indexes are transformed into scale of gold measure using the model.

numLevels = 4
validation_size = 400
# Study A
measureError_A = 5
studyData_A = createStudyData(raw_data = raw_data, measurement_error = measureError_A, number_of_indices = numLevels)
validation_data_A = createValidationData(val_size = validation_size, measurement_error = measureError_A, number_of_indices = numLevels )
# Study B
measureError_B = 10
studyData_B = createStudyData(raw_data = raw_data, measurement_error = measureError_B, number_of_indices = numLevels)
validation_data_B = createValidationData(val_size = validation_size, measurement_error = measureError_B, number_of_indices = numLevels )
# Study C
measureError_C = 60
studyData_C = createStudyData(raw_data = raw_data, measurement_error = measureError_C, number_of_indices = numLevels)
validation_data_C = createValidationData(val_size = validation_size, measurement_error = measureError_C, number_of_indices = numLevels )

## Transformation of Indexes to Gold Scale
index_transformation = as.formula(gold ~ index)
# Study A
lambda_lm_A <- lm(formula=index_transformation, data=validation_data_A)
studyData_A$index = as.numeric(studyData_A$index)
new = data.frame(index = studyData_A$index)
studyData_A$index_scaled = predict(lambda_lm_A, newdata = new)
# Study B
lambda_lm_B <- lm(formula=index_transformation, data=validation_data_B)
new = data.frame(index = studyData_B$index)
studyData_B$index_scaled = predict(lambda_lm_B, newdata = new)
# Study C
lambda_lm_C <- lm(formula=index_transformation, data=validation_data_C)
new = data.frame(index = studyData_C$index)
studyData_C$index_scaled = predict(lambda_lm_C, newdata = new)


# Study A
lm_A <- lm(formula=outcome~index, data=studyData_A)
estimate_A = lm_A$coefficients["index.L"]
stdError_A = summary(lm_A)$coefficients["index.L","Std. Error"]

# Study B
lm_B <- lm(formula=outcome~index, data=studyData_B)
estimate_B = lm_B$coefficients["index.L"]
stdError_B = summary(lm_B)$coefficients["index.L","Std. Error"]

# Study C
lm_C <- lm(formula=outcome~index, data=studyData_C)
estimate_C = lm_C$coefficients["index.L"]
stdError_C = summary(lm_C)$coefficients["index.L","Std. Error"]

























##########################################################################################
##########################################################################################
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
