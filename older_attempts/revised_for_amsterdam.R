## Title: Harmonisation of PAEE

## Papers to find comparison results:
## 1. Physical activity reduces the risk of incident type 2 diabetes in 
## general and in abdominally lean and obese men and women: the EPIC-InterAct 
## Study
## 2. Validity of a short questionnaire to assess physical activity in 10 European
## countries 
## Author: Paul Scherer
## Date: 10/11/2016
## University of Cambridge

###############################################################################
############################### SET UP ########################################
###############################################################################
# Working Directory
setwd("V:/Studies/InterConnect/Internal/Latent variable harmonisation")

# Libraries
library("survival")
library("graphics")
library("metafor")

# fitting of distributions via MLE
library("fitdistrplus")


###############################################################################
###############################################################################
###############################################################################

#  _   _         _  _      _         _    _                ______        _         
# | | | |       | |(_)    | |       | |  (_)               |  _  \      | |        
# | | | |  __ _ | | _   __| |  __ _ | |_  _   ___   _ __   | | | | __ _ | |_  __ _ 
# | | | | / _` || || | / _` | / _` || __|| | / _ \ | '_ \  | | | |/ _` || __|/ _` |
# \ \_/ /| (_| || || || (_| || (_| || |_ | || (_) || | | | | |/ /| (_| || |_| (_| |
#  \___/  \__,_||_||_| \__,_| \__,_| \__||_| \___/ |_| |_| |___/  \__,_| \__|\__,_|

# Read the actiheart and EPIC study data (sweden UMEA is missing in 
# main file, needs to be loaded seperately)
actiheart_summary <- read.csv("PHIA0000112014_IA85_12Mar/PAQIA0000112014_actiheart_summary.csv", header=T)
epic <- read.csv("PHIA0000112014_IA85_12Mar/PAQIA0000112014_epic.csv", header=T)

#UMEA file is separate, see details later
epic_umea <- read.csv("V:/Studies/InterConnect/Internal/Latent variable harmonisation/epic_umea.csv", header=T)

# actiheart - 0 = women, 1 = men, based on mean weights, assume men higher:
# aggregate(actiheart_summary$weight, by=list(sex=actiheart_summary$sex), mean, na.rm = TRUE)

# Create a list where actiheart results are split into each participant
act_sorted <- actiheart_summary[order(actiheart_summary$universal_id, actiheart_summary$visit), ]
act_split <- split(x = act_sorted, f = act_sorted$universal_id)

# create the TRUE PAEE values column
## Create a new data frame using the act_split list with PAEE values
output <- do.call(rbind,(lapply(act_split,function(x){
  # PAEE logic is to use Branch4 and then Branch7
  # Only use values where they have worn it for more than 3 days
  # Weight PAEE in proportion to wear duration
  universal_id <- x$universal_id[1]
  country <- x$country[1]
  age <- mean(x$age, na.rm=TRUE)
  height <- mean(x$height, na.rm=TRUE)
  weight <- mean(x$weight, na.rm=TRUE)
  sex <- x$sex[1]
  
  if (is.na(x$PAEE_Branch4[1])) {
    if (is.na(x$PAEE_Branch7[1])){
      PAEE_1 <- 0
    }
    else {
      PAEE_1 <- x$PAEE_Branch7[1]
    }
  }
  else {
    PAEE_1 <- x$PAEE_Branch4[1]
  }
  
  if (is.na(x$PAEE_Branch4[2])) {
    if (is.na(x$PAEE_Branch7[2])){
      PAEE_2 <- 0
    }
    else {
      PAEE_2 <- x$PAEE_Branch7[2]
    }
  }
  else {
    PAEE_2 <- x$PAEE_Branch4[2]
  }
  
  Pwear_1 <- x$Pwear[1]
  Pwear_2 <- x$Pwear[2]
  
  if (!is.na(Pwear_1) && !is.na(Pwear_2)) {
    if (Pwear_1 >= 96 && Pwear_2 >= 96){
      
      PAEE <- (PAEE_1  + PAEE_2) / 2
      
    } else if (Pwear_1 >= 24 & Pwear_2 >= 24) {
      
      PAEE <- (PAEE_1 * Pwear_1 + PAEE_2 * Pwear_2)/(Pwear_1 + Pwear_2)
      
    } else if (Pwear_1 >= 24 & Pwear_2 < 24) {
      
      PAEE <- PAEE_1
      
    } else if (Pwear_2 >= 24 & Pwear_1 < 24) {
      
      PAEE <- PAEE_2
    } else {
      PAEE <- NA
    }
    
  } else if (!is.na(Pwear_1)) {
    
    if (Pwear_1 >= 24) {
      PAEE <- PAEE_1
    } else {
      PAEE <- NA
    }
    
  } else if (!is.na(Pwear_2)) {
    
    if (Pwear_2 >= 24) {
      PAEE <- PAEE_2
    } else {
      PAEE <- NA
    }
    
  } else {
    PAEE <- NA
  }
  
  
  out <- data.frame(universal_id, country, height, weight, sex, age, PAEE, Pwear_1, Pwear_2, PAEE_1, PAEE_2)
  return(out)
})
))
row.names(output)<-NULL

###############################################################################
###############################################################################
###############################################################################
## In this section we regenerate cam_index column ie check logic - seems to match up

# fix employment status

epic$pa_workini <- as.vector(do.call(rbind,lapply(X = epic$pa_workini, FUN = function(x){
  if ( x == 9){
    ind = 5
  }
  else { 
    ind = x
  }
  return(ind)
})
))

# camMets_ind for cam_matrix
# totalMets_index calculation

epic$camMets_ind <- as.vector(do.call(rbind, lapply(X=epic$pa_lt, FUN = function(x){
  if(is.na(x)){
    ind = NA
  }
  else if (x == 0){
    ind = 1
  }
  else if (x <= 3.5) {
    ind = 2
  }
  else if (x <= 7) {
    ind = 3
  }
  else if (x > 7) {
    ind = 4
  } 
  else {
    ind = NA
  }
  return(ind)
})
))

#unemployed or missing not set to sedentary to match TP paper
cam_matrix = matrix(byrow = TRUE,
                     c(1, 2, 3, 4, 
                       2, 3, 4, 4, 
                       3, 4, 4, 4, 
                       4, 4, 4, 4,
                       1, 2, 3, 4), 
                     nrow=5, 
                     ncol=4)

epic$cam_index2 <- apply(X = epic[,c('pa_workini', 'camMets_ind')], MARGIN = 1, FUN = function(x){
  output <- cam_matrix[x[1], x[2]]
  return(output)
})

#if test has zero length then we are happy with the original value of cam_index!
test <- epic[epic$cam_index!=epic$cam_index2,]

# so we trim down the dataset ready to merge with sweden
epic_trim <- data.frame(universal_id = epic$universal_id, pa_workini = epic$pa_workini, pa_lt = epic$pa_lt)

                        

###############################################################################
###############to do: Add in UMEA        #####################################
###############################################################################
# add Umea , if possible (not sure if all indices can be generated)
# 
# the csv was generated from:
# V:\Studies\InterAct_WP_2.6\StudyData\HANDOVER TO DM\links\umea_idmerge.dta
# V:\Studies\InterAct_WP_2.6\UMEÅ_EPIC_PAQ\UMEA_EPIC_final.xls
# with logic from V:\Studies\InterAct_WP_2.6\UMEÅ_EPIC_PAQ\Description_Sweden_index_TP.doc
# to understand the meaning of the g2 and g6 variables

epic_umea_new <- data.frame(universal_id = epic_umea$universal_id)
#epic_umea_new$country <- as.factor('SWEDEN')
epic_umea_new$pa_workini <- apply(X = epic_umea[,c('g2a', 'g2b', 'g2c', 'g2d', 'g2e')], MARGIN = 1, FUN = function(x){                                    
  if (!is.na(x[1])){
    x_ind = 1
  }
  else if (!is.na(x[2])){
    x_ind = 2
  }
  else if (!is.na(x[3])){
    x_ind = 2
  }
  else if (!is.na(x[4])){
    x_ind = 3
  }
  else if (!is.na(x[5])){
    x_ind = 4
  }
  else {
    x_ind = 1
  }
  return(x_ind)
})

# coding from original interact analysis - doesn't look right

# epic_umea_new$pa_lt <- sapply(X = epic_umea[,'g6'],FUN = function(x){       
#   if (is.na(x)){
#     x_ind = 0
#   }
#   else if (x == 1){
#     x_ind = 3.5
#   }
#   else if (x == 2){
#     x_ind = 3.5
#   }
#   else if (x == 3){
#     x_ind = 5
#   }
#   else if (x == 4){
#     x_ind = 8
#   }
#   else if (x == 5){
#     x_ind = 4
#   }
#   else {
#     x_ind = 0
#   }
#   return(x_ind)
# })

epic_umea_new$pa_lt <- sapply(X = epic_umea[,'g6'],FUN = function(x){
  if (is.na(x)){
    x_ind = 0
  }
  else if (x == 1){
    x_ind = 0
  }
  else if (x == 2){
    x_ind = 0
  }
  else if (x == 3){
    x_ind = 3.5
  }
  else if (x == 4){
    x_ind = 5
  }
  else if (x == 5){
    x_ind = 8
  }
  else {
    x_ind = 0
  }
  return(x_ind)
})

# join Sweden to other epic data
epic_final <- rbind(epic_umea_new, epic_trim)

### validation_data creation
# merge or perform a "natural join" on the trimmed epic table and 
# cleaned actiheart table by universal id
validation_data <- merge(output, epic_final, by = 'universal_id')


#calculate cam_index
validation_data$camMets_ind <- as.vector(do.call(rbind, lapply(X=validation_data$pa_lt, FUN = function(x){
  if(is.na(x)){
    ind = NA
  }
  else if (x == 0){
    ind = 1
  }
  else if (x <= 3.5) {
    ind = 2
  }
  else if (x <= 7) {
    ind = 3
  }
  else if (x > 7) {
    ind = 4
  } 
  else {
    ind = NA
  }
  return(ind)
})
))

#unemployed or missing all not set to sedentary to match TP paper
cam_matrix = matrix( byrow = TRUE,
                     c(1, 2, 3, 4, 
                       2, 3, 4, 4, 
                       3, 4, 4, 4, 
                       4, 4, 4, 4,
                       1, 2, 3, 4), 
                     nrow=5, 
                     ncol=4)

validation_data$cam_index <- apply(X = validation_data[,c('pa_workini', 'camMets_ind')], MARGIN = 1, FUN = function(x){
  output <- cam_matrix[x[1], x[2]]
  return(output)
})


# calculate the BMI using the weight and height and include this in the table to make
# a model for PAEE
# bmi_calc <- function(weight, height){
#   bmi = (weight/height)/height
#   return(bmi)
# }
# validation_data$bmi <- mapply(FUN=bmi_calc, validation_data$weight, validation_data$height)
# 



###############################################################################
###############################################################################
# Results 2, 3


#record means etc for the cam index and binary index
# Note that these match the Peters et al paper! Hurrah!
cam_index_means = aggregate(validation_data$PAEE, by=list(cam_index=validation_data$cam_index, sex=validation_data$sex), mean, na.rm = TRUE)
colnames(cam_index_means)[3] = 'cam_index_mean'
cam_index_medians = aggregate(validation_data$PAEE, by=list(cam_index=validation_data$cam_index, sex=validation_data$sex), median, na.rm = TRUE)
colnames(cam_index_medians)[3] = 'cam_index_median'
cam_index_counts =aggregate(validation_data$PAEE, by=list(cam_index=validation_data$cam_index, sex=validation_data$sex), length)
colnames(cam_index_counts)[3] = 'cam_index_count'

# now define other indices
# Cam index is Method A

# Method B - ie sedentary (1) (cam index 1) or active (2) (cam index 2,3,4)
validation_data$binary_index <- unlist(lapply(validation_data$cam_index, FUN=function(x){
  if (x[1] == 1) {output = 1 }
  else {output = 2}
  return(output)
}))

# record results
binary_index_means = aggregate(validation_data$PAEE, by=list(binary_index=validation_data$binary_index, sex=validation_data$sex), mean, na.rm = TRUE)
colnames(binary_index_means)[3] = 'binary_index_mean'
binary_index_medians = aggregate(validation_data$PAEE, by=list(binary_index=validation_data$binary_index, sex=validation_data$sex), median, na.rm = TRUE)
colnames(binary_index_medians)[3] = 'binary_index_median'
binary_index_counts =aggregate(validation_data$PAEE, by=list(binary_index=validation_data$binary_index, sex=validation_data$sex), length)
colnames(binary_index_counts)[3] = 'binary_index_count'

###############################################################################
###############################################################################
# Method C - regression based (but still essentially 16 categories) on components of cam index
#create factors for regression
validation_data$camMets_ind_fact <- as.factor(validation_data$camMets_ind)
validation_data$pa_workini_fact <- as.factor(validation_data$pa_workini)

validation_data_women <- subset(validation_data, sex==0)
validation_data_men <- subset(validation_data, sex==1)

reg_out_men <- glm(formula=PAEE ~ pa_workini_fact + camMets_ind_fact , data = validation_data_men, family='gaussian')
reg_out_women <- glm(formula=PAEE ~ pa_workini_fact + camMets_ind_fact , data = validation_data_women, family='gaussian')

#store coefficients
#coeffs <- summary(reg_out)$coefficients[, 1]
coeffs_men <- summary(reg_out_men)$coefficients[, 1]
coeffs_women <- summary(reg_out_women)$coefficients[, 1]

###############################################################################
###############################################################################
# Result 3
# calculate the lambdas (to be used much later!)

# Setting the 'observed' PAEE values by the mean

# cam index - do a 'vlookup'
temp <- merge(validation_data, cam_index_means, by=c('cam_index', 'sex'))[,c('universal_id', 'sex', 
        'country', 'PAEE', 'camMets_ind',  'cam_index', 'binary_index',
        'camMets_ind_fact', 'pa_workini_fact', 'cam_index_mean')]
new_data <- temp[order(temp$universal_id),]

# binary index - 'vlookup' again
temp <- merge(new_data, binary_index_means, by=c('binary_index', 'sex'))
new_data <- temp[order(temp$universal_id),]

# regression based
# convert work index into binaries to make use of regression coeffs
temp <- data.frame(sapply(levels(new_data$pa_workini_fact), function(x) as.integer(x == new_data$pa_workini_fact)))
colnames(temp) <- c('pa_workini_1','pa_workini_2','pa_workini_3','pa_workini_4','pa_workini_5')
new_data <- cbind(new_data, temp)

# convert ltpa index into binaries
temp <- data.frame(sapply(levels(new_data$camMets_ind_fact), function(x) as.integer(x == new_data$camMets_ind_fact)))
colnames(temp) <- c('camMets_ind_1','camMets_ind_2','camMets_ind_3','camMets_ind_4')
new_data <- cbind(new_data, temp)
rm(temp)

# use regression coefficients to estimate PAEE
final_output_women <- subset(new_data, sex==0)
final_output_men <- subset(new_data, sex==1)

final_output_women$reg_value <- coeffs_women['(Intercept)'] + final_output_women$pa_workini_2 * coeffs_women['pa_workini_fact2'] +
  final_output_women$pa_workini_3 * coeffs_women['pa_workini_fact3'] + final_output_women$pa_workini_4 * coeffs_women['pa_workini_fact4'] +
  final_output_women$pa_workini_5 * coeffs_women['pa_workini_fact5'] +
  final_output_women$camMets_ind_2 * coeffs_women['camMets_ind_fact2'] + final_output_women$camMets_ind_3 * coeffs_women['camMets_ind_fact3'] +
  final_output_women$camMets_ind_4 * coeffs_women['camMets_ind_fact4']

final_output_men$reg_value <- coeffs_men['(Intercept)'] + final_output_men$pa_workini_2 * coeffs_men['pa_workini_fact2'] +
  final_output_men$pa_workini_3 * coeffs_men['pa_workini_fact3'] + final_output_men$pa_workini_4 * coeffs_men['pa_workini_fact4'] +
  final_output_men$pa_workini_5 * coeffs_men['pa_workini_fact5'] +
  final_output_men$camMets_ind_2 * coeffs_men['camMets_ind_fact2'] + final_output_men$camMets_ind_3 * coeffs_men['camMets_ind_fact3'] +
  final_output_men$camMets_ind_4 * coeffs_men['camMets_ind_fact4']

# Now do the lambda calculations - these give values of 1!!! as we originally predicted...
rdr_regression_fit_binary_men <- lm(formula=PAEE~binary_index_mean, data=final_output_men)
rdr_regression_fit_binary_women <- lm(formula=PAEE~binary_index_mean, data=final_output_women)
rdr_regression_fit_cam_men <- lm(formula=PAEE~cam_index_mean, data=final_output_men)
rdr_regression_fit_cam_women <- lm(formula=PAEE~cam_index_mean, data=final_output_women)
rdr_regression_fit_reg_men <- lm(formula=PAEE~reg_value, data=final_output_men)
rdr_regression_fit_reg_women <- lm(formula=PAEE~reg_value, data=final_output_women)


png( "lambda_cam_men.png")
plot(x = final_output_men$cam_index_mean, y= final_output_men$PAEE, xlim=c(30,60))
abline(0,1)
dev.off()

png( "lambda_reg_men.png")
plot(x = final_output_men$reg_value, y= final_output_men$PAEE, xlim=c(30,60))
abline(0,1)
dev.off()

png( "lambda_bin_men.png")
plot(x = final_output_men$binary_index_mean, y= final_output_men$PAEE, xlim=c(30,60))
abline(0,1)
dev.off()

png( "lambda_cam_women.png")
plot(x = final_output_women$cam_index_mean, y= final_output_women$PAEE, xlim=c(30,60))
abline(0,1)
dev.off()

png( "lambda_reg_women.png")
plot(x = final_output_women$reg_value, y= final_output_women$PAEE, xlim=c(30,60))
abline(0,1)
dev.off()

png( "lambda_bin_women.png")
plot(x = final_output_women$binary_index_mean, y= final_output_women$PAEE, xlim=c(30,60))
abline(0,1)
dev.off()


# histograms
WDforPlots = "V:/Studies/InterConnect/Internal/Latent variable harmonisation/paee_harmonisation/plots"
setwd(WDforPlots)

counts = table(final_output_men$cam_index)
png( "bar_plot_men_cam.png")
barplot(counts, main="Cam Index Distributions for Men", xlab="Cambridge Index", col=c("darkred", "red", "orange", "green"), beside=TRUE)
dev.off()

counts = table(final_output_men$binary_index)
png( "bar_plot_men_bin.png")
barplot(counts, main="Binary Index Distributions for Men", xlab="Binary Index", col=c("darkred", "red"), beside=TRUE)
dev.off()


# barplot(counts, main="LTPA/OPA Distributions for Men", xlab="LTPA Category", col=c("darkred", "red", "orange", "green", "blue"), beside=TRUE)
# barplot(counts[,1], main="OPA Distributions for Men LTPA=1", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
# barplot(counts[,2], main="OPA Distributions for Men LTPA=2", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
# barplot(counts[,3], main="OPA Distributions for Men LTPA=3", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
# barplot(counts[,4], main="OPA Distributions for Men LTPA=4", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))

counts = table(final_output_men$pa_workini_fact,final_output_men$camMets_ind_fact)
png( "bar_plot_men_reg.png")
barplot(counts, main="LTPA/OPA Distributions for Men", xlab="LTPA Category", col=c("darkred", "red", "orange", "green", "blue"), beside=TRUE)
dev.off()
png( "bar_plot_men_reg_1.png")
barplot(counts[,1], main="OPA Distributions for Men LTPA=1", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()
png( "bar_plot_men_reg_2.png")
barplot(counts[,2], main="OPA Distributions for Men LTPA=2", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()
png( "bar_plot_men_reg_3.png")
barplot(counts[,3], main="OPA Distributions for Men LTPA=3", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()
png( "bar_plot_men_reg_4.png")
barplot(counts[,4], main="OPA Distributions for Men LTPA=4", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()

counts = table(final_output_women$cam_index)
png( "bar_plot_women_cam.png")
barplot(counts, main="Cam Index Distributions for Women", xlab="Cambridge Index", col=c("darkred", "red", "orange", "green"), beside=TRUE)
dev.off()

counts = table(final_output_women$binary_index)
png( "bar_plot_women_bin.png")
barplot(counts, main="Binary Index Distributions for Women", xlab="Binary Index", col=c("darkred", "red"), beside=TRUE)
dev.off()

# barplot(counts, main="LTPA/OPA Distributions for Women", xlab="LTPA Category", col=c("darkred", "red", "orange", "green", "blue"), beside=TRUE)
# barplot(counts[,1], main="OPA Distributions for women LTPA=1", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
# barplot(counts[,2], main="OPA Distributions for women LTPA=2", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
# barplot(counts[,3], main="OPA Distributions for women LTPA=3", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
# barplot(counts[,4], main="OPA Distributions for women LTPA=4", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
counts = table(final_output_women$pa_workini_fact,final_output_women$camMets_ind_fact)
png( "bar_plot_women_reg.png")
barplot(counts, main="LTPA/OPA Distributions for Women", xlab="LTPA Category", col=c("darkred", "red", "orange", "green", "blue"), beside=TRUE)
dev.off()
png( "bar_plot_women_reg_1.png")
barplot(counts[,1], main="OPA Distributions for women LTPA=1", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()
png( "bar_plot_women_reg_2.png")
barplot(counts[,2], main="OPA Distributions for women LTPA=2", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()
png( "bar_plot_women_reg_3.png")
barplot(counts[,3], main="OPA Distributions for women LTPA=3", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()
png( "bar_plot_women_reg_4.png")
barplot(counts[,4], main="OPA Distributions for women LTPA=4", xlab="OPA Category", col=c("darkred", "red", "orange", "green", "blue"))
dev.off()

# means bar graphs


setwd("V:/Studies/InterConnect/Internal/Latent variable harmonisation")