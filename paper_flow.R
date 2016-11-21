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
        'camMets_ind_fact',	'pa_workini_fact', 'cam_index_mean')]
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

plot(x = final_output_men$cam_index_mean, y= final_output_men$PAEE)
abline(0,1)
plot(x = final_output_men$reg_value, y= final_output_men$PAEE)
abline(0,1)
plot(x = final_output_men$binary_index_mean, y= final_output_men$PAEE)
abline(0,1)


#  _____  _               _         ______        _         
# /  ___|| |             | |        |  _  \      | |        
# \ `--. | |_  _   _   __| | _   _  | | | | __ _ | |_  __ _ 
#  `--. \| __|| | | | / _` || | | | | | | |/ _` || __|/ _` |
# /\__/ /| |_ | |_| || (_| || |_| | | |/ /| (_| || |_| (_| |
# \____/  \__| \__,_| \__,_| \__, | |___/  \__,_| \__|\__,_|
#                             __/ |                         
#                            |___/                          

### Setup the data necessary to calculate the cox regression on the prediction set.
study_data <- read.csv("PHIA0000232016_IA88_25Jul/PHIA0000232016.csv")

# these data - 1 = men, 2 = women, based on mean weights, assume men higher:
# aggregate(study_data$weight_adj, by=list(sex=study_data$sex), mean, na.rm = TRUE)
# this is opposite to the validation data

# Fix sex
study_data$sex <- unlist(lapply(study_data$sex, FUN=function(x){
  if (x==1){
    x = 1
  } else {
    x = 0
  }
}))


#fix work
study_data$pa_workini <- unlist(lapply(study_data$pa_work, FUN=function(x){
  if (x == 6){
    out = 5
  } else {
    out = x
  }
  return (out)
}))

# simplify columns
study_data <- study_data[c('MRCid_IAp_13', 'country', 'centre', 'dmstatus_ver_outc', 'age_recr_max', 'fup_time',
  'sex', 'l_school', 'bmi_adj', 'pa_workini', 'm_cycl', 'm_sport', 'alc_re', 'qe_energy', 'smoke_stat')]

study_data$pa_workini_fact <- factor(study_data$pa_work)

## Calculating X_t, with X_t0
tempfupdiff <- (study_data$fup_time/365.25)
study_data$ageEnd <- study_data$age_recr_max + tempfupdiff

# Calculating X_d
eventCens <- lapply(study_data$dmstatus_ver_outc, FUN=function(x){
  if (x==1 || x==2){
    x = 1
  } else {
    x = 0
  }
})

study_data$eventCens <- unlist(eventCens)

# prentice weighted age starts
study_data$age_recr_prentice <- study_data$age_recr_max
for (i in 1:length(study_data$age_recr_max)){
  if(study_data$dmstatus_ver_outc[i] == 2){
    study_data$age_recr_prentice[i] = study_data$ageEnd[i] - 0.00001 
  }
}


# smoke, lschool, country factorization and leveling
study_data$country <- factor(study_data$country, labels=c("FRANCE", "ITALY", "SPAIN", "UK", "NETHERLANDS", "GERMANY", "SWEDEN", "DENMARK"))
study_data$smoke_stat[study_data$smoke_stat == 4] <- NA
study_data$smoke_stat <- factor(study_data$smoke_stat, labels=c("NEVER", "FORMER", "SMOKER"))
study_data$l_school[study_data$l_school == 5] <- NA
study_data$l_school <- factor(study_data$l_school, labels=c("NONE", "PRIMARY", "TECHNICAL/PROFESSIONAL", "SECONDARY", "LONGER EDUCATION/UNI"))

## Calculating own Cambridge Index
study_data$cam_total <- c(rep(0,nrow(study_data)))
for (i in 1:nrow(study_data)){
  study_data$cam_total[i] = sum(study_data$m_cycl[i]/6, study_data$m_sport[i]/6, na.rm=TRUE)    
}

# camMets_ind for cam_matrix
# totalMets_index calculation
study_data$camMets_ind <- as.vector(do.call(rbind, lapply(X=study_data$cam_total, FUN = function(x){
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


study_data$cam_index <- apply(X = study_data[,c('pa_workini', 'camMets_ind')], MARGIN = 1, FUN = function(x){
  output <- cam_matrix[x[1], x[2]]
  return(output)
})

 study_data$camMets_fact <- as.factor(study_data$camMets_ind)
 study_data$cam_index_fact <- as.factor(study_data$cam_index)

# create binary index
study_data$binary_index <- unlist(lapply(study_data$cam_index, FUN=function(x){
  if (x[1] == 1) {output = 1 }
  else {output = 2}
  return(output)
}))

study_data$binary_index <- as.factor(study_data$binary_index)
og_men <- subset(study_data, sex==1)
og_women <- subset(study_data, sex==0)


# fix empty france level for men
og_men$country <- droplevels(og_men$country)

###############################################################################
###############################################################################
###############################################################################
#test the process
#  matches pa_methods in Stata for 'cont' version
test_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_fact + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
test_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_fact + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)

test_reg_cont_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
test_reg_cont_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)

#
# Setting the 'observed' PAEE values by the mean
# work on non-gender split data for convenience.
# cam index - do a 'vlookup'
rm(temp)
rm(new_study_data)
temp <- merge(study_data, cam_index_means, by=c('cam_index', 'sex'))[,c('MRCid_IAp_13', 'sex', 'age_recr_prentice', 'ageEnd','eventCens',
                                                                    'country', 'centre', 'camMets_ind', 'camMets_fact',
                                                                    'cam_index', 'cam_index_fact', 'binary_index',
                                                                    'pa_workini_fact', 'cam_index_mean', 'bmi_adj',
                                                                    'smoke_stat', 'l_school', 'alc_re', 'qe_energy')]
new_study_data <- temp[order(temp$MRCid_IAp_13),]

# binary index - 'vlookup' again
temp <- merge(new_study_data, binary_index_means, by=c('binary_index', 'sex'))
new_study_data <- temp[order(temp$MRCid_IAp_13),]

# regression based
# convert work index into binaries to make use of regression coeffs
temp <- data.frame(sapply(levels(new_study_data$pa_workini_fact), function(x) as.integer(x == new_study_data$pa_workini_fact)))
colnames(temp) <- c('pa_workini_1','pa_workini_2','pa_workini_3','pa_workini_4','pa_workini_5')
new_study_data <- cbind(new_study_data, temp)

# convert ltpa index into binaries
temp <- data.frame(sapply(levels(new_study_data$camMets_fact), function(x) as.integer(x == new_study_data$camMets_fact)))
colnames(temp) <- c('camMets_ind_1','camMets_ind_2','camMets_ind_3','camMets_ind_4')
new_study_data <- cbind(new_study_data, temp)

# use regression coefficients to estimate PAEE
new_study_data_men <- subset(new_study_data, sex==1)
new_study_data_women <- subset(new_study_data, sex==0)

# fix empty france level for men
new_study_data_men$country <- droplevels(new_study_data_men$country)

new_study_data_women$reg_value <- coeffs_women['(Intercept)'] + new_study_data_women$pa_workini_2 * coeffs_women['pa_workini_fact2'] +
  new_study_data_women$pa_workini_3 * coeffs_women['pa_workini_fact3'] + new_study_data_women$pa_workini_4 * coeffs_women['pa_workini_fact4'] +
  new_study_data_women$pa_workini_5 * coeffs_women['pa_workini_fact5'] +
  new_study_data_women$camMets_ind_2 * coeffs_women['camMets_ind_fact2'] + new_study_data_women$camMets_ind_3 * coeffs_women['camMets_ind_fact3'] +
  new_study_data_women$camMets_ind_4 * coeffs_women['camMets_ind_fact4']

new_study_data_men$reg_value <- coeffs_men['(Intercept)'] + new_study_data_men$pa_workini_2 * coeffs_men['pa_workini_fact2'] +
  new_study_data_men$pa_workini_3 * coeffs_men['pa_workini_fact3'] + new_study_data_men$pa_workini_4 * coeffs_men['pa_workini_fact4'] +
  new_study_data_men$pa_workini_5 * coeffs_men['pa_workini_fact5'] +
  new_study_data_men$camMets_ind_2 * coeffs_men['camMets_ind_fact2'] + new_study_data_men$camMets_ind_3 * coeffs_men['camMets_ind_fact3'] +
  new_study_data_men$camMets_ind_4 * coeffs_men['camMets_ind_fact4']


####Cox regressions

# we would expect the binary estimate to be the "worst" (ie hazard ratio closer to 1),
# but it is not!
# for women, cam index is the "worst"
# for men, regression is the "worst"

# this is equivalent to active and not active categories
binary_cox_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ binary_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_study_data_men, robust=TRUE)
binary_cox_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ binary_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_study_data_women, robust=TRUE)

cam_cox_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_mean + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_study_data_men, robust=TRUE)
cam_cox_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_mean + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_study_data_women, robust=TRUE)

reg_cox_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ reg_value + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_study_data_men, robust=TRUE)
reg_cox_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ reg_value + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_study_data_women, robust=TRUE)

# comparison of associations using binary index means and REMA
## MEN
countries = levels(new_study_data_men$country)

estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = new_study_data_men, subset = new_study_data_men$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ binary_index_mean + bmi_adj
    + smoke_stat + l_school + alc_re + qe_energy, 
    data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('binary_index_mean', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('binary_index_mean', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',.(round(res$QEp,3)),')')), 
  xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ', .(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Binary Index - Men', adj = c( 0, 1 ))


## WOMEN
countries = levels(new_study_data_women$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = new_study_data_women, subset = new_study_data_women$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ binary_index_mean + bmi_adj
    + smoke_stat + l_school + alc_re + qe_energy, data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('binary_index_mean', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('binary_index_mean', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ', .(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ', .(round(res$pval,3)))), atransf = exp)

usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Binary Index - women', adj = c( 0, 1 ))


# comparison of associations using cam index means and REMA
#MEN
countries = levels(new_study_data_men$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = new_study_data_men, subset = new_study_data_men$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_mean + bmi_adj 
    + smoke_stat + l_school + alc_re + qe_energy, data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('cam_index_mean', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('cam_index_mean', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)
forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ', .(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ',.(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Cam Index - Men', adj = c( 0, 1 ))


#WOMEN
countries = levels(new_study_data_women$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = new_study_data_women, subset = new_study_data_women$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_mean + bmi_adj + 
    smoke_stat + l_school + alc_re + qe_energy, data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('cam_index_mean', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('cam_index_mean', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',.(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ',.(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Cam Index - women', adj = c( 0, 1 ))

# comparison of associations using regression values and REMA
#MEN
countries = levels(new_study_data_men$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = new_study_data_men, subset = new_study_data_men$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ reg_value + bmi_adj + smoke_stat + l_school + alc_re + qe_energy, 
    data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('reg_value', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('reg_value', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',.(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ',.(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Regression values - Men', adj = c( 0, 1 ))


#WOMEN
countries = levels(new_study_data_women$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = new_study_data_women, subset = new_study_data_women$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ reg_value + bmi_adj + 
    smoke_stat + l_school + alc_re + qe_energy, data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('reg_value', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('reg_value', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',.(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ',.(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Regression values - women', adj = c( 0, 1 ))

#  _____ _                 _       _   _                 
# /  ___(_)               | |     | | (_)                
# \ `--. _ _ __ ___  _   _| | __ _| |_ _  ___  _ __  ___ 
#  `--. \ | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_ \/ __|
# /\__/ / | | | | | | |_| | | (_| | |_| | (_) | | | \__ \
# \____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|___/


# ______         _            _   _             _     _   
# | ___ \       | |          | | (_)           (_)   | |  
# | |_/ /___  __| |_   _  ___| |_ _  ___  _ __  _ ___| |_ 
# |    // _ \/ _` | | | |/ __| __| |/ _ \| '_ \| / __| __|
# | |\ \  __/ (_| | |_| | (__| |_| | (_) | | | | \__ \ |_ 
# \_| \_\___|\__,_|\__,_|\___|\__|_|\___/|_| |_|_|___/\__|
                                                        

# MEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS - cam
# GERMANY SWEDEN DENMARK - reg

# Simulation
############################################################################
# Reductionist Harmonisation: Make everything binary dependent on thresholds
############################################################################

reductionist_men <- new_study_data_men[,c('country','bmi_adj', 'smoke_stat', 
                                          'l_school', 'alc_re', 'qe_energy',
                                          'age_recr_prentice', 'ageEnd', 'eventCens',
                                          'binary_index', 'cam_index', 'reg_value')]

# red_harm can take a value of 1 meaning inactive, or 2 meaning active

reductionist_men$red_harm <- apply(X = reductionist_men[,c('country','bmi_adj', 'smoke_stat', 
                                                   'l_school', 'alc_re', 'qe_energy',
                                                   'age_recr_prentice', 'ageEnd', 'eventCens',
                                                   'binary_index', 'cam_index', 'reg_value')], MARGIN = 1, 
                             FUN = function(x){
                               if(x[1] == 'ITALY' | x[1] == 'SPAIN'){
                                 output = x[10]
                               }
                               # in this conversion of cam index, we will use a different rule
                               # to that which was used to define binary index
                               if(x[1] == 'UK' | x[1] == 'NETHERLANDS'){

                                 if(x[11] == 1){
                                   output = 1
                                 }
                                 else {
                                   output = 2
                                 }
                               }
                               # arbitrarily choose 45 as cut point between active and non active
                               if(x[1] == 'GERMANY' | x[1] == 'SWEDEN' | x[1] == 'DENMARK'){
                                 if (x[12] < 45){
                                  output = 1
                                 }
                                 else {
                                   output = 2
                                 }
                               }
                               return(output)
                             }
)

countries = levels(reductionist_men$country)

  estimates = vector()
  s_errors = vector()

for (i in 1:length(countries)){
    temp = droplevels(subset(x = reductionist_men, subset = reductionist_men$country == countries[i]))
    model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ red_harm + bmi_adj
      + smoke_stat + l_school + alc_re + qe_energy, data = temp, robust=TRUE))
    model_coeffs <- model$coefficients[,c('coef', 'robust se')]
    estimates = rbind(estimates,model_coeffs[grep('red_harm', row.names(model_coeffs)),"coef"])
    s_errors = rbind(s_errors,model_coeffs[grep('red_harm', row.names(model_coeffs)),"robust se"])
}
  
res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',.(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ',.(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Reductionist - Men', adj = c( 0, 1 ))


# WOMEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS & FRANCE - cam
# GERMANY SWEDEN DENMARK - reg

# Simulation
############################################################################
# Reductionist Harmonisation: Make everything binary dependent on thresholds
############################################################################

reductionist_women <- new_study_data_women[,c('country','bmi_adj', 'smoke_stat', 
                                          'l_school', 'alc_re', 'qe_energy',
                                          'age_recr_prentice', 'ageEnd', 'eventCens',
                                          'binary_index', 'cam_index', 'reg_value')]

# red_harm can take a value of 1 meaning inactive, or 2 meaning active
reductionist_women$red_harm <- apply(X = reductionist_women[,c('country','bmi_adj', 'smoke_stat', 
                                                           'l_school', 'alc_re', 'qe_energy',
                                                           'age_recr_prentice', 'ageEnd', 'eventCens',
                                                           'binary_index', 'cam_index', 'reg_value')], MARGIN = 1, 
                                   FUN = function(x){
                                     if(x[1] == 'ITALY' | x[1] == 'SPAIN'){
                                       output = x[10]
                                     }
                                     # in this conversion of cam index, we will use a different rule
                                     # to that which was used to define binary index
                                     if(x[1] == 'UK' | x[1] == 'NETHERLANDS' |x[1] == 'FRANCE'){
                                       
                                       if(x[11] == 1){
                                         output = 1
                                       }
                                       else {
                                         output = 2
                                       }
                                     }
                                     # arbitrarily choose 45 as cut point between active and non active
                                     if(x[1] == 'GERMANY' | x[1] == 'SWEDEN' | x[1] == 'DENMARK'){
                                       if (x[12] < 45){
                                         output = 1
                                       }
                                       else {
                                         output = 2
                                       }
                                     }
                                     return(output)
                                   }
)



countries = levels(reductionist_women$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = reductionist_women, subset = reductionist_women$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ red_harm + bmi_adj
    + smoke_stat + l_school + alc_re + qe_energy, 
    data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('red_harm', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('red_harm', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',
                                        .(round(res$QEp,3)),')')),
       xlab=bquote(paste('Test of H'[0]*': true hazard ratio = 1, p = ',
                         .(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Reductionist - women', adj = c( 0, 1 ))

#  _   _       _ _     _       _   _             
# | | | |     | (_)   | |     | | (_)            
# | | | | __ _| |_  __| | __ _| |_ _  ___  _ __  
# | | | |/ _` | | |/ _` |/ _` | __| |/ _ \| '_ \ 
# \ \_/ / (_| | | | (_| | (_| | |_| | (_) | | | |
#  \___/ \__,_|_|_|\__,_|\__,_|\__|_|\___/|_| |_|

# ___  ___ _____ _____ _   _ ___________   __  
# |  \/  ||  ___|_   _| | | |  _  |  _  \ /  | 
# | .  . || |__   | | | |_| | | | | | | | `| | 
# | |\/| ||  __|  | | |  _  | | | | | | |  | | 
# | |  | || |___  | | | | | \ \_/ / |/ /  _| |_
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \___/
                                               
# Simulation
############################################################################
# Validation Harmonisation: Convert to PAEE using means
############################################################################
# MEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS - cam
# GERMANY SWEDEN DENMARK - reg

validation_men <- new_study_data_men[,c('country','bmi_adj', 'smoke_stat', 
                                          'l_school', 'alc_re', 'qe_energy',
                                          'age_recr_prentice', 'ageEnd', 'eventCens',
                                          'binary_index_mean', 'cam_index_mean', 'reg_value')]

# select the appropriate index mean
validation_men$val_harm <- apply(X = validation_men[,c('country','bmi_adj', 'smoke_stat', 
                                                           'l_school', 'alc_re', 'qe_energy',
                                                           'age_recr_prentice', 'ageEnd', 'eventCens',
                                                           'binary_index_mean', 'cam_index_mean', 'reg_value')], MARGIN = 1, 
                                   FUN = function(x){
                                     if(x[1] == 'ITALY' | x[1] == 'SPAIN'){
                                       # choose binary index mean PAEE
                                       output = x[10]
                                     }
                                     
                                     if(x[1] == 'UK' | x[1] == 'NETHERLANDS'){
                                       # choose cam index mean
                                       output = x[11]
                                     }
                                     
                                     if(x[1] == 'GERMANY' | x[1] == 'SWEDEN' | x[1] == 'DENMARK'){
                                       # choose reg mean
                                       output = x[12]
                                     }
                                     return(as.numeric(output))
                                   }
)

countries = levels(validation_men$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = validation_men, subset = validation_men$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ val_harm + bmi_adj
    + smoke_stat + l_school + alc_re + qe_energy, data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('val_harm', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('val_harm', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',.(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true relative risk = 1, p = ',.(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Validation - Men', adj = c( 0, 1 ))


# WOMEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS & FRANCE - cam
# GERMANY SWEDEN DENMARK - reg

# validation harmonisation method. Convert all to PAEE value
validation_women <- new_study_data_women[,c('country','bmi_adj', 'smoke_stat', 
                                        'l_school', 'alc_re', 'qe_energy',
                                        'age_recr_prentice', 'ageEnd', 'eventCens',
                                        'binary_index_mean', 'cam_index_mean', 'reg_value')]

# select the appropriate index mean
validation_women$val_harm <- apply(X = validation_women[,c('country','bmi_adj', 'smoke_stat', 
                                                       'l_school', 'alc_re', 'qe_energy',
                                                       'age_recr_prentice', 'ageEnd', 'eventCens',
                                                       'binary_index_mean', 'cam_index_mean', 'reg_value')], MARGIN = 1, 
                                 FUN = function(x){
                                   if(x[1] == 'ITALY' | x[1] == 'SPAIN'){
                                     # choose binary index mean PAEE
                                     output = x[10]
                                   }
                                   
                                   if(x[1] == 'UK' | x[1] == 'NETHERLANDS' |  x[1] == 'FRANCE'){
                                     # choose cam index mean
                                     output = x[11]
                                   }
                                   
                                   if(x[1] == 'GERMANY' | x[1] == 'SWEDEN' | x[1] == 'DENMARK'){
                                     # choose reg mean
                                     output = x[12]
                                   }
                                   return(as.numeric(output))
                                 }
)

countries = levels(validation_women$country)
estimates = vector()
s_errors = vector()

for (i in 1:length(countries)){
  temp = droplevels(subset(x = validation_women, subset = validation_women$country == countries[i]))
  model = summary(coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ val_harm + bmi_adj
    + smoke_stat + l_school + alc_re + qe_energy, data = temp, robust=TRUE))
  model_coeffs <- model$coefficients[,c('coef', 'robust se')]
  estimates = rbind(estimates,model_coeffs[grep('val_harm', row.names(model_coeffs)),"coef"])
  s_errors = rbind(s_errors,model_coeffs[grep('val_harm', row.names(model_coeffs)),"robust se"])
}

res <- rma(yi = estimates, sei = s_errors, method='DL', slab = countries)

forest(res, digits=3, mlab=bquote(paste('Overall (I'^2*' = ', .(round(res$I2)),'%, p = ',.(round(res$QEp,3)),')')),
  xlab=bquote(paste('Test of H'[0]*': true relative risk = 1, p = ',.(round(res$pval,3)))), atransf = exp)
usr <- par("usr")
text(usr[2], usr[4], "Hazard ratio [95% CI]", adj = c(1, 4))
text(usr[1], usr[4], 'Validation - women', adj = c( 0, 1 ))


# ___  ___ _____ _____ _   _ ___________   _____   ___  
# |  \/  ||  ___|_   _| | | |  _  |  _  \ / __  \ / _ \ 
# | .  . || |__   | | | |_| | | | | | | | `' / /'/ /_\ \
# | |\/| ||  __|  | | |  _  | | | | | | |   / /  |  _  |
# | |  | || |___  | | | | | \ \_/ / |/ /  ./ /___| | | |
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \_____/\_| |_/           
# Simulation
############################################################################
# Validation approach 2a: Sampling Directly from data
############################################################################
# We have to first make 3 new columns in both the validation and study data set for sampled values
# 
# MEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS - cam
# GERMANY SWEDEN DENMARK - reg


# WOMEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS & FRANCE - cam
# GERMANY SWEDEN DENMARK - reg

# ___  ___ _____ _____ _   _ ___________   _____ ______ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \ / __  \| ___ \
# | .  . || |__   | | | |_| | | | | | | | `' / /'| |_/ /
# | |\/| ||  __|  | | |  _  | | | | | | |   / /  | ___ \
# | |  | || |___  | | | | | \ \_/ / |/ /  ./ /___| |_/ /
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \_____/\____/ 
# Simulation
############################################################################
# Validation approach 2b: Sampling Directly from data
############################################################################
# We have to first make 3 new columns in both the validation and study data set for sampled values
# from the subsets of the indices (sex stratified). Then calculate lambdas and make the simulations as well.

# binary_sampled (2 groups * 2 genders)
# cam_sampled (4 groups * 2 genders)
# reg_sampled (16 groups * 2 genders)



# MEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS - cam
# GERMANY SWEDEN DENMARK - reg


# WOMEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS & FRANCE - cam
# GERMANY SWEDEN DENMARK - reg

# ___  ___ _____ _____ _   _ ___________   _____ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \ |____ |
# | .  . || |__   | | | |_| | | | | | | |     / /
# | |\/| ||  __|  | | |  _  | | | | | | |     \ \
# | |  | || |___  | | | | | \ \_/ / |/ /  .___/ /
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \____/     
# We fit distributions to the validation data and then generate samples from these for the study data
# Simulation
############################################################################
# Validation approach 3: Sampling from Fitted Distribution
############################################################################


# binary_fitted (2 groups * 2 genders)
# cam_fitted (4 groups * 2 genders)
# reg_fitted (16 groups * 2 genders)


# MEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS - cam
# GERMANY SWEDEN DENMARK - reg


# WOMEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS & FRANCE - cam
# GERMANY SWEDEN DENMARK - reg



# ___  ___ _____ _____ _   _ ___________     ___ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \   /   |
# | .  . || |__   | | | |_| | | | | | | |  / /| |
# | |\/| ||  __|  | | |  _  | | | | | | | / /_| |
# | |  | || |___  | | | | | \ \_/ / |/ /  \___  |
# \_|  |_/\____/  \_/ \_| |_/\___/|___/       |_/
#                                                
# We estimate the density of the pdf using kernel based estimation methods. Then sample directly using the kernels.
# Simulation
############################################################################
# Validation approach 3: Sampling from Kernel Density Estimated Distribution
############################################################################

# binary_kden (2 groups * 2 genders)
# cam_kden (4 groups * 2 genders)
# reg_kden (16 groups * 2 genders)

# MEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS - cam
# GERMANY SWEDEN DENMARK - reg


# WOMEN
# ITALY & SPAIN - binary
# UK & NETHERLANDS & FRANCE - cam
# GERMANY SWEDEN DENMARK - reg
