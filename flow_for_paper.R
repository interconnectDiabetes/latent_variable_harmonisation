##########################################################################################################
##########################################################################################################
############################### Lambda Measurement Error Calculation #####################################
##########################################################################################################
##########################################################################################################

## Paper(s):
## 1. Physical activity reduces the risk of incident type 2 diabetes in 
## general and in abdominally lean and obese men and women: the EPIC-InterAct 
## Study
## 2. Validity of a short questionnaire to assess physical activity in 10 European
## countries 
## Author: Paul Scherer
## Date: 12/08/2016

###############################################################################
########################### DATA AND SETTINGS #################################
###############################################################################
## WARNING THIS SCRIPT WILL ONLY RUN IF WORKING DIRECTORY IS SET TO THE RIGHT LOCATION
setwd("V:/Studies/InterConnect/Internal/Latent variable harmonisation")
library("survival")

# Read the actiheart and EPIC study data (sweden UMEA is missing in 
# main file, needs to be loaded seperately)
actiheart_summary <- read.csv("PHIA0000112014_IA85_12Mar/PAQIA0000112014_actiheart_summary.csv", header=T)
epic <- read.csv("PHIA0000112014_IA85_12Mar/PAQIA0000112014_epic.csv", header=T)

# actiheart - 0 = women, 1 = men, based on mean weights, assume men higher:
# aggregate(actiheart_summary$weight, by=list(sex=actiheart_summary$sex), mean, na.rm = TRUE)

###############################################################################
###############################################################################
###############################################################################

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
###############to do: Add in UMEA        #####################################
###############################################################################
# add Umea , if possible (not sure if all indices can be generated)
# 
# epic_umea <- read.csv("V:/Studies/InterConnect/Internal/Latent variable harmonisation/epic_umea.csv", header=T)
# the csv was generated from:
# V:\Studies\InterAct_WP_2.6\StudyData\HANDOVER TO DM\links\umea_idmerge.dta
# V:\Studies\InterAct_WP_2.6\UMEÅ_EPIC_PAQ\UMEA_EPIC_final.xls


### merged_output creation
# merge or perform a "natural join" on the trimmed epic table and 
# cleaned actiheart table by universal id
merged_output <- merge(output, epic, by = 'universal_id')

# fix employment status

merged_output$pa_workini <- as.vector(do.call(rbind,lapply(X = merged_output$pa_workini, FUN = function(x){
  
  if ( x == 9){
    ind = 5
  }
  else { 
    ind = x
  }
  return(ind)
})
))

# calculate the BMI using the weight and height and include this in the table to make
# a model for PAEE
bmi_calc <- function(weight, height){
  bmi = (weight/height)/height
  return(bmi)
}
merged_output$bmi <- mapply(FUN=bmi_calc, merged_output$weight, merged_output$height)


###############################################################################
###############################################################################
###############################################################################
## In this section we regenerate cam_index column ie check logic - seems to match up

# camMets_ind for cam_matrix
# totalMets_index calculation

merged_output$camMets_ind <- as.vector(do.call(rbind, lapply(X=merged_output$pa_lt, FUN = function(x){
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

#unemployed or missing all sedentary

cam_matrix = matrix( byrow = TRUE,
                     c(1, 2, 3, 4, 
                       2, 3, 4, 4, 
                       3, 4, 4, 4, 
                       4, 4, 4, 4,
                       1, 1, 1, 1), 
                     nrow=5, 
                     ncol=4)

merged_output$cam_index2 <- apply(X = merged_output[,c('pa_workini', 'camMets_ind')], MARGIN = 1,
                                  FUN = function(x){
                                   
                                    output <- cam_matrix[x[1], x[2]]
                                    return(output)
                                  })

#if test has zero length then we are happy with the original value of cam_index!

test <- merged_output[merged_output$cam_index!=merged_output$cam_index2,]

###############################################################################
###############################################################################
###############################################################################
#record means etc for the cam index and binary index

cam_index_means = aggregate(merged_output$PAEE, by=list(cam_index=merged_output$cam_index, sex=merged_output$sex), mean, na.rm = TRUE)
colnames(cam_index_means)[3] = 'cam_index_mean'
cam_index_medians = aggregate(merged_output$PAEE, by=list(cam_index=merged_output$cam_index, sex=merged_output$sex), median, na.rm = TRUE)
colnames(cam_index_medians)[3] = 'cam_index_median'
cam_index_counts =aggregate(merged_output$PAEE, by=list(cam_index=merged_output$cam_index, sex=merged_output$sex), length)
colnames(cam_index_counts)[3] = 'cam_index_count'

# now define other indices
# Cam index is Method A

# Method B - ie sedentary (1) (cam index 1) or active (2) (cam index 2,3,4)

merged_output$binary_index <- unlist(lapply(merged_output$cam_index, FUN=function(x){
  if (x[1] == 1) {output = 1 }
  else {output = 2}
  
  return(output)
}))

# record results

binary_index_means = aggregate(merged_output$PAEE, by=list(binary_index=merged_output$binary_index, sex=merged_output$sex), mean, na.rm = TRUE)
colnames(binary_index_means)[3] = 'binary_index_mean'
binary_index_medians = aggregate(merged_output$PAEE, by=list(binary_index=merged_output$binary_index, sex=merged_output$sex), median, na.rm = TRUE)
colnames(binary_index_medians)[3] = 'binary_index_median'
binary_index_counts =aggregate(merged_output$PAEE, by=list(binary_index=merged_output$binary_index, sex=merged_output$sex), length)
colnames(binary_index_counts)[3] = 'binary_index_count'

###############################################################################
###############################################################################
###############################################################################
# Method C - regression based (but still essentially 16 categories) on components of cam index



#create factors for regression
merged_output$camMets_ind_fact <- as.factor(merged_output$camMets_ind)
merged_output$pa_work_ind <- as.factor(merged_output$pa_workini)

merged_output_women <- subset(merged_output, sex==0)
merged_output_men <- subset(merged_output, sex==1)

reg_out_men <- glm(formula=PAEE ~ pa_work_ind + camMets_ind_fact , data = merged_output_men, family='gaussian')
reg_out_women <- glm(formula=PAEE ~ pa_work_ind + camMets_ind_fact , data = merged_output_women, family='gaussian')

#store coefficients

#coeffs <- summary(reg_out)$coefficients[, 1]

coeffs_men <- summary(reg_out_men)$coefficients[, 1]
coeffs_women <- summary(reg_out_women)$coefficients[, 1]

###############################################################################
###############################################################################
###############################################################################
# calculate the lambdas (to be used much later!)

# Setting the 'observed' PAEE values by the mean

# cam index - do a 'vlookup'

temp <- merge(merged_output, cam_index_means, by=c('cam_index', 'sex'))[,c('universal_id', 'sex', 
        'country.x', 'PAEE', 'center', 'camMets_ind',  'cam_index', 'binary_index',
        'camMets_ind_fact',	'pa_work_ind', 'pa_total', 'bmi', 'cam_index_mean')]

new_data <- temp[order(temp$universal_id),]

# binary index - 'vlookup' again

temp <- merge(new_data, binary_index_means, by=c('binary_index', 'sex'))

new_data <- temp[order(temp$universal_id),]

# regression based

# convert work index into binaries to make use of regression coeffs

temp <- data.frame(sapply(levels(new_data$pa_work_ind), function(x) as.integer(x == new_data$pa_work_ind)))
colnames(temp) <- c('pa_work_ind_1','pa_work_ind_2','pa_work_ind_3','pa_work_ind_4','pa_work_ind_5')
new_data <- cbind(new_data, temp)

# convert ltpa index into binaries

temp <- data.frame(sapply(levels(new_data$camMets_ind_fact), function(x) as.integer(x == new_data$camMets_ind_fact)))
colnames(temp) <- c('camMets_ind_1','camMets_ind_2','camMets_ind_3','camMets_ind_4')
new_data <- cbind(new_data, temp)

# use regression coefficients to estimate PAEE

final_output_women <- subset(new_data, sex==0)
final_output_men <- subset(new_data, sex==1)

final_output_women$reg_value <- coeffs_women['(Intercept)'] + final_output_women$pa_work_ind_2 * coeffs_women['pa_work_ind2'] +
  final_output_women$pa_work_ind_3 * coeffs_women['pa_work_ind3'] + final_output_women$pa_work_ind_4 * coeffs_women['pa_work_ind4'] +
  final_output_women$pa_work_ind_5 * coeffs_women['pa_work_ind5'] +
  final_output_women$camMets_ind_2 * coeffs_women['camMets_ind_fact2'] + final_output_women$camMets_ind_3 * coeffs_women['camMets_ind_fact3'] +
  final_output_women$camMets_ind_4 * coeffs_women['camMets_ind_fact4']

final_output_men$reg_value <- coeffs_men['(Intercept)'] + final_output_men$pa_work_ind_2 * coeffs_men['pa_work_ind2'] +
  final_output_men$pa_work_ind_3 * coeffs_men['pa_work_ind3'] + final_output_men$pa_work_ind_4 * coeffs_men['pa_work_ind4'] +
  final_output_men$pa_work_ind_5 * coeffs_men['pa_work_ind5'] +
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


###############################################################################
###############################################################################
###############################################################################
###
### Cox Regression bit
###

### Setup the data necessary to calculate the cox regression on the prediction set.
og_data <- read.csv("PHIA0000232016_IA88_25Jul/PHIA0000232016.csv")


# these data - 1 = men, 2 = women, based on mean weights, assume men higher:
# aggregate(og_data$weight_adj, by=list(sex=og_data$sex), mean, na.rm = TRUE)
# this is opposite to the validation data

# Fix sex
og_data$sex <- unlist(lapply(og_data$sex, FUN=function(x){
  if (x==1){
    x = 1
  } else {
    x = 0
  }
}))


# simplify columns

og_data <- og_data[c('MRCid_IAp_13', 'country', 'centre', 'dmstatus_ver_outc', 'age_recr_max', 'fup_time',
                     'sex', 'l_school', 'bmi_adj', 'pa_work', 'm_cycl', 'm_sport',
                     'alc_re', 'qe_energy', 'smoke_stat')]




og_data$pa_work <- unlist(lapply(og_data$pa_work, FUN=function(x){
  if (x == 6){
    out = 5
  } else {
    out = x
  }
  return (out)
}))

og_data$pa_workini <- factor(og_data$pa_work)

## Calculating X_t, with X_t0
tempfupdiff <- (og_data$fup_time/365.25)
og_data$ageEnd <- og_data$age_recr_max + tempfupdiff

# Calculating X_d
eventCens <- lapply(og_data$dmstatus_ver_outc, FUN=function(x){
  if (x==1 || x==2){
    x = 1
  } else {
    x = 0
  }
})

og_data$eventCens <- unlist(eventCens)

# prentice weighted age starts
og_data$age_recr_prentice <- og_data$age_recr_max
for (i in 1:length(og_data$age_recr_max)){
  if(og_data$dmstatus_ver_outc[i] == 2){
    og_data$age_recr_prentice[i] = og_data$ageEnd[i] - 0.00001 
  }
}


# smoke, lschool, country factorization and leveling
og_data$country <- factor(og_data$country, labels=c("FRANCE", "ITALY", 
                                                    "SPAIN", "UK", "NETHERLANDS", "GERMANY", "SWEDEN", "DENMARK"))

og_data$smoke_stat <- factor(og_data$smoke_stat, labels=c("NEVER", "FORMER", 
                                                          "SMOKER", "UNKOWN"))

og_data$l_school <- factor(og_data$l_school, labels=c("NONE", "PRIMARY", 
                                                      "TECHNICAL/PROFESSIONAL", "SECONDARY", "LONGER EDUCATION/UNI", "NOT SPECIFIED"))


## Calculating own Cambridge Index
og_data$cam_total <- c(rep(0,nrow(og_data)))
for (i in 1:nrow(og_data)){
  og_data$cam_total[i] = sum(og_data$m_cycl[i]/6, og_data$m_sport[i]/6, na.rm=TRUE)    
}

# camMets_ind for cam_matrix
# totalMets_index calculation
og_data$camMets_ind <- as.vector(do.call(rbind, lapply(X=og_data$cam_total, FUN = function(x){
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


og_data$cam_index <- apply(X = og_data[,c('pa_work', 'camMets_ind')], MARGIN = 1, 
                           FUN = function(x){
                             output <- cam_matrix[x[1], x[2]]
                             return(output)
                           }
)

 og_data$camMets_fact <- as.factor(og_data$camMets_ind)
 og_data$cam_index_fact <- as.factor(og_data$cam_index)

# create binary index

og_data$binary_index <- unlist(lapply(og_data$cam_index, FUN=function(x){
  if (x[1] == 1) {output = 1 }
  else {output = 2}
  
  return(output)
}))

og_data$binary_index <- as.factor(og_data$binary_index)

og_men <- subset(og_data, sex==1)
og_women <- subset(og_data, sex==0)


# fix empty france level for men
og_men$country <- droplevels(og_men$country)



###############################################################################
###############################################################################
###############################################################################
#test the process
# vaguely matches pa_methods in Stata

test_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_fact + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
test_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_fact + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)

test_reg_cont_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
test_reg_cont_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)


#
# Setting the 'observed' PAEE values by the mean
# work on non-gender split data for convenience.

# cam index - do a 'vlookup'

rm(temp)
rm(new_og_data)

temp <- merge(og_data, cam_index_means, by=c('cam_index', 'sex'))[,c('MRCid_IAp_13', 'sex', 'age_recr_prentice', 'ageEnd',
                                                                     'eventCens',
                                                                    'country', 'centre', 'camMets_ind', 'camMets_fact',
                                                                    'cam_index', 'cam_index_fact', 'binary_index',
                                                                    'pa_workini', 'cam_index_mean', 'bmi_adj',
                                                                     'smoke_stat', 'l_school', 'alc_re', 'qe_energy')]



new_og_data <- temp[order(temp$MRCid_IAp_13),]

# binary index - 'vlookup' again

temp <- merge(new_og_data, binary_index_means, by=c('binary_index', 'sex'))

new_og_data <- temp[order(temp$MRCid_IAp_13),]

# regression based

# convert work index into binaries to make use of regression coeffs

temp <- data.frame(sapply(levels(new_og_data$pa_workini), function(x) as.integer(x == new_og_data$pa_workini)))
colnames(temp) <- c('pa_workini_1','pa_workini_2','pa_workini_3','pa_workini_4','pa_workini_5')
new_og_data <- cbind(new_og_data, temp)

# convert ltpa index into binaries

temp <- data.frame(sapply(levels(new_og_data$camMets_fact), function(x) as.integer(x == new_og_data$camMets_fact)))
colnames(temp) <- c('camMets_ind_1','camMets_ind_2','camMets_ind_3','camMets_ind_4')
new_og_data <- cbind(new_og_data, temp)

# use regression coefficients to estimate PAEE


new_data_men <- subset(new_og_data, sex==1)
new_data_women <- subset(new_og_data, sex==0)

# fix empty france level for men
new_data_men$country <- droplevels(new_data_men$country)

new_data_women$reg_value <- coeffs_women['(Intercept)'] + new_data_women$pa_workini_2 * coeffs_women['pa_work_ind2'] +
  new_data_women$pa_workini_3 * coeffs_women['pa_work_ind3'] + new_data_women$pa_workini_4 * coeffs_women['pa_work_ind4'] +
  new_data_women$pa_workini_5 * coeffs_women['pa_work_ind5'] +
  new_data_women$camMets_ind_2 * coeffs_women['camMets_ind_fact2'] + new_data_women$camMets_ind_3 * coeffs_women['camMets_ind_fact3'] +
  new_data_women$camMets_ind_4 * coeffs_women['camMets_ind_fact4']

new_data_men$reg_value <- coeffs_men['(Intercept)'] + new_data_men$pa_workini_2 * coeffs_men['pa_work_ind2'] +
  new_data_men$pa_workini_3 * coeffs_men['pa_work_ind3'] + new_data_men$pa_workini_4 * coeffs_men['pa_work_ind4'] +
  new_data_men$pa_workini_5 * coeffs_men['pa_work_ind5'] +
  new_data_men$camMets_ind_2 * coeffs_men['camMets_ind_fact2'] + new_data_men$camMets_ind_3 * coeffs_men['camMets_ind_fact3'] +
  new_data_men$camMets_ind_4 * coeffs_men['camMets_ind_fact4']


####Cox regressions

# As a baseline, do the regression assuming we are harmonising to the lowest common denominator
# which in this case is binary index
# similarly, there is a singularity error on men
# also this gives 'better' results for men ie smaller hazard ratio, when they should be 'worse'


binary_cox_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ binary_index_mean + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_data_men, robust=TRUE)
binary_cox_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ binary_index_mean + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_data_women, robust=TRUE)


cam_cox_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_mean + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_data_men, robust=TRUE)
cam_cox_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_mean + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_data_women, robust=TRUE)

reg_cox_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ reg_value + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_data_men, robust=TRUE)
reg_cox_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ reg_value + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = new_data_women, robust=TRUE)

