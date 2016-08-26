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
## setwd("V:/Studies/InterConnect/Internal/Latent variable harmonisation")
library("survival")

# Read the actiheart and EPIC study data (sweden UMEA is missing in 
# main file, needs to be loaded seperately)
actiheart_summary <- read.csv("PHIA0000112014_IA85_12Mar/PAQIA0000112014_actiheart_summary.csv", header=T)
epic <- read.csv("PHIA0000112014_IA85_12Mar/PAQIA0000112014_epic.csv", header=T)


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

### merged_output creation
# merge or perform a "natural join" on the trimmed epic table and 
# cleaned actiheart table by universal id
merged_output <- merge(output, epic, by = 'universal_id')

# calculate the BMI using the weight and height and include this in the table to make
# a model for PAEE
bmi_calc <- function(weight, height){
  bmi = (weight/height)/height
  return(bmi)
}
merged_output$bmi <- mapply(FUN=bmi_calc, merged_output$weight, merged_output$height)

# Make sure pa_workini is seen as a categorical variable and not a double
merged_output$pa_workini <- as.factor(merged_output$pa_workini)
merged_output$new_ltpa <- mapply(FUN=sum, merged_output$m_walk, merged_output$m_floors, 
                                 merged_output$m_cycl, merged_output$m_sport, merged_output$m_houswk, merged_output$m_vigpa, 
                                 merged_output$m_gard, merged_output$m_diy)

## Calculating own Cambridge Index
merged_output$cam_total <- c(rep(0,nrow(merged_output)))
for (i in 1:nrow(merged_output)){
  merged_output$cam_total[i] = sum(merged_output$m_cycl[i]/6, merged_output$m_sport[i]/6, na.rm=TRUE)    
}

# camMets_ind for cam_matrix
# totalMets_index calculation
merged_output$camMets_ind <- as.vector(do.call(rbind, lapply(X=merged_output$cam_total, FUN = function(x){
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


### 
### Setting the 'observed' PAEE values by the mean
###
## Cam index means
merged_output$cam_index_means <- unlist(mapply(merged_output$cam_index, merged_output$sex, SIMPLIFY = FALSE, FUN=function(x,y){
  if (y == 0) {
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = 35.6
    } else if (x == 2) {
      output = 43.7
    } else if (x == 3) {
      output = 49.0
    } else if (x == 4) {
      output = 56.2
    } else {
      output = NA
    }
  } else if (y == 1){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = 36.5
    } else if (x == 2) {
      output = 39.8
    } else if (x == 3) {
      output = 43.6
    } else if (x == 4) {
      output = 48.2
    } else {
      output = NA
    } 
  } else {
    output = NA
  }
  return(output)
}
))

###############################################################################################
##################################### FIRST RUN THROUGH #######################################
###############################################################################################

###
### Regression and splitting of the genders
###
merged_output_men <- subset(merged_output, sex==0)
merged_output_women <- subset(merged_output, sex==1)

rdr_regression_fit_men <- lm(formula=PAEE~cam_index_means, data=merged_output_men)
rdr_regression_fit_women <- lm(formula=PAEE~cam_index_means, data=merged_output_women)

lambda_men <- rdr_regression_fit_men$coefficients[2]
lambda_women <- rdr_regression_fit_women$coefficients[2]

lambda_men_var <- (summary(rdr_regression_fit_men)$coefficients[4])^2
lambda_women_var <- (summary(rdr_regression_fit_women)$coefficients[4])^2

# Store the lambda in a list for comparison
lambda_list_men <- c(lambda_men)
lambda_list_women <- c(lambda_women)


###
### Cox Regression bit
###

### Setup the data necessary to calculate the cox regression on the prediction set.
og_data <- read.csv("PHIA0000232016_IA88_25Jul/PHIA0000232016.csv")
og_data$new_ltpa <- mapply(FUN=sum, og_data$m_walk, og_data$m_floors, 
                           og_data$m_cycl, og_data$m_sport, og_data$m_houswrk, og_data$m_vigpa,
                           og_data$m_gard, og_data$m_diy, na.rm=TRUE)
og_data$bmi <- og_data$bmi_adj
og_data$sex <- unlist(lapply(og_data$sex, FUN=function(x){
  if (x == 1) {
    out = 0
  } else {
    out = 1
  }
  return(out)
}))

og_data$pa_workini <- unlist(lapply(og_data$pa_work, FUN=function(x){
  if (x == 6){
    out = 9
  } else {
    out = x
  }
  return (out)
}))

og_data$pa_workini <- factor(og_data$pa_workini)

## Calculating X_t, with X_t0
tempfupdiff <- (og_data$fup_time/365.25)
og_data$ageEnd <- og_data$age_recr_max + tempfupdiff

og_data$age <- mapply(og_data$age_recr_max, og_data$ageEnd, FUN=function(x,y){
  out = (x+y)/2
  return(out)
})

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

## smoke, lschool, country factorization and leveling
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


# Create Cambridge index matrix
cam_matrix = matrix(byrow = TRUE,
                    c(1, 2, 3, 4, 
                      2, 3, 4, 4, 
                      3, 4, 4, 4, 
                      4, 4, 4, 4,
                      1, 2, 3, 4), 
                    nrow=5, 
                    ncol=4)

og_data$cam_index <- apply(X = og_data[,c('pa_work', 'camMets_ind')], MARGIN = 1, 
                           FUN = function(x){
                             if(x[1]==6){
                               x_ind = 5
                             } else {
                               x_ind = x[1]
                             }
                             output <- cam_matrix[x_ind, x[2]]
                             return(output)
                           }
)

og_data$camMets_ind <- as.factor(og_data$camMets_ind)
og_data$cam_index <- as.factor(og_data$cam_index)

## Cam index means
og_data$cam_index_means <- unlist(mapply(og_data$cam_index, og_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
  if (y == 0) {
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = 35.6
    } else if (x == 2) {
      output = 43.7
    } else if (x == 3) {
      output = 49.0
    } else if (x == 4) {
      output = 56.2
    } else {
      output = NA
    }
  } else if (y == 1){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = 36.5
    } else if (x == 2) {
      output = 39.8
    } else if (x == 3) {
      output = 43.6
    } else if (x == 4) {
      output = 48.2
    } else {
      output = NA
    } 
  } else {
    output = NA
  }
  return(output)
}
))


og_men <- subset(og_data, sex==0)
og_women <- subset(og_data, sex==1)

cox_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = og_men, robust=TRUE)
cox_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = og_women)

beta_men <- cox_regression_men$coefficients
beta_se_men <- summary(cox_regression_men)$coefficients[4]

beta_women <- cox_regression_women$coefficients
beta_se_women <- summary(cox_regression_women)$coefficients[4]

upper_ci_men <- summary(cox_regression_men)$conf.int[4]
lower_ci_men <- summary(cox_regression_men)$conf.int[3]
upper_ci_women <- summary(cox_regression_women)$conf.int[4]
lower_ci_women <- summary(cox_regression_women)$conf.int[3]

upper_ci_list_men <- c(upper_ci_men)
lower_ci_list_men <- c(lower_ci_men)
upper_ci_list_women <- c(upper_ci_women)
lower_ci_list_women <- c(lower_ci_women)

# Store the lambda in a list for comparison
beta_list_men <- c(beta_men)
beta_list_women <- c(beta_women)

# calculate the calibrated beta and store in list
cali_beta_men <- beta_men/lambda_men
cali_beta_women <- beta_women/lambda_women

cali_beta_list_men <- c(cali_beta_men)
cali_beta_list_women <- c(cali_beta_women)

# Calibrated Confidence Interval calc helper
cal_beta_sd_men <- (beta_se_men^2)/(lambda_men^2) + ((beta_men/lambda_men^2)^2)*lambda_men_var # formula 12
cal_beta_sd_women <- (beta_se_women^2)/(lambda_women^2) + ((beta_women/lambda_women^2)^2)*lambda_women_var # formula 12

cal_confidence_interval_help_men <- 1.96*sqrt(cal_beta_sd_men)
cal_confidence_interval_help_women <- 1.96*sqrt(cali_beta_women)

cal_upper_ci_men <- cali_beta_men + cal_confidence_interval_help_men
cal_lower_ci_men <- cali_beta_men - cal_confidence_interval_help_men
cal_upper_ci_women <- cali_beta_women + cal_confidence_interval_help_women
cal_lower_ci_women <- cali_beta_women - cal_confidence_interval_help_women

cal_upper_ci_list_men <- c(cal_upper_ci_men)
cal_lower_ci_list_men <- c(cal_lower_ci_men)
cal_upper_ci_list_women <- c(cal_upper_ci_women)
cal_lower_ci_list_women <- c(cal_lower_ci_women)


###############################################################################################
####################################  SAMPLING RUN THROUGH  ###################################
###############################################################################################
# Set seed for random pull. 
set.seed(999)

## Set lists of values for each PA categorisation that will be used to draw random numbers for
## the validation set's observed PAEE values

cat1_men <- mapply(merged_output_men$PAEE, merged_output_men$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat1_men <- cat1_men[!sapply(cat1_men,is.na)]


cat2_men <- mapply(merged_output_men$PAEE, merged_output_men$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat2_men <- cat2_men[!sapply(cat2_men,is.na)]

cat3_men <- mapply(merged_output_men$PAEE, merged_output_men$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat3_men <- cat3_men[!sapply(cat3_men,is.na)]

cat4_men <- mapply(merged_output_men$PAEE, merged_output_men$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
})
cat4_men <- cat4_men[!sapply(cat4_men,is.na)]

cat1_women <- mapply(merged_output_women$PAEE, merged_output_women$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output)
})
cat1_women <- cat1_women[!sapply(cat1_women,is.na)]

cat2_women <- mapply(merged_output_women$PAEE, merged_output_women$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output)
}) 
cat2_women <- cat2_women[!sapply(cat2_women,is.na)]

cat3_women <- mapply(merged_output_women$PAEE, merged_output_women$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output)
}) 
cat3_women <- cat3_women[!sapply(cat3_women,is.na)]

cat4_women <- mapply(merged_output_women$PAEE, merged_output_women$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
}) 
cat4_women <- cat4_women[!sapply(cat4_women,is.na)]

## First run through we set the mean to be the first thing we sample per category
## and calculate the lm many times

for (i in 1:1000) {
  cat1_men_choice <- sample(cat1_men,1,replace=TRUE)
  cat2_men_choice <- sample(cat2_men,1,replace=TRUE)
  cat3_men_choice <- sample(cat3_men,1,replace=TRUE)
  cat4_men_choice <- sample(cat4_men,1,replace=TRUE)
  cat1_women_choice <- sample(cat1_women,1,replace=TRUE)
  cat2_women_choice <- sample(cat2_women,1,replace=TRUE)
  cat3_women_choice <- sample(cat3_women,1,replace=TRUE)
  cat4_women_choice <- sample(cat4_women,1,replace=TRUE)

  merged_output$cam_index_means <- unlist(mapply(merged_output$cam_index, merged_output$sex, SIMPLIFY = FALSE, FUN=function(x,y){
  if (y == 0) {
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = cat1_men_choice
    } else if (x == 2) {
      output = cat2_men_choice
    } else if (x == 3) {
      output = cat3_men_choice
    } else if (x == 4) {
      output = cat4_men_choice
    } else {
      output = NA
    }
  } else if (y == 1){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = cat1_women_choice
    } else if (x == 2) {
      output = cat2_women_choice
    } else if (x == 3) {
      output = cat3_women_choice
    } else if (x == 4) {
      output = cat4_women_choice
    } else {
      output = NA
    } 
  } else {
    output = NA
  }
  return(output)
  }
  ))

  ###
  ### Regression and splitting of the genders
  ###
  merged_output_men <- subset(merged_output, sex==0)
  merged_output_women <- subset(merged_output, sex==1)

  rdr_regression_fit_men <- lm(formula=PAEE~cam_index_means, data=merged_output_men)
  rdr_regression_fit_women <- lm(formula=PAEE~cam_index_means, data=merged_output_women)

  lambda_men <- rdr_regression_fit_men$coefficients[2]
  lambda_women <- rdr_regression_fit_women$coefficients[2]

  lambda_men_var <- (summary(rdr_regression_fit_men)$coefficients[4])^2
  lambda_women_var <- (summary(rdr_regression_fit_women)$coefficients[4])^2

  # Store the lambda in a list for comparison
  lambda_list_men <- c(lambda_list_men, lambda_men)
  lambda_list_women <- c(lambda_list_women, lambda_women)

  ## Cam index means
  og_data$cam_index_means <- unlist(mapply(og_data$cam_index, og_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
  if (y == 0) {
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = cat1_men_choice
    } else if (x == 2) {
      output = cat2_men_choice
    } else if (x == 3) {
      output = cat3_men_choice
    } else if (x == 4) {
      output = cat4_men_choice
    } else {
      output = NA
    }
  } else if (y == 1){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = cat1_women_choice
    } else if (x == 2) {
      output = cat2_women_choice
    } else if (x == 3) {
      output = cat3_women_choice
    } else if (x == 4) {
      output = cat4_women_choice
    } else {
      output = NA
    } 
  } else {
    output = NA
  }
  return(output)
  }
  ))

  # Predictions per model

  og_men <- subset(og_data, sex==0)
  og_women <- subset(og_data, sex==1)

  cox_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = og_men, robust=TRUE)
  cox_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = og_women)

  beta_men <- cox_regression_men$coefficients
  beta_se_men <- summary(cox_regression_men)$coefficients[4]

  beta_women <- cox_regression_women$coefficients
  beta_se_women <- summary(cox_regression_women)$coefficients[4]

  upper_ci_men <- summary(cox_regression_men)$conf.int[4]
  lower_ci_men <- summary(cox_regression_men)$conf.int[3]
  upper_ci_women <- summary(cox_regression_women)$conf.int[4]
  lower_ci_women <- summary(cox_regression_women)$conf.int[3]

  beta_list_men <- c(beta_list_men, beta_men)
  beta_list_women <- c(beta_list_women, beta_women)
  upper_ci_list_men <- c(upper_ci_list_men, upper_ci_men)
  lower_ci_list_men <- c(lower_ci_list_men, lower_ci_men)
  upper_ci_list_women <- c(upper_ci_list_women, upper_ci_women)
  lower_ci_list_women <- c(lower_ci_list_women, lower_ci_women)

  # calculate the calibrated beta and store in list
  cali_beta_men <- beta_men/lambda_men
  cali_beta_women <- beta_women/lambda_women

  cali_beta_list_men <- c(cali_beta_list_men, cali_beta_men)
  cali_beta_list_women <- c(cali_beta_list_women, cali_beta_women)

  # Calibrated Confidence Interval calc helper
  cal_beta_sd_men <- (beta_se_men^2)/(lambda_men^2) + ((beta_men/lambda_men^2)^2)*lambda_men_var # formula 12
  cal_beta_sd_women <- (beta_se_women^2)/(lambda_women^2) + ((beta_women/lambda_women^2)^2)*lambda_women_var # formula 12

  cal_confidence_interval_help_men <- 1.96*sqrt(cal_beta_sd_men)
  cal_confidence_interval_help_women <- 1.96*sqrt(cal_beta_sd_women)

  cal_upper_ci_men <- cali_beta_men + cal_confidence_interval_help_men
  cal_lower_ci_men <- cali_beta_men - cal_confidence_interval_help_men
  cal_upper_ci_women <- cali_beta_women + cal_confidence_interval_help_women
  cal_lower_ci_women <- cali_beta_women - cal_confidence_interval_help_women

  cal_upper_ci_list_men <- c(cal_upper_ci_list_men,cal_upper_ci_men)
  cal_lower_ci_list_men <- c(cal_lower_ci_list_men,cal_lower_ci_men)
  cal_upper_ci_list_women <- c(cal_upper_ci_list_women,cal_upper_ci_women)
  cal_lower_ci_list_women <- c(cal_lower_ci_list_women,cal_lower_ci_women)
}

# get rid of the silly names r puts on everything
lambda_list_men <- unname(lambda_list_men)
lambda_list_women <- unname(lambda_list_women)

beta_list_men <- unname(beta_list_men)
beta_list_women <- unname(beta_list_women)
upper_ci_list_men <- unname(upper_ci_list_men)
lower_ci_list_men <- unname(lower_ci_list_men)
upper_ci_list_women <- unname(upper_ci_list_women)
lower_ci_list_women <- unname(lower_ci_list_women)

cali_beta_list_men <- unname(cali_beta_list_men)
cali_beta_list_women <- unname(cali_beta_list_women)
cal_upper_ci_list_men <- unname(cal_upper_ci_list_men)
cal_lower_ci_list_men <- unname(cal_lower_ci_list_men)
cal_upper_ci_list_women <- unname(cal_upper_ci_list_women)
cal_lower_ci_list_women <- unname(cal_lower_ci_list_women)

### Put everything in a nice data frame
# men
men_dataframe <- data.frame(lambda = lambda_list_men, beta=beta_list_men, lower95=lower_ci_list_men, 
  upper95=upper_ci_list_men, calibratedBeta=cali_beta_list_men, calBetaLower95=cal_lower_ci_list_men,
  calBetaUpper95=cal_upper_ci_list_men)

# women
women_dataframe <- data.frame(lambda = lambda_list_women, beta=beta_list_women, lower95=lower_ci_list_women, 
  upper95=upper_ci_list_women, calibratedBeta=cali_beta_list_women, calBetaLower95=cal_lower_ci_list_women,
  calBetaUpper95=cal_upper_ci_list_women)

