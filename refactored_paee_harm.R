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

###############################################################################
############################### FUNCTIONS #####################################
###############################################################################
bmi_calc <- function(weight, height){
  bmi = (weight/height)/height
  return(bmi)
}

gaussian_index_sample <- function(x,gender){
  # returns a datapoint sampled from the index distribution
  # This can be seen as an exposure to be used in the data_generator
  # :param: index = cambridge index
  if (gender == 0){
    if (x == 1){
      index_mean <- index_mean1_men
      index_stdev <- index_stdev1_men
    } else if (x == 2){
      index_mean <- index_mean2_men
      index_stdev <- index_stdev2_men
    } else if (x == 3){
      index_mean <- index_mean3_men
      index_stdev <- index_stdev3_men
    } else {
      index_mean <- index_mean4_men
      index_stdev <- index_stdev4_men
    }
  } else {
    if (x == 1){
      index_mean <- index_mean1_women
      index_stdev <- index_stdev1_women
    } else if (x == 2){
      index_mean <- index_mean2_women
      index_stdev <- index_stdev2_women
    } else if (x == 3){
      index_mean <- index_mean3_women
      index_stdev <- index_stdev3_women
    } else {
      index_mean <- index_mean4_women
      index_stdev <- index_stdev4_women
    }
  }
  data_point <- rnorm(1, index_mean, index_stdev)
  return (data_point)
}


###############################################################################
################################# DATA ########################################
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
act_sorted <- actiheart_summary[order(actiheart_summary$universal_id, actiheart_summary$visit), ]
act_split <- split(x = act_sorted, f = act_sorted$universal_id)
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

# Merging of actiheart and epic set by natural join to create the validation dataset
val_data <- merge(output, epic, by = 'universal_id')
remove(list = c('epic', 'act_sorted', 'actiheart_summary', 'output')) #remove unnecessary dataframes
val_data$bmi <- mapply(FUN=bmi_calc, val_data$weight, val_data$height)
val_data$sex <- unlist(lapply(val_data$sex, FUN=function(x){
  if (x == 1) {
    out = 0
  } else {
    out = 1
  }
  return(out)
}))

# Make sure pa_workini is seen as a categorical variable and not a double
val_data$pa_workini <- as.factor(val_data$pa_workini)
val_data$new_ltpa <- mapply(FUN=sum, val_data$m_walk, val_data$m_floors, 
                                 val_data$m_cycl, val_data$m_sport, val_data$m_houswk, val_data$m_vigpa, 
                                 val_data$m_gard, val_data$m_diy)

## Calculating own Cambridge Index
val_data$cam_total <- c(rep(0,nrow(val_data)))
for (i in 1:nrow(val_data)){
  val_data$cam_total[i] = sum(val_data$m_cycl[i]/6, val_data$m_sport[i]/6, na.rm=TRUE)    
}

# In House Calculation of Cambridge Index (camMets_ind for cam_matrix)
val_data$camMets_ind <- as.vector(do.call(rbind, lapply(X=val_data$cam_total, FUN = function(x){
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

#  _____  _               _         ______        _         
# /  ___|| |             | |        |  _  \      | |        
# \ `--. | |_  _   _   __| | _   _  | | | | __ _ | |_  __ _ 
#  `--. \| __|| | | | / _` || | | | | | | |/ _` || __|/ _` |
# /\__/ /| |_ | |_| || (_| || |_| | | |/ /| (_| || |_| (_| |
# \____/  \__| \__,_| \__,_| \__, | |___/  \__,_| \__|\__,_|
#                             __/ |                         
#                            |___/                          
study_data <- read.csv("PHIA0000232016_IA88_25Jul/PHIA0000232016.csv")
study_data$new_ltpa <- mapply(FUN=sum, study_data$m_walk, study_data$m_floors, 
                           study_data$m_cycl, study_data$m_sport, study_data$m_houswrk, study_data$m_vigpa,
                           study_data$m_gard, study_data$m_diy, na.rm=TRUE)
study_data$bmi <- study_data$bmi_adj

study_data$sex <- unlist(lapply(study_data$sex, FUN=function(x){
  if (x == 1) {
    out = 0
  } else {
    out = 1
  }
  return(out)
}))

study_data$pa_workini <- unlist(lapply(study_data$pa_work, FUN=function(x){
  if (x == 6){
    out = 9
  } else {
    out = x
  }
  return (out)
}))

study_data$pa_workini <- factor(study_data$pa_workini)

## Calculating X_t, with X_t0
tempfupdiff <- (study_data$fup_time/365.25)
study_data$ageEnd <- study_data$age_recr_max + tempfupdiff

study_data$age <- mapply(study_data$age_recr_max, study_data$ageEnd, FUN=function(x,y){
  out = (x+y)/2
  return(out)
})

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

## smoke, lschool, country factorization and leveling
study_data$country <- factor(study_data$country, labels=c("FRANCE", "ITALY", 
                                                    "SPAIN", "UK", "NETHERLANDS", "GERMANY", "SWEDEN", "DENMARK"))
study_data$smoke_stat <- factor(study_data$smoke_stat, labels=c("NEVER", "FORMER", 
                                                          "SMOKER", "UNKOWN"))
study_data$l_school <- factor(study_data$l_school, labels=c("NONE", "PRIMARY", 
                                                      "TECHNICAL/PROFESSIONAL", "SECONDARY", "LONGER EDUCATION/UNI", "NOT SPECIFIED"))

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

# Create Cambridge index matrix
cam_matrix = matrix(byrow = TRUE,
                    c(1, 2, 3, 4, 
                      2, 3, 4, 4, 
                      3, 4, 4, 4, 
                      4, 4, 4, 4,
                      1, 2, 3, 4), 
                    nrow=5, 
                    ncol=4)

study_data$cam_index <- apply(X = study_data[,c('pa_work', 'camMets_ind')], MARGIN = 1, 
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

study_data$camMets_ind <- as.factor(study_data$camMets_ind)
study_data$cam_index <- as.factor(study_data$cam_index)

###############################################################################
################################# METHODS #####################################
###############################################################################

# ___  ___ _____ _____ _   _ ___________   __  
# |  \/  ||  ___|_   _| | | |  _  |  _  \ /  | 
# | .  . || |__   | | | |_| | | | | | | | `| | 
# | |\/| ||  __|  | | |  _  | | | | | | |  | | 
# | |  | || |___  | | | | | \ \_/ / |/ /  _| |_
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \___/

# Harmonisation using the paee means from each of the cambridge indices.

# Setting the means for the validation data
val_data$cam_index_means <- unlist(mapply(val_data$cam_index, val_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
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

# val_data$cam_index_means <- unlist(mapply(val_data$cam_index, val_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
#   if (y == 0) {
#     if (is.na(x)) {
#       output = NA
#     } else if (x == 1){
#       output = 32.49748
#     } else if (x == 2) {
#       output = 42.040295
#     } else if (x == 3) {
#       output = 45.94756
#     } else if (x == 4) {
#       output = 55.79739
#     } else {
#       output = NA
#     }
#   } else if (y == 1){
#     if (is.na(x)) {
#       output = NA
#     } else if (x == 1){
#       output = 35.41059
#     } else if (x == 2) {
#       output = 38.77947
#     } else if (x == 3) {
#       output = 42.77706
#     } else if (x == 4) {
#       output = 45.80403
#     } else {
#       output = NA
#     }
#   } else {
#     output = NA
#   }
#   return(output)
# }
# ))

# Set the means for the study data
study_data$cam_index_means <- unlist(mapply(study_data$cam_index, study_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
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

# study_data$cam_index_means <- unlist(mapply(study_data$cam_index, study_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
#   if (y == 0) {
#     if (is.na(x)) {
#       output = NA
#     } else if (x == 1){
#       output = 32.49748
#     } else if (x == 2) {
#       output = 42.040295
#     } else if (x == 3) {
#       output = 45.94756
#     } else if (x == 4) {
#       output = 55.79739
#     } else {
#       output = NA
#     }
#   } else if (y == 1){
#     if (is.na(x)) {
#       output = NA
#     } else if (x == 1){
#       output = 35.41059
#     } else if (x == 2) {
#       output = 38.77947
#     } else if (x == 3) {
#       output = 42.77706
#     } else if (x == 4) {
#       output = 45.80403
#     } else {
#       output = NA
#     }
#   } else {
#     output = NA
#   }
#   return(output)
# }
# ))

## Lambda Calculation
val_data_men <- subset(val_data, sex==0)
val_data_women <- subset(val_data, sex==1)

rdr_regression_fit_men <- lm(formula=PAEE~cam_index_means, data=val_data_men)
rdr_regression_fit_women <- lm(formula=PAEE~cam_index_means, data=val_data_women)

lambda_men <- rdr_regression_fit_men$coefficients["cam_index_means"]
lambda_women <- rdr_regression_fit_women$coefficients["cam_index_means"]

lambda_men_var <- (summary(rdr_regression_fit_men)$coefficients["cam_index_means","Std. Error"])^2
lambda_women_var <- (summary(rdr_regression_fit_women)$coefficients["cam_index_means","Std. Error"])^2

# Store the lambda in a list for comparison
lambda_list_men <- c(lambda_men)
lambda_list_women <- c(lambda_women)

## Beta Calculation
study_data_men <- subset(study_data, sex==0)
study_data_women <- subset(study_data, sex==1)

cox_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = study_data_men, robust=TRUE)
cox_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = study_data_women, robust=TRUE)

beta_men <- cox_regression_men$coefficients
beta_se_men <- summary(cox_regression_men)$coefficients[,"robust se"]

beta_women <- cox_regression_women$coefficients
beta_se_women <- summary(cox_regression_women)$coefficients[,"robust se"]

upper_ci_men <- log(summary(cox_regression_men)$conf.int[,"upper .95"])
lower_ci_men <- log(summary(cox_regression_men)$conf.int[,"lower .95"])
upper_ci_women <- log(summary(cox_regression_women)$conf.int[,"upper .95"])
lower_ci_women <- log(summary(cox_regression_women)$conf.int[,"lower .95"])

upper_ci_list_men <- c(upper_ci_men)
lower_ci_list_men <- c(lower_ci_men)
upper_ci_list_women <- c(upper_ci_women)
lower_ci_list_women <- c(lower_ci_women)

# Store the beta in a list for comparison
beta_list_men <- c(beta_men)
beta_list_women <- c(beta_women)

# calculate the calibrated beta and store in list
cal_beta_men <- beta_men/lambda_men
cal_beta_women <- beta_women/lambda_women

cal_beta_list_men <- c(cal_beta_men)
cal_beta_list_women <- c(cal_beta_women)

# Calibrated Confidence Interval calc helper
cal_beta_sd_men <- (beta_se_men^2)/(lambda_men^2) + ((beta_men/lambda_men^2)^2)*lambda_men_var # formula 12
cal_beta_sd_women <- (beta_se_women^2)/(lambda_women^2) + ((beta_women/lambda_women^2)^2)*lambda_women_var # formula 12

cal_confidence_interval_help_men <- 1.96*sqrt(cal_beta_sd_men)
cal_confidence_interval_help_women <- 1.96*sqrt(cal_beta_sd_women)

cal_upper_ci_men <- cal_beta_men + cal_confidence_interval_help_men
cal_lower_ci_men <- cal_beta_men - cal_confidence_interval_help_men
cal_upper_ci_women <- cal_beta_women + cal_confidence_interval_help_women
cal_lower_ci_women <- cal_beta_women - cal_confidence_interval_help_women

cal_upper_ci_list_men <- c(cal_upper_ci_men)
cal_lower_ci_list_men <- c(cal_lower_ci_men)
cal_upper_ci_list_women <- c(cal_upper_ci_women)
cal_lower_ci_list_women <- c(cal_lower_ci_women)




# ___  ___ _____ _____ _   _ ___________   _____   ___  
# |  \/  ||  ___|_   _| | | |  _  |  _  \ / __  \ / _ \ 
# | .  . || |__   | | | |_| | | | | | | | `' / /'/ /_\ \
# | |\/| ||  __|  | | |  _  | | | | | | |   / /  |  _  |
# | |  | || |___  | | | | | \ \_/ / |/ /  ./ /___| | | |
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \_____/\_| |_/                                     
###############################################################################################################
################################## WITH SAMPLED PAEE VALS #####################################################
# Set seed for pseudo random selection. 
set.seed(999)

## Set lists of values for each PA categorisation that will be used to draw random numbers for
## the validation set's observed PAEE values
cat1_men <- mapply(val_data_men$PAEE, val_data_men$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat1_men <- cat1_men[!sapply(cat1_men,is.na)]


cat2_men <- mapply(val_data_men$PAEE, val_data_men$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat2_men <- cat2_men[!sapply(cat2_men,is.na)]

cat3_men <- mapply(val_data_men$PAEE, val_data_men$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output) 
}) 
cat3_men <- cat3_men[!sapply(cat3_men,is.na)]

cat4_men <- mapply(val_data_men$PAEE, val_data_men$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
})
cat4_men <- cat4_men[!sapply(cat4_men,is.na)]

cat1_women <- mapply(val_data_women$PAEE, val_data_women$cam_index, FUN=function(x,y){
  if (y == 1){
    output = x
  } else {
    output = NA
  }
  return(output)
})
cat1_women <- cat1_women[!sapply(cat1_women,is.na)]

cat2_women <- mapply(val_data_women$PAEE, val_data_women$cam_index, FUN=function(x,y){
  if (y == 2){
    output = x
  } else {
    output = NA
  }
  return(output)
}) 
cat2_women <- cat2_women[!sapply(cat2_women,is.na)]

cat3_women <- mapply(val_data_women$PAEE, val_data_women$cam_index, FUN=function(x,y){
  if (y == 3){
    output = x
  } else {
    output = NA
  }
  return(output)
}) 
cat3_women <- cat3_women[!sapply(cat3_women,is.na)]

cat4_women <- mapply(val_data_women$PAEE, val_data_women$cam_index, FUN=function(x,y){
  if (y == 4){
    output = x
  } else {
    output = NA
  }
  return(output)
}) 
cat4_women <- cat4_women[!sapply(cat4_women,is.na)]

##
## First run through we set the mean to be the first thing we sample per category
## and calculate the lm many times
##
for (i in 1:10) {
  cat1_men_choice <- sample(cat1_men,1,replace=TRUE)
  cat2_men_choice <- sample(cat2_men,1,replace=TRUE)
  cat3_men_choice <- sample(cat3_men,1,replace=TRUE)
  cat4_men_choice <- sample(cat4_men,1,replace=TRUE)
  cat1_women_choice <- sample(cat1_women,1,replace=TRUE)
  cat2_women_choice <- sample(cat2_women,1,replace=TRUE)
  cat3_women_choice <- sample(cat3_women,1,replace=TRUE)
  cat4_women_choice <- sample(cat4_women,1,replace=TRUE)

  val_data$cam_index_means <- unlist(mapply(val_data$cam_index, val_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
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
  val_data_men <- subset(val_data, sex==0)
  val_data_women <- subset(val_data, sex==1)

  rdr_regression_fit_men <- lm(formula=PAEE~cam_index_means, data=val_data_men)
  rdr_regression_fit_women <- lm(formula=PAEE~cam_index_means, data=val_data_women)

  lambda_men <- rdr_regression_fit_men$coefficients["cam_index_means"]
  lambda_women <- rdr_regression_fit_women$coefficients["cam_index_means"]

  lambda_men_var <- (summary(rdr_regression_fit_men)$coefficients["cam_index_means","Std. Error"])^2
  lambda_women_var <- (summary(rdr_regression_fit_women)$coefficients["cam_index_means","Std. Error"])^2

  # Store the lambda in a list for comparison
  lambda_list_men <- c(lambda_list_men, lambda_men)
  lambda_list_women <- c(lambda_list_women, lambda_women)

  ## Cam index means
  study_data$cam_index_means <- unlist(mapply(study_data$cam_index, study_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
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
  study_data_men <- subset(study_data, sex==0)
  study_data_women <- subset(study_data, sex==1)

  cox_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = study_data_men, robust=TRUE)
  cox_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = study_data_women, robust=TRUE)

  beta_men <- cox_regression_men$coefficients
  beta_se_men <- summary(cox_regression_men)$coefficients[,"robust se"]

  beta_women <- cox_regression_women$coefficients
  beta_se_women <- summary(cox_regression_women)$coefficients[,"robust se"]

  upper_ci_men <- log(summary(cox_regression_men)$conf.int[,"upper .95"])
  lower_ci_men <- log(summary(cox_regression_men)$conf.int[,"lower .95"])
  upper_ci_women <- log(summary(cox_regression_women)$conf.int[,"upper .95"])
  lower_ci_women <- log(summary(cox_regression_women)$conf.int[,"lower .95"])

  beta_list_men <- c(beta_list_men, beta_men)
  beta_list_women <- c(beta_list_women, beta_women)
  upper_ci_list_men <- c(upper_ci_list_men, upper_ci_men)
  lower_ci_list_men <- c(lower_ci_list_men, lower_ci_men)
  upper_ci_list_women <- c(upper_ci_list_women, upper_ci_women)
  lower_ci_list_women <- c(lower_ci_list_women, lower_ci_women)

  # calculate the calibrated beta and store in list
  cal_beta_men <- beta_men/lambda_men
  cal_beta_women <- beta_women/lambda_women

  cal_beta_list_men <- c(cal_beta_list_men, cal_beta_men)
  cal_beta_list_women <- c(cal_beta_list_women, cal_beta_women)

  # Calibrated Confidence Interval calc helper
  cal_beta_sd_men <- (beta_se_men^2)/(lambda_men^2) + ((beta_men/lambda_men^2)^2)*lambda_men_var # formula 12
  cal_beta_sd_women <- (beta_se_women^2)/(lambda_women^2) + ((beta_women/lambda_women^2)^2)*lambda_women_var # formula 12

  cal_confidence_interval_help_men <- 1.96*sqrt(cal_beta_sd_men)
  cal_confidence_interval_help_women <- 1.96*sqrt(cal_beta_sd_women)

  cal_upper_ci_men <- cal_beta_men + cal_confidence_interval_help_men
  cal_lower_ci_men <- cal_beta_men - cal_confidence_interval_help_men
  cal_upper_ci_women <- cal_beta_women + cal_confidence_interval_help_women
  cal_lower_ci_women <- cal_beta_women - cal_confidence_interval_help_women

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

cal_beta_list_men <- unname(cal_beta_list_men)
cal_beta_list_women <- unname(cal_beta_list_women)
cal_upper_ci_list_men <- unname(cal_upper_ci_list_men)
cal_lower_ci_list_men <- unname(cal_lower_ci_list_men)
cal_upper_ci_list_women <- unname(cal_upper_ci_list_women)
cal_lower_ci_list_women <- unname(cal_lower_ci_list_women)


### Put everything in a data frame
# men
men_dataframe <- data.frame(lambda = lambda_list_men, beta=beta_list_men, lower95=lower_ci_list_men, 
  upper95=upper_ci_list_men, calibratedBeta=cal_beta_list_men, calBetaLower95=cal_lower_ci_list_men,
  calBetaUpper95=cal_upper_ci_list_men)

# women
women_dataframe <- data.frame(lambda = lambda_list_women, beta=beta_list_women, lower95=lower_ci_list_women, 
  upper95=upper_ci_list_women, calibratedBeta=cal_beta_list_women, calBetaLower95=cal_lower_ci_list_women,
  calBetaUpper95=cal_upper_ci_list_women)


# ___  ___ _____ _____ _   _ ___________   _____ ______ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \ / __  \| ___ \
# | .  . || |__   | | | |_| | | | | | | | `' / /'| |_/ /
# | |\/| ||  __|  | | |  _  | | | | | | |   / /  | ___ \
# | |  | || |___  | | | | | \ \_/ / |/ /  ./ /___| |_/ /
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \_____/\____/ 

for (i in 1:10) {
  ###
  ### Lambda Calculation
  ###
  val_data_men <- subset(val_data, sex==0)
  val_data_women <- subset(val_data, sex==1)

  rdr_regression_fit_men <- lm(formula=PAEE~cam_index_means, data=val_data_men)
  rdr_regression_fit_women <- lm(formula=PAEE~cam_index_means, data=val_data_women)

  lambda_men <- rdr_regression_fit_men$coefficients["cam_index_means"]
  lambda_women <- rdr_regression_fit_women$coefficients["cam_index_means"]

  lambda_men_var <- (summary(rdr_regression_fit_men)$coefficients["cam_index_means","Std. Error"])^2
  lambda_women_var <- (summary(rdr_regression_fit_women)$coefficients["cam_index_means","Std. Error"])^2

  # Store the lambda in a list for comparison
  lambda_list_men <- c(lambda_list_men, lambda_men)
  lambda_list_women <- c(lambda_list_women, lambda_women)

  ##
  ## Beta Calculation
  ##
  study_data$cam_index_means <- unlist(mapply(study_data$cam_index, study_data$sex, SIMPLIFY = FALSE, FUN=function(x,y){
  if (y == 0) {
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = sample(cat1_men,1,replace=TRUE)
    } else if (x == 2) {
      output = sample(cat2_men,1,replace=TRUE)
    } else if (x == 3) {
      output = sample(cat3_men,1,replace=TRUE)
    } else if (x == 4) {
      output = sample(cat4_men,1,replace=TRUE)
    } else {
      output = NA
    }
  } else if (y == 1){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      output = sample(cat1_women,1,replace=TRUE)
    } else if (x == 2) {
      output = sample(cat2_women,1,replace=TRUE)
    } else if (x == 3) {
      output = sample(cat3_women,1,replace=TRUE)
    } else if (x == 4) {
      output = sample(cat4_women,1,replace=TRUE)
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
  study_data_men <- subset(study_data, sex==0)
  study_data_women <- subset(study_data, sex==1)

  cox_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = study_data_men, robust=TRUE)
  cox_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = study_data_women, robust=TRUE)

  beta_men <- cox_regression_men$coefficients
  beta_se_men <- summary(cox_regression_men)$coefficients[,"robust se"]

  beta_women <- cox_regression_women$coefficients
  beta_se_women <- summary(cox_regression_women)$coefficients[,"robust se"]

  upper_ci_men <- log(summary(cox_regression_men)$conf.int[,"upper .95"])
  lower_ci_men <- log(summary(cox_regression_men)$conf.int[,"lower .95"])
  upper_ci_women <- log(summary(cox_regression_women)$conf.int[,"upper .95"])
  lower_ci_women <- log(summary(cox_regression_women)$conf.int[,"lower .95"])

  beta_list_men <- c(beta_list_men, beta_men)
  beta_list_women <- c(beta_list_women, beta_women)
  upper_ci_list_men <- c(upper_ci_list_men, upper_ci_men)
  lower_ci_list_men <- c(lower_ci_list_men, lower_ci_men)
  upper_ci_list_women <- c(upper_ci_list_women, upper_ci_women)
  lower_ci_list_women <- c(lower_ci_list_women, lower_ci_women)

  # calculate the calibrated beta and store in list
  cal_beta_men <- beta_men/lambda_men
  cal_beta_women <- beta_women/lambda_women

  cal_beta_list_men <- c(cal_beta_list_men, cal_beta_men)
  cal_beta_list_women <- c(cal_beta_list_women, cal_beta_women)

  # Calibrated Confidence Interval calc helper
  cal_beta_sd_men <- (beta_se_men^2)/(lambda_men^2) + ((beta_men/lambda_men^2)^2)*lambda_men_var # formula 12
  cal_beta_sd_women <- (beta_se_women^2)/(lambda_women^2) + ((beta_women/lambda_women^2)^2)*lambda_women_var # formula 12

  cal_confidence_interval_help_men <- 1.96*sqrt(cal_beta_sd_men)
  cal_confidence_interval_help_women <- 1.96*sqrt(cal_beta_sd_women)

  cal_upper_ci_men <- cal_beta_men + cal_confidence_interval_help_men
  cal_lower_ci_men <- cal_beta_men - cal_confidence_interval_help_men
  cal_upper_ci_women <- cal_beta_women + cal_confidence_interval_help_women
  cal_lower_ci_women <- cal_beta_women - cal_confidence_interval_help_women

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

cal_beta_list_men <- unname(cal_beta_list_men)
cal_beta_list_women <- unname(cal_beta_list_women)
cal_upper_ci_list_men <- unname(cal_upper_ci_list_men)
cal_lower_ci_list_men <- unname(cal_lower_ci_list_men)
cal_upper_ci_list_women <- unname(cal_upper_ci_list_women)
cal_lower_ci_list_women <- unname(cal_lower_ci_list_women)


### Put everything in a data frame
# men
men_dataframe <- data.frame(lambda = lambda_list_men, beta=beta_list_men, lower95=lower_ci_list_men, 
  upper95=upper_ci_list_men, calibratedBeta=cal_beta_list_men, calBetaLower95=cal_lower_ci_list_men,
  calBetaUpper95=cal_upper_ci_list_men)

# women
women_dataframe <- data.frame(lambda = lambda_list_women, beta=beta_list_women, lower95=lower_ci_list_women, 
  upper95=upper_ci_list_women, calibratedBeta=cal_beta_list_women, calBetaLower95=cal_lower_ci_list_women,
  calBetaUpper95=cal_upper_ci_list_women)

# ___  ___ _____ _____ _   _ ___________   _____ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \ |____ |
# | .  . || |__   | | | |_| | | | | | | |     / /
# | |\/| ||  __|  | | |  _  | | | | | | |     \ \
# | |  | || |___  | | | | | \ \_/ / |/ /  .___/ /
# \_|  |_/\____/  \_/ \_| |_/\___/|___/   \____/     
# We fit distributions to the validation data and then generate samples from these for the study data

# fitting of distributions via MLE
library(fitdistrplus)
# normal distributions
dist1_men <- fitdist(cat1_men, "norm", method='mle')
dist2_men <- fitdist(cat2_men, "norm", method='mle')
dist3_men <- fitdist(cat3_men, "norm", method='mle')
dist4_men <- fitdist(cat4_men, "norm", method='mle')

dist1_women <- fitdist(cat1_women, "norm", method='mle')
dist2_women <- fitdist(cat2_women, "norm", method='mle')
dist3_women <- fitdist(cat3_women, "norm", method='mle')
dist4_women <- fitdist(cat4_women, "norm", method='mle')

dist_list = list(dist1_men,dist2_men,dist3_men,dist4_men,dist1_women,dist2_women,dist3_women,dist4_women)

for (i in 1:8){
  plot(dist_list[[i]])
}

index_mean1_men <- dist1_men[[1]][1]
index_stdev1_men <- dist1_men[[1]][2]
index_mean2_men <- dist2_men[[1]][1]
index_stdev2_men <- dist2_men[[1]][2]
index_mean3_men <- dist3_men[[1]][1]
index_stdev3_men <- dist3_men[[1]][2]
index_mean4_men <- dist4_men[[1]][1]
index_stdev4_men <- dist4_men[[1]][2]
mean_paee_difference_men <- abs(mean(c(index_mean1_men-index_mean2_men, index_mean2_men-index_mean3_men, index_mean3_men-index_mean4_men)))

index_mean1_women <- dist1_women[[1]][1]
index_stdev1_women <- dist1_women[[1]][2]
index_mean2_women <- dist2_women[[1]][1]
index_stdev2_women <- dist2_women[[1]][2]
index_mean3_women <- dist3_women[[1]][1]
index_stdev3_women <- dist3_women[[1]][2]
index_mean4_women <- dist4_women[[1]][1]
index_stdev4_women <- dist4_women[[1]][2]
mean_paee_difference_women <- abs(mean(c(index_mean1_women-index_mean2_women, index_mean2_women-index_mean3_women, index_mean3_women-index_mean4_women)))

val_data_men$cam_index_dfit <- unlist(mapply(val_data_men$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else {
      output = gaussian_index_sample(x,gender=0)
    }
    return(output)
  }
))

study_data_men$cam_index_dfit <- unlist(mapply(study_data_men$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else {
      output = gaussian_index_sample(x,gender=0)
    }
    return(output)
  }
))

val_data_women$cam_index_dfit <- unlist(mapply(val_data_women$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else {
      output = gaussian_index_sample(x,gender=1)
    }
    return(output)
  }
))

study_data_women$cam_index_dfit <- unlist(mapply(study_data_women$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else {
      output = gaussian_index_sample(x,gender=1)
    }
    return(output)
  }
))

### Men
## coef = beta
pc_study_data_dfit_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_dfit, data = study_data_men, robust=TRUE)
cc_study_data_dfit_men <- pc_study_data_dfit_men$coefficients["cam_index_dfit"]
# stdError
stdError_cc_study_data_dfit_men <- (summary(pc_study_data_dfit_men)$coefficients[,"robust se"])[-1]
# confidence intervals
lci_beta_hat_star_dfit_men <- cc_study_data_dfit_men - 1.96*stdError_cc_study_data_dfit_men
uci_beta_hat_star_dfit_men <- cc_study_data_dfit_men + 1.96*stdError_cc_study_data_dfit_men

## coef = lambda
pc_val_data_dfit_men <- lm(formula=PAEE~cam_index_dfit, data=val_data_men)
cc_val_data_dfit_men <- pc_val_data_dfit_men$coefficients["cam_index_dfit"]
# stdError
stdError_cc_val_data_dfit_men <- (summary(pc_val_data_dfit_men)$coefficients[,"Std. Error"])[-1]
# confidence intervals
lci_lambda_hat_dfit_men <- cc_val_data_dfit_men - 1.96*(stdError_cc_val_data_dfit_men)
uci_lambda_hat_dfit_men <- cc_val_data_dfit_men + 1.96*(stdError_cc_val_data_dfit_men)

# RC 
beta_hat_star_dfit_men <- cc_study_data_dfit_men
lambda_hat_dfit_men <- cc_val_data_dfit_men
beta_hat_dfit_men <- beta_hat_star_dfit_men/lambda_hat_dfit_men
var_beta_dfit_men <- stdError_cc_study_data_dfit_men/(lambda_hat_dfit_men^2) + (beta_hat_star_dfit_men/lambda_hat_dfit_men^2)^2*stdError_cc_val_data_dfit_men
# confidence intervals
uci_beta_hat_dfit_men <- beta_hat_dfit_men - 1.96*sqrt(var_beta_dfit_men)
lci_beta_hat_dfit_men <- beta_hat_dfit_men + 1.96*sqrt(var_beta_dfit_men)

### Women
## coef = beta
pc_study_data_dfit_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_dfit, data = study_data_women, robust=TRUE)
cc_study_data_dfit_women <- pc_study_data_dfit_women$coefficients["cam_index_dfit"]
# stdError
stdError_cc_study_data_dfit_women <- (summary(pc_study_data_dfit_women)$coefficients[,"robust se"])[-1]
# confidence intervals
lci_beta_hat_star_dfit_women <- cc_study_data_dfit_women - 1.96*stdError_cc_study_data_dfit_women
uci_beta_hat_star_dfit_women <- cc_study_data_dfit_women + 1.96*stdError_cc_study_data_dfit_women

## coef = lambda
pc_val_data_dfit_women <- lm(formula=PAEE~cam_index_dfit, data=val_data_women)
cc_val_data_dfit_women <- pc_val_data_dfit_women$coefficients["cam_index_dfit"]
# stdError
stdError_cc_val_data_dfit_women <- (summary(pc_val_data_dfit_women)$coefficients[,"Std. Error"])[-1]
# confidence intervals
lci_lambda_hat_dfit_women <- cc_val_data_dfit_women - 1.96*(stdError_cc_val_data_dfit_women)
uci_lambda_hat_dfit_women <- cc_val_data_dfit_women + 1.96*(stdError_cc_val_data_dfit_women)

# RC 
beta_hat_star_dfit_women <- cc_study_data_dfit_women
lambda_hat_dfit_women <- cc_val_data_dfit_women
beta_hat_dfit_women <- beta_hat_star_dfit_women/lambda_hat_dfit_women
var_beta_dfit_women <- stdError_cc_study_data_dfit_women/(lambda_hat_dfit_women^2) + (beta_hat_star_dfit_women/lambda_hat_dfit_women^2)^2*stdError_cc_val_data_dfit_women
# confidence intervals
uci_beta_hat_dfit_women <- beta_hat_dfit_women - 1.96*sqrt(var_beta_dfit_women)
lci_beta_hat_dfit_women <- beta_hat_dfit_women + 1.96*sqrt(var_beta_dfit_women)

# ___  ___ _____ _____ _   _ ___________     ___ 
# |  \/  ||  ___|_   _| | | |  _  |  _  \   /   |
# | .  . || |__   | | | |_| | | | | | | |  / /| |
# | |\/| ||  __|  | | |  _  | | | | | | | / /_| |
# | |  | || |___  | | | | | \ \_/ / |/ /  \___  |
# \_|  |_/\____/  \_/ \_| |_/\___/|___/       |_/
#                                                
# We estimate the density of the pdf using kernel based estimation methods. Then sample directly using the kernels.
# fitting of distributions via MLE
library(stats)
# normal distributions
density1_men <- density(cat1_men)
density2_men <- density(cat2_men)
density3_men <- density(cat3_men)
density4_men <- density(cat4_men)

density1_women <- density(cat1_women)
density2_women <- density(cat2_women)
density3_women <- density(cat3_women)
density4_women <- density(cat4_women)

density_list = list(density1_men,density2_men,density3_men,density4_men,density1_women,density2_women,density3_women,density4_women)

for (i in 1:8){
  plot(density_list[[i]])
}

# Bandwidths necessary for sampling
bw1_men <- density1_men$bw
bw2_men <- density2_men$bw
bw3_men <- density3_men$bw
bw4_men <- density4_men$bw

bw1_women <- density1_women$bw
bw2_women <- density2_women$bw
bw3_women <- density3_women$bw
bw4_women <- density4_women$bw

# Sampling for each participant
study_data_men$cam_kernelEst <- unlist(mapply(study_data_men$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      km_1_men <- sample(cat1_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_1_men, sd = bw1_men)
    } else if (x == 2) {
      km_2_men <- sample(cat2_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_2_men, sd = bw2_men)
    } else if (x == 3) {
      km_3_men <- sample(cat3_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_3_men, sd = bw3_men)
    } else if (x == 4) {
      km_4_men <- sample(cat4_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_4_men, sd = bw4_men)
    } else {
      output = NA
    }
    return(output)
  }
))

val_data_men$cam_kernelEst <- unlist(mapply(val_data_men$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      km_1_men <- sample(cat1_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_1_men, sd = bw1_men)
    } else if (x == 2) {
      km_2_men <- sample(cat2_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_2_men, sd = bw2_men)
    } else if (x == 3) {
      km_3_men <- sample(cat3_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_3_men, sd = bw3_men)
    } else if (x == 4) {
      km_4_men <- sample(cat4_men, 1, replace=TRUE)
      output = rnorm(1, mean = km_4_men, sd = bw4_men)
    } else {
      output = NA
    }
    return(output)
  }
))

study_data_women$cam_kernelEst <- unlist(mapply(study_data_women$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      km_1_women <- sample(cat1_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_1_women, sd = bw1_women)
    } else if (x == 2) {
      km_2_women <- sample(cat2_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_2_women, sd = bw2_women)
    } else if (x == 3) {
      km_3_women <- sample(cat3_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_3_women, sd = bw3_women)
    } else if (x == 4) {
      km_4_women <- sample(cat4_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_4_women, sd = bw4_women)
    } else {
      output = NA
    }
    return(output)
  }
))

val_data_women$cam_kernelEst <- unlist(mapply(val_data_women$cam_index, SIMPLIFY = FALSE, FUN=function(x){
    if (is.na(x)) {
      output = NA
    } else if (x == 1){
      km_1_women <- sample(cat1_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_1_women, sd = bw1_women)
    } else if (x == 2) {
      km_2_women <- sample(cat2_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_2_women, sd = bw2_women)
    } else if (x == 3) {
      km_3_women <- sample(cat3_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_3_women, sd = bw3_women)
    } else if (x == 4) {
      km_4_women <- sample(cat4_women, 1, replace=TRUE)
      output = rnorm(1, mean = km_4_women, sd = bw4_women)
    } else {
      output = NA
    }
    return(output)
  }
))


############################################################################################
############################ Pretty Printing Results #######################################
############################################################################################

##
## By index (by category) (without calibration) - we want the ci (no lambda)
##
rdr_regression_fit_men <- lm(formula=PAEE~cam_index, data=val_data_men)
rdr_regression_fit_women <- lm(formula=PAEE~cam_index, data=val_data_women)
summary(rdr_regression_fit_men)
summary(rdr_regression_fit_women)

# Beta
cox_regression_men_index <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index, data = study_data_men, robust=TRUE)
cox_regression_women_index <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index, data = study_data_women, robust=TRUE)
beta_inc_men <- sum(summary(cox_regression_men_index)$coef[,1])/4
beta_inc_women <- sum(summary(cox_regression_women_index)$coef[,1])/4

l95_inc_men <- sum(log(summary(cox_regression_men_index)$conf.int[,3]))/4
l95_inc_women <- sum(log(summary(cox_regression_women_index)$conf.int[,3]))/4

u95_inc_men <- sum(log(summary(cox_regression_men_index)$conf.int[,4]))/4
u95_inc_women <- sum(log(summary(cox_regression_women_index)$conf.int[,4]))/4

ex_beta_inc_men <- exp(sum(summary(cox_regression_men_index)$coef[,1])/4)
ex_beta_inc_women <- exp(sum(summary(cox_regression_women_index)$coef[,1])/4)

ex_l95_inc_men <-exp( sum(log(summary(cox_regression_men_index)$conf.int[,3]))/4)
ex_l95_inc_women <- exp(sum(log(summary(cox_regression_women_index)$conf.int[,3]))/4)

ex_u95_inc_men <-exp( sum(log(summary(cox_regression_men_index)$conf.int[,4]))/4)
ex_u95_inc_women <- exp(sum(log(summary(cox_regression_women_index)$conf.int[,4]))/4)

##
## Method 1
##

# then by means (with calibration)
# this beta is pre calculated as the first row of the men/women_dataframe
# The result for uncalibrated and calibrated are present in this dataframe 
beta_by_means_table_men <- (as.numeric(men_dataframe[1,]))
beta_by_means_table_women <- as.numeric(women_dataframe[1,])
beta_by_means_table <- data.frame(rbind(beta_by_means_table_men,beta_by_means_table_women))
rownames(beta_by_means_table) <- c("male", "female")
colnames(beta_by_means_table) <- c("lambda", "beta", "betaLower95", "betaUpper95", "calibratedBeta", "calBetaLower95", "calBetaUpper95")
exponentiated_by_means_table <- exp(beta_by_means_table)^6.8 

##
## Method 2 (A or B)
##
# Because R makes summary dataframes into multinested char arrays (non ideal) we deal with each one seperately
summary(men_dataframe$lambda)["Median"]
summary(men_dataframe$beta)["Median"]  
summary(men_dataframe$lower95)["Median"]
summary(men_dataframe$upper95)["Median"]
summary(men_dataframe$calibratedBeta)["Median"]
summary(men_dataframe$calBetaLower95)["Median"]
summary(men_dataframe$calBetaUpper95)["Median"]

summary(women_dataframe$lambda)["Median"]
summary(women_dataframe$beta)["Median"]  
summary(women_dataframe$lower95)["Median"]
summary(women_dataframe$upper95)["Median"]
summary(women_dataframe$calibratedBeta)["Median"]
summary(women_dataframe$calBetaLower95)["Median"]
summary(women_dataframe$calBetaUpper95)["Median"]

exp(summary(men_dataframe$lambda)["Median"])^6.8
exp(summary(men_dataframe$beta)["Median"])^6.8  
exp(summary(men_dataframe$lower95)["Median"])^6.8
exp(summary(men_dataframe$upper95)["Median"])^6.8
exp(summary(men_dataframe$calibratedBeta)["Median"])^6.8
exp(summary(men_dataframe$calBetaLower95)["Median"])^6.8
exp(summary(men_dataframe$calBetaUpper95)["Median"])^6.8

exp(summary(women_dataframe$lambda)["Median"])^6.8
exp(summary(women_dataframe$beta)["Median"])^6.8  
exp(summary(women_dataframe$lower95)["Median"])^6.8
exp(summary(women_dataframe$upper95)["Median"])^6.8
exp(summary(women_dataframe$calibratedBeta)["Median"])^6.8
exp(summary(women_dataframe$calBetaLower95)["Median"])^6.8
exp(summary(women_dataframe$calBetaUpper95)["Median"])^6.8

##
## Method 3
##

##
## Method 4
##

