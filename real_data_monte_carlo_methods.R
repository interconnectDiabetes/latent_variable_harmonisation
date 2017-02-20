## Title: Harmonisation of PAEE

# Preceded by initialFlow.R.
# This script runs univariate linear harmonisation methods on different categorical representations of physical activity
# using 4 different methods to generate samples.

## Papers to find comparison results:
## 1. Physical activity reduces the risk of incident type 2 diabetes in 
## general and in abdominally lean and obese men and women: the EPIC-InterAct 
## Study
## 2. Validity of a short questionnaire to assess physical activity in 10 European
## countries 
## Author: Paul Scherer
## Date: 10/11/2016
## University of Cambridge


#TODO: Some serious refactoring of the final parts to make it less verbose. 

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
  if (gender == 1){
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
    if (x == 0){
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

gaussian_index_sample_bin <- function(x,gender){
  # returns a datapoint sampled from the index distribution
  # This can be seen as an exposure to be used in the data_generator
  if (gender == 1){
    if (x == 1){
      index_mean <- mean_bin1_men
      index_stdev <- stdev_bin1_men
    } else {
      index_mean <- mean_bin2_men
      index_stdev <- stdev_bin2_men 
    }
  } else {
    if (x == 1){
      index_mean <- mean_bin1_women
      index_stdev <- stdev_bin1_women
    } else {
      index_mean <- mean_bin2_women
      index_stdev <- stdev_bin2_women
    } 
  }
  data_point <- rnorm(1, index_mean, index_stdev)
  return (data_point)
}

gaussian_index_sample_reg <- function(x,y,gender){
  # returns a datapoint sampled from the index distribution
  # This can be seen as an exposure to be used in the data_generator
  if (gender == 1){
    if (x == 1){
      if (y == 1){
        index_mean <- mean_reg1_1_men
        index_stdev <- stdev_reg1_1_men
      }
      else if (y==2){
        index_mean <- mean_reg1_2_men
        index_stdev <- stdev_reg1_2_men
      }
      else if (y==3){
        index_mean <- mean_reg1_3_men
        index_stdev <- stdev_reg1_3_men
      }
      else if (y==4){
        index_mean <- mean_reg1_4_men
        index_stdev <- stdev_reg1_4_men
      }
      else if (y==5){
        index_mean <- mean_reg1_5_men
        index_stdev <- stdev_reg1_5_men
      }
    } else if (x == 2){
      if (y == 1){
        index_mean <- mean_reg2_1_men
        index_stdev <- stdev_reg2_1_men
      }
      else if (y==2){
        index_mean <- mean_reg2_2_men
        index_stdev <- stdev_reg2_2_men
      }
      else if (y==3){
        index_mean <- mean_reg2_3_men
        index_stdev <- stdev_reg2_3_men
      }
      else if (y==4){
        index_mean <- mean_reg2_4_men
        index_stdev <- stdev_reg2_4_men
      }
      else if (y==5){
        index_mean <- mean_reg2_5_men
        index_stdev <- stdev_reg2_5_men
      }
    } else if (x == 3){
      if (y == 1){
        index_mean <- mean_reg3_1_men
        index_stdev <- stdev_reg3_1_men
      }
      else if (y==2){
        index_mean <- mean_reg3_2_men
        index_stdev <- stdev_reg3_2_men
      }
      else if (y==3){
        index_mean <- mean_reg3_3_men
        index_stdev <- stdev_reg3_3_men
      }
      else if (y==4){
        index_mean <- mean_reg3_4_men
        index_stdev <- stdev_reg3_4_men
      }
      else if (y==5){
        index_mean <- mean_reg3_5_men
        index_stdev <- stdev_reg3_5_men
      }
    } else {
      if (y == 1){
        index_mean <- mean_reg4_1_men
        index_stdev <- stdev_reg4_1_men
      }
      else if (y==2){
        index_mean <- mean_reg4_2_men
        index_stdev <- stdev_reg4_2_men
      }
      else if (y==3){
        index_mean <- mean_reg4_3_men
        index_stdev <- stdev_reg4_3_men
      }
      else if (y==4){
        index_mean <- mean_reg4_4_men
        index_stdev <- stdev_reg4_4_men
      }
      else if (y==5){
        index_mean <- mean_reg4_5_men
        index_stdev <- stdev_reg4_5_men
      }
    }
  } else {
    if (x == 0){
      if (y == 1){
        index_mean <- mean_reg1_1_women
        index_stdev <- stdev_reg1_1_women
      }
      else if (y==2){
        index_mean <- mean_reg1_2_women
        index_stdev <- stdev_reg1_2_women
      }
      else if (y==3){
        index_mean <- mean_reg1_3_women
        index_stdev <- stdev_reg1_3_women
      }
      else if (y==4){
        index_mean <- mean_reg1_4_women
        index_stdev <- stdev_reg1_4_women
      }
      else if (y==5){
        index_mean <- mean_reg1_5_women
        index_stdev <- stdev_reg1_5_women
      }
    } else if (x == 2){
      if (y == 1){
        index_mean <- mean_reg2_1_women
        index_stdev <- stdev_reg2_1_women
      }
      else if (y==2){
        index_mean <- mean_reg2_2_women
        index_stdev <- stdev_reg2_2_women
      }
      else if (y==3){
        index_mean <- mean_reg2_3_women
        index_stdev <- stdev_reg2_3_women
      }
      else if (y==4){
        index_mean <- mean_reg2_4_women
        index_stdev <- stdev_reg2_4_women
      }
      else if (y==5){
        index_mean <- mean_reg2_5_women
        index_stdev <- stdev_reg2_5_women
      }
    } else if (x == 3){
      if (y == 1){
        index_mean <- mean_reg3_1_women
        index_stdev <- stdev_reg3_1_women
      }
      else if (y==2){
        index_mean <- mean_reg3_2_women
        index_stdev <- stdev_reg3_2_women
      }
      else if (y==3){
        index_mean <- mean_reg3_3_women
        index_stdev <- stdev_reg3_3_women
      }
      else if (y==4){
        index_mean <- mean_reg3_4_women
        index_stdev <- stdev_reg3_4_women
      }
      else if (y==5){
        index_mean <- mean_reg3_5_women
        index_stdev <- stdev_reg3_5_women
      }
    } else {
      if (y == 1){
        index_mean <- mean_reg4_1_women
        index_stdev <- stdev_reg4_1_women
      }
      else if (y==2){
        index_mean <- mean_reg4_2_women
        index_stdev <- stdev_reg4_2_women
      }
      else if (y==3){
        index_mean <- mean_reg4_3_women
        index_stdev <- stdev_reg4_3_women
      }
      else if (y==4){
        index_mean <- mean_reg4_4_women
        index_stdev <- stdev_reg4_4_women
      }
      else if (y==5){
        index_mean <- mean_reg4_5_women
        index_stdev <- stdev_reg4_5_women
      }
    }
  }
  data_point <- rnorm(1, index_mean, index_stdev)
  return (data_point)
}


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
                                                                             'camMets_ind_fact',  'pa_workini_fact', 'cam_index_mean')]
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

plot(x = final_output_women$cam_index_mean, y= final_output_women$PAEE)
abline(0,1)
plot(x = final_output_women$reg_value, y= final_output_women$PAEE)
abline(0,1)
plot(x = final_output_women$binary_index_mean, y= final_output_women$PAEE)
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

og_men_lengths = lapply(X = split(x = og_men$cam_index, f = og_men$cam_index_fact), FUN = length)

index_list_men = split(x=final_output_men$PAEE, f= as.factor(final_output_men$cam_index))

comb_list_men = mapply(FUN = list, og_men_lengths, index_list_men, SIMPLIFY=FALSE)


og_women_lengths = lapply(X = split(x = og_women$cam_index, f = og_women$cam_index_fact), FUN = length)

index_list_women = split(x=final_output_women$PAEE, f= as.factor(final_output_women$cam_index))

comb_list_women = mapply(FUN = list, og_women_lengths, index_list_women, SIMPLIFY=FALSE)

#number of draws
numtrials <- 200

results_df <- data.frame()

betas_per_men <- vector("numeric")
std_errs_per_men <- vector("numeric")

betas_ind_men <- vector("numeric")
std_errs_ind_men <- vector("numeric")


for (j in 1:numtrials){

  ptm <- proc.time()
  og_men$paee_sample_per =  unlist(unname(lapply(X = comb_list_men,
                                  FUN = function(input){
                                    
                                    output = sample(x = input[[2]], size = input[[1]], replace = TRUE)
                                    return(output)
                                  }
                       
                       )))
  time1 = proc.time() - ptm
  
  ptm <- proc.time()
  og_men$paee_sample_ind =  unlist(unname(lapply(X = comb_list_men,
                                FUN = function(input){
                                  
                                  output = rep(sample(x = input[[2]], size = 1, replace = TRUE), times = input[[1]])
                                  return(output)
                                }
                                
  )))
  time2 = proc.time() - ptm
  ptm <- proc.time()
  # calculate and store the coefficients per person
  reg_out_per <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ paee_sample_per + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
  reg_coeff_per <- reg_out_per$coefficients["paee_sample_per"]
  #reg_std_per <- (summary(reg_out_per)$coefficients[,"Std. Error"])["paee_sample_per"]
  betas_per_men <- c(betas_per_men, reg_coeff_per)
  #std_errs_per <- c(std_errs_per, reg_std_per)
  
  time3 = proc.time() - ptm
  
  ptm <- proc.time()
  # calculate and store the coefficients per index
  reg_out_ind <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ paee_sample_ind + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
  reg_coeff_ind <- reg_out_ind$coefficients["paee_sample_ind"]
  #reg_std_ind <- (summary(reg_out_ind)$coefficients[,"Std. Error"])["paee_sample_ind"]
  betas_ind_men <- c(betas_ind_men, reg_coeff_ind)
  #std_errs_ind <- c(std_errs_ind, reg_std_ind)
  time4 = proc.time() - ptm
}


betas_per_women <- vector("numeric")
std_errs_per_women <- vector("numeric")

betas_ind_women <- vector("numeric")
std_errs_ind_women <- vector("numeric")


for (j in 1:numtrials){
  
  ptm <- proc.time()
  og_women$paee_sample_per =  unlist(unname(lapply(X = comb_list_women,
                                                 FUN = function(input){
                                                   
                                                   output = sample(x = input[[2]], size = input[[1]], replace = TRUE)
                                                   return(output)
                                                 }
                                                 
  )))
  time1 = proc.time() - ptm
  
  ptm <- proc.time()
  og_women$paee_sample_ind =  unlist(unname(lapply(X = comb_list_women,
                                                 FUN = function(input){
                                                   
                                                   output = rep(sample(x = input[[2]], size = 1, replace = TRUE), times = input[[1]])
                                                   return(output)
                                                 }
                                                 
  )))
  time2 = proc.time() - ptm
  ptm <- proc.time()
  # calculate and store the coefficients per person
  reg_out_per <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ paee_sample_per + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)
  reg_coeff_per <- reg_out_per$coefficients["paee_sample_per"]
  #reg_std_per <- (summary(reg_out_per)$coefficients[,"Std. Error"])["paee_sample_per"]
  betas_per_women <- c(betas_per_women, reg_coeff_per)
  #std_errs_per <- c(std_errs_per, reg_std_per)
  
  time3 = proc.time() - ptm
  
  ptm <- proc.time()
  # calculate and store the coefficients per index
  reg_out_ind <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ paee_sample_ind + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)
  reg_coeff_ind <- reg_out_ind$coefficients["paee_sample_ind"]
  #reg_std_ind <- (summary(reg_out_ind)$coefficients[,"Std. Error"])["paee_sample_ind"]
  betas_ind_women <- c(betas_ind_women, reg_coeff_ind)
  #std_errs_ind <- c(std_errs_ind, reg_std_ind)
  time4 = proc.time() - ptm
}

library(ggplot2)
#for visualisation, plot the last per person set of values
print(ggplot(og_men,aes(x=paee_sample_per,group=cam_index,fill=factor(cam_index)))+
        geom_histogram(position="identity",alpha=0.5,binwidth=1)+theme_bw()+
        ggtitle(label = "PAEE person men")+
        scale_fill_manual(values = c("red", "grey", "seagreen3","blue")))

#for visualisation, plot the last per index set of values
print(ggplot(og_men,aes(x=paee_sample_ind,group=cam_index,fill=factor(cam_index)))+
        geom_histogram(position="identity",alpha=0.5,binwidth=1)+theme_bw()+
        ggtitle(label = "PAEE index men")+
        scale_fill_manual(values = c("red", "grey", "seagreen3","blue")))

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
