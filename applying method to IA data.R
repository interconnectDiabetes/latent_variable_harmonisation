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

validation_data$index16 = with(validation_data, interaction(as.factor(validation_data$camMets_ind),  as.factor(validation_data$pa_workini)), drop=TRUE)

#fix levels in validation data to use same names as main study
levels(validation_data$country) = c("DENMARK", "FRANCE", "GERMANY", "GREECE", "ITALY", "NETHERLANDS", 
                            "NORWAY", "SPAIN", "UK", "SWEDEN")

# add in BMI
validation_data$bmi_adj = validation_data$weight/validation_data$height^2


###############################################################################
###############################################################################

#record means etc for the cam index and binary index
# Note that these match the Peters et al paper! Hurrah!
cam_index_means = aggregate(validation_data$PAEE, by=list(cam_index=validation_data$cam_index, sex=validation_data$sex), mean, na.rm = TRUE)
colnames(cam_index_means)[3] = 'cam_index_mean'
cam_index_medians = aggregate(validation_data$PAEE, by=list(cam_index=validation_data$cam_index, sex=validation_data$sex), median, na.rm = TRUE)
colnames(cam_index_medians)[3] = 'cam_index_median'
cam_index_counts =aggregate(validation_data$PAEE, by=list(cam_index=validation_data$cam_index, sex=validation_data$sex), length)
colnames(cam_index_counts)[3] = 'cam_index_count'



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
study_data$cam_index <- as.factor(study_data$cam_index)
study_data$cam_index_cont <- as.numeric(study_data$cam_index)

# create binary index
study_data$binary_index <- unlist(lapply(study_data$cam_index, FUN=function(x){
  if (x[1] == 1) {output = 1 }
  else {output = 2}
  return(output)
}))
study_data$binary_index <- as.factor(study_data$binary_index)

#16 level index

study_data$index16 = with(study_data, interaction(as.factor(study_data$camMets_ind),  as.factor(study_data$pa_workini)), drop=TRUE)


og_men <- subset(study_data, sex==1)
og_women <- subset(study_data, sex==0)


# fix empty france level for men
og_men$country <- droplevels(og_men$country)

###############################################################################
###############################################################################
###############################################################################
#test the process
#  matches pa_methods in Stata for 'cont' version
test_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
test_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)

test_reg_cont_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_cont + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
test_reg_cont_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_cont + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)

################################################
### Now try applying our method to the IA data
###############################################

IA_df_men = data.frame()
IA_df_women = data.frame()

for (i in 1:20){

validation_data$cam_index = as.factor(validation_data$cam_index)

validation_data = validation_data[sample(nrow(validation_data),nrow(validation_data)),]

validation_data = validation_data[order(validation_data$cam_index),]

validation_data_men = validation_data[which(validation_data$sex == 1),]
validation_data_women = validation_data[which(validation_data$sex == 0),]

validation_data_men$split = unlist(lapply(X = split(x = validation_data_men, validation_data_men$cam_index),FUN = function(X){
  y1 = rep(x=1,times=ceiling(length(X$cam_index)/2))
  y2 = rep(x=2,times=length(X$cam_index)-ceiling(length(X$cam_index)/2))
  y = c(y1,y2)
  return(y)
}))

validation_data_women$split = unlist(lapply(X = split(x = validation_data_women, validation_data_women$cam_index),FUN = function(X){
  y1 = rep(x=1,times=ceiling(length(X$cam_index)/2))
  y2 = rep(x=2,times=length(X$cam_index)-ceiling(length(X$cam_index)/2))
  y = c(y1,y2)
  return(y)
}))

validation_data_men_A = validation_data_men[validation_data_men$split == 1,]
validation_data_men_B = validation_data_men[validation_data_men$split == 2,]

validation_data_women_A = validation_data_women[validation_data_women$split == 1,]
validation_data_women_B = validation_data_women[validation_data_women$split == 2,]

# Do the mapping using data from validation sample A, and apply to A, B and study

mapping_formula_men = as.formula(PAEE ~ cam_index + 0) # 0 forces no intercept
mapping_model_men <- lm(formula=mapping_formula_men, data=validation_data_men_A)
validation_data_men_A$index_mapped = predict(mapping_model_men, newdata = validation_data_men_A)
validation_data_men_B$index_mapped = predict(mapping_model_men, newdata = validation_data_men_B)
og_men$index_mapped = predict(mapping_model_men, newdata = og_men)

mapping_formula_women = as.formula(PAEE ~ cam_index + 0) # 0 forces no intercept
mapping_model_women <- lm(formula=mapping_formula_women, data=validation_data_women_A)
validation_data_women_A$index_mapped = predict(mapping_model_women, newdata = validation_data_women_A)
validation_data_women_B$index_mapped = predict(mapping_model_women, newdata = validation_data_women_B)
og_women$index_mapped = predict(mapping_model_women, newdata = og_women)

#get mapped exposure to outcome relationship (naive model)

print('got here')

mapped_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
estimate_mapped_men = mapped_model_men$coefficients["index_mapped"]
std_error_mapped_men = summary(mapped_model_men)$coefficients["index_mapped","robust se"]

mapped_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)
estimate_mapped_women = mapped_model_women$coefficients["index_mapped"]
std_error_mapped_women = summary(mapped_model_women)$coefficients["index_mapped","robust se"]

#error in mapping (no calibration, "lambda" is 1) - use validation A

mapping_error_formula_men = as.formula(PAEE ~ index_mapped)
mapping_error_model_men <- lm(formula=mapping_error_formula_men, data=validation_data_men_A)
estimate_mapping_men = mapping_error_model_men$coefficients["index_mapped"] # should be 1
std_error_mapping_men = summary(mapping_error_model_men)$coefficients["index_mapped","Std. Error"]

mapping_error_formula_women = as.formula(PAEE ~ index_mapped)
mapping_error_model_women <- lm(formula=mapping_error_formula_women, data=validation_data_women_A)
estimate_mapping_women = mapping_error_model_women$coefficients["index_mapped"] # should be 1
std_error_mapping_women = summary(mapping_error_model_women)$coefficients["index_mapped","Std. Error"]


#calibration & associated error - use validation B and apply calibration to study

calib_formula_men = as.formula(PAEE ~ index_mapped)
calib_model_men <- lm(formula=calib_formula_men, data=validation_data_men_B)
og_men$exposure_corrected = predict(calib_model_men, newdata = og_men)
estimate_calib_men = calib_model_men$coefficients["index_mapped"]
std_error_calib_men = summary(calib_model_men)$coefficients["index_mapped","Std. Error"]

calib_formula_women = as.formula(PAEE ~ index_mapped)
calib_model_women <- lm(formula=calib_formula_women, data=validation_data_women_B)
og_women$exposure_corrected = predict(calib_model_women, newdata = og_women)
estimate_calib_women = calib_model_women$coefficients["index_mapped"]
std_error_calib_women = summary(calib_model_women)$coefficients["index_mapped","Std. Error"]

# Use mapped and calibrated values to generate a corrected model - estimate beta using the corrected exposure and outcome
# store the estimate and its standard error

corrected_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
estimate_corrected_men = corrected_model_men$coefficients["exposure_corrected"]
std_error_corrected_men = summary(corrected_model_men)$coefficients["exposure_corrected","robust se"]

corrected_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)
estimate_corrected_women = corrected_model_women$coefficients["exposure_corrected"]
std_error_corrected_women = summary(corrected_model_women)$coefficients["exposure_corrected","robust se"]

# refine the standard error of the main model estimate using delta function

variance_beta_men = std_error_mapped_men^2
beta_lambda_div_sq_men = (unname(unlist(estimate_mapped_men/(estimate_calib_men)^2)))^2
var_lambda_men = (std_error_calib_men)^2
delta_variance_men = (variance_beta_men / (estimate_calib_men)^2) + (beta_lambda_div_sq_men * var_lambda_men)
delta_std_error_corrected_men = sqrt(delta_variance_men)

variance_beta_women = std_error_mapped_women^2
beta_lambda_div_sq_women = (unname(unlist(estimate_mapped_women/(estimate_calib_women)^2)))^2
var_lambda_women = (std_error_calib_women)^2
delta_variance_women = (variance_beta_women / (estimate_calib_women)^2) + (beta_lambda_div_sq_women * var_lambda_women)
delta_std_error_corrected_women = sqrt(delta_variance_women)

# Store estimate and corrected std error

estimates_men = estimate_corrected_men
std_errors_men = delta_std_error_corrected_men

estimates_women = estimate_corrected_women
std_errors_women = delta_std_error_corrected_women

thing_to_store_men = data.frame(naive_estimate_men = estimate_mapped_men, naive_se_men = std_error_mapped_men,
                                corrected_estimate_men = estimates_men, corrected_se_men = delta_std_error_corrected_men,
                                lambda_men = estimate_calib_men, lambda_se_men = std_error_calib_men)

thing_to_store_women = data.frame(naive_estimate_women = estimate_mapped_women, naive_se_women = std_error_mapped_women,
                                corrected_estimate_women = estimates_women, corrected_se_women = delta_std_error_corrected_women,
                                lambda_women = estimate_calib_women, lambda_se_women = std_error_calib_women)

IA_df_men = rbind(IA_df_men, thing_to_store_men)
IA_df_women = rbind(IA_df_women, thing_to_store_women)

}

#############################################################
############ 16 level VERSION
############################################################


validation_data$index16 = as.factor(validation_data$index16)

validation_data = validation_data[sample(nrow(validation_data),nrow(validation_data)),]

validation_data = validation_data[order(validation_data$index16),]

validation_data_men = validation_data[which(validation_data$sex == 1),]
validation_data_women = validation_data[which(validation_data$sex == 0),]

validation_data_men$split = unlist(lapply(X = split(x = validation_data_men, validation_data_men$index16),FUN = function(X){
  y1 = rep(x=1,times=ceiling(length(X$index16)/2))
  y2 = rep(x=2,times=length(X$index16)-ceiling(length(X$index16)/2))
  y = c(y1,y2)
  return(y)
}))

validation_data_women$split = unlist(lapply(X = split(x = validation_data_women, validation_data_women$index16),FUN = function(X){
  y1 = rep(x=1,times=ceiling(length(X$index16)/2))
  y2 = rep(x=2,times=length(X$index16)-ceiling(length(X$index16)/2))
  y = c(y1,y2)
  return(y)
}))

validation_data_men_A = validation_data_men[validation_data_men$split == 1,]
validation_data_men_B = validation_data_men[validation_data_men$split == 2,]

validation_data_women_A = validation_data_women[validation_data_women$split == 1,]
validation_data_women_B = validation_data_women[validation_data_women$split == 2,]

# Do the mapping using data from validation sample A, and apply to A, B and study

mapping_formula_men = as.formula(PAEE ~ index16 + 0) # 0 forces no intercept
mapping_model_men <- lm(formula=mapping_formula_men, data=validation_data_men_A)
validation_data_men_A$index_mapped = predict(mapping_model_men, newdata = validation_data_men_A)
validation_data_men_B$index_mapped = predict(mapping_model_men, newdata = validation_data_men_B)
og_men$index_mapped = predict(mapping_model_men, newdata = og_men)

mapping_formula_women = as.formula(PAEE ~ index16 + 0) # 0 forces no intercept
mapping_model_women <- lm(formula=mapping_formula_women, data=validation_data_women_A)
validation_data_women_A$index_mapped = predict(mapping_model_women, newdata = validation_data_women_A)
validation_data_women_B$index_mapped = predict(mapping_model_women, newdata = validation_data_women_B)
og_women$index_mapped = predict(mapping_model_women, newdata = og_women)

#get mapped exposure to outcome relationship (naive model)

mapped_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
estimate_mapped_men = mapped_model_men$coefficients["index_mapped"]
std_error_mapped_men = summary(mapped_model_men)$coefficients["index_mapped","robust se"]

mapped_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)
estimate_mapped_women = mapped_model_women$coefficients["index_mapped"]
std_error_mapped_women = summary(mapped_model_women)$coefficients["index_mapped","robust se"]

#error in mapping (no calibration, "lambda" is 1) - use validation A

mapping_error_formula_men = as.formula(PAEE ~ index_mapped)
mapping_error_model_men <- lm(formula=mapping_error_formula_men, data=validation_data_men_A)
estimate_mapping_men = mapping_error_model_men$coefficients["index_mapped"] # should be 1
std_error_mapping_men = summary(mapping_error_model_men)$coefficients["index_mapped","Std. Error"]

mapping_error_formula_women = as.formula(PAEE ~ index_mapped)
mapping_error_model_women <- lm(formula=mapping_error_formula_women, data=validation_data_women_A)
estimate_mapping_women = mapping_error_model_women$coefficients["index_mapped"] # should be 1
std_error_mapping_women = summary(mapping_error_model_women)$coefficients["index_mapped","Std. Error"]


#calibration & associated error - use validation B and apply calibration to study

calib_formula_men = as.formula(PAEE ~ index_mapped)
calib_model_men <- lm(formula=calib_formula_men, data=validation_data_men_B)
og_men$exposure_corrected = predict(calib_model_men, newdata = og_men)
estimate_calib_men = calib_model_men$coefficients["index_mapped"]
std_error_calib_men = summary(calib_model_men)$coefficients["index_mapped","Std. Error"]

calib_formula_women = as.formula(PAEE ~ index_mapped)
calib_model_women <- lm(formula=calib_formula_women, data=validation_data_women_B)
og_women$exposure_corrected = predict(calib_model_women, newdata = og_women)
estimate_calib_women = calib_model_women$coefficients["index_mapped"]
std_error_calib_women = summary(calib_model_women)$coefficients["index_mapped","Std. Error"]

# Use mapped and calibrated values to generate a corrected model - estimate beta using the corrected exposure and outcome
# store the estimate and its standard error

corrected_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_men, robust=TRUE)
estimate_corrected_men = corrected_model_men$coefficients["exposure_corrected"]
std_error_corrected_men = summary(corrected_model_men)$coefficients["exposure_corrected","robust se"]

corrected_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + smoke_stat + l_school + alc_re + qe_energy + country, data = og_women, robust=TRUE)
estimate_corrected_women = corrected_model_women$coefficients["exposure_corrected"]
std_error_corrected_women = summary(corrected_model_women)$coefficients["exposure_corrected","robust se"]

# refine the standard error of the main model estimate using delta function

variance_beta_men = std_error_mapped_men^2
beta_lambda_div_sq_men = (unname(unlist(estimate_mapped_men/(estimate_calib_men)^2)))^2
var_lambda_men = (std_error_calib_men)^2
delta_variance_men = (variance_beta_men / (estimate_calib_men)^2) + (beta_lambda_div_sq_men * var_lambda_men)
delta_std_error_corrected_men = sqrt(delta_variance_men)

variance_beta_women = std_error_mapped_women^2
beta_lambda_div_sq_women = (unname(unlist(estimate_mapped_women/(estimate_calib_women)^2)))^2
var_lambda_women = (std_error_calib_women)^2
delta_variance_women = (variance_beta_women / (estimate_calib_women)^2) + (beta_lambda_div_sq_women * var_lambda_women)
delta_std_error_corrected_women = sqrt(delta_variance_women)

# Store estimate and corrected std error

estimates_men = estimate_corrected_men
std_errors_men = delta_std_error_corrected_men

estimates_women = estimate_corrected_women
std_errors_women = delta_std_error_corrected_women


################################################
### Simple, matching model
###############################################


test_reg_men_simple <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + country, data = og_men, robust=TRUE)
# cam_index3 0.8645 se 0.062316 (on log scale)
# cam_index3 0.7497   0.068050 
# cam_index4 0.6977  0.069207
# all differ by about 0.87 => linear assumption ok


test_reg_women_simple <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index + bmi_adj + country, data = og_women, robust=TRUE)
# cam_index2            0.8344     0.055931
# cam_index3            0.8529     0.062997
# cam_index4            0.7278     0.068548
# linear assumption more dubious?


test_reg_cont_men_simple <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_cont + bmi_adj + country, data = og_men, robust=TRUE)
# -0.123712 0.021898

test_reg_cont_women_simple <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_cont + bmi_adj + country, data = og_women, robust=TRUE)
# -0.094726   0.021738

lin_test_simple_men = lm(PAEE~as.numeric(cam_index), data=validation_data_men)
# men relationship 6.5
lin_test_simple_women = lm(PAEE~as.numeric(cam_index), data=validation_data_women)
# women relationship 3.7

# men baseline
# -0.1237/6.5 = -0.0190 kj/kg/day SE: 0.0219/6.5 = 0.00337

# women baseline
# -0.0947/3.7 = -0.0256 kj/kg/day SE: 0.0217/3.7 = 0.00588




IA_df_men = data.frame()
IA_df_women = data.frame()


for (i in 1:1001){

  set.seed(i)
  
validation_data$cam_index = as.factor(validation_data$cam_index)

# validation_data = validation_data[sample(nrow(validation_data),nrow(validation_data)),]
# 
# validation_data = validation_data[order(validation_data$cam_index),]

validation_data_men = validation_data[which(validation_data$sex == 1),]
validation_data_women = validation_data[which(validation_data$sex == 0),]

# validation_data_men$split = unlist(lapply(X = split(x = validation_data_men, validation_data_men$cam_index),FUN = function(X){
#   y1 = rep(x=1,times=ceiling(length(X$cam_index)/2))
#   y2 = rep(x=2,times=length(X$cam_index)-ceiling(length(X$cam_index)/2))
#   y = c(y1,y2)
#   return(y)
# }))
# 
# validation_data_women$split = unlist(lapply(X = split(x = validation_data_women, validation_data_women$cam_index),FUN = function(X){
#   y1 = rep(x=1,times=ceiling(length(X$cam_index)/2))
#   y2 = rep(x=2,times=length(X$cam_index)-ceiling(length(X$cam_index)/2))
#   y = c(y1,y2)
#   return(y)
# }))

# split the validation study in half randomly
set.seed(i)
validation_data_men$split = unlist(lapply(X = split(x = validation_data_men, validation_data_men$index),FUN = function(X){
  y1 = rep(x=1,times=ceiling(length(X$index)/2))
  y2 = rep(x=2,times=length(X$index)-ceiling(length(X$index)/2))
  y = c(y1,y2)
  y = y[sample(1:length(y))]
  return(y)
}))
set.seed(i)
validation_data_women$split = unlist(lapply(X = split(x = validation_data_women, validation_data_women$index),FUN = function(X){
  y1 = rep(x=1,times=ceiling(length(X$index)/2))
  y2 = rep(x=2,times=length(X$index)-ceiling(length(X$index)/2))
  y = c(y1,y2)
  y = y[sample(1:length(y))]
  return(y)
}))


validation_data_men_A = validation_data_men[validation_data_men$split == 1,]
validation_data_men_B = validation_data_men[validation_data_men$split == 2,]

validation_data_women_A = validation_data_women[validation_data_women$split == 1,]
validation_data_women_B = validation_data_women[validation_data_women$split == 2,]

# Do the mapping using data from validation sample A, and apply to A, B and study

mapping_formula_men = as.formula(PAEE ~ cam_index + bmi_adj + country + 0) # 0 forces no intercept
mapping_model_men <- lm(formula=mapping_formula_men, data=validation_data_men_A)
validation_data_men_A$index_mapped = predict(mapping_model_men, newdata = validation_data_men_A)
validation_data_men_B$index_mapped = predict(mapping_model_men, newdata = validation_data_men_B)
og_men$index_mapped = predict(mapping_model_men, newdata = og_men)

mapping_formula_women = as.formula(PAEE ~ cam_index + bmi_adj + country + 0) # 0 forces no intercept
mapping_model_women <- lm(formula=mapping_formula_women, data=validation_data_women_A)
validation_data_women_A$index_mapped = predict(mapping_model_women, newdata = validation_data_women_A)
validation_data_women_B$index_mapped = predict(mapping_model_women, newdata = validation_data_women_B)
og_women$index_mapped = predict(mapping_model_women, newdata = og_women)

#get mapped exposure to outcome relationship (naive model)

print('gothere')

mapped_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + country + 0, data = og_men, robust=TRUE)
estimate_mapped_men = mapped_model_men$coefficients["index_mapped"]
std_error_mapped_men = summary(mapped_model_men)$coefficients["index_mapped","robust se"]

mapped_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + country + 0, data = og_women, robust=TRUE)
estimate_mapped_women = mapped_model_women$coefficients["index_mapped"]
std_error_mapped_women = summary(mapped_model_women)$coefficients["index_mapped","robust se"]

#error in mapping (no calibration, "lambda" is 1) - use validation A

mapping_error_formula_men = as.formula(PAEE ~ index_mapped + bmi_adj + country + 0)
mapping_error_model_men <- lm(formula=mapping_error_formula_men, data=validation_data_men_A)
estimate_mapping_men = mapping_error_model_men$coefficients["index_mapped"] # should be 1
std_error_mapping_men = summary(mapping_error_model_men)$coefficients["index_mapped","Std. Error"]

mapping_error_formula_women = as.formula(PAEE ~ index_mapped + bmi_adj + country + 0)
mapping_error_model_women <- lm(formula=mapping_error_formula_women, data=validation_data_women_A)
estimate_mapping_women = mapping_error_model_women$coefficients["index_mapped"] # should be 1
std_error_mapping_women = summary(mapping_error_model_women)$coefficients["index_mapped","Std. Error"]


#calibration & associated error - use validation B and apply calibration to study

calib_formula_men = as.formula(PAEE ~ index_mapped + bmi_adj + country)
calib_model_men <- lm(formula=calib_formula_men, data=validation_data_men_B)
og_men$exposure_corrected = predict(calib_model_men, newdata = og_men)
estimate_calib_men = calib_model_men$coefficients["index_mapped"]
std_error_calib_men = summary(calib_model_men)$coefficients["index_mapped","Std. Error"]

calib_formula_women = as.formula(PAEE ~ index_mapped + bmi_adj + country)
calib_model_women <- lm(formula=calib_formula_women, data=validation_data_women_B)
og_women$exposure_corrected = predict(calib_model_women, newdata = og_women)
estimate_calib_women = calib_model_women$coefficients["index_mapped"]
std_error_calib_women = summary(calib_model_women)$coefficients["index_mapped","Std. Error"]

# Use mapped and calibrated values to generate a corrected model - estimate beta using the corrected exposure and outcome
# store the estimate and its standard error

corrected_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + country, data = og_men, robust=TRUE)
estimate_corrected_men = corrected_model_men$coefficients["exposure_corrected"]
std_error_corrected_men = summary(corrected_model_men)$coefficients["exposure_corrected","robust se"]

corrected_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + country, data = og_women, robust=TRUE)
estimate_corrected_women = corrected_model_women$coefficients["exposure_corrected"]
std_error_corrected_women = summary(corrected_model_women)$coefficients["exposure_corrected","robust se"]

# refine the standard error of the main model estimate using delta function

variance_beta_men = std_error_mapped_men^2
beta_lambda_div_sq_men = (unname(unlist(estimate_mapped_men/(estimate_calib_men)^2)))^2
var_lambda_men = (std_error_calib_men)^2
delta_variance_men = (variance_beta_men / (estimate_calib_men)^2) + (beta_lambda_div_sq_men * var_lambda_men)
delta_std_error_corrected_men = sqrt(delta_variance_men)

variance_beta_women = std_error_mapped_women^2
beta_lambda_div_sq_women = (unname(unlist(estimate_mapped_women/(estimate_calib_women)^2)))^2
var_lambda_women = (std_error_calib_women)^2
delta_variance_women = (variance_beta_women / (estimate_calib_women)^2) + (beta_lambda_div_sq_women * var_lambda_women)
delta_std_error_corrected_women = sqrt(delta_variance_women)

# Store estimate and corrected std error

estimates_men = estimate_corrected_men
std_errors_men = delta_std_error_corrected_men

estimates_women = estimate_corrected_women
std_errors_women = delta_std_error_corrected_women

thing_to_store_men = data.frame(naive_estimate_men = estimate_mapped_men, naive_se_men = std_error_mapped_men,
                                corrected_estimate_men = estimates_men, corrected_se_men = delta_std_error_corrected_men,
                                lambda_men = estimate_calib_men, lambda_se_men = std_error_calib_men)

thing_to_store_women = data.frame(naive_estimate_women = estimate_mapped_women, naive_se_women = std_error_mapped_women,
                                  corrected_estimate_women = estimates_women, corrected_se_women = delta_std_error_corrected_women,
                                  lambda_women = estimate_calib_women, lambda_se_women = std_error_calib_women)

IA_df_men = rbind(IA_df_men, thing_to_store_men)
print(i)
IA_df_women = rbind(IA_df_women, thing_to_store_women)
print(i)

}


IA_df_men[which((IA_df_men$corrected_estimate_men) == median(IA_df_men$corrected_estimate_men)),]
IA_df_women[which((IA_df_women$corrected_estimate_women) == median(IA_df_women$corrected_estimate_women)),]

# In previous loop, do it 100 times to account for split variation
# Then take median
#naive_estimate_men naive_se_men corrected_estimate_men corrected_se_men lambda_men
#index_mapped34        -0.01845768  0.003349215            -0.01954376      0.004985866  0.9444284
#lambda_se_men
#index_mapped34     0.1693589

#naive_estimate_women naive_se_women corrected_estimate_women corrected_se_women
#index_mapped29          -0.02287536    0.005595674               -0.0258762        0.007435163
#lambda_women lambda_se_women
#index_mapped29    0.8840309       0.1332662



##
# Investigating effects of sampling

set.seed(123456)

menA = data.frame()
menB= data.frame()
womenA = data.frame()
womenB = data.frame()

for (i in 1:10){
  
  validation_data$cam_index = as.factor(validation_data$cam_index)
  
  validation_data = validation_data[sample(nrow(validation_data),nrow(validation_data)),]
  
  validation_data = validation_data[order(validation_data$cam_index),]
  
  validation_data_men = validation_data[which(validation_data$sex == 1),]
  validation_data_women = validation_data[which(validation_data$sex == 0),]
  
  validation_data_men$split = unlist(lapply(X = split(x = validation_data_men, validation_data_men$cam_index),FUN = function(X){
    y1 = rep(x=1,times=ceiling(length(X$cam_index)/2))
    y2 = rep(x=2,times=length(X$cam_index)-ceiling(length(X$cam_index)/2))
    y = c(y1,y2)
    return(y)
  }))
  
  validation_data_women$split = unlist(lapply(X = split(x = validation_data_women, validation_data_women$cam_index),FUN = function(X){
    y1 = rep(x=1,times=ceiling(length(X$cam_index)/2))
    y2 = rep(x=2,times=length(X$cam_index)-ceiling(length(X$cam_index)/2))
    y = c(y1,y2)
    return(y)
  }))
  
  validation_data_men_A = validation_data_men[validation_data_men$split == 1,]
  validation_data_men_B = validation_data_men[validation_data_men$split == 2,]
  
  validation_data_women_A = validation_data_women[validation_data_women$split == 1,]
  validation_data_women_B = validation_data_women[validation_data_women$split == 2,]
  
  temp_men_A = aggregate(validation_data_men_A$PAEE, by=list(validation_data_men_A$cam_index), mean, na.rm = TRUE)
  temp_men_B = aggregate(validation_data_men_B$PAEE, by=list(validation_data_men_B$cam_index), mean, na.rm = TRUE)
  temp_women_A = aggregate(validation_data_women_A$PAEE, by=list(validation_data_women_A$cam_index), mean, na.rm = TRUE)
  temp_women_B = aggregate(validation_data_women_B$PAEE, by=list(validation_data_women_B$cam_index), mean, na.rm = TRUE)
  
  menA = rbind(menA, temp_men_A[,2])
  menB = rbind(menB, temp_men_B[,2])
  womenA = rbind(womenA, temp_women_A[,2])
  womenB = rbind(womenB, temp_women_B[,2])
  
}

##Harmonisation and uncertatinty only

validation_data$cam_index = as.factor(validation_data$cam_index)

validation_data_men = validation_data[which(validation_data$sex == 1),]
validation_data_women = validation_data[which(validation_data$sex == 0),]

# Do the mapping using data from validation sample, and apply to study

mapping_formula_men = as.formula(PAEE ~ cam_index + bmi_adj + country + 0) # 0 forces no intercept
mapping_model_men <- lm(formula=mapping_formula_men, data=validation_data_men)
validation_data_men$index_mapped = predict(mapping_model_men, newdata = validation_data_men)
og_men$index_mapped = predict(mapping_model_men, newdata = og_men)

mapping_formula_women = as.formula(PAEE ~ cam_index + bmi_adj + country + 0) # 0 forces no intercept
mapping_model_women <- lm(formula=mapping_formula_women, data=validation_data_women)
validation_data_women$index_mapped = predict(mapping_model_women, newdata = validation_data_women)
og_women$index_mapped = predict(mapping_model_women, newdata = og_women)

#get mapped exposure to outcome relationship (naive model)

print('gothere')

mapped_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + country + 0, data = og_men, robust=TRUE)
estimate_mapped_men = mapped_model_men$coefficients["index_mapped"]
std_error_mapped_men = summary(mapped_model_men)$coefficients["index_mapped","robust se"]

mapped_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ index_mapped + bmi_adj + country + 0, data = og_women, robust=TRUE)
estimate_mapped_women = mapped_model_women$coefficients["index_mapped"]
std_error_mapped_women = summary(mapped_model_women)$coefficients["index_mapped","robust se"]

#error in mapping (no calibration, "lambda" is 1) - use validation 

mapping_error_formula_men = as.formula(PAEE ~ index_mapped + bmi_adj + country + 0)
mapping_error_model_men <- lm(formula=mapping_error_formula_men, data=validation_data_men)
estimate_mapping_men = mapping_error_model_men$coefficients["index_mapped"] # should be 1
std_error_mapping_men = summary(mapping_error_model_men)$coefficients["index_mapped","Std. Error"]

mapping_error_formula_women = as.formula(PAEE ~ index_mapped + bmi_adj + country + 0)
mapping_error_model_women <- lm(formula=mapping_error_formula_women, data=validation_data_women)
estimate_mapping_women = mapping_error_model_women$coefficients["index_mapped"] # should be 1
std_error_mapping_women = summary(mapping_error_model_women)$coefficients["index_mapped","Std. Error"]


#calibration & associated error - use validation B and apply calibration to study

calib_formula_men = as.formula(PAEE ~ index_mapped + bmi_adj + country)
calib_model_men <- lm(formula=calib_formula_men, data=validation_data_men)
og_men$exposure_corrected = predict(calib_model_men, newdata = og_men)
estimate_calib_men = calib_model_men$coefficients["index_mapped"]
std_error_calib_men = summary(calib_model_men)$coefficients["index_mapped","Std. Error"]

calib_formula_women = as.formula(PAEE ~ index_mapped + bmi_adj + country)
calib_model_women <- lm(formula=calib_formula_women, data=validation_data_women)
og_women$exposure_corrected = predict(calib_model_women, newdata = og_women)
estimate_calib_women = calib_model_women$coefficients["index_mapped"]
std_error_calib_women = summary(calib_model_women)$coefficients["index_mapped","Std. Error"]

# Use mapped and calibrated values to generate a corrected model - estimate beta using the corrected exposure and outcome
# store the estimate and its standard error

corrected_model_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + country, data = og_men, robust=TRUE)
estimate_corrected_men = corrected_model_men$coefficients["exposure_corrected"]
std_error_corrected_men = summary(corrected_model_men)$coefficients["exposure_corrected","robust se"]

corrected_model_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ exposure_corrected + bmi_adj + country, data = og_women, robust=TRUE)
estimate_corrected_women = corrected_model_women$coefficients["exposure_corrected"]
std_error_corrected_women = summary(corrected_model_women)$coefficients["exposure_corrected","robust se"]

# refine the standard error of the main model estimate using delta function

variance_beta_men = std_error_mapped_men^2
beta_lambda_div_sq_men = (unname(unlist(estimate_mapped_men/(estimate_calib_men)^2)))^2
var_lambda_men = (std_error_calib_men)^2
delta_variance_men = (variance_beta_men / (estimate_calib_men)^2) + (beta_lambda_div_sq_men * var_lambda_men)
delta_std_error_corrected_men = sqrt(delta_variance_men)

variance_beta_women = std_error_mapped_women^2
beta_lambda_div_sq_women = (unname(unlist(estimate_mapped_women/(estimate_calib_women)^2)))^2
var_lambda_women = (std_error_calib_women)^2
delta_variance_women = (variance_beta_women / (estimate_calib_women)^2) + (beta_lambda_div_sq_women * var_lambda_women)
delta_std_error_corrected_women = sqrt(delta_variance_women)

# Store estimate and corrected std error

estimates_men = estimate_corrected_men
std_errors_men = delta_std_error_corrected_men

estimates_women = estimate_corrected_women
std_errors_women = delta_std_error_corrected_women

thing_to_store_men = data.frame(naive_estimate_men = estimate_mapped_men, naive_se_men = std_error_mapped_men,
                                corrected_estimate_men = estimates_men, corrected_se_men = delta_std_error_corrected_men,
                                lambda_men = estimate_calib_men, lambda_se_men = std_error_calib_men)

thing_to_store_women = data.frame(naive_estimate_women = estimate_mapped_women, naive_se_women = std_error_mapped_women,
                                  corrected_estimate_women = estimates_women, corrected_se_women = delta_std_error_corrected_women,
                                  lambda_women = estimate_calib_women, lambda_se_women = std_error_calib_women)

## Forest plot

## Plot  
#set up plot

svg(filename = 'InterAct_results.svg', width = 8, height = 4.5)
par(mar=c(4.1,4.1,2.1,2.1))

labels = c('Baseline, per category (Men)',
           'Baseline, per category (Women)',
           'Baseline, per unit (Men)',
           'Baseline, per unit (Women)',
           'Harmonised, uncorrected, per unit (Men)',
           'Harmonised, uncorrected, per unit (Women)',
           'Harmonised, corrected, per unit (Men)',
           'Harmonised, corrected, per unit (Women)')

estimates = c(-0.124,-0.0947,-0.0190,-0.0256,-0.0190,-0.0253,-0.0190,-0.0253)
std_error = c(0.0219, 0.0217, 0.00337, 0.00588, 0.00345, 0.00623, 0.00403,0.00676)

forest(x = estimates, sei = std_error,
       xlab='Hazard Ratio',digits=3,
       refline = 1.0, slab = labels, transf=exp, alim=c(0.8,1.0), xlim = c(0.6,1.2),
       psize = c(1,1,1,1,1,1,1,1))


usr <- par("usr")
text(usr[2], usr[4], "Hazard Ratio [95% CI]", adj = c(1, 6),cex=1)
text(usr[1], usr[4], "Model", adj = c( 0, 6 ),cex=1)
dev.off()
