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

# Predictions per model

og_men <- subset(og_data, sex==0)
og_women <- subset(og_data, sex==1)

## LINEAR MODEL
cox_regression_men <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = og_men, robust=TRUE)
cox_regression_women <- coxph(Surv(age_recr_prentice,ageEnd,eventCens) ~ cam_index_means, data = og_women)

beta_men <- cox_regression_men$coefficients
beta_women <- cox_regression_women$coefficients


###############################################################################################
####################################  SAMPLING RUN THROUGH  ###################################
###############################################################################################



