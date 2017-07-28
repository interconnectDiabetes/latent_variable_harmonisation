# Succeeded by paper_flow.R
# Initial outline for paper where regression based harmonisation methods are applied to varying categorical independent variables.


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

# to do: Add in UMEA, if possible (not sure if all indices can be generated)
# 
# epic_umea <- read.csv("V:/Studies/InterConnect/Internal/Latent variable harmonisation/epic_umea.csv", header=T)
# the csv was generated from:
# V:\Studies\InterAct_WP_2.6\StudyData\HANDOVER TO DM\links\umea_idmerge.dta
# V:\Studies\InterAct_WP_2.6\UMEÅ_EPIC_PAQ\UMEA_EPIC_final.xls


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

## In this section we regenerate cam_index column ie check logic

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
                                    if (x[1] == 9){
                                      x_ind = 5
                                    } else {
                                      x_ind = x[1]
                                    }
                                    output <- cam_matrix[x_ind, x[2]]
                                    return(output)
                                  })

#if test has zero length then we are happy with the definition of cam_index!

test <- merged_output[merged_output$cam_index!=merged_output$cam_index2,]

#record results

cam_index_means = aggregate(merged_output$PAEE, by=list(cam_index=merged_output$cam_index, sex=merged_output$sex), mean, na.rm = TRUE)
cam_index_medians = aggregate(merged_output$PAEE, by=list(cam_index=merged_output$cam_index, sex=merged_output$sex), median, na.rm = TRUE)
cam_index_counts =aggregate(merged_output$PAEE, by=list(cam_index=merged_output$cam_index, sex=merged_output$sex), length)

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
binary_index_medians = aggregate(merged_output$PAEE, by=list(binary_index=merged_output$binary_index, sex=merged_output$sex), median, na.rm = TRUE)
binary_index_counts =aggregate(merged_output$PAEE, by=list(binary_index=merged_output$binary_index, sex=merged_output$sex), length)


# Method C - regression based (but still essentially 16 categories) on components of cam index

merged_output$pa_work_ind <- as.factor(as.vector(do.call(rbind,lapply(X = merged_output$pa_workini, FUN = function(x){
  
  if (x == 5 || x == 9){
    ind = 1
  }
  else { 
    ind = x
  }
  return(ind)
})
)))

#create factors for regression
merged_output$camMets_ind_fact <- as.factor(merged_output$camMets_ind)
merged_output$sex_fact <- as.factor(merged_output$sex)

reg_out <- glm(formula=PAEE ~ pa_work_ind + camMets_ind_fact + sex_fact, data = merged_output, family='gaussian')




# calculate the lambdas (to be used much later!)

# Setting the 'observed' PAEE values by the mean

# cam index

temp <- merge(merged_output, cam_index_means, by=c('cam_index', 'sex'))[,c('universal_id', 'sex', 'sex_fact', 
        'country.x', 'height', 'weight','age', 'PAEE', 'center', 'camMets_ind',  'cam_index',
        'binary_index',	'camMets_ind_fact',	'pa_work_ind', 'pa_total', 'bmi')]



new_data <- temp[order(temp$universal_id),]

# binary index

temp <- merge(new_data, binary_index_means, by=c('binary_index', 'sex'))

new_data <- temp[order(temp$universal_id),]

# regression based

# convert work index into binaries to make use of regression coeffs

temp <- data.frame(sapply(levels(new_data$pa_work_ind), function(x) as.integer(x == new_data$pa_work_ind)))
colnames(temp) <- c('pa_work_ind_1','pa_work_ind_2','pa_work_ind_3','pa_work_ind_4')
new_data <- cbind(new_data, temp)

# convert ltpa index into binaries

temp <- data.frame(sapply(levels(new_data$camMets_ind_fact), function(x) as.integer(x == new_data$camMets_ind_fact)))
colnames(temp) <- c('camMets_ind_1','camMets_ind_2','camMets_ind_3','camMets_ind_4')
new_data <- cbind(new_data, temp)

new_data$reg_value <- 

