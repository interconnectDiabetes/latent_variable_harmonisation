# Some code to run simulations using the real validation data

## Author: Tom Bishop
## Date: 23.01.2017

# Working Directory
setwd("V:/Studies/InterConnect/Internal/Latent variable harmonisation")


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




#####################some plots############
# These plots are to show the overlap in the distributions

cam_index_counts =aggregate(validation_data$PAEE, by=list(cam_index=validation_data$cam_index, sex=validation_data$sex), length)
colnames(cam_index_counts)[3] = 'cam_index_count'

validation_data_men <- na.omit(validation_data[validation_data$sex == 1,])
paee_1 <-validation_data_men$PAEE[validation_data_men$cam_index == 1]
paee_2 <-validation_data_men$PAEE[validation_data_men$cam_index == 2]
paee_3 <-validation_data_men$PAEE[validation_data_men$cam_index == 3]
paee_4 <-validation_data_men$PAEE[validation_data_men$cam_index == 4]

p1 <- hist(paee_1,breaks=c(seq(0,200, 5)))
p2 <- hist(paee_2,breaks=c(seq(0,200, 5)))
p3 <- hist(paee_3,breaks=c(seq(0,200, 5)))
p4 <- hist(paee_4,breaks=c(seq(0,200, 5)))

plot( p1, col=rgb(0,0,1,1/4), xlim = c(0,200), ylim = c(0,40), main = 'Men')
plot( p2, col=rgb(1,0,0,1/4), add=T)
plot( p3, col=rgb(1,1,0,1/4), add=T)
plot( p4, col=rgb(1,1,1,1/4), add=T)

validation_data_women <- na.omit(validation_data[validation_data$sex == 0,])
paee_1 <-validation_data_women$PAEE[validation_data_women$cam_index == 1]
paee_2 <-validation_data_women$PAEE[validation_data_women$cam_index == 2]
paee_3 <-validation_data_women$PAEE[validation_data_women$cam_index == 3]
paee_4 <-validation_data_women$PAEE[validation_data_women$cam_index == 4]

p1 <- hist(paee_1,breaks=c(seq(0,200, 5)))
p2 <- hist(paee_2,breaks=c(seq(0,200, 5)))
p3 <- hist(paee_3,breaks=c(seq(0,200, 5)))
p4 <- hist(paee_4,breaks=c(seq(0,200, 5)))

plot( p1, col=rgb(0,0,1,1/4), xlim = c(0,200), ylim = c(0,100), main = 'women')
plot( p2, col=rgb(1,0,0,1/4), add=T)
plot( p3, col=rgb(1,1,0,1/4), add=T)
plot( p4, col=rgb(1,1,1,1/4), add=T)

############################################################################
# Now run the simulations - set up

## Seed
set.seed(66)

## Parameters
# The linear correlation coefficient between exposure and outcome.
set_beta <- 0.5
constant <- 20

# Study data set of 20000 points for which we calculate the betas, (test)

paee = split(x = validation_data_women$PAEE, f = as.factor(validation_data_women$cam_index))

study_index_size = 5000
study_data = data.frame(PAEE = unlist(unname(lapply(X = paee, FUN = sample, size = study_index_size,replace=TRUE))))
study_data$cam_index = rep(x = as.numeric(names(paee)),each = study_index_size )
study_data$foo <-   rnorm(20000,(set_beta*study_data$PAEE) + constant,2.5)

####################################################################################################
################################## Regression Fitting ##############################################
####################################################################################################

# Initialise empty lists to store coeffcients
per_betas <- vector("numeric")
per_std_err <- vector("numeric")

ind_betas <- vector("numeric")
ind_std_err <- vector("numeric")


numtrials <- 100
for (i in 1:numtrials){
  # we assign random values from each category list of the validation set for each participant in the
  # study data set and then find the regression equation using this data. 
  # store the regression coefficient into list to create a dataframe afterwards
  
  study_data$paee_sample_per <- unlist(unname(lapply(X = paee, FUN = sample, size = study_index_size,replace=TRUE)))
  
  study_data$paee_sample_ind <- unlist(unname(lapply(X = paee,
                                                     FUN = function(paee_vals){
                                                       output = rep(x = sample(x = paee_vals, size = 1), times = study_index_size)
                                                       return(output)
                                                       
                                                     })))
  

  
  # calculate and store the coefficients 
  per_reg <- lm(formula=foo~paee_sample_per, data=study_data)
  per_coeff <- per_reg$coefficients["paee_sample_per"]
  per_std <- (summary(per_reg)$coefficients[,"Std. Error"])["paee_sample_per"]
  per_betas <- c(per_betas, per_coeff)
  per_std_err <- c(per_std_err, per_std)
  
  # calculate and store the coefficients 
  ind_reg <- lm(formula=foo~paee_sample_ind, data=study_data)
  ind_coeff <- ind_reg$coefficients["paee_sample_ind"]
  ind_std <- (summary(ind_reg)$coefficients[,"Std. Error"])["paee_sample_ind"]
  ind_betas <- c(ind_betas, ind_coeff)
  ind_std_err <- c(ind_std_err, ind_std)
  
  
}

study_data$paee_sample_ind_mean <- unlist(unname(lapply(X = paee,
                                                        FUN = function(paee_vals){
                                                          output = rep(x = mean(paee_vals), times = study_index_size)
                                                          return(output)
                                                          
                                                        })))   

reg_out_ind_mean <- lm(formula=foo~paee_sample_ind_mean, data=study_data)
reg_coeff_ind_mean <- reg_out_ind_mean$coefficients["paee_sample_ind_mean"]
reg_std_ind_mean <- (summary(reg_out_ind_mean)$coefficients[,"Std. Error"])["paee_sample_ind_mean"]



# Store values into dataframe
# Store values into dataframe
results <- as.data.frame(c(1:numtrials))
colnames(results) <- c("Trial")
results$reg_coeff_per <- per_betas
results$reg_stdError_per <- per_std_err
results$reg_coeff_ind <- ind_betas
results$reg_stdError_ind <- ind_std_err
results$reg_coeff_per_mean <- reg_coeff_ind_mean
results$reg_stdError_per_mean <- reg_std_ind_mean

lapply(X = results[,2:7],FUN =  quantile, probs=c(0.05,0.5,0.95), names=TRUE)

