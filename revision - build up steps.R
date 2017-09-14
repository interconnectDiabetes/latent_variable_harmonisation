library(ggplot2)
library(reshape2)
library(parallel)


# Properties of ground truth relationship between exposure and outcome
set_beta <- 0.5
constant <- 20

# Properties of the study
study_size <- 20000
mean_exp <- 45
std_dev_exp = 15
number_of_indices = 4

# Properties of the instrument for measuring exposure
exp_error = 15

#Properties of the validation study
val_size = 400

# set up study data

study_data = data.frame(exp =  rnorm(n = study_size, mean = mean_exp, sd = std_dev_exp))
study_data$outcome = set_beta*study_data$exp
exp_errors = rnorm(n = study_size, mean = 0, sd = exp_error)
study_data$exp_error = study_data$exp + exp_errors
study_data$index = cut_number(x = study_data$exp_error,n = number_of_indices, labels =FALSE)
study_data = study_data[with(study_data, order(index)),]

# set up validation data

validation_data = data.frame(exp =  rnorm(n = val_size, mean = mean_exp, sd = std_dev_exp))
exp_errors = rnorm(n = val_size, mean = 0, sd = exp_error)
validation_data$exp_error = validation_data$exp + exp_errors
validation_data$index = cut_number(x = validation_data$exp_error,n = number_of_indices, labels =FALSE)

# calculate validation study per index mean

index_means_list = vector(mode="list", length = number_of_indices)
for (i in 1:number_of_indices) {
  index_means_list[i] = mean(unname(unlist((split(x=validation_data$exp, f= as.factor(validation_data$index)))[i])))
}

study_means_list = vector(mode="list", length = number_of_indices)
for (i in 1:number_of_indices) {
  study_means_list[i] = mean(unname(unlist((split(x=study_data$exp, f= as.factor(study_data$index)))[i])))
}

#apply this index mean to the study data

study_data$exp_ind_mean <- unlist(lapply(X=study_data$index, FUN=function(index_val){
  output =  index_means_list[index_val]
}))

validation_data$exp_ind_mean <- unlist(lapply(X=validation_data$index, FUN=function(index_val){
  output =  index_means_list[index_val]
}))




# print(ggplot(study_data,aes(x=exp,group=index,fill=factor(index)))+
#         geom_histogram(position="identity",alpha=0.5,binwidth=1)+
#         ggtitle(label = "Study data"))
# 
# print(ggplot(validation_data,aes(x=exp,group=index,fill=factor(index)))+
#         geom_histogram(position="identity",alpha=0.5,binwidth=1)+
#         ggtitle(label = "Validation data"))


# if I knew the true paee
#true_result = lm(formula=outcome~exp, data=study_data)
#summary(true_result)

# of course, it is close to 0.5

# if I have paee with random error
#error_result = lm(formula=outcome~exp_error, data=study_data)
#summary(error_result)

# There has been attenuation due to the random error added, get
# an answer of about 0.35 instead

# The classic using the index mean

#plot(formula=outcome~exp_ind_mean, data=study_data)

index_mean_result = lm(formula=outcome~exp_ind_mean, data=study_data)
summary(index_mean_result)

# This is very close to using the study error values, both in terms of estimate and
# standard error. This feels too good to be true. Either something wrong with the 
# simulation or we need to penalise the estimate with a larger SE

# regression calibration - not really sure if this is relevant here.

validation_data$exp_means <- unlist(lapply(X=validation_data$index, FUN=function(index_val){
  output =  index_means_list[index_val]
}))
reg_lambda <- lm(formula =exp_error~exp_means, data=validation_data)
summary(reg_lambda)

# of course, lambda is 1

# try and give some error to the mean by doing the bootstraps

validation_index_size <- round(val_size/number_of_indices)

results <- data.frame()
for (i in 1:100){
  bootstrap_validation <- data.frame(index = rep(x = c(1:number_of_indices), each= validation_index_size))
  bootstrap_validation$exp <- unlist(unname(lapply(X = split(x=validation_data$exp, f= as.factor(validation_data$index)), 
                                                    FUN = sample, size = validation_index_size, replace=TRUE)))
  
  means_boots_list = vector(mode="list", length = number_of_indices)
  for (i in 1:number_of_indices) {
    means_boots_list[i] = mean(unname(unlist((split(x=bootstrap_validation$exp, f= as.factor(bootstrap_validation$index)))[i])))
  }
  
  study_data$exp_ind_mean <- unlist(lapply(X=study_data$index, FUN=function(index_val){
    output =  means_boots_list[index_val]
  }))
  
  reg_out_ind_mean <- lm(formula=outcome~exp_ind_mean, data=study_data)
  summary(reg_out_ind_mean)
  results = rbind(results, reg_out_ind_mean$coefficients["exp_ind_mean"])
}


######################################################################


# Properties of ground truth relationship between exposure and outcome
set_beta <- 0.5
constant <- 20

# Properties of the study
study_size <- 20000
mean_exp <- 45
std_dev_exp = 15
number_of_indices = 4
bins = c(-Inf, 20,45,65, Inf)

# Properties of the instrument for measuring exposure
exp_error = 15

#Properties of the validation study
val_size = 400

# set up study data

study_data = data.frame(exp =  rnorm(n = study_size, mean = mean_exp, sd = std_dev_exp))
study_data$outcome = set_beta*study_data$exp
exp_errors = rnorm(n = study_size, mean = 0, sd = exp_error)
study_data$exp_error = study_data$exp + exp_errors

study_data$index = cut(x = study_data$exp_error, breaks = bins, labels = FALSE)

#study_data$index = cut_number(x = study_data$exp_error,n = number_of_indices, labels =FALSE)
study_data = study_data[with(study_data, order(index)),]

# set up validation data

validation_data = data.frame(exp =  rnorm(n = val_size, mean = mean_exp, sd = std_dev_exp))
exp_errors = rnorm(n = val_size, mean = 0, sd = exp_error)
validation_data$exp_error = validation_data$exp + exp_errors
validation_data$index = cut(x = validation_data$exp_error, breaks = bins, labels = FALSE)
                            
# calculate validation study per index mean

index_means_list = vector(mode="list", length = number_of_indices)
for (i in 1:number_of_indices) {
  index_means_list[i] = mean(unname(unlist((split(x=validation_data$exp, f= as.factor(validation_data$index)))[i])))
}

study_means_list = vector(mode="list", length = number_of_indices)
for (i in 1:number_of_indices) {
  study_means_list[i] = mean(unname(unlist((split(x=study_data$exp, f= as.factor(study_data$index)))[i])))
}

#apply this index mean to the study data

study_data$exp_ind_mean <- unlist(lapply(X=study_data$index, FUN=function(index_val){
  output =  index_means_list[index_val]
}))

validation_data$exp_ind_mean <- unlist(lapply(X=validation_data$index, FUN=function(index_val){
  output =  index_means_list[index_val]
}))

index_mean_result = lm(formula=outcome~exp_ind_mean, data=study_data)
summary(index_mean_result)                
