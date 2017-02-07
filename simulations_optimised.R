# Using size_investigation as a template, this refactored code attempts to
# reduce time to run code and allow investigation of factors such as:
# 1. size of validation study,
# 2. amount of overlap,
# 3. number of index levels,
# 4. Type of draw - per person or per index etc.


## Author: Tom Bishop
## Date: 23.01.2017

###############################################################################
########################### DATA AND SETTINGS #################################
###############################################################################
  

library(ggplot2)

  ## Seed
  set.seed(66)
  
  ## Parameters

  # The linear correlation coefficient between exposure and outcome.
  set_beta <- 0.5
  constant <- 20


  # Index properties (Assumption that they are Gaussian)
  
  indices <- data.frame (mean = c(30, 40, 50, 60)
                        ,std_dev = c(5, 5, 5, 5)
                         ,index_indicator = c(1,2,3,4))
  
  # Defines runs to do with different sizes of validation study
 # validation_index_size <- c(seq(5,95, 10), seq(100, 1000, 100))
  validation_index_size <- c(5,100, 1000)  
#validation_index_size <- 5

  study_index_size = 5000
  
  #number of draws
  numtrials <- 100
  
  results_df <- data.frame()
  
#get our study data
study_data = data.frame(index = rep(x = indices$index_indicator,each = study_index_size ))

# for each index, pick all the values at once
study_data$paee = unlist(unname(lapply(X = split(x = indices, f = as.factor(indices$index_indicator)), 
                                       FUN = function(x){
                                         
                                         output = rnorm(n = study_index_size, mean = x$mean, sd = x$std_dev)
                                         return (output)
                                       })))

study_data$foo <-   rnorm(length(study_data$paee),(set_beta*study_data$paee) + constant,10)

# plot it


print(ggplot(study_data,aes(x=paee,group=index,fill=factor(index)))+
        geom_histogram(position="identity",alpha=0.5,binwidth=1)+
        ggtitle(label = "Study data")+scale_fill_manual(values = c("red", "grey", "seagreen3","blue")))


# this loop is for the different sizes of validation study
 for (i in 1:length(validation_index_size)){
   
   #create index column
   validation_study = data.frame(index = rep(x = indices$index_indicator,each = validation_index_size[i] ))
   
   #for each index, pick all the values at once from the appropriate distribution
   validation_study$paee = unlist(unname(lapply(X = split(x = indices, f = as.factor(indices$index_indicator)), 
                                                FUN = function(x){
     
                               output = rnorm(n = validation_index_size[i], 
                                              mean = x$mean, sd = x$std_dev)
                                return (output)
                              })))
   
   print(ggplot(validation_study,aes(x=paee,group=index,fill=factor(index)))+
     geom_histogram(position="identity",alpha=0.5,binwidth=1)+theme_bw()+
     ggtitle(label = paste("validation study size = ",validation_index_size[i]))+
       scale_fill_manual(values = c("red", "grey", "seagreen3","blue")))


   # Initialise some structures to store coeffcients
   
   betas_per <- vector("numeric")
   std_errs_per <- vector("numeric")
   
   betas_ind <- vector("numeric")
   std_errs_ind <- vector("numeric")
   
   # plot.new()
  # plot(x = study_data$paee, y= study_data$foo)
   
     for (j in 1:numtrials){
       # we assign random values from each category list of the validation set for each participant in the
       # study data set and then find the regression equation using this data. 
       # store the regression coefficient into list to create a dataframe afterwards
     
       #this is a draw per person, although it is done on a per index basis
       # for example if there are 100 people in the index it draws 100
       # values from the appropriate validation study index
       study_data$paee_sample_per <- unlist(unname(lapply(X = split(x=validation_study$paee, f= as.factor(validation_study$index)),
                                                      FUN = sample, size = study_index_size,replace=TRUE)))
       
       #this is a draw per index
       # for example if there are 100 people in the index it draws one
       # value from the appropriate validation study index
       # and replicates this 100 times
       study_data$paee_sample_ind <- unlist(unname(lapply(X = split(x=validation_study$paee, f= as.factor(validation_study$index)),
                                                      FUN = function(paee_vals){
                                                          output = rep(x = sample(x = paee_vals, size = 1), times = study_index_size)
                                                          return(output)
                                                        
                                                       })))
    
       # calculate and store the coefficients per person
       reg_out_per <- lm(formula=foo~paee_sample_per, data=study_data)
       reg_coeff_per <- reg_out_per$coefficients["paee_sample_per"]
       reg_std_per <- (summary(reg_out_per)$coefficients[,"Std. Error"])["paee_sample_per"]
       betas_per <- c(betas_per, reg_coeff_per)
       std_errs_per <- c(std_errs_per, reg_std_per)
       
       # calculate and store the coefficients per index
       reg_out_ind <- lm(formula=foo~paee_sample_ind, data=study_data)
       reg_coeff_ind <- reg_out_ind$coefficients["paee_sample_ind"]
       reg_std_ind <- (summary(reg_out_ind)$coefficients[,"Std. Error"])["paee_sample_ind"]
       betas_ind <- c(betas_ind, reg_coeff_ind)
       std_errs_ind <- c(std_errs_ind, reg_std_ind)
  
       
       
       #blah_per = lm(formula=foo~paee_sample_per, data=study_data)
       #abline(reg_out_per)
  
       
     }
  
  print(ggplot(study_data,aes(x=paee_sample_per,group=index,fill=factor(index)))+
          geom_histogram(position="identity",alpha=0.5,binwidth=1)+theme_bw()+
          ggtitle(label = paste("PAEE person - validation study size = ",validation_index_size[i]))+
          scale_fill_manual(values = c("red", "grey", "seagreen3","blue")))
  
  
  print(ggplot(study_data,aes(x=paee_sample_ind,group=index,fill=factor(index)))+
          geom_histogram(position="identity",alpha=0.5,binwidth=1)+theme_bw()+
          ggtitle(label = paste(" PAEE index - validation study size = ",validation_index_size[i]))+
          scale_fill_manual(values = c("red", "grey", "seagreen3","blue")))
  
   
   # Do the regression using validation study means (for comparison)
   
   study_data$paee_sample_ind_mean <- unlist(unname(lapply(X = split(x=validation_study$paee, f= as.factor(validation_study$index)),
                                                      FUN = function(paee_vals){
                                                        output = rep(x = mean(paee_vals), times = study_index_size)
                                                        return(output)
                                                        
                                                      })))   
   
   reg_out_ind_mean <- lm(formula=foo~paee_sample_ind_mean, data=study_data)
   reg_coeff_ind_mean <- reg_out_ind_mean$coefficients["paee_sample_ind_mean"]
   reg_std_ind_mean <- (summary(reg_out_ind_mean)$coefficients[,"Std. Error"])["paee_sample_ind_mean"]
   
   # Store values into dataframe
   results <- as.data.frame(c(1:numtrials))
   colnames(results) <- c("Trial")
   results$valid_size <- rep(x=validation_index_size[i], times = numtrials)
   results$reg_coeff_per <- betas_per
   results$reg_stdError_per <- std_errs_per
   results$reg_coeff_ind <- betas_ind
   results$reg_stdError_ind <- std_errs_ind
   results$reg_coeff_per_mean <- reg_coeff_ind_mean
   results$reg_stdError_per_mean <- reg_std_ind_mean
   
   results_df <- rbind(results_df, results)
 }

  # this is code to get the results data frame into a summary
  temp_output <- aggregate(results_df[,3:8], by=list(valid_size = results_df$valid_size), quantile, probs=c(0.05,0.5,0.95), names=TRUE)
  final_output = data.frame(trials = temp_output[,1])
  for (k in 2:ncol(temp_output)){
   
    temp = as.data.frame(temp_output[,k])
    colnames(temp) <- paste(colnames(temp_output)[k], colnames(temp), sep = "_")
    final_output = cbind(final_output, temp)
  }
  
  
# some code for showing some plots - not really finished yet

 
#   plot(x = study_data$paee_sample_per, y= study_data$foo)
#   blah_per = lm(formula=foo~paee_sample_per, data=study_data)
#   abline(blah_per)
#   
#   plot(x = study_data$paee_sample_ind, y= study_data$foo)
#   blah_ind = lm(formula=foo~paee_sample_ind, data=study_data)
#   abline(blah_ind)