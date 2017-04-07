# BootStrap the Bootstrap with a more realistic study generator.
# Author: Paul Scherer and Tom Bishop
# Date: 20.03.2017

###############################################################################
###################### R Environment Settings #################################
###############################################################################
setwd('V:/Studies/InterConnect/Internal/Latent variable harmonisation/plots')
library(ggplot2)
library(reshape2)
library(parallel)

###############################################################################
########################### DATA AND SETTINGS #################################
###############################################################################
# Calculate the number of cores and initiate cluster
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

# Base Data set properties
set_beta <- 0.5
constant <- 20

study_size <- 20000
mean_paee <- 45

# Test Properties
num_trials <- 10

std_dev_paee = 15
study_data = data.frame(paee =  rnorm(n = study_size, mean = mean_paee, sd = std_dev_paee))

###############################################################################
########################### Functions #########################################
###############################################################################
updateStudyData <- function(coh_base, study_data, number_of_indices){
	paee_errors = rnorm(n = study_size, mean = 0, sd = coh_base$std_dev[1])
	study_data$paee_error = study_data$paee + paee_errors
	study_data$index = cut_number(x = study_data$paee_error,n = number_of_indices, labels =FALSE)
	study_data$foo = rnorm(length(study_data$paee), (set_beta*study_data$paee) + constant, 10)
	study_data = study_data[with(study_data, order(index)),]
	return(study_data)
}

createValidationData <- function(coh_base, val_size) {
	number_of_indices = length(coh_base$indices)

	validation_data = data.frame(paee =  rnorm(n = val_size, mean = mean_paee, sd = std_dev_paee))
	paee_errors = rnorm(n = val_size, mean = 0, sd = coh_base$std_dev[1])
	validation_data$paee_error = validation_data$paee + paee_errors
	validation_data$index = cut_number(x = validation_data$paee_error,n = number_of_indices, labels =FALSE)

	return (validation_data)
}

bootstrapRun <- function(coh_base, study_size, val_size, study_data) {
	number_of_indices = length(coh_base$indices)
	validation_index_size <- round(val_size/number_of_indices)
	study_index_size <- round(study_size/number_of_indices)
	# spawn the study data
	study_data_cp = updateStudyData(coh_base, study_data, number_of_indices)
	# spawn the validation data
	validation_data = createValidationData(coh_base, val_size)

	# per cohort base result vector
	betas <- vector("numeric")
	std_errs <- vector("numeric")
	# overall results dataframe
	results_df <- data.frame()

	results <- as.data.frame(c(1:num_trials))
	colnames(results) <- c("NumTrial")
	results$valid_size <- rep(x=validation_index_size, times = num_trials)
	results$standard_deviation <- rep(x=coh_base$std_dev[1], times = num_trials)
	
	chickens  <- parLapply(cl, X=1:num_trials, fun=function(x){
		# Create the bootstrap, sampling from the validation data
		bootstrap_validation <- data.frame(index = rep(x = coh_base$indices, each= validation_index_size))
		bootstrap_validation$paee <- unlist(unname(lapply(X = split(x=validation_data$paee, f= as.factor(validation_data$index)), 
			FUN = sample, size = validation_index_size,replace=TRUE)))

		# # Regression and storing of regression coefficients and stderrs
		# study_data_cp$paee_sample_ind_mean <- unlist(unname(lapply(X = split(x=bootstrap_validation$paee, f= as.factor(bootstrap_validation$index)),
		#   FUN = function(paee_vals){
		#     output = rep(x = mean(paee_vals), times = study_index_size)
		#     return(output)
		#   })))
		# reg_out_ind_mean <- lm(formula=foo~paee_sample_ind_mean, data=study_data_cp)

		# new bootstrap study data as well through sampling rows
		means_of_boots <- (split(x=bootstrap_validation$paee, f= as.factor(bootstrap_validation$index)))
		means_boots_list = vector(mode="list", length = length(means_of_boots))
		for (i in 1:length(means_of_boots)) {
			means_boots_list[i] = mean(unname(unlist((split(x=bootstrap_validation$paee, f= as.factor(bootstrap_validation$index)))[i])))
		}
		study_data_cp_bt = study_data_cp[sample(nrow(study_data_cp), 20000, replace=TRUE), ]
		study_data_cp_bt$paee_sample_ind_mean <- unlist(lapply(X=study_data_cp_bt$index, FUN=function(index_val){
			output =  means_boots_list[index_val]
			}))
		reg_out_ind_mean <- lm(formula=foo~paee_sample_ind_mean, data=study_data_cp_bt)

		# Create the validation data copy (because of parallelization) for the lambda calculation		
		validation_data_cp <- validation_data
		validation_data_cp = validation_data_cp[with(validation_data_cp, order(index)),]
		validation_data_cp$paee_means_bt <- unlist(unname(lapply(X = split(x=bootstrap_validation$paee, f= as.factor(bootstrap_validation$index)),
		  FUN = function(paee_vals){
		    output = rep(x = mean(paee_vals), times = validation_index_size)
		    return(output)
		  })))
		lambda <- lm(formula =paee~paee_means_bt, data=validation_data_cp)

		# Estimate the standard error of the corrected estimate as currently it doesnt take the 2nd order variability of lambda
		# USING DELTA METHOD
		# variance = stdError * sqrt(numberofpeople) all squared
		var_lambda = (sqrt(val_size) * summary(reg_out_ind_mean)$coefficients["paee_sample_ind_mean","Std. Error"])^2
		var_Beta = (sqrt(val_size) * summary(lambda)$coefficients["paee_means_bt","Std. Error"])^2
		lambda_pure = unlist(unname(lambda$coefficients["paee_means_bt"]))
		beta_lambda_div_sq = (unname(unlist(reg_out_ind_mean$coefficients["paee_sample_ind_mean"]/(lambda$coefficients["paee_means_bt"])^2)))^2
		delta_variance = (var_Beta / (lambda_pure)^2) + beta_lambda_div_sq * var_lambda
		delta_stdError = sqrt(delta_variance)/sqrt(val_size) 

		# Return dual output into a list
		output = data.frame(c(reg_out_ind_mean$coefficients["paee_sample_ind_mean"], lambda$coefficients["paee_means_bt"], delta_stdError))
		return (output)
	})

	results$reg_cor_per_mean <- unlist(lapply(X = chickens, FUN = function(x){output = x[[1]][1]}))
	results$lambda <- unlist(lapply(X = chickens, FUN = function(x){output = x[[1]][2]}))
	results$corrected_cor <- results$reg_cor_per_mean/results$lambda
	results$delta_errors <- unlist(lapply(X = chickens, function(x){output = 2*1.96*(x[[1]][3])}))
	results_df <- rbind(results_df, results)
	
	# Summarizing the results dataframe
	temp_output <- aggregate(results_df[,4:7], by=list(valid_size = results_df$valid_size), quantile, probs=c(0.025,0.5,0.975), names=TRUE)
	final_output = data.frame(validation_size = temp_output[,1])

	for (k in 2:ncol(temp_output)){
	temp = as.data.frame(temp_output[,k])
	colnames(temp) <- paste(colnames(temp_output)[k], colnames(temp), sep = "_")
	final_output = cbind(final_output, temp)
	}
	final_output$val_size <- val_size
	return (final_output)
}

absDiff <- function(x,y){
	return (abs(x-y))
}



###############################################################################
########################### Simulation Section ################################
###############################################################################


# for later summation of values in the results dataframe
numSeeds <- 5
minmax_list = vector(mode="list", length = numSeeds)
accuracy_list = vector(mode="list", length = numSeeds)
minmax_list_cor = vector(mode="list", length = numSeeds)
accuracy_list_cor = vector(mode="list", length = numSeeds)

delta_errors_list = vector(mode="list", length = numSeeds)


number_of_indices <- 4

#progress bar for seeds completed
pbt <- txtProgressBar(min = 1, max = numSeeds, style = 3)
for (seeds in 1:numSeeds){
  set.seed(seeds)
	results = data.frame()
	for (val_size in seq(from=100, to=400, by=20)){
		for (standard_dev in 5:15){
			# Defining a base generator for one cohort (which spawns validation, study data)
			coh_base = data.frame(indices = c(1:number_of_indices), std_dev = standard_dev)
			bootRunResult = bootstrapRun(coh_base, study_size, val_size, study_data)
			bootRunResult = cbind(bootRunResult, standard_dev)
			results = rbind(results, bootRunResult)
		}
	}

	# Create a heatmap of 'accuracy' through absolute difference in the 2.5% and 97.5% tiles
	results$minMaxDiff <- unlist(unname(mapply(FUN=absDiff, results$`reg_cor_per_mean_2.5%`, results$`reg_cor_per_mean_97.5%`)))
	results$minMaxDiff_cor <- unlist(unname(mapply(FUN=absDiff, results$`corrected_cor_2.5%`, results$`corrected_cor_97.5%`)))
	# Create a heatmap of 'accuracy' through absolute difference of the true 0.5 and the reported median value
	results$accuracyDiff <- unlist(unname(mapply(FUN=absDiff, set_beta, results$`reg_cor_per_mean_50%`)))
	results$accuracyDiff_cor <- unlist(unname(mapply(FUN=absDiff, set_beta, results$`corrected_cor_50%`)))

	minmax_list[[seeds]] <- (results$minMaxDiff)
	minmax_list_cor[[seeds]] <- (results$minMaxDiff_cor)
	accuracy_list[[seeds]] <- (results$accuracyDiff)
	accuracy_list_cor[[seeds]] <- (results$accuracyDiff_cor)

	delta_errors_list[[seeds]] <- (results$'delta_errors_50%')

	setTxtProgressBar(pbt, seeds)
}
close(pbt)

# Then put the sum of the minmaxes together into the results dataframe
results$minMaxDiff <- Reduce("+", x = minmax_list)
results$accuracyDiff <- Reduce("+", x = accuracy_list)

# Divide by number of seed for scaling
results$minMaxDiff <- unlist(unname(lapply(X = results$minMaxDiff, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results$accuracyDiff <- unlist(unname(lapply(X = results$accuracyDiff, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

# heatmap confidence
ggplot(results, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = minMaxDiff), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference of 95% Interval")

# heatmap accuracy
ggplot(results, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = accuracyDiff), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference from Truth")

# Then put the sum of the minmaxes together into the results dataframe
results$minMaxDiff_cor <- Reduce("+", x = minmax_list_cor)
results$accuracyDiff_cor <- Reduce("+", x = accuracy_list_cor)

results$delta_error <- Reduce("+", x = delta_errors_list)

# Divide by number of seed for scaling
results$minMaxDiff_cor <- unlist(unname(lapply(X = results$minMaxDiff_cor, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results$accuracyDiff_cor <- unlist(unname(lapply(X = results$accuracyDiff_cor, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results$delta_error <- unlist(unname(lapply(X = results$delta_error, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))



# heatmap confidence corrected
ggplot(results, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = minMaxDiff_cor), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference of 95% Interval")

# heatmap accuracy corrected
ggplot(results, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = accuracyDiff_cor), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference from Truth")

# heatmap accuracy corrected
ggplot(results, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = delta_error), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Delta Interval Calculation")


































# for later summation of values in the results dataframe
numSeeds <- 5
minmax_list = vector(mode="list", length = numSeeds)
accuracy_list = vector(mode="list", length = numSeeds)
minmax_list_cor = vector(mode="list", length = numSeeds)
accuracy_list_cor = vector(mode="list", length = numSeeds)

delta_errors_list = vector(mode="list", length = numSeeds)


number_of_indices <- 2

#progress bar for seeds completed
pbt <- txtProgressBar(min = 1, max = numSeeds, style = 3)
for (seeds in 1:numSeeds){
  set.seed(seeds)
	results_2 = data.frame()
	for (val_size in seq(from=100, to=400, by=20)){
		for (standard_dev in 5:15){
			# Defining a base generator for one cohort (which spawns validation, study data)
			coh_base = data.frame(indices = c(1:number_of_indices), std_dev = standard_dev)
			bootRunResult = bootstrapRun(coh_base, study_size, val_size, study_data)
			bootRunResult = cbind(bootRunResult, standard_dev)
			results_2 = rbind(results_2, bootRunResult)
		}
	}

	# Create a heatmap of 'accuracy' through absolute difference in the 2.5% and 97.5% tiles
	results_2$minMaxDiff <- unlist(unname(mapply(FUN=absDiff, results_2$`reg_cor_per_mean_2.5%`, results_2$`reg_cor_per_mean_97.5%`)))
	results_2$minMaxDiff_cor <- unlist(unname(mapply(FUN=absDiff, results_2$`corrected_cor_2.5%`, results_2$`corrected_cor_97.5%`)))
	# Create a heatmap of 'accuracy' through absolute difference of the true 0.5 and the reported median value
	results_2$accuracyDiff <- unlist(unname(mapply(FUN=absDiff, set_beta, results_2$`reg_cor_per_mean_50%`)))
	results_2$accuracyDiff_cor <- unlist(unname(mapply(FUN=absDiff, set_beta, results_2$`corrected_cor_50%`)))

	minmax_list[[seeds]] <- (results_2$minMaxDiff)
	minmax_list_cor[[seeds]] <- (results_2$minMaxDiff_cor)
	accuracy_list[[seeds]] <- (results_2$accuracyDiff)
	accuracy_list_cor[[seeds]] <- (results_2$accuracyDiff_cor)

	delta_errors_list[[seeds]] <- (results_2$'delta_errors_50%')

	setTxtProgressBar(pbt, seeds)
}
close(pbt)

# Then put the sum of the minmaxes together into the results dataframe
results_2$minMaxDiff <- Reduce("+", x = minmax_list)
results_2$accuracyDiff <- Reduce("+", x = accuracy_list)

# Divide by number of seed for scaling
results_2$minMaxDiff <- unlist(unname(lapply(X = results_2$minMaxDiff, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results_2$accuracyDiff <- unlist(unname(lapply(X = results_2$accuracyDiff, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

# heatmap confidence
ggplot(results_2, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = minMaxDiff), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference of 95% Interval")

# heatmap accuracy
ggplot(results_2, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = accuracyDiff), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference from Truth")

# Then put the sum of the minmaxes together into the results dataframe
results_2$minMaxDiff_cor <- Reduce("+", x = minmax_list_cor)
results_2$accuracyDiff_cor <- Reduce("+", x = accuracy_list_cor)

results_2$delta_error <- Reduce("+", x = delta_errors_list)

# Divide by number of seed for scaling
results_2$minMaxDiff_cor <- unlist(unname(lapply(X = results_2$minMaxDiff_cor, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results_2$accuracyDiff_cor <- unlist(unname(lapply(X = results_2$accuracyDiff_cor, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results_2$delta_error <- unlist(unname(lapply(X = results_2$delta_error, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))



# heatmap confidence corrected
ggplot(results_2, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = minMaxDiff_cor), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference of 95% Interval")

# heatmap accuracy corrected
ggplot(results_2, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = accuracyDiff_cor), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference from Truth")

# heatmap accuracy corrected
ggplot(results_2, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = delta_error), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Delta Interval Calculation")

































# for later summation of values in the results dataframe
numSeeds <- 5
minmax_list = vector(mode="list", length = numSeeds)
accuracy_list = vector(mode="list", length = numSeeds)
minmax_list_cor = vector(mode="list", length = numSeeds)
accuracy_list_cor = vector(mode="list", length = numSeeds)

delta_errors_list = vector(mode="list", length = numSeeds)


number_of_indices <- 16

#progress bar for seeds completed
pbt <- txtProgressBar(min = 1, max = numSeeds, style = 3)
for (seeds in 1:numSeeds){
  set.seed(seeds)
	results_16 = data.frame()
	for (val_size in seq(from=100, to=400, by=20)){
		for (standard_dev in 5:15){
			# Defining a base generator for one cohort (which spawns validation, study data)
			coh_base = data.frame(indices = c(1:number_of_indices), std_dev = standard_dev)
			bootRunResult = bootstrapRun(coh_base, study_size, val_size, study_data)
			bootRunResult = cbind(bootRunResult, standard_dev)
			results_16 = rbind(results_16, bootRunResult)
		}
	}

	# Create a heatmap of 'accuracy' through absolute difference in the 2.5% and 97.5% tiles
	results_16$minMaxDiff <- unlist(unname(mapply(FUN=absDiff, results_16$`reg_cor_per_mean_2.5%`, results_16$`reg_cor_per_mean_97.5%`)))
	results_16$minMaxDiff_cor <- unlist(unname(mapply(FUN=absDiff, results_16$`corrected_cor_2.5%`, results_16$`corrected_cor_97.5%`)))
	# Create a heatmap of 'accuracy' through absolute difference of the true 0.5 and the reported median value
	results_16$accuracyDiff <- unlist(unname(mapply(FUN=absDiff, set_beta, results_16$`reg_cor_per_mean_50%`)))
	results_16$accuracyDiff_cor <- unlist(unname(mapply(FUN=absDiff, set_beta, results_16$`corrected_cor_50%`)))

	minmax_list[[seeds]] <- (results_16$minMaxDiff)
	minmax_list_cor[[seeds]] <- (results_16$minMaxDiff_cor)
	accuracy_list[[seeds]] <- (results_16$accuracyDiff)
	accuracy_list_cor[[seeds]] <- (results_16$accuracyDiff_cor)

	delta_errors_list[[seeds]] <- (results_16$'delta_errors_50%')

	setTxtProgressBar(pbt, seeds)
}
close(pbt)

# Then put the sum of the minmaxes together into the results dataframe
results_16$minMaxDiff <- Reduce("+", x = minmax_list)
results_16$accuracyDiff <- Reduce("+", x = accuracy_list)

# Divide by number of seed for scaling
results_16$minMaxDiff <- unlist(unname(lapply(X = results_16$minMaxDiff, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results_16$accuracyDiff <- unlist(unname(lapply(X = results_16$accuracyDiff, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

# heatmap confidence
ggplot(results_16, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = minMaxDiff), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference of 95% Interval")

# heatmap accuracy
ggplot(results_16, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = accuracyDiff), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference from Truth")

# Then put the sum of the minmaxes together into the results dataframe
results_16$minMaxDiff_cor <- Reduce("+", x = minmax_list_cor)
results_16$accuracyDiff_cor <- Reduce("+", x = accuracy_list_cor)

results_16$delta_error <- Reduce("+", x = delta_errors_list)

# Divide by number of seed for scaling
results_16$minMaxDiff_cor <- unlist(unname(lapply(X = results_16$minMaxDiff_cor, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results_16$accuracyDiff_cor <- unlist(unname(lapply(X = results_16$accuracyDiff_cor, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))

results_16$delta_error <- unlist(unname(lapply(X = results_16$delta_error, FUN = function(x){
	output = x/numSeeds
	return(output)
	})))



# heatmap confidence corrected
ggplot(results_16, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = minMaxDiff_cor), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference of 95% Interval")

# heatmap accuracy corrected
ggplot(results_16, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = accuracyDiff_cor), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Absolute Difference from Truth")

# heatmap accuracy corrected
ggplot(results_16, aes(val_size, standard_dev )) +
  geom_tile(aes(fill = delta_error), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Standard Deviation") +
  xlab("Validation Size") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Delta Interval Calculation")



# stopcluster
stopCluster(cl)



