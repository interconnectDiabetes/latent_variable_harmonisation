# BootStrap the Bootstrap
# Author: Paul Scherer and Tom Bishop
# Date: 03.02.2017

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
no_cores <- 3
cl <- makeCluster(no_cores)

# Base Data set properties
set_beta <- 0.5
constant <- 20

# Study Data properties
study_index_size <- 5000

# Test Properties
num_trials <- 250
###############################################################################
########################### Functions #########################################
###############################################################################
createStudyData <- function(coh_base, study_index_size){
	# given a cohort base, ie means, index, stdev it returns a study data dataframe
	# at this stage it can only accomodate the same index size for each index
	study_data = data.frame(index = rep(x = coh_base$indices, each=study_index_size))
	study_data$paee = unlist(unname(lapply(X = split(x = coh_base, f = as.factor(coh_base$indices)), 
		FUN = function(x){
			output = rnorm(n =study_index_size, mean = x$means, sd = x$std_dev)
			return (output)
		})))
	study_data$foo <- rnorm(length(study_data$paee), (set_beta*study_data$paee) + constant, 10) #note 0 noise
	return(study_data)
}

createValidationData <- function(coh_base, validation_index_size) {
	# given a cohort base it generates the validation data
	validation_data = data.frame(index = rep(x = coh_base$indices, each = validation_index_size))
	validation_data$paee = unlist(unname(lapply(X = split(x = coh_base, f = as.factor(coh_base$indices)), 
		FUN=function(x){
			output = rnorm(n = validation_index_size, mean=x$mean, sd=x$std_dev)
		})))
	return (validation_data)
}

bootstrapRun <- function(coh_base, study_index_size, val_size) {
	validation_index_size <- round(val_size/length(coh_base$means))
	# spawn the study data
	study_data = createStudyData(coh_base, study_index_size)
	# spawn the validation data
	validation_data = createValidationData(coh_base, validation_index_size)

	# per cohort base result vector
	betas <- vector("numeric")
	std_errs <- vector("numeric")
	# overall results dataframe
	results_df <- data.frame()

	results <- as.data.frame(c(1:num_trials))
	colnames(results) <- c("NumTrial")
	results$valid_size <- rep(x=validation_index_size, times = num_trials)
	results$standard_deviation <- rep(x=coh_base$std_dev[1], times = num_trials)
	results$reg_coeff_per_mean  <- parLapply(cl, X=1:num_trials, fun=function(x){
		bootstrap_validation <- data.frame(index = rep(x = coh_base$indices, each= validation_index_size))
		bootstrap_validation$paee <- unlist(unname(lapply(X = split(x=validation_data$paee, f= as.factor(validation_data$index)), 
			FUN = sample, size = validation_index_size,replace=TRUE)))

		# Regression and storing of regression coefficients and stderrs
		study_data$paee_sample_ind_mean <- unlist(unname(lapply(X = split(x=bootstrap_validation$paee, f= as.factor(bootstrap_validation$index)),
		  FUN = function(paee_vals){
		    output = rep(x = mean(paee_vals), times = study_index_size)
		    return(output)
		  })))
		reg_out_ind_mean <- lm(formula=foo~paee_sample_ind_mean, data=study_data)
		return (reg_out_ind_mean$coefficients["paee_sample_ind_mean"])
	})
	results$reg_coeff_per_mean <- unlist(results$reg_coeff_per_mean)
	results_df <- rbind(results_df, results)

	# Summarizing the results dataframe
	temp_output <- aggregate(results_df[,4], by=list(valid_size = results_df$valid_size), quantile, probs=c(0.025,0.5,0.975), names=TRUE)
	final_output = data.frame(validation_size = temp_output[,1])

	# # Summarizing the results dataframe
	# temp_output <- aggregate(results_df[,3:10], by=list(valid_size = results_df$valid_size), quantile, probs=c(0.5), names=TRUE)
	# final_output = data.frame(validation_size = temp_output[,1])

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
# 16
# for later summation of values in the results dataframe
numSeeds <- 25
minmax_list = vector(mode="list", length = numSeeds)
accuracy_list = vector(mode="list", length = numSeeds)

# Study Data properties
study_index_size <- 1250

#progress bar for seeds completed
pbt <- txtProgressBar(min = 1, max = numSeeds, style = 3)
for (seeds in 1:numSeeds){
  set.seed(seeds)
	results = data.frame()
	for (val_size in seq(from=100, to=400, by=20)){
		for (standard_dev in 5:15){
			# Defining a base generator for one cohort (which spawns validation, study data)
			coh_base = data.frame(means = c(30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60), indices = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), std_dev = standard_dev) # change stddev and index size to variables later
			bootRunResult = bootstrapRun(coh_base, study_index_size, val_size)
			bootRunResult = cbind(bootRunResult, standard_dev)
			results = rbind(results, bootRunResult)
		}
	}
	# Create a heatmap of 'accuracy' through absolute difference in the 2.5% and 97.5% tiles
	results$minMaxDiff <- unlist(unname(mapply(FUN=absDiff, results$`x_2.5%`, results$`x_97.5%`)))
	# Create a heatmap of 'accuracy' through absolute difference of the true 0.5 and the reported median value
	results$accuracyDiff <- unlist(unname(mapply(FUN=absDiff, set_beta, results$`x_50%`)))

	minmax_list[[seeds]] <- (results$minMaxDiff)
	accuracy_list[[seeds]] <- (results$accuracyDiff)
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


####### 32 ########

# for later summation of values in the results dataframe
numSeeds <- 25
minmax_list = vector(mode="list", length = numSeeds)
accuracy_list = vector(mode="list", length = numSeeds)

# Study Data properties
study_index_size <- 625

#progress bar for seeds completed
pbt <- txtProgressBar(min = 1, max = numSeeds, style = 3)
for (seeds in 1:numSeeds){
  set.seed(seeds)
	results_2 = data.frame()
	for (val_size in seq(from=100, to=400, by=20)){
		for (standard_dev in 5:15){
			# Defining a base generator for one cohort (which spawns validation, study data)
			coh_base = data.frame(means = c(30, 30.96774194, 31.93548387, 32.90322581, 33.87096774, 34.83870968, 
			35.80645161, 36.77419355, 37.74193548,  38.70967742,  39.67741935,  40.64516129, 41.61290323,  42.58064516, 
			43.5483871 ,  44.51612903, 45.48387097,  46.4516129 ,  47.41935484,  48.38709677, 49.35483871,  50.32258065, 
			51.29032258,  52.25806452, 53.22580645,  54.19354839,  55.16129032,  56.12903226, 57.09677419,  58.06451613, 
			59.03225806, 60), indices = c(1:32), std_dev = standard_dev) # change stddev and index size to variables later
			bootRunResult = bootstrapRun(coh_base, study_index_size, val_size)
			bootRunResult = cbind(bootRunResult, standard_dev)
			results_2 = rbind(results_2, bootRunResult)
		}
	}
	# Create a heatmap of 'accuracy' through absolute difference in the 2.5% and 97.5% tiles
	results_2$minMaxDiff <- unlist(unname(mapply(FUN=absDiff, results_2$`x_2.5%`, results_2$`x_97.5%`)))
	# Create a heatmap of 'accuracy' through absolute difference of the true 0.5 and the reported median value
	results_2$accuracyDiff <- unlist(unname(mapply(FUN=absDiff, set_beta, results_2$`x_50%`)))

	minmax_list[[seeds]] <- (results_2$minMaxDiff)
	accuracy_list[[seeds]] <- (results_2$accuracyDiff)
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


# #stop our cluster
stopCluster(cl)