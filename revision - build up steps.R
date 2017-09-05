# Base Data set properties
set_beta <- 0.5
constant <- 20

study_size <- 20000
mean_paee <- 45
std_dev_paee = 15
study_data = data.frame(paee =  rnorm(n = study_size, mean = mean_paee, sd = std_dev_paee))
study_data$foo = rnorm(length(study_data$paee), (set_beta*study_data$paee) + constant, 10)
