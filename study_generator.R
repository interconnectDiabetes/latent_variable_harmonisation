study_size = 20000
mean_paee = 45
std_dev_paee = 15

# generate the study data only once and use this for all experiments
# This is the main advantage over previous methods where the study data itself changes
# as the validation study instrument was changing
# Just a normal distribution

study_data = data.frame(paee =  rnorm(n = study_size, mean = mean_paee, sd = std_dev_paee))

hist(x=study_data$paee)

# assign people into an index level, the instrument error defines the amount of overlap in the levels
# what we do it add an error to each person's PAEE
# then create the equally sized index levels with cut_number
# the instrument error is essentially the amount of overlap / standard deviation 
# when we plot the indices, we see the overlap

instrument_error = 5
number_of_indices = 4

paee_errors = rnorm(n = study_size, mean = 0, sd = instrument_error)

study_data$paee_error = study_data$paee + paee_errors

study_data$index = cut_number(x = study_data$paee_error,n = number_of_indices, labels =FALSE)

print(ggplot(study_data,aes(x=paee,group=index,fill=factor(index)))+
        geom_histogram(position="identity",alpha=0.5,binwidth=1)+
        ggtitle(label = "Study data")+scale_fill_manual(values = c("red", "grey", "seagreen3","blue","cyan","pink","yellow","brown")))

aggregate(x = study_data$paee, list(study_data$index), mean)
aggregate(x = study_data$paee, list(study_data$index), sd)

# create validation study in the same way, assume population has same properties as actual study, just smaller
# again we plot to see the overlap

Validation_study_size = 100

validation_study_data = data.frame(paee =  rnorm(n = Validation_study_size, mean = mean_paee, sd = std_dev_paee))

paee_errors = rnorm(n = Validation_study_size, mean = 0, sd = instrument_error)

validation_study_data$paee_error = validation_study_data$paee + paee_errors

validation_study_data$index = cut_number(x = validation_study_data$paee_error,n = number_of_indices, labels =FALSE)

print(ggplot(validation_study_data,aes(x=paee,group=index,fill=factor(index)))+
        geom_histogram(position="identity",alpha=0.5,binwidth=1)+
        ggtitle(label = "Validation study data")+scale_fill_manual(values = c("red", "grey", "seagreen3","blue","cyan","pink","yellow","brown")))

