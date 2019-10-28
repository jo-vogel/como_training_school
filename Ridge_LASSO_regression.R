##########################################################################################
###########               Ridge and LASSO regression                           ###########
###########   on unstandardised and standardised crop model data               ###########
###########                                                                    ###########
###########       Author: Pauline Rivoire                                      ###########
##########################################################################################



#############################
##### Standardised data #####
#############################

##### Initialisation, librairies, data #####

library(glmnet);library(InformationValue)

# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

# Seasonal & monthly variables
load('C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Seasonal_monthly_variables.RData')
Season_month_variables_stand <- apply(Data, MARGIN = 2, FUN = scale) # standardisation

Yields_stand <- Season_month_variables_stand[,1]

Season_variables_stand <- cbind(Yields_stand, Season_month_variables_stand[,-1][,(1:(dim(Season_month_variables_stand)[2]-1))%%13 %in% c(11,12,0)])
Month_variables_stand <- cbind(Yields_stand, Season_month_variables_stand[,-1][,as.logical(1-(1:(dim(Season_month_variables_stand)[2]-1))%%13 %in% c(11,12,0))])


# Choose seasonal or monthly analysis

# Variables <- Season_variables_stand
Variables <- Month_variables_stand

#Percentile wanted, in c(0.025, 0.05, 0.1)
percentile <- 0.05
bad_yield_stand_threshold <- quantile(Yields_stand, percentile)


###### Ridge regression #####

set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Variables[training_indices,]
Testing_Data <- Variables[testing_indices,]



T1 <- Sys.time()
CV.Ridge_regression <- cv.glmnet(x = Training_Data[,-1],
                                 y = Training_Data[,"Yields_stand"]>bad_yield_stand_threshold,
                                 family = "binomial", alpha = 0)
T2 <- Sys.time()
difftime(T2,T1)

plot(CV.Ridge_regression, xvar = "lambda", label = TRUE, main="Ridge regression")

lambda_VALS <- c("lambda.min", "lambda.1se")
lambda_val <- lambda_VALS[2]


coeff_ridge <- coef(CV.Ridge_regression, s=lambda_val)

sorted_indices <- sort(abs(as.numeric(coeff_ridge)), decreasing = T, index.return = T)$ix
cbind(coeff_ridge@Dimnames[[1]][sorted_indices], round(coeff_ridge[sorted_indices], digits = 5))


predCV <- predict(CV.Ridge_regression, newx = Testing_Data[,-1], 
                  s = lambda_val,
                  type = "response")



#### Performence of the model ####


#Misclassification error
misClassError(Testing_Data[,1]>bad_yield_stand_threshold, predCV)

#ROC Curve
plotROC(Testing_Data[,1]>bad_yield_stand_threshold, predCV)


# Confusion matrix from package InformationValue
InformationValue::confusionMatrix(actuals = Testing_Data[,1]>bad_yield_stand_threshold,
                                  predictedScores = predCV, threshold = 0.5)


##### Play on the segregation threshold for optimizing the specificity #####
#cannot apply sensitivity function on vector: do it manually :/
Thresh_segreg <- seq(0.5,0.99, by = 0.01)
my_speci <- numeric()
my_sensi <- numeric()
for (index in 1:length(Thresh_segreg)) {
  my_sensi[index] <- InformationValue::sensitivity(actuals = Testing_Data[,1]>bad_yield_stand_threshold,
                                                   predictedScores = predCV, threshold = Thresh_segreg[index])
  my_speci[index] <- InformationValue::specificity(actuals = Testing_Data[,1]>bad_yield_stand_threshold,
                                                   predictedScores = predCV, threshold = Thresh_segreg[index])
  
}#end for index

plot(Thresh_segreg, my_sensi, ylim = c(0,1), type = "l",
     xlab = "segregation threshold",
     ylab = "sensitivity and specifity",
     main = paste("Ridge regression, percentile=", percentile))
text(x=0.55, y=0.9, "sensitivity")
lines(Thresh_segreg, my_speci, col="red")
text(x=0.55, y=0.6, "specificity", col="red")

# Sensitivity and Specificity using package InformationValue
InformationValue::sensitivity(Testing_Data[,1]>bad_yield_stand_threshold, predCV)
InformationValue::specificity(Testing_Data[,1]>bad_yield_stand_threshold, predCV)



# # Get AIC and BIC
# ridge.fit <- glmnet(x = Training_Data[,-1], y = Training_Data[,"Yields_stand"]>bad_yield_stand_threshold,
#                     family = "binomial", alpha = 0, lambda = CV.Ridge_regression$lambda.1se)
# 
# tLL <- ridge.fit$nulldev - deviance(ridge.fit) # tLL is 2 times the likelihood of the model
# k <- ridge.fit$df                              # number of nonzero coefficients
# n <- ridge.fit$nobs                            # number of observations
# 
# 
# AIC <- -tLL+2*k+2*k*(k+1)/(n-k-1)             #AIC corrected = AIC with penalty for number of parameters
# print(paste("AIC", AIC))
# 
# BIC<-log(n)*k - tLL
# print(paste("BIC", BIC))


##### scatterplot yield against 1 meteo variable


#extract the variables with largest coefficient
nb_important_indices <- 10 + 1

important_indices <- sorted_indices[1:nb_important_indices]
important_indices <- important_indices[-which(important_indices==1)]

# print these coefficients
cbind(coeff_ridge@Dimnames[[1]][important_indices], round(coeff_ridge[important_indices], digits = 5))


# scatterplot with yield for these variables

par(mar=c(4,4,1,1))
for (index in 1:((nb_important_indices-1)/2)) {
  plot(Yields_stand, Variables[,coeff_ridge@Dimnames[[1]][important_indices[index]]],
       ylab = coeff_ridge@Dimnames[[1]][important_indices[index]])
  points(Yields_stand[Yields_stand<bad_yield_stand_threshold],
         Variables[,coeff_ridge@Dimnames[[1]][important_indices[index]]][Yields_stand<bad_yield_stand_threshold],
         col="#d95f02")
}





# # Plot proba of bad yield against important variable
# # Mean value of the standardized data = 0, easier to get proba bad yield of the model
# 
# proba_bad_yield <- function(variable_name, stand_var_values,
#                             additionnal_var_name, additionnal_stand_var_value){
#   intercept <- coeff_ridge@x[coeff_ridge@Dimnames[[1]]=="(Intercept)"]
#   coeff <- coeff_ridge@x[coeff_ridge@Dimnames[[1]]==variable_name]
#   coeff_additional <- coeff_ridge@x[coeff_ridge@Dimnames[[1]]==additionnal_var_name]
#   return(1/(1+exp(intercept + coeff*stand_var_values + coeff_additional*additionnal_stand_var_value)))
# }#end for proba_bad_yield function
# 
# 
# #first variable with 3 values for the 3rd FOR 5th PERCENTILE
# range(Season_month_variables_stand[,"dps_July"])
# 
# par(mar=c(4,4,1,1))
# plot(seq(-3,3, by = 0.1),
#      proba_bad_yield(variable_name = "dps_July", stand_var_values = seq(-3,3, by = 0.1),
#                      additionnal_var_name="pr_July",
#                      additionnal_stand_var_value = 0),
#      ylab="proba bad yield", xlab = "Standardised dps_July")
# 
# points(seq(-3,3, by = 0.1),
#        proba_bad_yield(variable_name = "dps_July", stand_var_values = seq(-3,3, by = 0.1),
#                        additionnal_var_name="pr_July",
#                        additionnal_stand_var_value = min(Season_month_variables_stand[,"pr_July"])),
#        col="brown")
# 
# points(seq(-3,3, by = 0.1),
#        proba_bad_yield(variable_name = "dps_July", stand_var_values = seq(-3,3, by = 0.1),
#                        additionnal_var_name="pr_July",
#                        additionnal_stand_var_value = max(Season_month_variables_stand[,"pr_July"])),
#        col="blue")
# 
# legend("topright", title="Standardised pr_July" ,legend=c("0", "min over all years", "max over all years"),
#        col=c("black", "brown", "blue"), pch=1)
# 
# 
# 
# #second variable with 3 values for the 3rd FOR 5th PERCENTILE
# range(Season_month_variables_stand[,"tasmax_July"])
# 
# par(mar=c(4,4,1,1))
# plot(seq(-3,4, by = 0.1),
#      proba_bad_yield(variable_name = "tasmax_July", stand_var_values = seq(-3,4, by = 0.1),
#                      additionnal_var_name="pr_July",
#                      additionnal_stand_var_value = 0),
#      ylab="proba bad yield", xlab = "Standardised tasmax_July")
# 
# points(seq(-3,4, by = 0.1),
#        proba_bad_yield(variable_name = "tasmax_July", stand_var_values = seq(-3,4, by = 0.1),
#                        additionnal_var_name="pr_July",
#                        additionnal_stand_var_value = min(Season_month_variables_stand[,"pr_July"])),
#        col="brown")
# 
# points(seq(-3,4, by = 0.1),
#        proba_bad_yield(variable_name = "tasmax_July", stand_var_values = seq(-3,4, by = 0.1),
#                        additionnal_var_name="pr_July",
#                        additionnal_stand_var_value = max(Season_month_variables_stand[,"pr_July"])),
#        col="blue")
# 
# legend("topleft", title="Standardised pr_July" ,legend=c("0", "min over all years", "max over all years"),
#        col=c("black", "brown", "blue"), pch=1)

















###### Lasso regression #####
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Variables[training_indices,]
Testing_Data <- Variables[testing_indices,]


T1 <- Sys.time()
CV.Lasso_regression <- cv.glmnet(x = Training_Data[,-1],
                                 y = Training_Data[,"Yields_stand"]>bad_yield_stand_threshold,
                                 family = "binomial", alpha = 1, nfolds = 10)
T2 <- Sys.time()
difftime(T2,T1)

plot(CV.Lasso_regression, xvar = "lambda", label = TRUE, main="Lasso regression")

lambda_VALS <- c("lambda.min", "lambda.1se")

lambda_val <- lambda_VALS[2]


coeff_lasso <- coef(CV.Lasso_regression, s=lambda_val)

indices_sorted <- sort(abs(as.numeric(coeff_lasso)), decreasing = T, index.return = T)$ix

cbind(coeff_lasso@Dimnames[[1]][indices_sorted],
      round(coeff_lasso[indices_sorted],5))

predCV_lasso <- predict(CV.Lasso_regression, newx = Testing_Data[,-1], 
                        s = lambda_val,
                        type = "response")



#### Performence of the model ####


#Misclassification error
misClassError(Testing_Data[,1]>bad_yield_stand_threshold, predCV_lasso)

#ROC Curve
plotROC(Testing_Data[,1]>bad_yield_stand_threshold, predCV_lasso)


# Confusion matrix from package InformationValue
InformationValue::confusionMatrix(Testing_Data[,1]>bad_yield_stand_threshold, predCV_lasso)

# Sensitivity and Specificity using package InformationValue
InformationValue::sensitivity(Testing_Data[,1]>bad_yield_stand_threshold, predCV_lasso)
InformationValue::specificity(Testing_Data[,1]>bad_yield_stand_threshold, predCV_lasso)




# # Get AIC and BIC with the selected lambda
# ridge.fit <- glmnet(x = Training_Data[,-1], y = Training_Data[,"Yields_stand"]>bad_yield_stand_threshold,
#                     family = "binomial", alpha = 0, lambda = CV.Lasso_regression$lambda.1se)
# 
# tLL <- ridge.fit$nulldev - deviance(ridge.fit) # tLL is 2 times the likelihood of the model
# k <- ridge.fit$df                              # number of nonzero coefficients
# n <- ridge.fit$nobs                            # number of observations
# 
# AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)             #AIC corrected = AIC with penalty for number of parameters
# print(paste("AIC corrected", AIC))
# 
# BIC <- log(n)*k - tLL
# print(paste("BIC", BIC))




##### scatterplot yield against 1 meteo variable


#extract the largest coefficient
nb_important_indices <- 5 + 1

important_indices <- indices_sorted[1:nb_important_indices]
important_indices <- important_indices[-which(important_indices==1)]



cbind(coeff_lasso@Dimnames[[1]][important_indices], round(coeff_lasso[important_indices], digits = 5))
# scatterplot with yield for these variables

par(mar=c(4,4,1,1))
for (index in 1:(nb_important_indices-1)) {
  plot(Yields_stand, Variables[,coeff_lasso@Dimnames[[1]][important_indices[index]]],
       ylab = coeff_lasso@Dimnames[[1]][important_indices[index]])
  points(Yields_stand[Yields_stand<bad_yield_stand_threshold],
         Variables[,coeff_lasso@Dimnames[[1]][important_indices[index]]][Yields_stand<bad_yield_stand_threshold],
         col="#d95f02")
}