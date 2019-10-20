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

# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

library(glmnet);library(InformationValue)

# Seasonal & monthly variables
load('C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Seasonal_monthly_variables.RData')
Season_month_variables_stand <- apply(Data, MARGIN = 2, FUN = scale) # standardisation

Yields_stand <- Season_month_variables_stand[,1]


#Percentile wanted
percentile <- 0.1
bad_yield_stand_threshold <- quantile(Yields_stand, percentile)


###### Ridge regression #####

set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Season_month_variables_stand[training_indices,]
Testing_Data <- Season_month_variables_stand[testing_indices,]



T1 <- Sys.time()
CV.Ridge_regression <- cv.glmnet(x = Training_Data[,-1],
                                 y = Training_Data[,"Yields"]>bad_yield_stand_threshold,
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

#Misclassification error
misClassError(Testing_Data[,1]>bad_yield_stand_threshold, predCV)

#Lambda kept
CV.Ridge_regression$lambda.1se

#ROC Curve
plotROC(Testing_Data[,1]>bad_yield_stand_threshold, predCV)

# Get AIC and BIC with the selected lambda
ridge.fit <- glmnet(x = Training_Data[,-1], y = Training_Data[,"Yields"]>bad_yield_stand_threshold,
                    family = "binomial", alpha = 0, lambda = CV.Ridge_regression$lambda.1se)

tLL <- ridge.fit$nulldev - deviance(ridge.fit) # tLL is 2 times the likelihood of the model
k <- ridge.fit$df                              # number of nonzero coefficients
n <- ridge.fit$nobs                            # number of observations


AIC <- -tLL+2*k+2*k*(k+1)/(n-k-1)             #AIC corrected = AIC with penalty for number of parameters
print(paste("AIC", AIC))

BIC<-log(n)*k - tLL
print(paste("BIC", BIC))


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
  plot(Yields_stand, Season_month_variables_stand[,coeff_ridge@Dimnames[[1]][important_indices[index]]],
       ylab = coeff_ridge@Dimnames[[1]][important_indices[index]])
  points(Yields_stand[Yields_stand<bad_yield_stand_threshold],
         Season_month_variables_stand[,coeff_ridge@Dimnames[[1]][important_indices[index]]][Yields_stand<bad_yield_stand_threshold],
         col="#d95f02")
}







###### Lasso regression #####
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Season_month_variables_stand[training_indices,]
Testing_Data <- Season_month_variables_stand[testing_indices,]


T1 <- Sys.time()
CV.Lasso_regression <- cv.glmnet(x = Training_Data[,-1],
                                 y = Training_Data[,"Yields"]>bad_yield_stand_threshold,
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

#Misclassification error
misClassError(Testing_Data[,1]>bad_yield_stand_threshold, predCV_lasso)

#Lambda kept
CV.Lasso_regression$lambda.1se

#ROC Curve
plotROC(Testing_Data[,1]>bad_yield_stand_threshold, predCV_lasso)

# Get AIC and BIC with the selected lambda
ridge.fit <- glmnet(x = Training_Data[,-1], y = Training_Data[,"Yields"]>bad_yield_stand_threshold,
                    family = "binomial", alpha = 0, lambda = CV.Lasso_regression$lambda.1se)

tLL <- ridge.fit$nulldev - deviance(ridge.fit) # tLL is 2 times the likelihood of the model
k <- ridge.fit$df                              # number of nonzero coefficients
n <- ridge.fit$nobs                            # number of observations

AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)             #AIC corrected = AIC with penalty for number of parameters
print(paste("AIC corrected", AIC))

BIC <- log(n)*k - tLL
print(paste("BIC", BIC))




##### scatterplot yield against 1 meteo variable


#extract the largest coefficient
nb_important_indices <- 5 + 1

important_indices <- indices_sorted[1:nb_important_indices]
important_indices <- important_indices[-which(important_indices==1)]



cbind(coeff_lasso@Dimnames[[1]][important_indices], round(coeff_lasso[important_indices], digits = 5))
# scatterplot with yield for these variables

par(mar=c(4,4,1,1))
for (index in 1:(nb_important_indices-1)) {
  plot(Yields_stand, Season_month_variables_stand[,coeff_lasso@Dimnames[[1]][important_indices[index]]],
       ylab = coeff_lasso@Dimnames[[1]][important_indices[index]])
  points(Yields_stand[Yields_stand<bad_yield_stand_threshold],
         Season_month_variables_stand[,coeff_lasso@Dimnames[[1]][important_indices[index]]][Yields_stand<bad_yield_stand_threshold],
         col="#d95f02")
}













###############################
##### Unstandardised data #####
###############################

##### Initialisation, librairies and data #####

# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

library(glmnet);library(InformationValue)

# Seasonal & monthly variables
load('C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Seasonal_monthly_variables.RData')


Yields <- Data[,1]

#Percentile wanted
percentile <- 0.025
bad_yield_threshold <- quantile(Yields, percentile)


###### Ridge regression #####

set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Data[training_indices,]
Testing_Data <- Data[testing_indices,]



T1 <- Sys.time()
CV.Ridge_regression <- cv.glmnet(x = Training_Data[,-1],
                                 y = Training_Data[,"Yields"]>bad_yield_threshold,
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
# Misclassification error
misClassError(Testing_Data[,1]>bad_yield_threshold, predCV)

# Lambda kept
CV.Ridge_regression$lambda.1se

# ROC Curve
plotROC(Testing_Data[,1]>bad_yield_threshold, predCV)



##### scatterplot yield against 1 meteo variable


# extract the largest coefficient
nb_important_indices <- 5 + 1

important_indices <- sorted_indices[1:nb_important_indices]
important_indices <- important_indices[-which(important_indices==1)]

# print them
cbind(coeff_ridge@Dimnames[[1]][important_indices], round(coeff_ridge[important_indices], digits = 5))

# scatterplot with yield for these variables

par(mar=c(4,4,1,1))
for (index in 1:(nb_important_indices-1)) {
  plot(Yields, Data[,coeff_ridge@Dimnames[[1]][important_indices[index]]],
       ylab = coeff_ridge@Dimnames[[1]][important_indices[index]])
  points(Yields[Yields<bad_yield_threshold],
         Data[,coeff_ridge@Dimnames[[1]][important_indices[index]]][Yields<bad_yield_threshold],col="#d95f02")
}







###### Lasso regression #####
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Data[training_indices,]
Testing_Data <- Data[testing_indices,]


T1 <- Sys.time()
CV.Lasso_regression <- cv.glmnet(x = Training_Data[,-1],
                                 y = Training_Data[,"Yields"]>bad_yield_threshold,
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

#Misclassification error
misClassError(Testing_Data[,1]>bad_yield_threshold, predCV_lasso)

#Lambda kept
CV.Lasso_regression$lambda.1se

#ROC Curve
plotROC(Testing_Data[,1]>bad_yield_threshold, predCV_lasso)




##### scatterplot yield against 1 meteo variable


# extract the largest coefficient
nb_important_indices <- 5 + 1

important_indices <- indices_sorted[1:nb_important_indices]
important_indices <- important_indices[-which(important_indices==1)]

# print them
cbind(coeff_lasso@Dimnames[[1]][important_indices], round(coeff_lasso[important_indices], digits = 5))

# scatterplot with yield for these variables

par(mar=c(4,4,1,1))
for (index in 1:(nb_important_indices-1)) {
  plot(Yields, Data[,coeff_lasso@Dimnames[[1]][important_indices[index]]],
       ylab = coeff_lasso@Dimnames[[1]][important_indices[index]])
  points(Yields[Yields<bad_yield_threshold],
         Data[,coeff_lasso@Dimnames[[1]][important_indices[index]]][Yields<bad_yield_threshold],col="#d95f02")
}
