# Author: Pauline

rm(list=ls(all=TRUE))
library(glmnet)

# Input features should be as matrix and not data.frame
# Syntax is glmnet(x = our input matrix, y = our response,
#                  family = the distribution to be chosen in
#                  c("gaussian","binomial","poisson","multinomial","cox","mgaussian")
#                  , alpha = 0)
# (aplha = 0 for ridge, 1 for LASSO)


load('C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Data.RData')
load('./Processed data/Seasonal_monthly_variables.RData') # Matrix with all pretictands

Yields <- Data[,1]

yield_5pctl <- quantile(Yields, 0.05)



###### Ridge regression #####

set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Data[training_indices,]
Testing_Data <- Data[testing_indices,]



CV.Ridge_regression <- cv.glmnet(x = Training_Data[,-1],
                                y = Training_Data[,"Yields"]>yield_5pctl,
                                family = "binomial", alpha = 0)

plot(CV.Ridge_regression, xvar = "lambda", label = TRUE, main="Ridge regression")

lambda_VALS <- c("lambda.min", "lambda.1se")
lambda_val <- lambda_VALS[2]


coeff_ridge <- coef(CV.Ridge_regression, s=lambda_val)
coeff_ridge@Dimnames[[1]][which(coeff_ridge>0)]
coeff_ridge@Dimnames[[1]][which(coeff_ridge<0)]



cbind(c(coeff_ridge@Dimnames[[1]][which(coeff_ridge>0)],
        coeff_ridge@Dimnames[[1]][which(coeff_ridge<0)]),
      round(c(coeff_ridge[which(coeff_ridge>0)],
              coeff_ridge[which(coeff_ridge<0)]),5))



predCV <- predict(CV.Ridge_regression, newx = Testing_Data[,-1], 
                  s = lambda_val,
                  type = "response")

misClassError(Testing_Data[,1]>yield_5pctl, predCV)

plotROC(Testing_Data[,1]>yield_5pctl, predCV)

# Ridge:
# lambda=       min     1se
# misclassEr    0.025   0.0266
# AUROC         0.9905  0.9909

# text(paste("Ridge regression, lambda", lambda_val), side=1, line=0, at=0.8)




###### Lasso regression #####
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Data[training_indices,]
Testing_Data <- Data[testing_indices,]


T1 <- Sys.time()
CV.Lasso_regression <- cv.glmnet(x = Training_Data[,-1],
                                 y = Training_Data[,"Yields"]>yield_5pctl,
                                 family = "binomial", alpha = 1, nfolds = 10)
T2 <- Sys.time()
difftime(T2,T1)

plot(CV.Lasso_regression, xvar = "lambda", label = TRUE, main="Lasso regression")

lambda_VALS <- c("lambda.min", "lambda.1se")

lambda_val <- lambda_VALS[1]


coeff_lasso <- coef(CV.Lasso_regression, s=lambda_val)
coeff_lasso@Dimnames[[1]][which(coeff_lasso>0)]
coeff_lasso@Dimnames[[1]][which(coeff_lasso<0)]

cbind(c(coeff_lasso@Dimnames[[1]][which(coeff_lasso>0)],
        coeff_lasso@Dimnames[[1]][which(coeff_lasso<0)]),
      round(c(coeff_lasso[which(coeff_lasso>0)],
        coeff_lasso[which(coeff_lasso<0)]),5))[-1,]


library(InformationValue)
predCV_lasso <- predict(CV.Lasso_regression, newx = Testing_Data[,-1], 
                        s = lambda_val,
                        type = "response")

misClassError(Testing_Data[,1]>yield_5pctl, predCV_lasso)

plotROC(Testing_Data[,1]>yield_5pctl, predCV_lasso)




# LASSO:
# lambda=       min     1se
# misclassEr    0.025   0.0281
# AUROC         0.9923  0.9931



# Ridge:
# lambda=       min     1se
# misclassEr    0.025   0.0266
# AUROC         0.9905  0.9909

library(glinternet)


# text(paste("Lasso regression, lambda", lambda_val), side=1, line=0, at=0.8)
