# GLM with cross-validation and lasso regression 


# Data loading and processing ####
# setwd("C:/Users/vogel/boxup/Promotion/Damocles COST Action/Training school/Group1_Project")
# load('./Processed data/Data.RData') # Matrix with all pretictands
load('./Processed data/Seasonal_monthly_variables.RData') # Matrix with all pretictands
load('./Code/Workspaces/Binomial logit regression model.RData')
# load('./Code/Workspaces/GLM_Lasso.RData') 

library(glmnet)
library(leaps)
library(bestglm)
library(car)
library(glinternet)
library(dplyr)
library(InformationValue)
library(ROCR)
library(tictoc)

# Season_month_variables_stand <- apply(Data,2,scale) # z-score standardisation

# Determine low yield based on percentile threshold
thresholds <- c(0.025,0.05,0.1)
low_yield <- quantile(cy_gsl_FR[,1],thresholds[3])
cy <- ifelse(cy_gsl_FR[,1]<low_yield,0,1)
Data[,1] <- cy # Replace actual crop yield by binary info on fail/success

# Split data into training and testing data set
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Data[training_indices,]
Testing_Data <- Data[testing_indices,]

Dpd <- sapply(1:13, function(x){Data[,x+1]-Data[,x+53]}) # calculate dew point depression
colnames(Dpd) <- c("Dpd_November","Dpd_December","Dpd_January","Dpd_February","Dpd_March","Dpd_April","Dpd_May","Dpd_June","Dpd_July","Dpd_August","Dpd_Winter","Dpd_Spring","Dpd_Summer")
# # Xy <- cbind(data.frame(Data[,2:16]),Data[,1]) # only dew point 
# # Xy <- cbind(Dpd[,11:13],data.frame(Data[,c(25:27,77:79,1)])) # dps, prcp, tasmin
# # Xy <- cbind(Dpd[,8:10],data.frame(Data[,c(22:24,1)])) # dps, prcp
# Xy <- data.frame(Data[,c(12:14,25:27,64:66,77:79,1)]) # dps, prcp, tmin, tmax for all seasons
vec <- c(11:13,24:26,37:39,50:52,63:65,76:78) # seasonal data
vec2 <- vec[-c(7:12)] # exclude wind and srds


AllData <- as.data.frame(Data)
AllTraining_Data <- as.data.frame(Training_Data)
AllTesting_Data <- as.data.frame(Testing_Data)

Xy <- data.frame(AllTraining_Data[,c(vec2+1,1)]) # Subset
Xy_test <- data.frame(AllTesting_Data[,c(vec2+1,1)])



# GLM using cross-validation for goodness-of-fit ####
#####################################################
tic()
the_best_glm <- bestglm(Xy = Xy, IC = "CV", CVArgs = list(Method = "HTF", K = 10, REP = 100), 
        family=binomial) # you cannot use more than 15 variables
toc()
# Normally K= 10 or K= 5 are used
# Which measure of fit is used here in the cross-validation? https://en.wikipedia.org/wiki/Cross-validation_(statistics)#Measures_of_fit
('message: use larger number of Rep: 100 or 1000')

the_best_glm_AIC <- bestglm(Xy = Xy, IC = "AIC", 
                        family=binomial)
# the_best_glm
print.bestglm(the_best_glm)
print.bestglm(the_best_glm_AIC)
# summary(the_best_glm)
summary.bestglm(the_best_glm)
summary.bestglm(the_best_glm_AIC)
the_best_glm$BestModel$aic # AIC value
the_best_glm_AIC$BestModel$aic # AIC value


#The best models for each subset size
the_best_glm$Subsets
the_best_glm_AIC$Subsets

# mypred <- predict(the_best_glm,Xy_test,type="response")
# fitted.results_bestglm <- ifelse(mypred > 0.5,1,0)

# misClassError(AllTesting_Data[,1],fitted.results_bestglm)



# Interaction terms ####

#' calculate interactions manually for 4 best variables and put it in the model
#' (because you cannot include interactions directly in the model)

subsetXy <- Xy[,c(3,5,10,12,13)]
interactions <- cbind(subsetXy[,1]*subsetXy[,2],subsetXy[,1]*subsetXy[,3],subsetXy[,1]*subsetXy[,4],subsetXy[,2]*subsetXy[,3]
                      ,subsetXy[,2]*subsetXy[,4],subsetXy[,3]*subsetXy[,4])
subsetXy_interactions <- cbind(subsetXy[,c(1:4)],interactions,subsetXy[,5])
colnames(subsetXy_interactions) <- c("dps_Summer","pr_Spring","tasmin_Winter","tasmin_Summer","dps_Summer_pr_Spring","dps_Summer_tasmin_Winter","dps_Summer_tasmin_Summer","pr_Spring_tasmin_Winter","pr_Spring_tasmin_Summer","tasmin_Winter_tasmin_Summer","Yields")


the_best_glm_intact <- bestglm(Xy = subsetXy_interactions, IC = "CV", CVArgs = list(Method = "HTF", K = 10, REP = 1), 
                        family=binomial) # you cannot use more than 15 variables
summary(the_best_glm_intact)
summary.bestglm(the_best_glm_intact)
print.bestglm(the_best_glm_intact)

# does not work
# vif(the_best_glm_intact) # multicollinearity (https://rpubs.com/ranvirkumarsah/LR)
vifx(the_best_glm_intact)					 
# 1/vif(the_best_glm_intact)
# mean(vif(the_best_glm_intact))



# Lasso regression with interactions terms: glinternet-package:  ####
#####################################################################


# mynum <- 79 # 10 min for about 29 variables
# X1 <- AllTraining_Data[,vec+1]
X1_train <- AllTraining_Data[,2:79]
y1_train <- AllTraining_Data[,1]
# subset_Testing_Data <-  AllTesting_Data[,vec+1]
x1_test <-  AllTesting_Data[,2:79]
y1_test <- AllTesting_Data[,1]

# numLevels <- AllTraining_Data[,vec+1] %>% sapply(nlevels)
numLevels <- AllTraining_Data[,2:79] %>% sapply(nlevels)
numLevels[numLevels==0] <- 1 # set to 1 for continuous variables

# Fit model
cv_fit <- glinternet.cv(X1_train, y1_train, numLevels,family = "binomial")
plot(cv_fit)

i_1Std <- which(cv_fit$lambdaHat1Std == cv_fit$lambda) # the preferential lambda (tuning parameter)

coefs <- coef(cv_fit$glinternetFit)[[i_1Std]] 

coefs$mainEffects # model part without interactions
names(numLevels)[coefs$mainEffects$cont] # Main effect variables (without interactions)

coefs$interactions # model part with interactions pairs
names(numLevels)[coefs$interactions$contcont] # Main effect variables (with interactions)

# Assessing performance ####
sqrt(cv_fit$cvErr[[i_1Std]]) # root mean squared error (RMSE) on validation data


# predict.glinternet.cv(cv_fit,subset_Testing_Data,type="response") # does not work
# glinternet::predict.glinternet.cv(cv_fit,AllTesting_Data,type="response") # does not work
# predict.glinternet(cv_fit,AllTesting_Data,type="response") # does not work
mypred <- predict(cv_fit,x1_test,type="response")
fitted.results_bestglm <- ifelse(mypred > 0.5,1,0)


misClassError(y1_test,fitted.results_bestglm)

# AIC(cv_fit$) # does not work, needs log-likelihood


# ROC ####

pr <- prediction(mypred, y1_test)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
plotROC(actuals=y1_test,predictedScores=fitted.results_bestglm)

auc <- performance(pr, measure = "auc")
auc@y.values[[1]]





# Linear model for comparison ####
# X2 <- model.matrix(Price ~ . - 1, df) # not needed (no categorical predictors)
X1m <- as.matrix(X1_train)
cv_glmnet <- cv.glmnet(X1m,y1_train)
# cv_glmnet <- cv.glmnet(X,y)
sqrt(min(cv_glmnet$cvm))



save.image('./Code/Workspaces/GLM_Lasso.RData')

							

