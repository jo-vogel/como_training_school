# GLM with cross-validation and lasso regression 


# Data loading and processing ####
# setwd("D:/Promotion/Group project")
setwd("D:/user/vogelj/Group_project")
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
Data_stand <- apply(Data,2,scale) # z-score standardisation


# Determine low yield based on percentile threshold
thresholds <- c(0.025,0.05,0.1)
low_yield <- quantile(cy_gsl_FR[,1],thresholds[3])
cy <- ifelse(cy_gsl_FR[,1]<low_yield,0,1)
Data_stand[,1] <- cy # Replace actual crop yield by binary info on fail/success


Dpd <- sapply(1:13, function(x){Data_stand[,x+1]-Data_stand[,x+53]}) # calculate dew point depression
colnames(Dpd) <- c("Dpd_November","Dpd_December","Dpd_January","Dpd_February","Dpd_March","Dpd_April","Dpd_May","Dpd_June","Dpd_July","Dpd_August","Dpd_Winter","Dpd_Spring","Dpd_Summer")
# # Xy <- cbind(data.frame(Data[,2:16]),Data[,1]) # only dew point 
# # Xy <- cbind(Dpd[,11:13],data.frame(Data[,c(25:27,77:79,1)])) # dps, prcp, tasmin
# # Xy <- cbind(Dpd[,8:10],data.frame(Data[,c(22:24,1)])) # dps, prcp
# Xy <- data.frame(Data[,c(12:14,25:27,64:66,77:79,1)]) # dps, prcp, tmin, tmax for all seasons
vec <- c(11:13,24:26,37:39,50:52,63:65,76:78)+1 # seasonal data
vec2 <- vec[-c(7:12)] # exclude wind and srds
Data_seas <- cbind(Data_stand[,c(1,vec)],Dpd[,11:13])
Data_month <- cbind(Data_stand[,-vec],Dpd[,1:10])


Data_stand <- Data_seas[,c(1,5:7,14:22)] # seasonal data (min/max temp, prec., dpd)
# Data_stand <- Data_month[,c(1,12:21,41:71)] # monthly data (min/max temp, prec., dpd)


# Split data into training and testing data set
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Data_stand[training_indices,]
Testing_Data <- Data_stand[testing_indices,]

AllData <- as.data.frame(Data_stand)
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
X1_train <- AllTraining_Data[,2:length(AllTraining_Data)] # predictors
y1_train <- AllTraining_Data[,1] # predictand
# subset_Testing_Data <-  AllTesting_Data[,vec+1]
x1_test <-  AllTesting_Data[,2:length(AllTesting_Data)] # predictors
y1_test <- AllTesting_Data[,1] # predictand

# numLevels <- AllTraining_Data[,vec+1] %>% sapply(nlevels)
numLevels <- AllTraining_Data[,2:length(AllTraining_Data)] %>% sapply(nlevels)
numLevels[numLevels==0] <- 1 # set to 1 for continuous variables

# Fit model
cv_fit <- glinternet.cv(X1_train, y1_train, numLevels,family = "binomial")
plot(cv_fit)

cv_fit$lambdaHat1Std
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

# Confusion matrix ####
# a) Confusion matrix manually calculated
obs_pred <- cbind(y1_test,fitted.results_bestglm)
tp <- sum(rowSums(obs_pred)==2)
tn <- sum(rowSums(obs_pred)==0)
fp <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
fn <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab1 <- matrix(c(tp,fn,fp,tn),nrow=2,ncol=2)
con_tab1b <- con_tab1
colnames(con_tab1) <- c('Actual TRUE','Actual FALSE');rownames(con_tab1) <- c('Predicted TRUE','Predicted FALSE')
con_tab1[,1] <- con_tab1[,1]/ sum(y1_test==1)
con_tab1[,2] <-con_tab1[,2]/sum(y1_test==0)
# b) Confusion matrix from package InformationValue
con_tab <- InformationValue::confusionMatrix(y1_test,fitted.results_bestglm)

# Sensitivity and Specificity
# a) Manually calculated
spec <- tn/(tn+fp) 
sens <- tp/(tp+fn) 
# b) Using package caret
caret::sensitivity(data=as.factor(fitted.results_bestglm),reference=as.factor(y1_test),positive="1",negative="0")
caret::specificity(data=as.factor(fitted.results_bestglm),reference=as.factor(y1_test),positive="1",negative="0")
# c) using package InformationValue
InformationValue::sensitivity(y1_test,fitted.results_bestglm)
InformationValue::specificity(y1_test,fitted.results_bestglm)


# ROC ####
pr <- prediction(mypred, y1_test)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
plotROC(actuals=y1_test,predictedScores=fitted.results_bestglm)
auc2 <- auc(y1_test,fitted.results_bestglm)
auc <- performance(pr, measure = "auc")
auc@y.values[[1]]

# Calculate Sensitivity and Specificity using performance command
spec2 <- performance(pr, measure = "spec") 
# spec2b <- performance(pr, measure = "tnr") # equivalent
plot(spec2) # you can see that at cutoff 0.5 it is equal to my result, however I don't know how to extract the value
# sens2 <- performance(pr, measure = "sens") # does now work
sens2 <- performance(pr, measure = "tpr")
sens_spec <- performance(pr, measure="sens", x.measure="spec")
plot(sens_spec) # inverted AUC



# AIC(cv_fit$) # does not work, needs log-likelihood


# Linear model for comparison ####
# # X2 <- model.matrix(Price ~ . - 1, df) # not needed (no categorical predictors)
# X1m <- as.matrix(X1_train)
# cv_glmnet <- cv.glmnet(X1m,y1_train)
# # cv_glmnet <- cv.glmnet(X,y)
# sqrt(min(cv_glmnet$cvm))



# save.image('./Code/Workspaces/GLM_Lasso.RData')

							

