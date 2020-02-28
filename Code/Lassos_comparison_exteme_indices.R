# Lasso Model comparison ###################
# Authors: Pauline, Johannes
# This file is meant to compare the performance of LASSO  from glinternet and from glmnet (monthly meteovar + extreme indices)
# It is structured in the following way:
# a) recreate the training and testing dataset (should be the same as the one used to fit our model, thanks to set.seed(1994))
# b) load model output for each model
# c) calculate performance metrics for each model
# d) use these metrics to create model comparison plots
######################################


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))


#threshold for bad yields in c(0.025,0.05,0.1)
threshold <- 0.05

#which cutoff level?
segreg_th <- 0.5

##### Initialisation, librairies, data #####

library(ncdf4);library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr);library(tictoc);library(ggplot2)



##### Load standardized Data and buid training and testing data list #####

#Pauline
load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/extremeindices_and_monthlymeteovar_rescaled.Rdata")

#Johannes
load("D:/user/vogelj/Data/Group project Como/extremeindices_and_monthlymeteovar_rescaled.Rdata")



yield_3dim <- array(Data_standardized$yield,dim=c(965,1,1600))
dtr_3dim <- array(Data_standardized$dtr,dim=c(965,1,1600))
frs_3dim <- array(Data_standardized$frs,dim=c(965,1,1600))
txx_3dim <- array(Data_standardized$txx,dim=c(965,1,1600))
tnn_3dim <- array(Data_standardized$tnn,dim=c(965,1,1600))
rx5_3dim <- array(Data_standardized$rx5,dim=c(965,1,1600))
tx90p_3dim <- array(Data_standardized$tx90p,dim=c(965,1,1600))
tn10p_3dim <- array(Data_standardized$tn10p,dim=c(965,1,1600))

Model_data <- abind(yield_3dim,dtr_3dim,frs_3dim,txx_3dim,tnn_3dim,rx5_3dim,tx90p_3dim,tn10p_3dim
                    ,Data_standardized$tasmax,Data_standardized$vpd,Data_standardized$pr,along=2)
colnames(Model_data) <- c("Yield", "dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p",
                          "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                          "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                          "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
                          "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                          "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                          "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2",
                          "pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                          "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                          "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")

pix_num <- dim(Model_data)[1]
Yield <- Data_standardized$yield
low_yield <- apply(Yield, MARGIN = 1, FUN=quantile, probs=threshold, na.rm=T)
cy <- t(sapply(1:pix_num,function(x){ifelse(Yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield

cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))

Model_data[,1,] <-cy_reshaped


# Exclude NA variable columns
na_col <- matrix(data=NA,nrow=pix_num,ncol=dim(Model_data)[2])
for (j in 1:pix_num){
  for (i in 1:dim(Model_data)[2]){
    na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)

na_time <- vector("list",length=pix_num) # for each pixel, the positions of NAs over time
for (j in 1:pix_num){
  na_time[[j]] <- which(is.na(Model_data[j,1,])) # locations of years with NA values
}


## Split data into training and testing data set
vec <- 1:1600

years_with_na <- vector("logical",length=pix_num)
for (i in 1:pix_num){
  years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
}

training_indices <- vector("list",length=pix_num)
testing_indices <- vector("list",length=pix_num)
set.seed(1994)
for (x in 1:pix_num) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((1600-length(na_time[[x]]))*0.6)))
    testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
  } else {
    training_indices[[x]] <- sort(sample(1:1600, size = floor(1600*0.6)))
    testing_indices[[x]] <- (1:1600)[-training_indices[[x]]]    
  }
}


Training_Data <- lapply(1:pix_num,function(x){Model_data[x,,training_indices[[x]]]})
Testing_Data <- lapply(1:pix_num,function(x){Model_data[x,,testing_indices[[x]]]})

pix_in <- 1:pix_num


x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand


var_num <- apply(non_na_col,1,sum)
numLevels_list <- sapply(1:pix_num, function(x){ rep(1,times=var_num[x])})
for (i in 1:pix_num){
  names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
}


##### Usefull functions ####
number_coeff_kept <- function(coeff_list){#give number of coeff !=0
  return(length(which(abs(coeff_list[-1,])>0)))
}

extreme_in_coeff <- function(coeff_list){ #function to check how many extreme indeices are kept as predictors
  extreme_indices <- c("dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p")
  if(max(abs(coeff_list[extreme_indices,]))==0){
    return(0)
  } else {
    return(length(which(abs(coeff_list[extreme_indices,])>0)))
  }
}#end func extreme_in_coeff






##### glmnet model performance #####

# On the Drive you can find data in:
# Models/LASSO-Ridge regression/regression_results_Global_wo_interactions/Lasso_lambda1se_month_xtrm_Lasso_threshbadyield005.RData

# Pauline's Laptop
load(paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambda1se_month_xtrm_Lasso_threshbadyield005.Rdata"))

# extract performance measures

lambda1se_glmnet <- numeric()
for (pix in 1:pix_num) {
  lambda1se_glmnet[pix] <- lasso_model_lambda1se[[pix]]$lambda
}#end for pix

coeff_glmnet  <-list()
speci_glmnet <- rep(NA, pix_num)
sensi_glmnet <- rep(NA, pix_num)
csi_glmnet <- rep(NA, pix_num)

for (pixel in 1:pix_num) {
    
    coeff_glmnet[[pixel]] <- coefficients(lasso_model_lambda1se[[pixel]])
    
    mypred <- predict(lasso_model_lambda1se[[pixel]], as.matrix(x1_test_list[[pixel]]),type="response")
    
    fitted_bad_yield <- ifelse(mypred > segreg_th,1,0)
    
    speci_glmnet[pixel] <- specificity(actuals = as.matrix(y1_test_list[[pixel]]),
                                predictedScores = fitted_bad_yield,
                                threshold = segreg_th)
    sensi_glmnet[pixel] <- sensitivity(actuals = as.matrix(y1_test_list[[pixel]]),
                                predictedScores = fitted_bad_yield,
                                threshold = segreg_th)
    
    con_tab <- confusionMatrix(actuals = as.matrix(y1_test_list[[pixel]]),
                               predictedScores = fitted_bad_yield,
                               threshold = segreg_th)
    csi_glmnet[pixel] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
    if(is.na(con_tab["0","0"])){
      csi_glmnet[pixel] <- 0
    }
}#end pixel



nb_extr_kept_glmnet <- numeric()
nb_coeff_kept_glmnet <- numeric()

for (pixel in 1:pix_num) {
    nb_extr_kept_glmnet[pixel] <- extreme_in_coeff(coeff_glmnet[[pixel]])
    nb_coeff_kept_glmnet[pixel] <- number_coeff_kept(coeff_glmnet[[pixel]])
}#end pixel




##### glinternet model performance #####

# On the Drive you can find data in:
# Models/Lasso(glinternet)/Lasso_without_interactions/cv_fit_monthly_without_int_incl_ext.Rdata

load(paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/OtherModels/cv_fit_monthly_without_int_incl_ext.Rdata"))