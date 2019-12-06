# Model comparison ###################
# Authors: Cristina, Johannes, Pauline
# This file is meant to compare the performance of all the models
# It is structured in the following way:
# a) recreate the training and testing dataset (should be the same as the one used to fit our model, thanks to set.seed(2019))
# b) load model output for each model
# c) calculate performance metrics for each model
# d) use these metrics to create model comparison plots
######################################


# Training and testing data ####
################################

##### Initialisation, librairies, data

library(ncdf4);library(RColorBrewer);library(viridis)
library(maps);library(mapdata);library(ggplot2)
library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr)
library(foreach);library(doParallel)
library(tictoc)



# recreate the training and testing dataset ####
################################################

source('./Code/Lasso_interact_global_preparation.R') # load necessary files
##### Standardised data
# threshold <- 0.05 #threshold for bad yield
# segreg_th <- 0.5 #If the fitted proba from the model is below this threshold then it's a bad year

# # Get the data
# 
# # Path for Pauline
# path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
# 
# nh_files <- list.files(path=path_to_NH_files,pattern="NH_yield*") # all files from northern hemisphere
# nh_data <- lapply(1:length(nh_files),
#                   FUN = function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
# yield <- ncvar_get(nh_data[[1]],"yield")
# tasmax <- ncvar_get(nh_data[[1]],"tasmax")
# vpd <- ncvar_get(nh_data[[1]],"vpd")
# pr <- ncvar_get(nh_data[[1]],"pr")
# lat_subset <- ncvar_get(nh_data[[1]],"lat")
# lon_subset <- ncvar_get(nh_data[[1]],"lon")
# yield_stand <- ncvar_get(nh_data[[2]],"yield")
# tasmax_stand <- ncvar_get(nh_data[[2]],"tasmax")
# vpd_stand <- ncvar_get(nh_data[[2]],"vpd")
# pr_stand <- ncvar_get(nh_data[[2]],"pr")
# lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})
# coord_subset <- cbind(lon_subset,lat_subset)
# 
# 
# # load all coordinates of northern hemisphere
# nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
# nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
# lat_all <- ncvar_get(nh_data[[1]],"lat")
# lon_all <- ncvar_get(nh_data[[1]],"lon")
# lati_all <- rep(lat_all,each=length(lon_all))
# long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
# coord_all <- cbind(long_all,lati_all)
# 
# lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})
# 
# 
# # Process data 
# 
# yields_3dim <- array(yield,dim=c(965,1,1600));yields_stand_3dim <- array(yield,dim=c(965,1,1600))
# Model_data <- abind(yields_3dim,tasmax,vpd,pr,along=2)
# Model_data_stand <- abind(yields_stand_3dim,tasmax_stand,vpd_stand,pr_stand,along=2)
# 
# 
# pix_num <- dim(Model_data)[1]
# low_yield <- sapply(1:pix_num,function(x) {quantile(yield[x,],threshold,na.rm=T)})
# cy <- t(sapply(1:pix_num,function(x){ifelse(yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield
# 
# cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))
# Model_data[,1,] <-cy_reshaped
# Model_data_stand[,1,] <-cy_reshaped
# 
# columnnames <- c("Yield",
#                  "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
#                  "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
#                  "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
#                  "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
#                  "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
#                  "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2",
#                  "pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
#                  "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
#                  "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")
# colnames(Model_data) <- columnnames
# colnames(Model_data_stand) <- columnnames
# 
# 
# 
# # Exclude NA variable columns
# na_col <- matrix(data=NA,nrow=pix_num,ncol=52)
# for (j in 1:pix_num){
#   for (i in 1:52){
#     na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
#   }
# }
# non_na_col <- !na_col # columns without NAs
# non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)
# 
# 
# na_time <- vector("list",length=pix_num) # for each pixel, the positions of NAs over time
# for (j in 1:pix_num){
#   na_time[[j]] <- which(is.na(Model_data[j,1,])) # locations of years with NA values
# }
# 
# 
# 
# # Split data into training and testing data set
# vec <- 1:1600
# 
# years_with_na <- vector("logical",length=pix_num)
# for (i in 1:pix_num){
#   years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
# }
# 
# training_indices <- vector("list",length=pix_num)
# testing_indices <- vector("list",length=pix_num)
# set.seed(1994)
# for (x in 1:pix_num) {
#   if (years_with_na[x]) {
#     training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((1600-length(na_time[[x]]))*0.6)))
#     testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
#   } else {
#     training_indices[[x]] <- sort(sample(1:1600, size = floor(1600*0.6)))
#     testing_indices[[x]] <- (1:1600)[-training_indices[[x]]]    
#   }
# }
# 
# 
# Training_Data <- lapply(1:pix_num,function(x){Model_data_stand[x,,training_indices[[x]]]})
# Testing_Data <- lapply(1:pix_num,function(x){Model_data_stand[x,,testing_indices[[x]]]})
# 
# pix_in <- 1:pix_num
# 
# x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
# y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
# x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
# y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand
# 
# var_num <- apply(non_na_col,1,sum)
# numLevels_list <- sapply(1:pix_num, function(x){ rep(1,times=var_num[x])})
# for (i in 1:pix_num){
#   names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
# }


# Load Model output for Ridge ####
####################################
# start with lambda.min

# On the Drive you can find my data in:
# Models/LASSO-Ridge regression/regression_results_Global_wo_interactions/Ridge_lambdamin_threshbadyield005.RData

load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Ridge_lambdamin_threshbadyield005.RData")



# Create specificity, CSI and EDI for Ridge ####
################################################
nb_pix_ridge <- length(ridge_model_lambdamin)

coefs_ridge <- lapply(1:nb_pix_ridge, function(x){coef(ridge_model_lambdamin[[x]])})



pred_ridge <- lapply(1:nb_pix_ridge, function(x){predict(ridge_model_lambdamin[[x]],
                                                         as.matrix(x1_test_list[[x]]),type="response")})



fitted_results_ridge <- lapply(1:nb_pix_ridge, function(x){ifelse(pred_ridge[[x]] > segreg_th,1,0)})

speci_ridge <- sapply(1:nb_pix_ridge, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                         fitted_results_ridge[[x]],
                                                                         threshold = segreg_th)})
con_tab_ridge <-  lapply(1:nb_pix_ridge, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                       fitted_results_ridge[[x]],
                                                                                       threshold = segreg_th)})

csi_ridge <- numeric()  #critical success index
F_ridge <- numeric()    #False alarm rate
H_ridge <- numeric()    #Hit rate

for(pix in 1:nb_pix_ridge){
  csi_ridge[pix] <- con_tab_ridge[[pix]]["0","0"]/(con_tab_ridge[[pix]]["0","0"] +
                                                     con_tab_ridge[[pix]]["1","0"] +
                                                     con_tab_ridge[[pix]]["0","1"])
  
  F_ridge[pix] <- con_tab_ridge[[pix]]["0","1"]/(con_tab_ridge[[pix]]["0","1"] +
                                                     con_tab_ridge[[pix]]["1","1"])
  
  H_ridge[pix] <- con_tab_ridge[[pix]]["0","0"]/(con_tab_ridge[[pix]]["0","0"] +
                                                   con_tab_ridge[[pix]]["1","0"])
  
  if(is.na(con_tab_ridge[[pix]]["0","0"])){ #no extreme event forecasted => no first line in contengency table
    csi_ridge[pix] <- 0
    H_ridge[pix] <- 0
    F_ridge[pix] <- 0
  }
  if(is.na(con_tab_ridge[[pix]]["1","0"])){ #No good year forecasted. Problematic pixels for this model
    csi_ridge[pix] <- NA
    H_ridge[pix] <- NA
    F_ridge[pix] <- NA
  }
}#end for pix

EDI_ridge <- (log(F_ridge)-log(H_ridge))/(log(F_ridge)+log(H_ridge))


Result_matrix_Ridge <- cbind(speci_ridge, csi_ridge, EDI_ridge, coord_subset[,1],coord_subset[,2])
colnames(Result_matrix_Ridge) = c("speci", "CSI", "EDI", "lon", "lat")


# Load Model output for Lasso w/o interactions ####
##################################################
# start with lambda.min

# On the Drive you can find my data in:
# Models/LASSO-Ridge regression/regression_results_Global_wo_interactions/Lasso_lambdamin_threshbadyield005.RData

load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambdamin_threshbadyield005.RData")


# Create specificity, CSI and EDI for Lasso w/o interactions ####
#################################################################
nb_pix_simplelasso <- length(lasso_model_lambdamin)

coefs_simplelasso <- lapply(1:nb_pix_simplelasso, function(x){coef(lasso_model_lambdamin[[x]])})



pred_simplelasso <- lapply(1:nb_pix_simplelasso, function(x){predict(lasso_model_lambdamin[[x]],
                                                         as.matrix(x1_test_list[[x]]),type="response")})



fitted_results_simplelasso <- lapply(1:nb_pix_simplelasso, function(x){ifelse(pred_simplelasso[[x]] > segreg_th,1,0)})

speci_simplelasso <- sapply(1:nb_pix_simplelasso, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                                fitted_results_simplelasso[[x]],
                                                                                threshold = segreg_th)})
con_tab_simplelasso <-  lapply(1:nb_pix_simplelasso, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                       fitted_results_simplelasso[[x]],
                                                                                       threshold = segreg_th)})

csi_simplelasso <- numeric()  #critical success index
F_simplelasso <- numeric()    #False alarm rate
H_simplelasso <- numeric()    #Hit rate

for(pix in 1:nb_pix_simplelasso){
  csi_simplelasso[pix] <- con_tab_simplelasso[[pix]]["0","0"]/(con_tab_simplelasso[[pix]]["0","0"] +
                                                     con_tab_simplelasso[[pix]]["1","0"] +
                                                     con_tab_simplelasso[[pix]]["0","1"])
  
  F_simplelasso[pix] <- con_tab_simplelasso[[pix]]["0","1"]/(con_tab_simplelasso[[pix]]["0","1"] +
                                                   con_tab_simplelasso[[pix]]["1","1"])
  
  H_simplelasso[pix] <- con_tab_simplelasso[[pix]]["0","0"]/(con_tab_simplelasso[[pix]]["0","0"] +
                                                   con_tab_simplelasso[[pix]]["1","0"])
  
  if(is.na(con_tab_simplelasso[[pix]]["0","0"])){ #no extreme event forecasted => no first line in contengency table
    csi_simplelasso[pix] <- 0
    H_simplelasso[pix] <- 0
    F_simplelasso[pix] <- 0
  }
  if(is.na(con_tab_simplelasso[[pix]]["1","0"])){ #No good year forecasted. Problematic pixels for this model
    csi_simplelasso[pix] <- NA
    H_simplelasso[pix] <- NA
    F_simplelasso[pix] <- NA
  }
}#end for pix

EDI_simplelasso <- (log(F_simplelasso)-log(H_simplelasso))/(log(F_simplelasso)+log(H_simplelasso))


Result_matrix_simplelasso <- cbind(speci_simplelasso, csi_simplelasso, EDI_simplelasso, coord_subset[,1],coord_subset[,2])
colnames(Result_matrix_simplelasso) = c("speci", "CSI", "EDI", "lon", "lat")




# Load Model output for bestglm ####
####################################

load(file="D:/PROJECTS/DAMOCLES/BestGLM_rep1000_worksp/BestGLm_complete.RData")


# Create specificity, CSI and EDI for bestglm ####
##################################################
pix_num<-965

test_length <- pix_num 

mypred <- lapply(Pixel_ok, function(x){predict(model_CV[[x]]$BestModel,Testing_Data[[x]][,1:(ncol(Testing_Data[[x]])-1)],type='response')}) 


work_pix<-Pixel_ok

segreg_th <- 0.5

fitted.results_model <- lapply(seq_along(work_pix), function(x){ifelse(mypred[[x]] > segreg_th,1,0)})


y1_test_list <- lapply(1:pix_num, function(x){Testing_Data[[x]][,ncol(Testing_Data[[x]])]}) # predictand

y1_test_list_red <- lapply(work_pix,function(work_pix){y1_test_list[[work_pix]]})

mis_clas_err <- rep(NA,965)
# mis_clas_err[work_pix] <- sapply(seq_along(work_pix), function(x){misClassError(y1_test_list_red[[x]],mypred[[x]])})
mis_clas_err[work_pix] <- sapply(seq_along(work_pix), function(x){misClassError(actuals = y1_test_list_red[[x]],
                                                                                predictedScores=mypred[[x]],
                                                                                threshold = segreg_th)})

con_tab <-  lapply(seq_along(work_pix), function(x){InformationValue::confusionMatrix(y1_test_list_red[[x]],fitted.results_model[[x]])})
sensi <- rep(NA,965)
speci <- rep(NA,965)
sensi[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::sensitivity(y1_test_list_red[[x]],fitted.results_model[[x]],threshold = segreg_th)})

speci[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::specificity(y1_test_list_red[[x]],fitted.results_model[[x]], threshold = segreg_th)})

obs_pred <- lapply(seq_along(work_pix), function(x){cbind(y1_test_list_red[[x]],fitted.results_model[[x]])})
tp <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==2)})
tn <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==0)})
fp <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==0 & obs_pred[[x]][,2]==1)})
fn <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==1 & obs_pred[[x]][,2]==0)})
con_tab1 <- sapply(seq_along(work_pix), function(x){matrix(c(tp[x],fn[x],fp[x],tn[x]),nrow=2,ncol=2)})
# spec <- tn/(tn+fp) 
# sens <- tp/(tp+fn) 
csi <- rep(NA,965)
csi[work_pix] <- sapply(seq_along(work_pix), function(x){tn[x]/(tn[x]+fp[x]+fn[x])})



# Load Model output for Lasso w interactions ####
###################################################

# source('./Code/Lasso_interact_global.R') # takes ca. 3 hours: load necessary files, calculate model
# source('./Code/Lasso_interact_global_preparation.R') # load necessary files
# On the Drive you can find my data in:
# Models/LASSO_with_interactions/cv_fit_complete.RData.RData
load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_complete.RData") # load model output

#location Pauline's Laptop
# load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/OtherModels/cv_fit_complete.Rdata")

# Identify pixels with failed runs
failed_pixels_lwi <- which(sapply(1:965, function(x) {is.character(cv_fit[[x]])})==1)
work_pix_lwi <- pix_in[-failed_pixels_lwi] # working pixels


# Create specificity, CSI and EDI for Lasso w interactions ####
################################################################

mypred_lwi <- lapply(work_pix_lwi, function(x){predict(cv_fit[[x]],x1_test_list[[x]],type="response")}) 
fitted.results_model_lwi <- lapply(seq_along(work_pix_lwi), function(x){ifelse(mypred_lwi[[x]] > segreg_th,1,0)})
y1_test_list_red_lwi <- lapply(work_pix_lwi,function(work_pix_lwi){y1_test_list[[work_pix_lwi]]})

obs_pred_lwi <- lapply(seq_along(work_pix_lwi), function(x){cbind(y1_test_list_red_lwi[[x]],fitted.results_model_lwi[[x]])})
# tp <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==2)}) # Correct rejections
tn_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(rowSums(obs_pred_lwi[[x]])==0)}) # Hits
fp_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(obs_pred_lwi[[x]][,1]==0 & obs_pred_lwi[[x]][,2]==1)}) # Miss
fn_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(obs_pred_lwi[[x]][,1]==1 & obs_pred_lwi[[x]][,2]==0)}) # False alarm
# con_tab <- sapply(seq_along(work_pix), function(x){matrix(c(tp[x],fn[x],fp[x],tn[x]),nrow=2,ncol=2)})


# Calculate specificity
speci_lwi <- rep(NA,965)
speci_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){InformationValue::specificity(y1_test_list_red_lwi[[x]],fitted.results_model_lwi[[x]],
                                                                                         threshold = segreg_th)})
# Calculate CSI 
csi_lwi <- rep(NA,965)
csi_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){tn_lwi[x]/(tn_lwi[x]+fn_lwi[x]+fp_lwi[x])})

# Calculate EDI 
EDI_lwi <- rep(NA,965)
EDI_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) (log(fn_lwi[x])-log(tn_lwi[x]))/(log(fn_lwi[x])+log(tn_lwi[x])))




# Create comparison maps ####
#############################



# Specificity map #####
#######################





# CSI map ####
##############





# EDI map ####
##############
