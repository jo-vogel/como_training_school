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

message('Adjust the path accordingly.')
# Pauline:
# path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
# output_path <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Images/Model_Comparison/"

# Johannes:
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
output_path <- "D:/user/vogelj/Group_project/Output/Plots"

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

# load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Ridge_lambdamin_threshbadyield005.RData")
load(file = "D:/user/vogelj/Group_project/Code/Workspaces/Ridge_lambdamin_threshbadyield005.RData")


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

# load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambdamin_threshbadyield005.RData")
load(file = "D:/user/vogelj/Group_project/Code/Workspaces/Lasso_lambdamin_threshbadyield005.RData")


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

# On the Drive you can find my data in:
# Models/BestGlm/BestGlm_complete.RData

# load(file="D:/PROJECTS/DAMOCLES/BestGLM_rep1000_worksp/BestGLm_complete.RData")
load(file = "D:/user/vogelj/Group_project/Code/Workspaces/BestGLm_complete.RData")


#location Pauline's Laptop
# load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/OtherModels/BestGLm_complete.RData")

# Create specificity, CSI and EDI for bestglm ####
##################################################
pix_num<-965

mypred <- lapply(Pixel_ok, function(x){predict(BestGlm_model[[x]]$BestModel,Testing_Data[[x]][,1:(ncol(Testing_Data[[x]])-1)],type='response')}) 


nb_pix_bestglm<-Pixel_ok

segreg_th <- 0.5

fitted.results_model_bestglm <- lapply(seq_along(nb_pix_bestglm), function(x){ifelse(mypred[[x]] > segreg_th,1,0)})


y1_test_list_bestglm <- lapply(1:pix_num, function(x){Testing_Data[[x]][,ncol(Testing_Data[[x]])]}) # predictand

y1_test_list_red <- lapply(nb_pix_bestglm,function(work_pix){y1_test_list_bestglm[[work_pix]]})

mis_clas_err_bestglm <- rep(NA,965)
# mis_clas_err[work_pix] <- sapply(seq_along(work_pix), function(x){misClassError(y1_test_list_red[[x]],mypred[[x]])})
mis_clas_err_bestglm[nb_pix_bestglm] <- sapply(seq_along(nb_pix_bestglm), function(x){misClassError(actuals = y1_test_list_red[[x]],
                                                                                                    predictedScores=mypred[[x]],
                                                                                                    threshold = segreg_th)})

sensi_bestglm <- rep(NA,965)
speci_bestglm <- rep(NA,965)
con_tab_bestglm<- rep(NA,965)

con_tab_bestglm[nb_pix_bestglm]  <-  lapply(seq_along(nb_pix_bestglm), function(x){InformationValue::confusionMatrix(y1_test_list_red[[x]],fitted.results_model_bestglm[[x]])})

sensi_bestglm[nb_pix_bestglm] <- sapply(seq_along(nb_pix_bestglm), function(x){InformationValue::sensitivity(y1_test_list_red[[x]],fitted.results_model_bestglm[[x]],threshold = segreg_th)})

speci_bestglm[nb_pix_bestglm] <- sapply(seq_along(nb_pix_bestglm), function(x){InformationValue::specificity(y1_test_list_red[[x]],fitted.results_model_bestglm[[x]], threshold = segreg_th)})

obs_pred <- lapply(seq_along(nb_pix_bestglm), function(x){cbind(y1_test_list_red[[x]],fitted.results_model_bestglm[[x]])})
tp <- sapply(seq_along(nb_pix_bestglm), function(x){sum(rowSums(obs_pred[[x]])==2)})
tn <- sapply(seq_along(nb_pix_bestglm), function(x){sum(rowSums(obs_pred[[x]])==0)})
fp <- sapply(seq_along(nb_pix_bestglm), function(x){sum(obs_pred[[x]][,1]==0 & obs_pred[[x]][,2]==1)})
fn <- sapply(seq_along(nb_pix_bestglm), function(x){sum(obs_pred[[x]][,1]==1 & obs_pred[[x]][,2]==0)})


csi_bestglm <- rep(NA,965)
csi_bestglm[nb_pix_bestglm] <- sapply(seq_along(nb_pix_bestglm), function(x){tn[x]/(tn[x]+fp[x]+fn[x])})

csi_bestglm <- rep(NA,965)  #critical success index
F_bestglm <- rep(NA,965)   #False alarm rate
H_bestglm <- rep(NA,965)   #Hit rate

for(pix in nb_pix_bestglm){
  
  csi_bestglm[pix] <- con_tab_bestglm[[pix]]["0","0"]/(con_tab_bestglm[[pix]]["0","0"] +
                                                         con_tab_bestglm[[pix]]["1","0"] +
                                                         con_tab_bestglm[[pix]]["0","1"])
  
  F_bestglm[pix] <- con_tab_bestglm[[pix]]["0","1"]/(con_tab_bestglm[[pix]]["0","1"] +
                                                       con_tab_bestglm[[pix]]["1","1"])
  
  H_bestglm[pix] <- con_tab_bestglm[[pix]]["0","0"]/(con_tab_bestglm[[pix]]["0","0"] +
                                                       con_tab_bestglm[[pix]]["1","0"])
  
  if(is.na(con_tab_bestglm[[pix]]["0","0"])){ #no extreme event forecasted => no first line in contengency table
    csi_bestglm[pix] <- 0
    H_bestglm[pix] <- 0
    F_bestglm[pix] <- 0
  }
  if(is.na(con_tab_bestglm[[pix]]["1","0"])){ #No good year forecasted. Problematic pixels for this model
    csi_bestglm[pix] <- NA
    H_bestglm[pix] <- NA
    F_bestglm[pix] <- NA
  }
}#end for pix


EDI_bestglm <- (log(F_bestglm)-log(H_bestglm))/(log(F_bestglm)+log(H_bestglm))


Result_matrix_bestglm<- cbind(speci_bestglm, csi_bestglm, EDI_bestglm, coord_subset[,1],coord_subset[,2])
colnames(Result_matrix_bestglm) = c("speci", "CSI", "EDI", "lon", "lat")



# Load Model output for Lasso w interactions ####
###################################################

# source('./Code/Lasso_interact_global.R') # takes ca. 3 hours: load necessary files, calculate model
# source('./Code/Lasso_interact_global_preparation.R') # load necessary files
# On the Drive you can find my data in:
# Models/LASSO_with_interactions/cv_fit_complete.RData.RData
load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_complete.RData") # load model output

#location on Pauline's Laptop
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
tp_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(rowSums(obs_pred_lwi[[x]])==2)}) # Correct rejections
tn_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(rowSums(obs_pred_lwi[[x]])==0)}) # Hits
fp_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(obs_pred_lwi[[x]][,1]==0 & obs_pred_lwi[[x]][,2]==1)}) # Misses
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
hr_lwi <- rep(NA,965)
hr_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) tn_lwi[x]/(tn_lwi[x]+fp_lwi[x])) # hit rate
far_lwi <- rep(NA,965)
far_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) fn_lwi[x]/(fn_lwi[x]+tp_lwi[x])) # false alarm rate
EDI_lwi <- rep(NA,965)
# EDI_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) (log(fn_lwi[x]/(fn_lwi[x]+tp_lwi[x]))-log(tn_lwi[x]/(tn_lwi[x]+fp_lwi[x])))/(log(fn_lwi[x]/(fn_lwi[x]+tp_lwi[x]))+log(tn_lwi[x]/(tn_lwi[x]+fp_lwi[x]))))
EDI_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) (log(far_lwi[x])-log(hr_lwi[x]))/(log(far_lwi[x])+log(hr_lwi[x])))



# Create comparison maps ####
#############################
world <- map_data("world")
substract_score_plot <- function(score_name, score_1, model1_name, score_2, model2_name){#score name in c("speci","csi", EDI")
  sub_score <- score_1 - score_2
  DF_sub <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sub_score = sub_score)
  
  ggplot(data = DF_sub, aes(x=DF_sub$lon, y=DF_sub$lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=DF_sub$sub_score)) +
    scale_color_gradient2(low = "blue", high = "red") +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color=paste("Diff in", score_name),
         title = paste("Difference",score_name,model1_name, "-", model2_name),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  
}


# Specificity map #####
#######################
pairs(cbind(speci_bestglm, speci_lwi, speci_ridge, speci_simplelasso))



possible_pairs <- c("Best GLM-Ridge",
                    "Lasso with interact-Ridge",
                    "Lasso simple-Ridge",
                    "Best GLM-Lasso with interact",
                    "Best GLM-Lasso simple",
                    "Lasso with intersact-Lasso simple")

score <- "Specificity"
model_1 <- "Best GLM"
score_1 <- speci_bestglm
model_2 <- "Ridge"
score_2 <- speci_ridge

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)

# CSI map ####
##############
pairs(cbind(csi_bestglm, csi_lwi, csi_ridge, csi_simplelasso))

score <- "CSI"
model_1 <- "Lasso with interact"
score_1 <- csi_lwi
model_2 <- "Lasso simple"
score_2 <- csi_simplelasso

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)



# EDI map ####
##############
pairs(cbind(EDI_bestglm, EDI_lwi, EDI_ridge, EDI_simplelasso))


possible_pairs <- c("Best GLM-Ridge",
                    "Lasso with interact-Ridge",
                    "Lasso simple-Ridge",
                    "Best GLM-Lasso with interact",
                    "Best GLM-Lasso simple",
                    "Lasso with intersact-Lasso simple")

score <- "EDI"
model_1 <- "Lasso with interact"
score_1 <- EDI_lwi
model_2 <- "Lasso simple"
score_2 <- EDI_simplelasso

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)
