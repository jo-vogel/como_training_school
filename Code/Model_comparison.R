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

##### Initialisation, libraries, data

library(ncdf4);library(RColorBrewer);library(viridis)
library(maps);library(mapdata);library(ggplot2)
library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr)
library(foreach);library(doParallel)
library(tictoc); library(pbapply)



# recreate the training and testing dataset ####
################################################

message('Adjust the path accordingly.')
# Pauline:
# path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
# output_path <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Images/Model_Comparison/"

# Johannes:
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
output_path <- "D:/user/vogelj/Group_project/Output/Plots/"

# source('./Code/Lasso_interact_global_preparation.R') # load necessary files
source('./Code/Lasso_interact_global_preparation_incl_ext_ind.R') # monthly data including extreme indices


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
} # end for pix

EDI_ridge <- (log(F_ridge)-log(H_ridge))/(log(F_ridge)+log(H_ridge))


Result_matrix_Ridge <- cbind(speci_ridge, csi_ridge, EDI_ridge, coord_subset[,1],coord_subset[,2])
colnames(Result_matrix_Ridge) = c("speci", "CSI", "EDI", "lon", "lat")


# Load Model output for Lasso w/o interactions ####
##################################################

# On the Drive you can find my data in:
# Models/LASSO-Ridge regression/regression_results_Global_wo_interactions/Lasso_lambda1se_threshbadyield005.RData

# load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambda1se_threshbadyield005.RData")
load(file = "D:/user/vogelj/Group_project/Code/Workspaces/Lasso_lambda1se_month_xtrm_Lasso_threshbadyield005.RData")


# Create specificity, CSI and EDI for Lasso w/o interactions ####
#################################################################
nb_pix_simplelasso <- length(lasso_model_lambda1se)

coefs_simplelasso <- lapply(1:nb_pix_simplelasso, function(x){coef(lasso_model_lambda1se[[x]])})

nb_coeff_simple_lasso <- numeric()

for (pix in pix_in) {
  nb_coeff_simple_lasso[pix]<- sum(coefs_simplelasso[[pix]]!=0)-1
}

pred_simplelasso <- lapply(1:nb_pix_simplelasso, function(x){predict(lasso_model_lambda1se[[x]],
                                                                     as.matrix(x1_test_list[[x]]),type="response")})



fitted_results_simplelasso <- lapply(1:nb_pix_simplelasso, function(x){ifelse(pred_simplelasso[[x]] > segreg_th,1,0)})


sensi_simplelasso <- sapply(1:nb_pix_simplelasso, function(x){InformationValue::sensitivity(as.matrix(y1_test_list[[x]]),
                                                                                            fitted_results_simplelasso[[x]],
                                                                                            threshold = segreg_th)})

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


# Create specificity, CSI and EDI for Lasso w/o interactions with cutoff level ####
###################################################################################
# Adjust cutoff level
# source("./Code/cutoff_adj_glmnet_lambda1se.R")
source("./Code/Simple_Lasso_Ridge_ElasticNet/cutoff_adj_glmnet_lambda1se.R")
y1_train_list_simple_lasso <- y1_train_list
x1_train_list_simple_lasso <- x1_train_list
cost_fp_simple_lasso <- 100 # Misses: this should be associated with a higher cost, as it is more detrimental
cost_fn_simple_lasso <- 100 # False alarms


# cutoff_simple_lasso <- adjust_cutoff(x1_train_list = x1_train_list_simple_lasso, y1_train_list = y1_train_list_simple_lasso,
#                                      work_pix = pix_in, cost_fp = cost_fp_simple_lasso, cost_fn= cost_fn_simple_lasso)
cutoff_simple_lasso <- adjust_cutoff(model_vector=lasso_model_lambda1se, x1_train_list = x1_train_list_simple_lasso, y1_train_list = y1_train_list_simple_lasso,
                                     work_pix = pix_in, cost_fp = cost_fp_simple_lasso, cost_fn= cost_fn_simple_lasso)
segreg_th_adj <- cutoff_simple_lasso # replace the default threshold = 0.5, by the calculated optimal cutoff

segreg_th_adj <- 0.666

fitted_results_simplelasso_adj <- lapply(1:nb_pix_simplelasso, function(x){ifelse(pred_simplelasso[[x]] > segreg_th_adj,1,0)})

speci_simplelasso_adj <- sapply(1:nb_pix_simplelasso, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                                                fitted_results_simplelasso_adj[[x]],
                                                                                                threshold = segreg_th_adj)})

sensi_simplelasso_adj <- sapply(1:nb_pix_simplelasso, function(x){InformationValue::sensitivity(as.matrix(y1_test_list[[x]]),
                                                                                                fitted_results_simplelasso_adj[[x]],
                                                                                                threshold = segreg_th_adj)})


con_tab_simplelasso_adj <-  lapply(1:nb_pix_simplelasso, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                                       fitted_results_simplelasso_adj[[x]],
                                                                                                       threshold = segreg_th_adj)})
csi_simplelasso_adj <- numeric()  #critical success index

for(pix in 1:nb_pix_simplelasso){
  csi_simplelasso_adj[pix] <- con_tab_simplelasso_adj[[pix]]["0","0"]/(con_tab_simplelasso_adj[[pix]]["0","0"] +
                                                                         con_tab_simplelasso_adj[[pix]]["1","0"] +
                                                                         con_tab_simplelasso_adj[[pix]]["0","1"])
  if(is.na(con_tab_simplelasso_adj[[pix]]["0","0"])){ #no extreme event forecasted => no first line in contengency table
    csi_simplelasso_adj[pix] <- 0
  }
  if(is.na(con_tab_simplelasso_adj[[pix]]["1","0"])){ #No good year forecasted. Problematic pixels for this model
    csi_simplelasso_adj[pix] <- NA
  }
}#end for pix







# Load Model output for ElasticNet 0.5 ####
##################################################
# start with lambda.min

# On the Drive you can find my data in:
# Models/LASSO-Ridge regression/regression_results_Global_wo_interactions/Elastic-net_lambdamin_alpha05_threshbadyield005.RData

load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Elastic-net_lambdamin_alpha05_threshbadyield005.RData")
elastic_model_lambdamin <- fit_pix_lambdamin
rm(fit_pix_lambdamin)

# Create specificity, CSI and EDI for ElasticNet 0.5 ####
#################################################################
nb_pix_elastic <- length(elastic_model_lambdamin)

coefs_elastic <- lapply(1:nb_pix_elastic, function(x){coef(elastic_model_lambdamin[[x]])})



pred_elastic <- lapply(1:nb_pix_elastic, function(x){predict(elastic_model_lambdamin[[x]],
                                                             as.matrix(x1_test_list[[x]]),type="response")})



fitted_results_elastic <- lapply(1:nb_pix_elastic, function(x){ifelse(pred_elastic[[x]] > segreg_th,1,0)})

speci_elastic <- sapply(1:nb_pix_elastic, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                                    fitted_results_elastic[[x]],
                                                                                    threshold = segreg_th)})
con_tab_elastic <-  lapply(1:nb_pix_elastic, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                           fitted_results_elastic[[x]],
                                                                                           threshold = segreg_th)})

csi_elastic <- numeric()  #critical success index
F_elastic <- numeric()    #False alarm rate
H_elastic <- numeric()    #Hit rate

for(pix in 1:nb_pix_elastic){
  csi_elastic[pix] <- con_tab_elastic[[pix]]["0","0"]/(con_tab_elastic[[pix]]["0","0"] +
                                                         con_tab_elastic[[pix]]["1","0"] +
                                                         con_tab_elastic[[pix]]["0","1"])
  
  F_elastic[pix] <- con_tab_elastic[[pix]]["0","1"]/(con_tab_elastic[[pix]]["0","1"] +
                                                       con_tab_elastic[[pix]]["1","1"])
  
  H_elastic[pix] <- con_tab_elastic[[pix]]["0","0"]/(con_tab_elastic[[pix]]["0","0"] +
                                                       con_tab_elastic[[pix]]["1","0"])
  
  if(is.na(con_tab_elastic[[pix]]["0","0"])){ #no extreme event forecasted => no first line in contengency table
    csi_elastic[pix] <- 0
    H_elastic[pix] <- 0
    F_elastic[pix] <- 0
  }
  if(is.na(con_tab_elastic[[pix]]["1","0"])){ #No good year forecasted. Problematic pixels for this model
    csi_elastic[pix] <- NA
    H_elastic[pix] <- NA
    F_elastic[pix] <- NA
  }
}#end for pix

EDI_elastic <- (log(F_elastic)-log(H_elastic))/(log(F_elastic)+log(H_elastic))


Result_matrix_elastic <- cbind(speci_elastic, csi_elastic, EDI_elastic, coord_subset[,1],coord_subset[,2])
colnames(Result_matrix_elastic) = c("speci", "CSI", "EDI", "lon", "lat")







# Load Model output for bestglm ####
####################################

# On the Drive you can find my data in:
# Models/BestGlm/BestGlm_complete.RData

 load(file="D:/PROJECTS/DAMOCLES/BestGLM_rep1000_worksp/BestGLm_complete.RData")
#load(file = "D:/user/vogelj/Group_project/Code/Workspaces/BestGLm_complete.RData")


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

x1_test_list_bestglm <- lapply(1:pix_num, function(x){if(ncol(Testing_Data[[x]])>1){Testing_Data[[x]][,1:(ncol(Testing_Data[[x]])-1)]}}) # predictors
y1_test_list_bestglm <- lapply(1:pix_num, function(x){Testing_Data[[x]][,ncol(Testing_Data[[x]])]}) # predictand

x1_train_list_bestglm <- lapply(1:pix_num, function(x){if(ncol(Training_Data[[x]])>1){Training_Data[[x]][,1:(ncol(Training_Data[[x]])-1)]}}) # predictors
y1_train_list_bestglm <- lapply(1:pix_num, function(x){Training_Data[[x]][,ncol(Training_Data[[x]])]}) # predictand




## Adjust cutoff level for BestGlm ##
source("./Code/cutoff_adj_glm.R")
#y1_train_list_bestglm <- y1_train_list
#x1_train_list_bestglm <- x1_train_list
cost_fp_bestglm <- 100 # Misses: this should be associated with a higher cost, as it is more detrimental
cost_fn_bestglm <- 100 # False alarms
work_pix_bestglm<-Pixel_ok


cutoff_bestglm <- adjust_cutoff_bestglm(x1_train_list = x1_train_list_bestglm, y1_train_list = y1_train_list_bestglm,
                                        work_pix = work_pix_bestglm, cost_fp = cost_fp_bestglm, cost_fn= cost_fn_bestglm)
segreg_th_adj <- cutoff_bestglm # replace the default threshold = 0.5, by the calculated optimal cutoff



fitted_results_bestglm_adj <- lapply(1:nb_pix_bestglm, function(x){ifelse(pred_bestglm[[x]] > segreg_th_adj,1,0)})

speci_bestglm_adj <- sapply(1:nb_pix_bestglm, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                                        fitted_results_bestglm_adj[[x]],
                                                                                        threshold = segreg_th_adj)})

sensi_bestglm_adj <- sapply(1:nb_pix_bestglm, function(x){InformationValue::sensitivity(as.matrix(y1_test_list[[x]]),
                                                                                        fitted_results_bestglm_adj[[x]],
                                                                                        threshold = segreg_th_adj)})


con_tab_bestglm_adj <-  lapply(1:nb_pix_bestglm, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                               fitted_results_bestglm_adj[[x]],
                                                                                               threshold = segreg_th_adj)})
csi_bestglm_adj <- numeric()  #critical success index

for(pix in 1:nb_pix_bestglm){
  csi_bestglm_adj[pix] <- con_tab_bestglm_adj[[pix]]["0","0"]/(con_tab_bestglm_adj[[pix]]["0","0"] +
                                                                 con_tab_bestglm_adj[[pix]]["1","0"] +
                                                                 con_tab_bestglm_adj[[pix]]["0","1"])
  if(is.na(con_tab_bestglm_adj[[pix]]["0","0"])){ #no extreme event forecasted => no first line in contengency table
    csi_bestglm_adj[pix] <- 0
  }
  if(is.na(con_tab_bestglm_adj[[pix]]["1","0"])){ #No good year forecasted. Problematic pixels for this model
    csi_bestglm_adj[pix] <- NA
  }
}#end for pix



# Load Model output for Lasso glinternet (with and without) interactions ####
###################################################

message('Indicate model data set to be loaded here')
# model_name <- cv_fit_monthly_no_int # monthly model without interactions
# model_name <- cv_fit_monthly_without_int_incl_ext # monthly model including extreme indices without interactions
model_names <- c("cv_fit_monthly_no_int","cv_fit_monthly_without_int_incl_ext")

# model_name <- "cv_fit_monthly_without_int_incl_ext"
for (model_name in model_names){

# On the Drive you can find my data in:
# with interactions
# folder: Models/Lasso (glinternet)/LASSO_with_interactions/cv_fit_complete.RData
if (model_name == "cv_fit_monthly_with_int") {
  source('./Code/Lasso_interact_global_preparation.R') # load necessary files
  load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_complete.RData") # load monthly model output
}
# load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_seasonal.RData") # load seasonal model output
# load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_monthly_with_int_incl_ext.RData") # monthly model including extreme indices with interactions

# without interactions
# folder: Models/Lasso (glinternet)/LASSO_without_interactions/cv_fit_no_int.RData
  if (model_name == "cv_fit_monthly_no_int") {
    source('./Code/Lasso_interact_global_preparation.R') # load necessary files
    load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_no_int.RData") # monthly model without interactions
  }
# load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_seasonal_no_int.RData") # seasonal model without interactions
  if (model_name == "cv_fit_monthly_without_int_incl_ext") {
    source('./Code/Lasso_interact_global_preparation_incl_ext_ind.R') # monthly data including extreme indices
    load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_monthly_without_int_incl_ext.RData") # monthly model including extreme indices without interactions
  } 

#location on Pauline's Laptop
# load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/OtherModels/cv_fit_complete.Rdata")
# load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/OtherModels/cv_fit_no_int.Rdata")


# Identify pixels with failed runs
failed_pixels_lwi <- which(sapply(1:965, function(x) {is.character(cv_fit[[x]])})==1)
# work_pix_lwi <- pix_in[-failed_pixels_lwi] # working pixels
work_pix_lwi <- if (length(failed_pixels_lwi)==0) {work_pix_lwi <- pix_in
} else  {work_pix_lwi <- pix_in[-failed_pixels_lwi]} # working pixels

# Adjust cutoff level
source("./Code/cutoff_adjustment.R")
y1_train_list_lwi <- y1_train_list
x1_train_list_lwi <- x1_train_list
cost_fp_lwi <- 100 # Misses: this should be associated with a higher cost, as it is more detrimental
cost_fn_lwi <- 100 # False alarms 
mypred_train_lwi <- lapply(work_pix_lwi, function(x){predict(cv_fit[[x]],x1_train_list_lwi[[x]],type="response", lambdaType="lambdaHat1Std")}) # recommended lambdaType
# mypred_train_lwi <- lapply(work_pix_lwi, function(x){predict(cv_fit[[x]],x1_train_list_lwi[[x]],type="response", lambdaType="lambdaHat")}) # default lambdaType

cutoff_lwi <- adjust_cutoff(x1_train_list = x1_train_list_lwi, y1_train_list = y1_train_list_lwi, mypred_train = mypred_train_lwi,
                            work_pix = work_pix_lwi, cost_fp = cost_fp_lwi, cost_fn= cost_fn_lwi)
cutoff_avg_lwi <- mean(cutoff_lwi)
segreg_th_adj <- cutoff_avg_lwi # replace the default threshold = 0.5, by the calculated optimal cutoff
segreg_th <- 0.5 # default

i_1Std <- sapply(work_pix_lwi, function(x){ which(cv_fit[[x]]$lambdaHat1Std == cv_fit[[x]]$lambda)}) # the preferential lambda (tuning parameter): lambdaHat1Std
# i_1Std <- sapply(work_pix_lwi, function(x){ which(cv_fit[[x]]$lambdaHat == cv_fit[[x]]$lambda)}) # the preferential lambda (tuning parameter): lambdaHat
i_1Std_all_pix <- rep(NA,965)
i_1Std_all_pix[work_pix_lwi] <- i_1Std # needed as a workaround (to have an object of lenght=965)

# Create specificity, CSI and EDI for Lasso w interactions ####
################################################################

mypred_lwi <- lapply(work_pix_lwi, function(x){predict(cv_fit[[x]],x1_test_list[[x]],type="response",lambdaType="lambdaHat1Std")}) # recommended lambdaType
# mypred_lwi <- lapply(work_pix_lwi, function(x){predict(cv_fit[[x]],x1_test_list[[x]],type="response",lambdaType="lambdaHat")}) # default lambdaType
fitted.results_model_lwi <- lapply(seq_along(work_pix_lwi), function(x){ifelse(mypred_lwi[[x]] > segreg_th,1,0)})
fitted.results_model_adj_lwi <- lapply(seq_along(work_pix_lwi), function(x){ifelse(mypred_lwi[[x]] > mean(segreg_th_adj),1,0)})
y1_test_list_red_lwi <- lapply(work_pix_lwi,function(work_pix_lwi){y1_test_list[[work_pix_lwi]]})

obs_pred_lwi <- lapply(seq_along(work_pix_lwi), function(x){cbind(y1_test_list_red_lwi[[x]],fitted.results_model_lwi[[x]])})
tp_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(rowSums(obs_pred_lwi[[x]])==2)}) # Correct rejections
tn_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(rowSums(obs_pred_lwi[[x]])==0)}) # Hits
fp_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(obs_pred_lwi[[x]][,1]==0 & obs_pred_lwi[[x]][,2]==1)}) # Misses
fn_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(obs_pred_lwi[[x]][,1]==1 & obs_pred_lwi[[x]][,2]==0)}) # False alarm
obs_pred_adj_lwi <- lapply(seq_along(work_pix_lwi), function(x){cbind(y1_test_list_red_lwi[[x]],fitted.results_model_adj_lwi[[x]])})
tp_adj_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(rowSums(obs_pred_adj_lwi[[x]])==2)}) # Correct rejections
tn_adj_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(rowSums(obs_pred_adj_lwi[[x]])==0)}) # Hits
fp_adj_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(obs_pred_adj_lwi[[x]][,1]==0 & obs_pred_adj_lwi[[x]][,2]==1)}) # Misses
fn_adj_lwi <- sapply(seq_along(work_pix_lwi), function(x){sum(obs_pred_adj_lwi[[x]][,1]==1 & obs_pred_adj_lwi[[x]][,2]==0)}) # False alarm
# con_tab <- sapply(seq_along(work_pix), function(x){matrix(c(tp[x],fn[x],fp[x],tn[x]),nrow=2,ncol=2)})


# Calculate specificity
speci_lwi <- rep(NA,965)
speci_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){InformationValue::specificity(y1_test_list_red_lwi[[x]],fitted.results_model_lwi[[x]],
                                                                                                     threshold = segreg_th)})

# Calculate sensitivity
sensi_lwi <- rep(NA,965)
sensi_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){InformationValue::sensitivity(y1_test_list_red_lwi[[x]],fitted.results_model_lwi[[x]],
                                                                                                     threshold = segreg_th)})

# Sensitivity and specificity with adjusted cutoff
sensi_adj_lwi <- rep(NA,965)
speci_adj_lwi <- rep(NA,965)
sensi_adj_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){InformationValue::sensitivity(y1_test_list_red_lwi[[x]],fitted.results_model_adj_lwi[[x]],
                                                                                                         threshold = mean(segreg_th_adj))}) # note: the threshold is actually not used here, it is already applied in fitted.results_model_adj 
speci_adj_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){InformationValue::specificity(y1_test_list_red_lwi[[x]],fitted.results_model_adj_lwi[[x]],
                                                                                                         threshold = mean(segreg_th_adj))})

# Calculate CSI 
csi_lwi <- rep(NA,965)
csi_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){tn_lwi[x]/(tn_lwi[x]+fn_lwi[x]+fp_lwi[x])})
# CSI with adjusted cutoff
csi_adj_lwi <- rep(NA,965)
csi_adj_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x){tn_adj_lwi[x]/(tn_adj_lwi[x]+fp_adj_lwi[x]+fn_adj_lwi[x])})





#If without interact
csi_lasso_wo_int_glmint <- csi_lwi
csi_lasso_wo_int_glmint_adj <- csi_adj_lwi
speci_lasso_wo_int_glmint <- speci_lwi
speci_lasso_wo_int_glmint_adj <- speci_adj_lwi



# Calculate EDI 
hr_lwi <- rep(NA,965)
hr_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) tn_lwi[x]/(tn_lwi[x]+fp_lwi[x])) # hit rate
far_lwi <- rep(NA,965)
far_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) fn_lwi[x]/(fn_lwi[x]+tp_lwi[x])) # false alarm rate
EDI_lwi <- rep(NA,965)
# EDI_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) (log(fn_lwi[x]/(fn_lwi[x]+tp_lwi[x]))-log(tn_lwi[x]/(tn_lwi[x]+fp_lwi[x])))/(log(fn_lwi[x]/(fn_lwi[x]+tp_lwi[x]))+log(tn_lwi[x]/(tn_lwi[x]+fp_lwi[x]))))
EDI_lwi[work_pix_lwi] <- sapply(seq_along(work_pix_lwi), function(x) (log(far_lwi[x])-log(hr_lwi[x]))/(log(far_lwi[x])+log(hr_lwi[x])))

# Calculate number of coefficients
coefs_lwi <- vector("list",length=965)
coefs_lwi[work_pix_lwi] <- pblapply(work_pix_lwi, function(x){coef(cv_fit[[x]]$glinternetFit)[[i_1Std_all_pix[[x]]]]})
coeff_kep_lwi <- numeric()
for (pix in 1:pix_num) {coeff_kep_lwi[pix] <- length(coefs_lwi[[pix]]$mainEffects$cont)}

# Calculate number of interactions
num_interact <- numeric()
for (pix in 1:pix_num) {
  if (is.null(dim(coefs_lwi[[pix]]$interactions$contcont)[1])){
    num_interact[pix] <- 0
  } else {
    num_interact[pix] <-  dim(coefs_lwi[[pix]]$interactions$contcont)[1]
  }
}


if (model_name == "cv_fit_monthly_no_int") {
  csi_adj_lwi_monthly_no_int <- csi_adj_lwi
  coefs_lwi_monthly_no_int <- coefs_lwi
  coeff_kep_lwi_monthly_no_int <- coeff_kep_lwi
}
if (model_name == "cv_fit_monthly_without_int_incl_ext"){
  csi_adj_lwi_monthly_without_int_incl_ext <- csi_adj_lwi
  coefs_lwi_monthly_without_int_incl_ext <- coefs_lwi
  coeff_kep_lwi_monthly_without_int_incl_ext <- coeff_kep_lwi
}


# print average performance
mean(csi_lwi,na.rm=T);mean(csi_adj_lwi,na.rm=T)
mean(speci_lwi,na.rm=T);mean(speci_adj_lwi,na.rm=T)
mean(sensi_lwi,na.rm=T);mean(sensi_adj_lwi,na.rm=T)
segreg_th_adj
mean(coeff_kep_lwi,na.rm=T)
mean(num_interact,na.rm=T)

# Calculaute preferential lambda
pref_lam <- rep(NA,965)
pref_lam[work_pix_lwi] <- sapply(work_pix_lwi, function(x) cv_fit[[x]]$lambdaHat1Std)
  # write.csv(pref_lam,file="glinternet_lasso_without_interactions_lambdaHat1Std.csv")
pref_lam_hat <- rep(NA,965)
pref_lam_hat[work_pix_lwi] <- sapply(work_pix_lwi, function(x) cv_fit[[x]]$lambdaHat)
  # write.csv(pref_lam_hat,file="glinternet_lasso_without_interactions_lambdaHat.csv")

source("./Code/Plots_performance_and_variables.r")
# source("./Code/Problematic_pixels.r")


}


# Create comparison maps ####
#############################
world <- map_data("world")
substract_score_plot <- function(score_name, score_1, model1_name, score_2, model2_name){#score name in c("speci","csi", EDI")
  sub_score <- score_1 - score_2
  DF_sub <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sub_score = sub_score)
  
  ggplot(data = DF_sub, aes(x=DF_sub$lon, y=DF_sub$lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_tile(aes(fill=DF_sub$sub_score)) +
    scale_fill_gradient2(low = "blue", high = "red") +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(fill=paste("Diff in", score_name),
         title = paste("Difference",score_name,model1_name, "-", model2_name),
         subtitle = paste("Adjusted cut-off", sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  
}


# Specificity map #####
#######################
pairs(cbind(speci_bestglm, speci_lwi, speci_ridge, speci_simplelasso, speci_elastic))



possible_pairs <- c("BestGLM-Ridge",
                    "Lasso_with_interact-Ridge",
                    "LassoSimple-Ridge",
                    "Best GLM-Lasso_with_interact",
                    "Best GLM-LassoSimple",
                    "Lasso_with_interact-LassoSimple",
                    "Elastic-BestGLM",
                    "Elastic-Ridge",
                    "Elastic-LassoSimple",
                    "Elastic-Lasso_with_interact")

score <- "Specificity"
model_1 <- "Lasso_wo_int_lambda1se_adj"
score_1 <- speci_simplelasso_adj
model_2 <- "Lasso_wo_int_glinternet_adj"
score_2 <- speci_lasso_wo_int_glmint_adj

pairs(cbind(speci_simplelasso, speci_simplelasso_adj, speci_lasso_wo_int_glmint,
            speci_lasso_wo_int_glmint_adj))

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)

# CSI map ####
##############
# pairs(cbind(csi_bestglm, csi_lwi, csi_ridge, csi_simplelasso, csi_elastic))

score <- "CSI"
# model_1 <- "Lasso_wo_int_lambda1se_adj"
# score_1 <- csi_simplelasso_adj
# model_2 <- "Lasso_wo_int_glinternet_adj"
# score_2 <- csi_lasso_wo_int_glmint_adj
# model_1 <- model_names[2]
# score_1 <- csi_adj_lwi_monthly_without_int_incl_ext
# model_2 <- model_names[1]
# score_2 <- csi_adj_lwi_monthly_no_int

# source("./Code/Random_forest.R")
model_1 <- "Lasso_glmnet"
score_1 <- csi_simplelasso_adj
model_2 <- "Random_forest"
score_2 <- csi_rf # taken from Random_forest.r

# pairs(cbind(csi_simplelasso, csi_simplelasso_adj, csi_lasso_wo_int_glmint, csi_lasso_wo_int_glmint_adj))

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)



# EDI map ####
##############
pairs(cbind(EDI_bestglm, EDI_lwi, EDI_ridge, EDI_simplelasso))


possible_pairs <- c("BestGLM-Ridge",
                    "Lasso_with_interact-Ridge",
                    "LassoSimple-Ridge",
                    "Best GLM-Lasso_with_interact",
                    "Best GLM-LassoSimple",
                    "Lasso_with_interact-LassoSimple")

score <- "EDI"
model_1 <- "BestGLM"
score_1 <- EDI_bestglm
model_2 <- "Lasso_with_interact"
score_2 <- EDI_lwi

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)


# Boxplots ####
###############
png(filename="D:/user/vogelj/Group_project/Output/Plots/Performance_boxplots.png",res=2000,units="cm",width=15,height=20)
box_names <- c("Bestglm","Lasso with \n interactions","Ridge","Lasso w/o interactions")
par(mfrow=c(3,1),mar=c(3, 4, 4, 2) + 0.1,oma=c(0,0,0,0))
boxplot(cbind(speci_bestglm, speci_lwi, speci_ridge, speci_simplelasso),main="Specificity",names=box_names, col='lightblue')
boxplot(cbind(csi_bestglm, csi_lwi, csi_ridge, csi_simplelasso),main="CSI",names=box_names, col='lightblue')
boxplot(cbind(EDI_bestglm, EDI_lwi, EDI_ridge, EDI_simplelasso),main="EDI",names=box_names, col='lightblue')
dev.off()



# Scatterplot performance VS mean yield

mean_yield <- apply(yield, MARGIN = 1, FUN = mean, na.rm=T)

par(mfrow=c(2,3))
plot(mean_yield, csi_simplelasso, xlab="Mean Yield", ylab="CSI", main="Lasso w/o interaction\nmonthly data")
plot(mean_yield, csi_elastic, xlab="Mean Yield", ylab="CSI", main="Elastic Net\nmonthly data")
plot(mean_yield, csi_ridge, xlab="Mean Yield", ylab="CSI", main="Ridge\nmonthly data")
plot(mean_yield, csi_bestglm, xlab="Mean Yield", ylab="CSI", main="Best GLM\nmonthly data")
plot(mean_yield, csi_lwi, xlab="Mean Yield", ylab="CSI", main="Lasso with interaction\nmonthly data")

par(mfrow=c(2,3))
plot(mean_yield, speci_simplelasso, xlab="Mean Yield", ylab="Specificity", main="Lasso w/o interaction\nmonthly data")
plot(mean_yield, speci_elastic, xlab="Mean Yield", ylab="Specificity", main="Elastic Net\nmonthly data")
plot(mean_yield, speci_ridge, xlab="Mean Yield", ylab="Specificity", main="Ridge\nmonthly data")
plot(mean_yield, speci_bestglm, xlab="Mean Yield", ylab="Specificity", main="Best GLM\nmonthly data")
plot(mean_yield, speci_lwi, xlab="Mean Yield", ylab="Specificity", main="Lasso with interaction\nmonthly data")




