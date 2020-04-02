
library(randomForest)
library(pbapply)
library(tidyverse)
library(viridis)

source('./Code/Lasso_interact_global_preparation_incl_ext_ind.R') # monthly data including extreme indices

# model_mtry_list <- vector("list",length=33)
# predValid_list <- vector("list",length=33)
#
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_3_8.RData")
# model_mtry_list[1:6] <- model_mtry_all
# predValid_list[1:6] <- predValid_all
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_9_14.RData")
# model_mtry_list[7:12] <- model_mtry_all
# predValid_list[7:12] <- predValid_all
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_15_20.RData")
# model_mtry_list[13:18] <- model_mtry_all
# predValid_list[13:18] <- predValid_all
#
# save(model_mtry_list,predValid_list,file="D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_3_20.RData")
# # restart R (for saving RAM)
#
# model_mtry_list <- vector("list",length=33)
# predValid_list <- vector("list",length=33)
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_21_25.RData")
# model_mtry_list[19:23] <- model_mtry_all
# predValid_list[19:23] <- predValid_all
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_26_30.RData")
# model_mtry_list[24:28] <- model_mtry_all
# predValid_list[24:28] <- predValid_all
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_31_35.RData")
# model_mtry_list[29:33] <- model_mtry_all
# predValid_list[29:33] <- predValid_all
# save(model_mtry_list,predValid_list,file="D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_21_35.RData")
# # restart R (for saving RAM)

# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_3_20.RData")
# model_mtry_list2 <- model_mtry_list
# predValid_list2 <- predValid_list
# rm(model_mtry_list);rm(predValid_list)
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_21_35.RData")
# model_mtry_list[1:18] <- model_mtry_list2[1:18]
# predValid_list[1:18] <- predValid_list2[1:18]
# save(model_mtry_list,predValid_list,file="D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_3_35.RData")

load("D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_3_20.RData")
# load("D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_21_35.RData"); predValid_list <- predValid_list[19:33]
# # load("D:/user/vogelj/Group_project/Code/Workspaces/rf_esssentials_3_35.RData")
# tree_num <- 18
tree_num <- 14
speci <- vector("numeric",length=tree_num)
sensi <- vector("numeric",length=tree_num)
obs_pred <- vector("list",length=tree_num)
tp <- vector("list",length=tree_num)
tn <- vector("list",length=tree_num)
fp <- vector("list",length=tree_num)
fn <- vector("list",length=tree_num)
csi_rf <- vector("list",length=tree_num)
work_pix <- 1:965
# for (i in 21:35) {
for (i in 3:20) {
  k <- i - 2
  # k <- i - 20
  speci[k] <- mean(sapply(1:length(x1_train_list), function(x){InformationValue::specificity(y1_test_list[[x]],(as.numeric(predValid_list[[k]][[x]])-1))}))
  sensi[k] <- mean(sapply(1:length(x1_train_list), function(x) {InformationValue::sensitivity(y1_test_list[[x]],(as.numeric(predValid_list[[k]][[x]])-1))}))
  obs_pred[[k]] <- lapply(seq_along(work_pix), function(x){cbind(y1_test_list[[x]],as.numeric(predValid_list[[k]][[x]])-1)})
  tp[[k]] <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[k]][[x]])==2)}) # Correct rejections
  tn[[k]] <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[k]][[x]])==0)}) # Hits
  fp[[k]] <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[k]][[x]][,1]==0 & obs_pred[[k]][[x]][,2]==1)}) # Misses
  fn[[k]] <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[k]][[x]][,1]==1 & obs_pred[[k]][[x]][,2]==0)}) # False alarm
  csi_rf[[k]] <- sapply(seq_along(work_pix), function(x){tn[[k]][x]/(tn[[k]][x]+fn[[k]][x]+fp[[k]][x])})
}

csi_mean <- sapply(csi_rf,mean)

plot(sensi)
plot(speci)
plot(csi_mean)
speci
csi_mean


dims <- sapply(1:965, function(x) dim(x1_train_list[[x]]))
plot(dims[1,])            
plot(dims[2,])
min(dims[2,]) # minimum variable number is 16: so I should not exceed that probably

# it does not get much higher with mtry > 16, until csi = 0.323 and speci = 0.381

