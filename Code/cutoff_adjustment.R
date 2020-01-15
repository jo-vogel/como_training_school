adjust_cutoff <- function(x1_train_list, y1_train_list, work_pix, cost_fp = 100, cost_fn = 100){
  


# Finding appropriate cutoff level ####
#######################################

# Procedure for whole training data set ####
source("./Code/unbalanced_funtions.R") 

y1_train_list_red <- lapply(work_pix,function(work_pix){y1_train_list[[work_pix]]})  
mypred_train <- lapply(work_pix, function(x){predict(cv_fit[[x]],x1_train_list[[x]],type="response")}) 

# Data set with actuals and predictions
data_train_all <- pblapply(1:length(work_pix), function(x){ data.frame("Actuals"=y1_train_list_red[[x]], "Predictions"=mypred_train[[x]])}) # train data
# Calculate confusion matrix
cm_info_all <-  pblapply(1:length(work_pix), function(x){ConfusionMatrixInfo( data = data_train_all[[x]], predict = "Predictions", 
                                                                              actual = "Actuals", cutoff = .5 )})
# Calculate ROC curve and cost function
roc_info_all <- pblapply(1:length(work_pix), function(x){ROCInfo( data = cm_info_all[[x]]$data, predict = "predict", 
                                                                  actual = "actual", cost.fp = cost_fp, cost.fn = cost_fn )}) # note: the cutoff from cm_info_all has no role here
cutoff_avg <- pbsapply(1:length(work_pix), function(x){roc_info_all[[x]]$cutoff}) # find the cutoff value
return(mean(cutoff_avg)) # calculate the average cutoff value
# boxplot(cutoff_avg)


# Assess performance with adjusted cutoff in test data set ####
# The performance i assessed in the test data set using the cutoff value determined with the training data set
# data_test_all <- pblapply(1:length(work_pix), function(x){ data.frame("Actuals"=y1_test_list_red[[x]], "Predictions"=mypred[[x]])})
# 
# cm_info_all_test<-  pblapply(1:length(work_pix), function(x){ConfusionMatrixInfo( data = data_test_all[[x]], predict = "Predictions", 
#                                                                                   actual = "Actuals", cutoff = 0.5 )})
# tp_pre <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_test[[x]]$data$type)[4]}) # True positive, Correct rejections
# tn_pre <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_test[[x]]$data$type)[3]}) # True negative, Hits
# fp_pre <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_test[[x]]$data$type)[2]}) # False positive, Misses
# fn_pre <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_test[[x]]$data$type)[1]}) # False negative, False alarm
# Adjust confusion matrix for test data using the suggested cutoff value
# cm_info_all_adj <-  pblapply(1:length(work_pix), function(x){ConfusionMatrixInfo( data = data_test_all[[x]], predict = "Predictions", 
#                                                                                   actual = "Actuals", cutoff = mean(cutoff_avg) )})
# 
# tp_adj <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_adj[[x]]$data$type)[4]}) # True positive, Correct rejections
# tn_adj <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_adj[[x]]$data$type)[3]}) # True negative, Hits
# fp_adj <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_adj[[x]]$data$type)[2]}) # False positive, Misses
# fn_adj <- sapply(seq_along(work_pix), function(x){summary(cm_info_all_adj[[x]]$data$type)[1]}) # False negative, False alarm

# Plotting of performance metrices before and after
# mean(tp_adj,na.rm=T);mean(tp_pre,na.rm=T) # less true positives after adjusting
# median(tp_adj,na.rm=T);median(tp_pre,na.rm=T) # less true positives after adjusting
# plot(tp_adj,col='green',ylim=c(100,950)); points(tp_pre,col='blue') # less true positives after adjusting
# mean(tn_adj,na.rm=T);mean(tn_pre,na.rm=T) # seems like less true negatives after adjusting, but this is due to outliers before
# median(tn_adj,na.rm=T);median(tn_pre,na.rm=T) # more true negatives after adjusting
# plot(tn_adj,col='green'); points(tn_pre,col='blue') # some weird outliers are produced
# plot(tn_adj,col='green',ylim=c(0,100)); points(tn_pre,col='blue') # more true negatives after adjusting
# mean(fp_adj,na.rm=T);mean(fp_pre,na.rm=T) # less false positives
# median(fp_adj,na.rm=T);median(fp_pre,na.rm=T) # less false positives
# plot(fp_adj,col='green'); points(fp_pre,col='blue') # some weird outliers are produced
# plot(fp_adj,col='green',ylim=c(0,100)); points(fp_pre,col='blue') # less false positives
# mean(fn_adj,na.rm=T);mean(fn_pre,na.rm=T) # more false negatives after adjusting
# median(fn_adj,na.rm = T);median(fn_pre,na.rm = T) # more false negatives after adjusting
# plot(fn_adj,col='green'); points(fn_pre,col='blue') # more false negatives after adjusting

# CSI with adjusted cutoff
# csi_adj <- rep(NA,965)
# csi_adj[work_pix] <- sapply(seq_along(work_pix), function(x){tn_adj[x]/(tn_adj[x]+fp_adj[x]+fn_adj[x])})

# Sensitivity and specificity with adjusted cutoff
# fitted.results_model_adj <- lapply(seq_along(work_pix), function(x){ifelse(mypred[[x]] > mean(cutoff_avg),1,0)})
# sensi_adj <- rep(NA,965)
# speci_adj <- rep(NA,965)
# sensi_adj[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::sensitivity(y1_test_list_red[[x]],fitted.results_model_adj[[x]],
#                                                                                              threshold = mean(cutoff_avg))})
# speci_adj[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::specificity(y1_test_list_red[[x]],fitted.results_model_adj[[x]],
#                                                                                              threshold = mean(cutoff_avg))})
# 
# # csi_diff <- csi_adj - csi
# plot(csi_diff);mean(csi_diff, na.rm=T)
# speci_diff <- speci_adj - speci
# plot(speci_diff);mean(speci_diff, na.rm=T)
# sensi_diff <- sensi_adj- sensi
# plot(sensi_diff);mean(sensi_diff, na.rm=T)
}