adjust_cutoff_bestglm <- function(x1_train_list, y1_train_list,  mypred_train, work_pix, cost_fp = 100, cost_fn = 100){
  
  
  
  # Finding appropriate cutoff level ####
  #######################################
  
  # Procedure for whole training data set ####
  source("./Code/unbalanced_funtions.R") 
  
  
  
  mypred_train <- lapply(Pixel_ok, function(x){predict(BestGlm_model[[x]]$BestModel,Training_Data[[x]][,1:(ncol(Training_Data[[x]])-1)],type='response')}) 
  y1_train_list_red <- lapply(work_pix,function(work_pix){y1_train_list[[work_pix]]})  
  # mypred_train <- lapply(work_pix, function(x){predict(cv_fit[[x]],x1_train_list[[x]],type="response")}) 
  
  # Data set with actuals and predictions
  data_train_all <- pblapply(1:length(work_pix), function(x){ data.frame("Actuals"= y1_train_list_red[[x]], "Predictions"=mypred_train[[x]])}) # train data
  # Calculate confusion matrix
  cm_info_all <-  pblapply(1:length(work_pix), function(x){ConfusionMatrixInfo( data = data_train_all[[x]], predict = "Predictions", 
                                                                                actual = "Actuals", cutoff = .5 )})
  # Calculate ROC curve and cost function
  roc_info_all <- pblapply(1:length(work_pix), function(x){ROCInfo( data = cm_info_all[[x]]$data, predict = "predict", 
                                                                    actual = "actual", cost.fp = cost_fp, cost.fn = cost_fn )}) # note: the cutoff from cm_info_all has no role here
  cutoff_avg <- pbsapply(1:length(work_pix), function(x){roc_info_all[[x]]$cutoff}) # find the cutoff value
  return(mean(cutoff_avg)) # calculate the average cutoff value
  # boxplot(cutoff_avg)
}
