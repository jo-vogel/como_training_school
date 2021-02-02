adjust_cutoff_local <- function(model_vector, x1_train, y1_train, cost_fp = 100, cost_fn = 100){
  
  # Finding appropriate cutoff level ####
  #######################################
  
  # Procedure for whole training data set ####
  source("./Code/unbalanced_funtions.R") 
  
  mypred_train <- predict(model_vector, x1_train, type="response")
  
  # Data set with actuals and predictions
  data_train_all <- data.frame("Actuals"=as.numeric(y1_train), "Predictions"=as.numeric(mypred_train)) # train data
  # Calculate confusion matrix
  cm_info_all <-  ConfusionMatrixInfo( data = data_train_all, predict = "Predictions", actual = "Actuals", cutoff = .5 )
  # Calculate ROC curve and cost function
  roc_info_all <- ROCInfo( data = cm_info_all$data, predict = "predict",
                           actual = "actual", cost.fp = cost_fp, cost.fn = cost_fn )
  # note: the cutoff from cm_info_all has no role here
  cutoff_avg <- roc_info_all$cutoff # find the cutoff value
  return(cutoff_avg) # calculate the average cutoff value
  
}