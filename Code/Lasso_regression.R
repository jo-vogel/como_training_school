# Fitting of the Lasso model
# Authors: Pauline Rivore, Johannes Vogel, Cristina Deidda

# Structure of the code
# a) Load the standardized data (meteorological variables and extremal indicators)
# b) Create training and testing dataset
# c) Run the cross validation to find lambda.min and lambda.1se
# d) Run the Lasso regression 


# Required libraries ####
library(glmnet);library(InformationValue);library(abind);library(stringr);library(ggplot2)

# Load standardized Data #####
path_data <- message("insert data directory here")
path_model <- message("insert model directory here")
# path_data <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/"
# path_model <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/"

load(paste0(path_data,"extremeindices_and_monthlymeteovar_rescaled_995pix.RData"))
pix_num <- length(Data_xtrm_standardized$longitudes)
nb_years <- dim(Data_xtrm_standardized$yield)[2]

# Create arrays of data #####
yield_3dim <- array(Data_xtrm_standardized$yield, dim = c(pix_num,1,nb_years))
dtr_3dim <- array(Data_xtrm_standardized$dtr, dim = c(pix_num,1,nb_years))
frs_3dim <- array(Data_xtrm_standardized$frs, dim = c(pix_num,1,nb_years))
txx_3dim <- array(Data_xtrm_standardized$txx, dim = c(pix_num,1,nb_years))
tnn_3dim <- array(Data_xtrm_standardized$tnn, dim = c(pix_num,1,nb_years))
rx5_3dim <- array(Data_xtrm_standardized$rx5, dim = c(pix_num,1,nb_years))
tx90p_3dim <- array(Data_xtrm_standardized$tx90p, dim = c(pix_num,1,nb_years))
tn10p_3dim <- array(Data_xtrm_standardized$tn10p, dim = c(pix_num,1,nb_years))

Model_data <- abind(yield_3dim, dtr_3dim, frs_3dim, txx_3dim, tnn_3dim, rx5_3dim, tx90p_3dim, tn10p_3dim
                    , Data_xtrm_standardized$tasmax, Data_xtrm_standardized$vpd, Data_xtrm_standardized$pr, along = 2)
colnames(Model_data) <- c("Yield", "dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p",
                          "tmax_Aug_Y1", "tmax_Sep_Y1", "tmax_Oct_Y1", "tmax_Nov_Y1", "tmax_Dec_Y1", "tmax_Jan_Y2",
                          "tmax_Feb_Y2", "tmax_Mar_Y2", "tmax_Apr_Y2", "tmax_May_Y2", "tmax_Jun_Y2", "tmax_Jul_Y2",
                          "tmax_Aug_Y2", "tmax_Sep_Y2", "tmax_Oct_Y2", "tmax_Nov_Y2", "tmax_Dec_Y2",
                          "vpd_Aug_Y1", "vpd_Sep_Y1", "vpd_Oct_Y1", "vpd_Nov_Y1", "vpd_Dec_Y1", "vpd_Jan_Y2",
                          "vpd_Feb_Y2", "vpd_Mar_Y2", "vpd_Apr_Y2", "vpd_May_Y2", "vpd_Jun_Y2", "vpd_Jul_Y2",
                          "vpd_Aug_Y2", "vpd_Sep_Y2", "vpd_Oct_Y2", "vpd_Nov_Y2", "vpd_Dec_Y2",
                          "pr_Aug_Y1", "pr_Sep_Y1", "pr_Oct_Y1", "pr_Nov_Y1", "pr_Dec_Y1", "pr_Jan_Y2",
                          "pr_Feb_Y2", "pr_Mar_Y2", "pr_Apr_Y2", "pr_May_Y2", "pr_Jun_Y2", "pr_Jul_Y2",
                          "pr_Aug_Y2", "pr_Sep_Y2", "pr_Oct_Y2", "pr_Nov_Y2", "pr_Dec_Y2")

Yield <- Data_xtrm_standardized$yield

threshold <- 0.05 #threshold for bad yields
low_yield <- apply(Yield, MARGIN = 1, FUN = quantile, probs = threshold, na.rm = T)
cy <- t(sapply(1:pix_num, function(x){ifelse(Yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield

cy_reshaped <- array(data = cy,dim = c(dim(cy)[1],1,nb_years))

Model_data[,1,] <-cy_reshaped


# Exclude NA variable columns
na_col <- matrix(data = NA, nrow = pix_num, ncol = dim(Model_data)[2])
for (j in 1:pix_num){
  for (i in 1:dim(Model_data)[2]){
    na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)

na_time <- vector("list", length = pix_num) # for each pixel, the positions of NAs over time
for (j in 1:pix_num){
  na_time[[j]] <- which(is.na(Model_data[j,1,])) # locations of years with NA values
}


# Split data into training and testing data set #####
years_with_na <- vector("logical", length = pix_num)
for (i in 1:pix_num){
  years_with_na[i] <- ifelse(length(na_time[[i]] ) == 0, F, T)
}

training_indices <- vector("list", length = pix_num)
testing_indices <- vector("list", length = pix_num)


# Parameters of the splitting
seed <- 1994
train_size <- 70

set.seed(seed)
for (x in 1:pix_num) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x = (1:nb_years)[-na_time[[x]]], size = floor((nb_years-length(na_time[[x]]))*(train_size/100))))
    testing_indices[[x]] <- (1:nb_years)[-c(na_time[[x]], training_indices[[x]])]
  } else {
    training_indices[[x]] <- sort(sample(1:nb_years, size = floor(nb_years*(train_size/100))))
    testing_indices[[x]] <- (1:nb_years)[-training_indices[[x]]]    
  }
}


Training_Data <- lapply(1:pix_num, function(x){Model_data[x,,training_indices[[x]]]})
Testing_Data <- lapply(1:pix_num, function(x){Model_data[x,,testing_indices[[x]]]})

pix_in <- 1:pix_num


x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand


var_num <- apply(non_na_col, 1, sum)
numLevels_list <- sapply(1:pix_num, function(x){ rep(1, times = var_num[x])})
for (i in 1:pix_num){
  names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
}



# Run the CrossValidation #####
# With the cross validation, we obtain lambda.min and lambda.1se
model_cv_fitting <- list()
nbyears_final_training_data <- numeric()

for (pixel in 1:pix_num) {
  
  var_pix <- as.matrix(x1_train_list[[pixel]])
  yield_pix <- as.matrix(y1_train_list[[pixel]])
  which_na_xtrm <- which(is.na(var_pix[,1]))
  nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))
  
  if (sum(yield_pix[-which_na_xtrm,])<=nbyears_final_training_data[pixel] &
      sum(yield_pix[-which_na_xtrm,])>=nbyears_final_training_data[pixel]-8){
    #For a correct fitting, we need at least 8 occurences of good years and bad years
    #If it is not the case, the model is not run
    model_cv_fitting[[pixel]] <- "Training years w/o na in extremes have less than 8 bad years"
  } else {
    if(length(which_na_xtrm)>0){ #remove potential missing data in extreme indices
      model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix[-which_na_xtrm,],
                                             y = yield_pix[-which_na_xtrm,],
                                             family = "binomial",
                                             alpha = 1, #alpha = 1 corresponds to the Lasso regression
                                             nfolds = 10) #number of folds for the cross validation
    } else {
      if(sum(yield_pix)>=nbyears_final_training_data[pixel]-8){
        #For a correct fitting, we need at least 8 occurences of good years and bad years
        #If it is not the case, the model is not run
        model_cv_fitting[[pixel]] <- "Training years have less than 8 bad years"
      } else {
        model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix,
                                               y = yield_pix,
                                               family = "binomial",
                                               alpha = 1,
                                               nfolds = 10)
      }#end if else not enough bad occurences
      
    }#end if exists na else
    
  }#end if all years kept are good else
  
  print(paste(pixel, "out of", pix_num))
}#end for pixel


save(model_cv_fitting, file = paste0(path_model, "cv_month_xtrm_LASSO_threshbadyield005_seed",
                                     seed, "_train", train_size,"_995pixels.RData"))


# Run the model with lambda1se and lambda min just obtained #####

load(file = paste0(path_model, "cv_month_xtrm_LASSO_threshbadyield005_seed",
                   seed, "_train", train_size, "_995pixels.RData"))


if (model_name == "Lasso"){
  lasso_model_lambdamin <- list()
  lasso_model_lambda1se <- list()
  
  for (pixel in 1:pix_num) {
    if(is.character(model_cv_fitting[[pixel]])){ #If the cross validation could not be run
      lasso_model_lambdamin[[pixel]] <- model_cv_fitting[[pixel]]
      lasso_model_lambda1se[[pixel]] <- model_cv_fitting[[pixel]]
    } else {
      var_pix <- as.matrix(x1_train_list[[pixel]])
      yield_pix <- as.matrix(y1_train_list[[pixel]])
      which_na_xtrm <- which(is.na(var_pix[,1]))
      nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))
      
      training_years_wo_na <- which(!is.na(x1_train_list[[pixel]]$dtr))
      
      if(length(which_na_xtrm)>0){
        
        lasso_model_lambdamin[[pixel]] <- glmnet(x = var_pix[-which_na_xtrm,],
                                                 y = yield_pix[-which_na_xtrm,],
                                                 family = "binomial",
                                                 alpha = 1,
                                                 lambda = model_cv_fitting[[pixel]]$lambda.min)
        
        lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix[-which_na_xtrm,],
                                                 y = yield_pix[-which_na_xtrm,],
                                                 family = "binomial",
                                                 alpha = 1,
                                                 lambda = model_cv_fitting[[pixel]]$lambda.1se)
      } else {
        
        lasso_model_lambdamin[[pixel]] <- glmnet(x = var_pix, y = yield_pix,
                                                 family = "binomial",
                                                 alpha = 1,
                                                 lambda = model_cv_fitting[[pixel]]$lambda.min)
        
        lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix, y = yield_pix,
                                                 family = "binomial",
                                                 alpha = 1,
                                                 lambda = model_cv_fitting[[pixel]]$lambda.1se)
      }#end if else there exists na in extremes
      
    } #end ifelse the cross validation could not be run
    print(paste(pixel, "out of", pix_num))
    
  }#end for pixel
  
  save(lasso_model_lambdamin, file = paste0(path_model, "Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
                                            seed, "_train", train_size,"_995pixels.RData"))
  save(lasso_model_lambda1se, file = paste0(path_model, "Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
                                            seed, "_train", train_size,"_995pixels.RData"))
}