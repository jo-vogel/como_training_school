##############################################################################
###########     Sensitivity analysis on number years available     ###########
###########              Lasso regression (glmnet)                 ###########
###########                                                        ###########
###########               Author: Pauline Rivoire                  ###########
##############################################################################
# It is structured in the following way:
# a) Load the data
# b) build training and testing dataset for seed=1994 training_size=70%, with the 1600 years
# b) exctrat several gripoints and create subsamples of different sizes from the training and testing data
# c) run the Lasso crossvalidation to get lambda min et lambda1se
# d) run Lasso glmnet on training data withlambdamin and lambda 1se to get the coefficients of the logistic regression
# e) Analyse results


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

# Model name and reference in the glmnet function
model_name <- "Lasso" ; no_model <- 1

#which lambda?
lambda_VALS <- c("lambda.min", "lambda.1se")

#threshold for bad yields
Thresholds_to_test <- c(0.03,0.05,0.1,0.15)


##### Initialisation, librairies, data #####

library(ncdf4);library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr);library(tictoc);library(ggplot2);library(viridis)




##### Load standardized Data #####
#On the drive in Data/Global_data/Final_Data
#Pauline's Laptop
load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/extremeindices_and_monthlymeteovar_rescaled_995pix.RData")

pix_num <- length(Data_xtrm_standardized$longitudes)
nb_years <- dim(Data_xtrm_standardized$yield)[2]

##### Process data #####
yield_3dim <- array(Data_xtrm_standardized$yield,dim=c(pix_num,1,nb_years))
dtr_3dim <- array(Data_xtrm_standardized$dtr,dim=c(pix_num,1,nb_years))
frs_3dim <- array(Data_xtrm_standardized$frs,dim=c(pix_num,1,nb_years))
txx_3dim <- array(Data_xtrm_standardized$txx,dim=c(pix_num,1,nb_years))
tnn_3dim <- array(Data_xtrm_standardized$tnn,dim=c(pix_num,1,nb_years))
rx5_3dim <- array(Data_xtrm_standardized$rx5,dim=c(pix_num,1,nb_years))
tx90p_3dim <- array(Data_xtrm_standardized$tx90p,dim=c(pix_num,1,nb_years))
tn10p_3dim <- array(Data_xtrm_standardized$tn10p,dim=c(pix_num,1,nb_years))

Model_data <- abind(yield_3dim,dtr_3dim,frs_3dim,txx_3dim,tnn_3dim,rx5_3dim,tx90p_3dim,tn10p_3dim
                    ,Data_xtrm_standardized$tasmax,Data_xtrm_standardized$vpd,Data_xtrm_standardized$pr,along=2)
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

Yield <- Data_xtrm_standardized$yield

segreg_th <- 0.6582418

# extract several the grippoints ####
load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/final_889pixels_coords.Rdata")


FR_pix <- which.min((final_pixels_coord$longitude-1.1)^2+(final_pixels_coord$latitude-47.7)^2)
US_pix <- which.min((final_pixels_coord$longitude+90)^2+(final_pixels_coord$latitude-44.3)^2)
CH_pix <- which.min((final_pixels_coord$longitude-118.1)^2+(final_pixels_coord$latitude-30.8)^2)


pix3_model_data <- Model_data[c(FR_pix, US_pix, CH_pix),,]

csi_FR <- matrix(data = NA, nrow = 50, ncol = length(Thresholds_to_test))
csi_US <- matrix(data = NA, nrow = 50, ncol = length(Thresholds_to_test))
csi_CH <- matrix(data = NA, nrow = 50, ncol = length(Thresholds_to_test))
nb_bad_y_test_FR <- matrix(data = NA, nrow = 50, ncol = length(Thresholds_to_test))
nb_bad_y_test_US <- matrix(data = NA, nrow = 50, ncol = length(Thresholds_to_test))
nb_bad_y_test_CH <- matrix(data = NA, nrow = 50, ncol = length(Thresholds_to_test))

for (thresh in 1:length(Thresholds_to_test)) {
  tic()
  low_yield <- apply(pix3_model_data[,"Yield",], MARGIN = 1, FUN=quantile, probs=Thresholds_to_test[thresh], na.rm=T)
  cy <- t(sapply(1:3,function(x){ifelse(pix3_model_data[x,"Yield",]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield
  
  cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,nb_years))
  pix3_model_data_tochange <- pix3_model_data
  pix3_model_data_tochange[,"Yield",] <-cy_reshaped
  # Exclude NA variable columns
  na_col <- matrix(data=NA ,nrow=3, ncol=dim(pix3_model_data_tochange)[2])
  for (j in 1:3){
    for (i in 1:dim(pix3_model_data_tochange)[2]){
      na_col[j,i] <- all(is.na(pix3_model_data_tochange[j,i,])) # TRUE if entire column is NA
    }
  }
  non_na_col <- !na_col # columns without NAs
  non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)
  
  na_time <- vector("list", length=3) # for each pixel, the positions of NAs over time
  for (j in 1:3){
    na_time[[j]] <- which(is.na(pix3_model_data_tochange[j,1,])) # locations of years with NA values
  }
  
  
  ##### Split data into training and testing data set #####
  vec <- 1:nb_years
  
  years_with_na <- vector("logical",length=pix_num)
  for (i in 1:3){
    years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
  }
  
  training_indices <- vector("list",length=3)
  testing_indices <- vector("list",length=3)
  
  # Parameters of the splitting
    train_size <- 70
  
  for (seed in c(1:49,1994)) {
    seed_nb <- which(c(1:49,1994)==seed)
    
    set.seed(seed)
    for (x in 1:3) {
      if (years_with_na[x]) {
        training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((nb_years-length(na_time[[x]]))*(train_size/100))))
        testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
      } else {
        training_indices[[x]] <- sort(sample(1:nb_years, size = floor(nb_years*(train_size/100))))
        testing_indices[[x]] <- (1:nb_years)[-training_indices[[x]]]    
      }
    }
    
    
    Training_Data <- lapply(1:3,function(x){pix3_model_data_tochange[x,,training_indices[[x]]]})
    Testing_Data <- lapply(1:3,function(x){pix3_model_data_tochange[x,,testing_indices[[x]]]})
    
    pix_in <- 1:3
    
    
    x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
    y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
    x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
    y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand
    
    
    var_num <- apply(non_na_col,1,sum)
    numLevels_list <- sapply(1:3, function(x){ rep(1,times=var_num[x])})
    for (i in 1:3){
      names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
    }
    
    
    
    #### Run the CrossValidation #####

    model_cv_fitting <- list()
    nbyears_final_training_data <- numeric()
    for (pixel in 1:3) {
      var_pix <- as.matrix(x1_train_list[[pixel]])
      yield_pix <- as.matrix(y1_train_list[[pixel]])
      which_na_xtrm <- which(is.na(var_pix[,1]))
      nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))

      if (sum(yield_pix[-which_na_xtrm,])<=nbyears_final_training_data[pixel] &
          sum(yield_pix[-which_na_xtrm,])>=nbyears_final_training_data[pixel]-8){
        #if all the year kept are actually good years, not prediction possible, or there are only 1 or 2 bad years : not possible
        model_cv_fitting[[pixel]] <- "Training years w/o na in extremes have less than 8 bad years"
      } else {
        if(length(which_na_xtrm)>0){
          model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix[-which_na_xtrm,],
                                                 y = yield_pix[-which_na_xtrm,],
                                                 family = "binomial", alpha = no_model, nfolds = 10)
        } else {
          if(sum(yield_pix)>=nbyears_final_training_data[pixel]-8){
            model_cv_fitting[[pixel]] <- "Training years have less than 8 bad years"
          } else {
            model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix,
                                                   y = yield_pix,
                                                   family = "binomial", alpha = no_model, nfolds = 10)
          }#end if else too half years == 0 kg/ha

        }#end if exists na else

      }#end if all years kept are good else


    }#end for pixel

    #30 sec

    save(model_cv_fitting, file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/SensitivityAnalysis/threshold_3GP/3GP_cv_month_xtrm_LASSO_threshbadyield",
                                         sprintf("%03d",100*Thresholds_to_test[thresh]),"_seed",seed,"_train70_995pixels.RData"))
    
    
    ##### Run the model with lambda1se and lambda min #####
    
    
    load(file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/SensitivityAnalysis/threshold_3GP/3GP_cv_month_xtrm_LASSO_threshbadyield",
                       sprintf("%03d",100*Thresholds_to_test[thresh]),"_seed",seed,"_train70_995pixels.RData"))
    
    
    lasso_model_lambdamin <- list()
    lasso_model_lambda1se <- list()
    
    for (pixel in 1:3) {
      if(is.character(model_cv_fitting[[pixel]])){
        lasso_model_lambdamin[[pixel]] <- model_cv_fitting[[pixel]]
        lasso_model_lambda1se[[pixel]] <- model_cv_fitting[[pixel]]
      } else {
        var_pix <- as.matrix(x1_train_list[[pixel]])
        yield_pix <- as.matrix(y1_train_list[[pixel]])
        which_na_xtrm <- which(is.na(var_pix[,1]))
        nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))
        
        training_years_wo_na <- which(!is.na(x1_train_list[[pixel]]$dtr))
        
        if(length(which_na_xtrm)>0){
          
          lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix[-which_na_xtrm,],
                                                   y = yield_pix[-which_na_xtrm,],
                                                   family = "binomial", alpha = no_model,
                                                   lambda = model_cv_fitting[[pixel]]$lambda.1se)
        } else {
          
          
          lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix, y = yield_pix,
                                                   family = "binomial", alpha = no_model,
                                                   lambda = model_cv_fitting[[pixel]]$lambda.1se)
        }#end if exists na else
        
      } #end ifelse
      
      
    }#end for pixel
    
    save(lasso_model_lambda1se, file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/SensitivityAnalysis/threshold_3GP/3GP_Lasso_lambdamin_month_xtrm_LASSO_threshbadyield",
                                              sprintf("%03d",100*Thresholds_to_test[thresh]),"_seed",seed,"_train70_995pixels.RData"))
    
    
    coeff  <-list()
    speci <- rep(NA, 3)
    sensi <- rep(NA, 3)
    
    mypred <- vector("list", 3)
    fitted_bad_yield <- vector("list", 3)
    nb_training_years <- rep(NA, 3)
    nb_testing_years <- rep(NA, 3)
    
    for (pixel in 1:3) {
      coeff[[pixel]] <- coefficients(lasso_model_lambda1se[[pixel]])
      
      mypred[[pixel]] <- predict(lasso_model_lambda1se[[pixel]], as.matrix(x1_test_list[[pixel]]),type="response")
      
      fitted_bad_yield[[pixel]] <- ifelse(mypred[[pixel]] > segreg_th,1,0)
      
      speci[pixel] <- InformationValue::specificity(actuals = as.matrix(y1_test_list[[pixel]]),
                                                    predictedScores = fitted_bad_yield[[pixel]],
                                                    threshold = segreg_th)
      sensi[pixel] <- InformationValue::sensitivity(actuals = as.matrix(y1_test_list[[pixel]]),
                                                    predictedScores = fitted_bad_yield[[pixel]],
                                                    threshold = segreg_th)
      
      con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(y1_test_list[[pixel]]),
                                                   predictedScores = fitted_bad_yield[[pixel]],
                                                   threshold = segreg_th)
      if(pixel==1){
        csi_FR[seed_nb,thresh] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
        if(is.na(con_tab["0","0"])){
          csi_FR[seed_nb,thresh] <- 0
        }
        nb_bad_y_test_FR[seed_nb,thresh] <- sum(1-y1_test_list[[pixel]])
      }#pix FR
      if(pixel==2){
        csi_US[seed_nb,thresh] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
        if(is.na(con_tab["0","0"])){
          csi_US[seed_nb,thresh] <- 0
        }
        nb_bad_y_test_US[seed_nb,thresh] <- sum(1-y1_test_list[[pixel]])
      }#pix US
      if(pixel==3){
        csi_CH[seed_nb,thresh] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
        if(is.na(con_tab["0","0"])){
          csi_CH[seed_nb,thresh] <- 0
        }
        nb_bad_y_test_CH[seed_nb,thresh] <- sum(1-y1_test_list[[pixel]])
      }#pix CH
      
      
      nb_training_years[pixel] <- length(training_indices[[pixel]])
      nb_testing_years[pixel] <- length(testing_indices[[pixel]])
      
    }#end pixel
    print(paste("seed ", seed_nb, "out of 50"))
  } #end for seed
  
  
  toc()
}#end for thresh

colnames(csi_FR) <- c("thresh 0.03", "thresh 0.05","thresh 0.10","thresh 0.15")
colnames(csi_US) <- c("thresh 0.03", "thresh 0.05","thresh 0.10","thresh 0.15")
colnames(csi_CH) <- c("thresh 0.03", "thresh 0.05","thresh 0.10","thresh 0.15")

RESU <- list(csi_FR=csi_FR, csi_US=csi_US, csi_CH=csi_CH,
             nb_bad_y_test_FR=nb_bad_y_test_FR,
             nb_bad_y_test_US=nb_bad_y_test_US,
             nb_bad_y_test_CH=nb_bad_y_test_CH)
save(RESU,
     file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/SensitivityAnalysis/threshold_3GP/sensi_3GP_csi_nbbadyearstest.RData")



load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/SensitivityAnalysis/threshold_3GP/sensi_3GP_csi_nbbadyearstest.RData")

boxplot(RESU$csi_FR, main="Grid point in France", ylab="CSI",
        xlab="Bad year threshold (percentile)", names=c("0.03", "0.05", "0.1", "0.15"))
boxplot(RESU$csi_US, main="Grid point in USA", ylab="CSI",
        xlab="Bad year threshold (percentile)", names=c("0.03", "0.05", "0.1", "0.15"))
boxplot(RESU$csi_CH, main="Grid point in China", ylab="CSI",
        xlab="Bad year threshold (percentile)", names=c("0.03", "0.05", "0.1", "0.15"))
plot(c(RESU$nb_bad_y_test_FR), c(RESU$csi_FR), ylab="CSI", xlab="Number of bad years in testing data",
     main="Grid point in France")
plot(c(RESU$nb_bad_y_test_US), c(RESU$csi_US), ylab="CSI", xlab="Number of bad years in testing data",
     main="Grid point in USA")
plot(c(RESU$nb_bad_y_test_CH), c(RESU$csi_CH), ylab="CSI", xlab="Number of bad years in testing data",
     main="Grid point in China")
