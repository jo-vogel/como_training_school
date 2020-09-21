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
threshold <- 0.05

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
low_yield <- apply(Yield, MARGIN = 1, FUN=quantile, probs=threshold, na.rm=T)
cy <- t(sapply(1:pix_num,function(x){ifelse(Yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield

cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,nb_years))

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


##### Split data into training and testing data set #####
vec <- 1:nb_years

years_with_na <- vector("logical",length=pix_num)
for (i in 1:pix_num){
  years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
}

training_indices <- vector("list",length=pix_num)
testing_indices <- vector("list",length=pix_num)





# Parameters of the splitting
seed=1994
train_size <- 70

set.seed(seed)
for (x in 1:pix_num) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((nb_years-length(na_time[[x]]))*(train_size/100))))
    testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
  } else {
    training_indices[[x]] <- sort(sample(1:nb_years, size = floor(nb_years*(train_size/100))))
    testing_indices[[x]] <- (1:nb_years)[-training_indices[[x]]]    
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

# extract several the grippoints ####
load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/final_889pixels_coords.Rdata")


FR_pix <- which.min((final_pixels_coord$longitude-1.1)^2+(final_pixels_coord$latitude-47.7)^2)
US_pix <- which.min((final_pixels_coord$longitude+90)^2+(final_pixels_coord$latitude-44.3)^2)
CH_pix <- which.min((final_pixels_coord$longitude-118.1)^2+(final_pixels_coord$latitude-30.8)^2)

final_pixels_coord[FR_pix,];final_pixels_coord[US_pix,];final_pixels_coord[CH_pix,]

pixels_number <- c(FR_pix, US_pix, CH_pix)

meteovar_train_FR_complete <- x1_train_list[[FR_pix]] ; meteovar_test_FR_complete <- x1_test_list[[FR_pix]]
yield_train_FR_complete <- y1_train_list[[FR_pix]] ; yield_test_FR_complete <- y1_test_list[[FR_pix]]

meteovar_train_US_complete <- x1_train_list[[US_pix]] ; meteovar_test_US_complete <- x1_test_list[[US_pix]]
yield_train_US_complete <- y1_train_list[[US_pix]] ; yield_test_US_complete <- y1_test_list[[US_pix]]

meteovar_train_CH_complete <- x1_train_list[[CH_pix]] ; meteovar_test_CH_complete <- x1_test_list[[CH_pix]]
yield_train_CH_complete <- y1_train_list[[CH_pix]] ; yield_test_CH_complete <- y1_test_list[[CH_pix]]


ratio_to_test <- c(1, 3/4, 2/3, 1/2, 1/3, 1/4) #1/4 is the limit to have enough bad years in testing data (we really don't have a lot of bad years with our definition of bad years)

RESULTS <- list()

for (GP in 1:3) { #GP <- 1
  
  pixel_number <- pixels_number[GP]
  
  meteovar_train_complete <- x1_train_list[[pixel_number]] ; meteovar_test_complete <- x1_test_list[[pixel_number]]
  yield_train_complete <- y1_train_list[[pixel_number]] ; yield_test_complete <- y1_test_list[[pixel_number]]
  
  for (ratio in 1:length(ratio_to_test)) { #ratio <- 2
    set.seed(1994)
    ind_sub_train <- sample(1:1120, size = floor(1120*ratio_to_test[ratio]), replace = F)
    ind_sub_test <- sample(1:420, size = floor(420*ratio_to_test[ratio]), replace = F)
    
    meteovar_train_sub <- meteovar_train_complete[ind_sub_train,] ; meteovar_test_sub <- meteovar_test_complete[ind_sub_test,]
    yield_train_sub <- yield_train_complete[ind_sub_train] ; yield_test_sub <- yield_test_complete[ind_sub_test]
    
    var_pix <- as.matrix(meteovar_train_sub)
    yield_pix <- as.matrix(yield_train_sub)
    which_na_xtrm <- which(is.na(var_pix[,1]))
    
    if(length(which_na_xtrm)>0){ #remove potential missing data in extreme indices
      model_cv_fitting <- cv.glmnet(x = var_pix[-which_na_xtrm,], y = yield_pix[-which_na_xtrm,],
                                    family = "binomial", alpha = 1, nfolds = 10)
    } else {
      model_cv_fitting <- cv.glmnet(x = var_pix, y = yield_pix,
                                    family = "binomial", alpha = 1, nfolds = 10)
    }#end if else na
    
    if(length(which_na_xtrm)>0){
      
      lasso_model_lambda1se <- glmnet(x = var_pix[-which_na_xtrm,], y = yield_pix[-which_na_xtrm,],
                                      family = "binomial", alpha = 1, lambda = model_cv_fitting$lambda.1se)
    } else {
      
      lasso_model_lambda1se <- glmnet(x = var_pix, y = yield_pix,
                                      family = "binomial", alpha = 1, lambda = model_cv_fitting$lambda.1se)
    }#end if else there exists na in extremes
    
    
    segreg_th <- 0.6582418
    
    
    coeff <- coefficients(lasso_model_lambda1se)
    
    mypred <- predict(lasso_model_lambda1se, as.matrix(meteovar_test_sub),type="response")
    
    fitted_bad_yield <- ifelse(mypred > segreg_th,1,0)
    
    speci <- InformationValue::specificity(actuals = as.matrix(yield_test_sub),
                                           predictedScores = fitted_bad_yield, threshold = segreg_th)
    sensi <- InformationValue::sensitivity(actuals = as.matrix(yield_test_sub),
                                           predictedScores = fitted_bad_yield, threshold = segreg_th)
    
    con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(yield_test_sub),
                                                 predictedScores = fitted_bad_yield, threshold = segreg_th)
    csi <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
    if(is.na(con_tab["0","0"])){
      csi <- 0
    }
    
    nb_training_years <- length(ind_sub_train)
    nb_testing_years <- length(ind_sub_test)
    
    nb_bad_years_train <- sum(abs(yield_train_sub-1))
    nb_bad_years_test <- sum(abs(yield_test_sub-1))
    
    results_GP <- list(coeff=coeff, speci=speci, sensi=sensi, csi=csi,
                       nb_bad_years_train=nb_bad_years_train, nb_bad_years_test=nb_bad_years_test)
    
    if(ratio==1){
      results_GP_1 <- results_GP
    }
    if(ratio==2){
      results_GP_34 <- results_GP
    }
    if(ratio==3){
      results_GP_23 <- results_GP
    }
    if(ratio==4){
      results_GP_2 <- results_GP
    }
    if(ratio==5){
      results_GP_3 <- results_GP
    }
    if(ratio==6){
      results_GP_4 <- results_GP
    }
    
  }#end for ratio
  
  RESULTS[[GP]] <- list(results_GP_1=results_GP_1, results_GP_34=results_GP_34, results_GP_23=results_GP_23,
                        results_GP_2=results_GP_2, results_GP_3=results_GP_3, results_GP_4=results_GP_4)
  
}#end for GP








##### Simulation to assess robustness of csi changing with length time serie #####

SIMU <- list()
ratio_to_test_2 <- c(3/4, 2/3, 1/2, 1/3, 1/4)
name_ratios <- c("3/4", "2/3", "1/2", "1/3", "1/4")
n_simu <- 100

set.seed(1994)

for (GP in 1:3) {
  
  
  matrix_csi <- matrix(NA, nrow = n_simu, ncol = length(ratio_to_test_2))
  colnames(matrix_csi) <- name_ratios
  matrix_nb_bad_years_test <- matrix(NA, nrow = n_simu, ncol = length(ratio_to_test_2))
  colnames(matrix_nb_bad_years_test) <- name_ratios
  
  for (simu in 1:n_simu) {
    pixel_number <- pixels_number[GP]
    
    meteovar_train_complete <- x1_train_list[[pixel_number]] ; meteovar_test_complete <- x1_test_list[[pixel_number]]
    yield_train_complete <- y1_train_list[[pixel_number]] ; yield_test_complete <- y1_test_list[[pixel_number]]
    
    for (ratio in 1:length(ratio_to_test_2)) { #ratio <- 2
      ind_sub_train <- sample(1:1120, size = floor(1120*ratio_to_test_2[ratio]), replace = F)
      ind_sub_test <- sample(1:420, size = floor(420*ratio_to_test_2[ratio]), replace = F)
      
      meteovar_train_sub <- meteovar_train_complete[ind_sub_train,] ; meteovar_test_sub <- meteovar_test_complete[ind_sub_test,]
      yield_train_sub <- yield_train_complete[ind_sub_train] ; yield_test_sub <- yield_test_complete[ind_sub_test]
      
      var_pix <- as.matrix(meteovar_train_sub)
      yield_pix <- as.matrix(yield_train_sub)
      which_na_xtrm <- which(is.na(var_pix[,1]))
      
      if(length(which_na_xtrm)>0){ #remove potential missing data in extreme indices
        model_cv_fitting <- cv.glmnet(x = var_pix[-which_na_xtrm,], y = yield_pix[-which_na_xtrm,],
                                      family = "binomial", alpha = 1, nfolds = 10)
      } else {
        model_cv_fitting <- cv.glmnet(x = var_pix, y = yield_pix,
                                      family = "binomial", alpha = 1, nfolds = 10)
      }#end if else na
      
      if(length(which_na_xtrm)>0){
        
        lasso_model_lambda1se <- glmnet(x = var_pix[-which_na_xtrm,], y = yield_pix[-which_na_xtrm,],
                                        family = "binomial", alpha = 1, lambda = model_cv_fitting$lambda.1se)
      } else {
        
        lasso_model_lambda1se <- glmnet(x = var_pix, y = yield_pix,
                                        family = "binomial", alpha = 1, lambda = model_cv_fitting$lambda.1se)
      }#end if else there exists na in extremes
      
      
      segreg_th <- 0.6582418
      
      
      coeff <- coefficients(lasso_model_lambda1se)
      
      mypred <- predict(lasso_model_lambda1se, as.matrix(meteovar_test_sub),type="response")
      
      fitted_bad_yield <- ifelse(mypred > segreg_th,1,0)
      
      speci <- InformationValue::specificity(actuals = as.matrix(yield_test_sub),
                                             predictedScores = fitted_bad_yield, threshold = segreg_th)
      sensi <- InformationValue::sensitivity(actuals = as.matrix(yield_test_sub),
                                             predictedScores = fitted_bad_yield, threshold = segreg_th)
      
      con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(yield_test_sub),
                                                   predictedScores = fitted_bad_yield, threshold = segreg_th)
      csi <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
      if(length(con_tab["0","0"])==0){
        csi <- 0
      }
      
      nb_training_years <- length(ind_sub_train)
      nb_testing_years <- length(ind_sub_test)
      
      nb_bad_years_train <- sum(abs(yield_train_sub-1))
      nb_bad_years_test <- sum(abs(yield_test_sub-1))
      
      
      matrix_csi[simu,ratio] <- csi
      matrix_nb_bad_years_test[simu,ratio] <- nb_bad_years_test
      
    }#end for ratio
  }#end for simu
  
  SIMU[[GP]] <- list(matrix_csi = matrix_csi, matrix_nb_bad_years_test)
}#end for GP

save(SIMU, file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/FR_US_CH_simulation_CSI.Rdata")


boxplot(SIMU[[1]]$matrix_csi, main = "CSI for 100 subsamples\ngridpoint in FR",
        xlab="nb on years availabale\n(relative size compared to 1600 years)", ylim=c(0,1))
abline(h=RESULTS[[1]]$results_GP_1$csi, col="chartreuse4", lwd=2)
text(x=4.6, y=RESULTS[[1]]$results_GP_1$csi+0.03, "1600 years", col="chartreuse4")
boxplot(SIMU[[2]]$matrix_csi, main = "CSI for 100 subsamples\ngridpoint in US", ylim=c(0,1),
        xlab="nb on years availabale\n(relative size compared to 1600 years)")
abline(h=RESULTS[[2]]$results_GP_1$csi, col="chartreuse4", lwd=2)
text(x=4.6, y=RESULTS[[2]]$results_GP_1$csi+0.03, "1600 years", col="chartreuse4")
boxplot(SIMU[[3]]$matrix_csi, main = "CSI for 100 subsamples\ngridpoint in CH",
        xlab="nb on years availabale\n(relative size compared to 1600 years)", ylim=c(0,1))
abline(h=RESULTS[[3]]$results_GP_1$csi, col="chartreuse4", lwd=2)
text(x=4.6, y=RESULTS[[3]]$results_GP_1$csi+0.03, "1600 years", col="chartreuse4")



load(file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/FR_US_CH_simulation_CSI.Rdata")


#plot csi as a function of nbof bad years in the testing data
pix_name <- c("FR", "US", "CH")
ratio_name <- c("3/4","2/3","1/2","1/3","1/4")
ratios <- c(3/4, 2/3, 1/2, 1/3, 1/4)
for (pixel in 1:3) {
  
  plot(c(SIMU[[pixel]][[2]]), c(SIMU[[pixel]]$matrix_csi),
       xlab="Nb bad years in testing data",
       ylab="CSI", main=paste(pix_name[pixel])
  )
  
  print(paste("Correlation test between CSI and nb years in testing data in", pix_name[pixel]))
  print(cor.test(c(SIMU[[pixel]][[2]]), c(SIMU[[pixel]]$matrix_csi))$p.value)
  # 
  # for (ratio in 1:5) {
  #   plot(SIMU[[pixel]][[2]][,ratio], SIMU[[pixel]]$matrix_csi[,ratio],
  #        xlab="Nb bad years in testing data",
  #        ylab="CSI", main=paste(pix_name[pixel], "\nSubsample size", floor(ratios[ratio]*1600))
  #        )
  # }#end for ratio
}# end for pixel

