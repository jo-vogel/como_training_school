###########################################
##### Various plot for one grid point #####
##### Ridge and Lasso                 #####
###########################################

# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

#Coordinate of my village
lon_maringes <- 4.3
lat_maringes <- 45.7


#threshold for bad yields in c(0.025,0.05,0.1)
threshold <- 0.05

##### Initialisation, librairies, data #####

library(ncdf4);library(rgdal);library(raster);library(RColorBrewer);library(viridis)
library(maps);library(mapdata);library(ggplot2)
library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr)
library(foreach);library(doParallel)
library(tictoc)


#####################################
##### Monthly standardised data #####
#####################################

# Pauline:
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
output_path <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Images/Model_Comparison/"

source('./Code/Lasso_interact_global_preparation.R') # load necessary files


# which method? model_name in c("Ridge", "Lasso", "ElasticNet")
model_name <- "ElasticNet"
stopifnot(model_name %in% c("Ridge", "Lasso", "ElasticNet"))


#which segregation threshold for the model?
segreg_th <- 0.5


# Load the fitted model ####
############################

studied_pix <- which.min((sqrt((coord_subset[,1]-lon_maringes)^2 + (coord_subset[,2]-lat_maringes)^2)))

if(model_name == "ElasticNet"){
  load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                    "Elastic-net_lambdamin_alpha05_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))
  model_mypix <- fit_pix_lambdamin[[studied_pix]]
} else {
  load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                    model_name,"_lambdamin_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))
  
  if(model_name == "Ridge"){
    model_mypix <- ridge_model_lambdamin[[studied_pix]]
  }
  if(model_name == "Lasso"){
    model_mypix <- lasso_model_lambdamin[[studied_pix]]
  }
  
} #end if else


#Meteo var for my pixel

Datastand_mypix <- Model_data_stand[studied_pix,,]
meanData_mypix <- apply(Datastand_mypix, MARGIN = 1, FUN = mean)
meanData_mypix <- meanData_mypix[!is.na(meanData_mypix)][-1]
Data_low5pctl_mypix <- apply(Datastand_mypix, MARGIN = 1, FUN = quantile, probs=0.05, na.rm=T)
Data_low5pctl_mypix <- Data_low5pctl_mypix[!is.na(Data_low5pctl_mypix)][-1]
Data_up95pctl_mypix <- apply(Datastand_mypix, MARGIN = 1, FUN = quantile, probs=0.95, na.rm=T)
Data_up95pctl_mypix <- Data_up95pctl_mypix[!is.na(Data_up95pctl_mypix)][-1]



#model perf for my pix

coefs <- coefficients(model_mypix)

my_pred <- predict(model_mypix, as.matrix(x1_test_list[[studied_pix]]),type="response")

misClassError(actuals = as.matrix(y1_test_list[[studied_pix]]),
              predictedScores=my_pred, threshold = segreg_th)


InformationValue::specificity(as.matrix(y1_test_list[[studied_pix]]),
                              my_pred, threshold = segreg_th)


InformationValue::confusionMatrix(as.matrix(y1_test_list[[studied_pix]]),
                                  my_pred, threshold = segreg_th)



#Plots


VAR_values <- matrix(rep(meanData_mypix, each=100), nrow = 100)
colnames(VAR_values) <- colnames(x1_test_list[[studied_pix]])

plot_2var <- function(name_var1, name_var2, xval=seq(-1,1, length.out = 100), lgd_pos = "bottomleft"){
  VAR_values_chge <- VAR_values
  VAR_values_chge[,name_var1] <- xval
  
  #Modify the variables you want to play with
  VAR_values_chge_lowvar2 <- VAR_values_chge
  VAR_values_chge_lowvar2[,name_var2] <- rep(Data_low5pctl_mypix[name_var2], 100)
  VAR_values_chge_highvar2 <- VAR_values_chge
  VAR_values_chge_highvar2[,name_var2] <- rep(Data_up95pctl_mypix[name_var2], 100)
  
  pred_highvar2 <- predict(model_mypix, VAR_values_chge_highvar2, type="response")
  pred_lowvar2 <- predict(model_mypix, VAR_values_chge_lowvar2, type="response")
  
  lim_y <- c(min(c(pred_highvar2,pred_lowvar2)),max(c(pred_highvar2,pred_lowvar2)))
  
  plot(xval,
       pred_highvar2,
       col="green",
       ylab="Proba good yield",
       xlab=paste("rescaled", name_var1),
       ylim=lim_y)
  points(xval,
         pred_lowvar2,
         col="orange")
  points(xval,
         predict(model_mypix, VAR_values_chge,
                 type="response"))
  
  legend(lgd_pos, pch = 1, col=c("green", "black", "orange"),
         legend = c(paste(name_var2, "95th perc"),
                    paste(name_var2, "mean"),
                    paste(name_var2, "5th perc")))
  abline(h=0.5, lty=2)
}#end plot_2var


source("./Code/get_first_coeff_function.R")
first4_coef <- get_firstcoeffs(coefs, nb_of_coeff = 4)

plot_2var(name_var1 = paste(first4_coef[1,1],"_",first4_coef[1,2], sep = ""),
          name_var2 = paste(first4_coef[3,1],"_",first4_coef[3,2], sep = ""))
