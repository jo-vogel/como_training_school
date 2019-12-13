##############################################################################
###########          Ridge and LASSO regression                    ###########
###########   on global standardised crop model data               ###########
##############################################################################


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

model_name <- "Elastic-net"


#threshold for bad yields in c(0.025,0.05,0.1)
threshold <- 0.05

#which segregation threshold for the model?
segreg_th <- 0.5

##### Initialisation, librairies, data #####

library(ncdf4);library(rgdal);library(raster);library(RColorBrewer);library(viridis)
library(maps);library(mapdata);library(ggplot2)
library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr)
library(foreach);library(iterators);library(parallel);library(doParallel)
library(tictoc)


#############################
##### Standardised data #####
#############################


# recreate the training and testing dataset ####
################################################

message('Adjust the path accordingly.')
# Laptop Pauline:
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
output_path <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Images/Model_Comparison/"

# Johannes:
# path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
# output_path <- "D:/user/vogelj/Group_project/Output/Plots/"

source('./Code/Lasso_interact_global_preparation.R') # load necessary files



##### Run the model ####

#Value for alpha, in (0;1)
alpha_elastic <- 0.5
set.seed(2019)

crossval <- list()

fit_pix_lambdamin <- list()
fit_pix_lambda1se <- list()

for (pixel in 1:dim(Model_data_stand)[1]) {
  print(paste("Pixel",pixel,"out of",dim(Model_data_stand)[1]))
  tic()
  crossval[[pixel]] <- cv.glmnet(x = as.matrix(x1_train_list[[pixel]]),
                                 y = as.matrix(y1_train_list[[pixel]]),
                                 family = "binomial", alpha = alpha_elastic, nfolds = 10)
  fit_pix_lambdamin[[pixel]] <- glmnet(x = as.matrix(x1_train_list[[pixel]]),
                                       y = as.matrix(y1_train_list[[pixel]]),
                                       family = "binomial", alpha = alpha_elastic,
                                       lambda = crossval[[pixel]]$lambda.min)
  fit_pix_lambda1se[[pixel]] <- glmnet(x = as.matrix(x1_train_list[[pixel]]),
                                       y = as.matrix(y1_train_list[[pixel]]),
                                       family = "binomial", alpha = alpha_elastic,
                                       lambda = crossval[[pixel]]$lambda.1se)
  toc()
}#end for pixel


save(crossval, file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/cv",
                            model_name,"_alpha",str_pad(alpha_elastic*10, 2, pad = "0"),"_threshbadyield",
                            str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))
save(fit_pix_lambdamin, file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                                     model_name,"_lambdamin_alpha",str_pad(alpha_elastic*10, 2, pad = "0"),"_threshbadyield",
                                     str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))

save(fit_pix_lambda1se, file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                                    model_name,"_lambda1se_alpha",str_pad(alpha_elastic*10, 2, pad = "0"),"_threshbadyield",
                                     str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))

# Cross validation to find the good alpha


# ELASTIC NET WITH 0 < ALPHA < 1


no_cores <- detectCores() / 2 - 1

# registerDoParallel(cores = no_cores)

cl<-makeCluster(no_cores)
clusterEvalQ(cl, {
  library(glmnet)
  library(dplyr)
}) # parallelisation has own environment, therefore some packages and variables need be loaded again
registerDoParallel(cl)

tic()


# cv_elastic_model_select <- foreach (pixel=1:dim(Model_data_stand)[1]) %dopar% {
set.seed(2019)
cv_elastic_model_select5<- foreach (pixel=757) %dopar% {
  
  a_to_test <- seq(0.1, 0.9, 0.05)
  search_list_1se <- matrix(data = NA, nrow = length(a_to_test), ncol = 3)
  search_list_min <- matrix(data = NA, nrow = length(a_to_test), ncol = 3)
  colnames(search_list_1se) <- c("cvm", "lambda.1se", "alpha")
  colnames(search_list_min) <- c("cvm", "lambda.min", "alpha")

  
  for (a in a_to_test) {
    crossval <- cv.glmnet(x = as.matrix(x1_train_list[[pixel]]),
                          y = as.matrix(y1_train_list[[pixel]]),
                          family = "binomial", alpha = a, nfolds = 10)
    search_list_1se[which(a_to_test==a),] <- c(crossval$cvm[crossval$lambda == crossval$lambda.1se],
                                               crossval$lambda.1se, a)
    search_list_min[which(a_to_test==a),] <- c(crossval$cvm[crossval$lambda == crossval$lambda.min],
                                               crossval$lambda.min, a)
  }#end for a
  
  min_cvm_lamda.1se <- search_list_1se[which.min(search_list_1se[,"cvm"]),]
  min_cvm_lamda.min <- search_list_min[which.min(search_list_min[,"cvm"]),]
  return(list(lambda1se = list(fit = glmnet(x = as.matrix(x1_train_list[[pixel]]),
                                            y = as.matrix(y1_train_list[[pixel]]),
                                            family = "binomial", lambda = min_cvm_lamda.1se["lambda.1se"],
                                            alpha = min_cvm_lamda.1se["alpha"]),
                               alph = min_cvm_lamda.1se["alpha"]),
              lambdamin = list(fit = glmnet(x = as.matrix(x1_train_list[[pixel]]),
                                            y = as.matrix(y1_train_list[[pixel]]),
                                            family = "binomial", lambda = min_cvm_lamda.min["lambda.min"],
                                            alpha = min_cvm_lamda.min["alpha"]),
                               alph = min_cvm_lamda.min["alpha"])))
}
cbind("lambda1se"=c(cv_elastic_model_select1[[1]]$lambda1se$alph,
                    cv_elastic_model_select2[[1]]$lambda1se$alph,
                    cv_elastic_model_select3[[1]]$lambda1se$alph,
                    cv_elastic_model_select4[[1]]$lambda1se$alph,
                    cv_elastic_model_select5[[1]]$lambda1se$alph),
      "lambdamin"=c(cv_elastic_model_select1[[1]]$lambdamin$alph,
                    cv_elastic_model_select2[[1]]$lambdamin$alph,
                    cv_elastic_model_select3[[1]]$lambdamin$alph,
                    cv_elastic_model_select4[[1]]$lambdamin$alph,
                    cv_elastic_model_select5[[1]]$lambdamin$alph))

toc()





# a_to_test <- seq(0.1, 0.9, 0.05)
# 
# 
# search_list <- matrix(data = NA, nrow = length(a_to_test), ncol = 3)
# colnames(search_list) <- c("cvm", "lambda.1se", "alpha")
# 
# tic()
# 
# for (a in a_to_test) {
#   crossval <- cv.glmnet(x = as.matrix(x1_train_list[[pixel]]),
#                         y = as.matrix(y1_train_list[[pixel]]),
#                         family = "binomial", alpha = a, nfolds = 10)
#   search_list[which(a_to_test==a),] <- c(crossval$cvm[crossval$lambda == crossval$lambda.1se],
#                                          crossval$lambda.1se, a)
# }#end for a
# toc()
# 
# 
# 
# min_cvm <- search_list[which.min(search_list[,"cvm"]),]
# result <- glmnet(x = as.matrix(x1_train_list[[pixel]]),
#                  y = as.matrix(y1_train_list[[pixel]]),
#                  family = "binomial", lambda = min_cvm["lambda.1se"], alpha = min_cvm["alpha"])
# coef(result)




##### Load the model ####
alpha_elastic <- 0.5

#which lambda?
lambda_VALS <- c("lambda.min", "lambda.1se")
lambda_val <- lambda_VALS[2]

if(lambda_val == "lambda.min"){
  load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                  model_name,"_lambdamin_alpha",str_pad(alpha_elastic*10, 2, pad = "0"),"_threshbadyield",
                  str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))
  model_fit <- fit_pix_lambdamin
  rm(fit_pix_lambdamin)
}

if(lambda_val == "lambda.1se"){
  load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                  model_name,"_lambda1se_alpha",str_pad(alpha_elastic*10, 2, pad = "0"),"_threshbadyield",
                  str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))
  model_fit <- fit_pix_lambda1se
  rm(fit_pix_lambda1se)
}



# Model performance ###
#######################
#Start with Lambda min

test_length <- length(model_fit)

coefs <- lapply(1:test_length, function(x){coef(model_fit[[x]])})



mypred <- lapply(1:test_length, function(x){predict(model_fit[[x]],
                                                    as.matrix(x1_test_list[[x]]),type="response")})

#which segregation threshold for the model?
segreg_th <- 0.5

fitted.results_model <- lapply(1:test_length, function(x){ifelse(mypred[[x]] > segreg_th,1,0)})

mis_clas_err <- sapply(1:test_length, function(x){misClassError(actuals = as.matrix(y1_test_list[[x]]),
                                                                predictedScores=mypred[[x]],
                                                                threshold = segreg_th)})

con_tab <-  lapply(1:test_length, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                fitted.results_model[[x]],
                                                                                threshold = segreg_th)})

speci <- sapply(1:test_length, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                         fitted.results_model[[x]],
                                                                         threshold = segreg_th)})

# CSI #
csi <- numeric()
for(pix in 1:length(coefs)){
  csi[pix] <- con_tab[[pix]]["0","0"]/(con_tab[[pix]]["0","0"] + con_tab[[pix]]["1","0"] + con_tab[[pix]]["0","1"])
  if(is.na(con_tab[[pix]]["0","0"])){
    csi[pix] <- 0
  }
}#end for pix

# Plots ####
world <- map_data("world")

# Plot specificity error ####

DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = speci)

ggplot(data = DF_speci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=speci)) +
  scale_color_gradient2(limits=c(0,1), midpoint = 0.5,
                        low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Specif.",
       title = paste("Monthly data: specificity, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


# plot CSI=(hits)/(hits + misses + false alarm) ###
DF_csi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], csi = csi)

ggplot(data = DF_csi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=csi)) +
  scale_color_gradient2(limits=c(0,1), midpoint = 0.5,
                        low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="csi",
       title = paste("Seasonal data: critical success index, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)



# Plot number of variables kept ###
coeff_kep <- numeric()

for (pix in 1:965) {
  # coeff_kep[pix] <- length(coefs[[pix]]$mainEffects$cont)
  coeff_kep[pix] <- sum(coefs[[pix]][row.names(coefs[[pix]])!="(Intercept)"]!=0)
  
}#end for pix



# Plot 

DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = coeff_kep)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=coeff_kep)) +
  scale_color_gradient(limits=c(1,39),
                       low = "pink", high = "darkblue") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of var.",
       title = paste("Number of variables kept, simple",model_name,"regression"),
       subtitle = lambda_val)+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
