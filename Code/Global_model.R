# Data processing for northern hemisphere
# Author: Johannes Vogel

library(ncdf4)
library(rgdal)
library(raster)
library(dplyr)
library(glinternet)
library(abind)
library(foreach) 
library(doParallel)
library(InformationValue)
library(ROCR)
library(glmnet)
library(tictoc)

# path_to_NH_files <- "D:/Local/Data/Group project Como"
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"

nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
nh_variables <- lapply(1:length(nh_files),function(x){ncvar_get(nh_data[[x]])})
lati <- ncvar_get(nh_data[[1]],"lat")
long <- ncvar_get(nh_data[[1]],"lon")
lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})  
pix_num <- dim(nh_variables[[1]])[1]*dim(nh_variables[[1]])[2] # number of pixels

# Calculate threshold
threshold <- 0.05
yields <- matrix(data=nh_variables[[3]],nrow=320*76,ncol=1600) # reshape
# identical(nh_variables[[3]][180,50,],yields[49*320+180,]) # reshaping worked
low_yield <- sapply(1:pix_num,function(x) {quantile(yields[x,],threshold,na.rm=T)})
cy <- t(sapply(1:pix_num,function(x){ifelse(yields[x,]<low_yield[x],0,1)}))
# sum(yields[49*320+180,]<low_yield[49*320+180]) # check if correctly assigned: yes, 80/1600 are below threshold
# sum(cy[49*320+180,]) # Number is accurate: 1520 entries (95) without high crop loss

threshold_exclusion <- 0.025 # the exclusion of pixels is always based on the lowest threhold
low_yield_exclusion <- sapply(1:pix_num,function(x) {quantile(yields[x,],threshold_exclusion,na.rm=T)})
cy_exclusion <- t(sapply(1:pix_num,function(x){ifelse(yields[x,]<low_yield_exclusion[x],0,1)}))
row_sums_cy <- rowSums(cy_exclusion)
# sum(!is.na(row_sums_cy)) # all pixels with data: 955
table(row_sums_cy) # works fine for most pixel, but not for a few it is problematic
# problem with ties, for some pixel (e.g. sum(yields[4352,]==0)), there are so many 0s that the treshold is at 0
# which(row_sums_cy==1600)
plot(yields[4352,]) # one of the problematic pixels: contains a lot of zeros
# 901/995 pixel are all right


# Plotting ####

message('does not work') 
lon <- rep(long,length(lati)) # coordinates rearranged
lat <- rep(lati,each=length(long))
coord <- cbind(lon,lat)
dat_coord <- cbind(coord,gs_length[,1])
# dat_ras <- rasterFromXYZ(dat_coord) # gives error


border <- readOGR('D:/user/vogelj/Data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')	
# download from https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-countries-2/



## Growing season #####
gs_start <- nh_variables[[2]][,,1] # Julian date of start of season
gs_length <- matrix(data=nh_variables[[2]],nrow=320*76,ncol=1600)
max_gs_length <- apply(gs_length,1,max) # maximum growing season length
# num_years <- 4
# nh_variables[[1]][180,50,] # Growing season length example for one pixel
# sum(!is.na(max_gs_length))/sum(is.na(max_gs_length)) # only 4% are not NAs and therefore contain wheat yield data
total_yield <- rowSums(yields)
max_gs_length <- matrix(data=max_gs_length,nrow=320,ncol=76) # reshape to 2dim
total_yield <- matrix(data=total_yield,nrow=320,ncol=76) # reshape to 2dim

max_temp_aug <- matrix(nh_variables[[8]][,,1,],nrow=pix_num,ncol=1600)
# identical(mean(max_temp_aug[49*320+180,]),mean(nh_variables[[8]][180,50,1,])) # correctly reshaped
mean_max_temp_aug <- rowMeans(max_temp_aug)
mean_max_temp_aug <- matrix(data=mean_max_temp_aug,nrow=320,ncol=76)

message('latitude needs to be inverted')
gs_length_ras <- raster(t(max_gs_length[,76:1]), xmn=min(long), xmx=max(long), ymn=min(lati), ymx=max(lati), crs=CRS(projection(border)))
gs_start_ras <- raster(t(gs_start[,76:1]), xmn=min(long), xmx=max(long), ymn=min(lati), ymx=max(lati), crs=CRS(projection(border)))
total_yield_ras <- raster(t(total_yield[,76:1]), xmn=min(long), xmx=max(long), ymn=min(lati), ymx=max(lati), crs=CRS(projection(border)))
mean_max_temp_aug_ras <- raster(t(mean_max_temp_aug[,76:1]), xmn=min(long), xmx=max(long), ymn=min(lati), ymx=max(lati), crs=CRS(projection(border)))

projection(border) # projections are equal
projection(gs_length_ras)

# png(file=paste0(getwd(),"/rasterplot.png"),height=8,width=10,unit="cm",
    # pointsize = 6, bg = "white", res = 1000, restoreConsole = TRUE, type = "windows")
x11()
plot(gs_length_ras,asp=1);plot(border,add=T)
plot(gs_start_ras,asp=1);plot(border,add=T)
plot(total_yield_ras,asp=1);plot(border,add=T)
plot(mean_max_temp_aug_ras,asp=1);plot(border,add=T)
# dev.off()




# Lasso model building ####
###########################


# make it for one pixel first
# then a few
# then parallise it and run globally
# you only need to run 955 pixels, the rest is NA

# Create data table ####
# monthly_data <- abind(nh_variables[[5]][,,1:12,],nh_variables[[8]][,,1:12,],nh_variables[[9]][,,1:12,],nh_variables[[10]][,,1:12,],along=3)
monthly_data <- abind(nh_variables[[5]][,,1:12,],nh_variables[[8]][,,1:12,],nh_variables[[10]][,,1:12,],along=3)
monthly_data3dim <- array(data=monthly_data,dim=c(pix_num,dim(monthly_data)[3],1600))
# frpix <- monthly_data3dim[13601,,] # presumably our example pixel from France

# Exclude irrelevant (without wheat yield) pixels
# Two exclusion criteria: a) NA values or b) many 0 values (so that quantile threshold becomes 0)
# for a) Exclusion of 23325 pixels
# for b): Exclusion of 65 problematic pixels with frequent yields of 0 (65 pixels in the case of lowest threshold 0.025)
# 930 pixels remain
# pix_exclude <- c(which(is.na(row_sums_cy)),which(row_sums_cy==1600)) # excluded pixels
pix_in <- which(row_sums_cy!=1600) # included pixels (NA pixels and pixels with many 0s are directly excluded with only 1 command)
pix_in <- pix_in[-10] # exclude erroneous pixel (see "Erroneous pixel.R" for details)
Model_data_wheat <- monthly_data3dim[pix_in,,]
cy_wheat <- cy[pix_in,]

# z-score standardisation
# Model_data_stand <- sapply(1:dim(Model_data_wheat)[1], function(x){apply(Model_data_wheat[x,,],2,scale)})
Model_data_stand <- array(data=NA,dim=dim(Model_data_wheat))
message("Very important: apply inverts the dimensions!")
for (i in 1:dim(Model_data_wheat)[1]){
  Model_data_stand[i,,] <- t(apply(Model_data_wheat[i,,],1,scale)) # reinvert dimensions (inverted by apply)
}

cy_reshaped <- array(data=cy_wheat,dim=c(dim(cy_wheat)[1],1,1600))
# Model_data <- abind(cy_reshaped,Model_data_stand,along=2)
Model_data <- abind(cy_reshaped,Model_data_wheat,along=2)
colnames(Model_data) <- c("Yield","pr_Aug","pr_Sep","pr_Oct","pr_Nov","pr_Dec","pr_Jan","pr_Feb","pr_Mar","pr_Apr","pr_May","pr_Jun","pr_Jul","tmax_Aug","tmax_Sep","tmax_Oct","tmax_Nov","tmax_Dec","tmax_Jan","tmax_Feb","tmax_Mar","tmax_Apr","tmax_May","tmax_Jun","tmax_Jul","vpd_Aug","vpd_Sep","vpd_Oct","vpd_Nov","vpd_Dec","vpd_Jan","vpd_Feb","vpd_Mar","vpd_Apr","vpd_May","vpd_Jun","vpd_Jul")

# Check if dimensions are correct
# # monthly_data works fine
# monthly_data[180,50,1,1]
# nh_variables[[5]][180,50,1,1]
# monthly_data[180,50,2,1]
# nh_variables[[5]][180,50,2,1]
# 
# monthly_data[180,50,,1]
# nh_variables[[5]][180,50,1:12,1] 
# 
# # monthly_data3dim works fine
# monthly_data3dim[49*320+180,1,1]
# nh_variables[[5]][180,50,1,1]
# monthly_data3dim[49*320+180,2,1]
# nh_variables[[5]][180,50,2,1]
# 
# monthly_data3dim[49*320+180,,1]
# nh_variables[[5]][180,50,1:12,1]
# 
# # Model_data works fine
# Model_data[49*320+180,2,1]
# nh_variables[[5]][180,50,1,1]
# Model_data[49*320+180,3,1]
# nh_variables[[5]][180,50,2,1]
# 
# Model_data[49*320+180,,1]
# nh_variables[[5]][180,50,1:12,1]

# Split data into training and testing data set
set.seed(1994)
training_indices <- sort(sample(1:1600, size = floor(1600*0.6)))
testing_indices <- (1:1600)[-training_indices]
Training_Data <- Model_data[,,training_indices]
Testing_Data <- Model_data[,,testing_indices]


cv_fit_list <- vector(mode='list',length=dim(Model_data)[1])
# count <- 1
# wheat_pixels <- which(!is.na(row_sums_cy)) # pixels with wheat yield data
numLevels <- rep(1,times=dim(Model_data)[2]-1)
names(numLevels) <- colnames(Model_data[,2:length(colnames(Model_data)),])

tic()
no_cores <- detectCores() / 2 - 1
cl<-makeCluster(no_cores)
clusterEvalQ(cl, {
  library(glinternet)
  library(dplyr)
}) # parallelisation has own environment, therefore some packages and variables need be loaded again
registerDoParallel(cl)

# cv_fit <- foreach (i=1:dim(Model_data)[1],.combine=list) %dopar% {
cv_fit <- foreach (i=1:10,.multicombine=TRUE) %dopar% {
# for (i in 1:dim(Model_data)[1][1:3]){
  AllData <- as.data.frame(t(Model_data[i,,]))
  AllTraining_Data <- as.data.frame(t(Training_Data[i,,]))
  AllTesting_Data <- as.data.frame(t(Testing_Data[i,,]))
  
  
  X1_train <- AllTraining_Data[,2:dim(AllTraining_Data)[2]] # predictors
  y1_train <- AllTraining_Data[,1] # predictand
  # subset_Testing_Data <-  AllTesting_Data[,vec+1]
  x1_test <-  AllTesting_Data[,2:dim(AllTesting_Data)[2]] # predictors
  y1_test <- AllTesting_Data[,1] # predictand
  

  # Fit model
  # cv_fit <- glinternet.cv(X1_train, y1_train, numLevels,family = "binomial")
  glinternet.cv(X1_train, y1_train, numLevels,family = "binomial")
  # cv_fit_list[[i]] <- cv_fit
  # count <- count + 1
}
stopCluster(cl)
toc()

plot(cv_fit)

i_1Std <- which(cv_fit$lambdaHat1Std == cv_fit$lambda) # the preferential lambda (tuning parameter)

coefs <- coef(cv_fit$glinternetFit)[[i_1Std]]

coefs$mainEffects # model part without interactions
names(numLevels)[coefs$mainEffects$cont] # Main effect variables (without interactions)

coefs$interactions # model part with interactions pairs
names(numLevels)[coefs$interactions$contcont] # Main effect variables (with interactions)


# Assessing performance ####
sqrt(cv_fit$cvErr[[i_1Std]]) # root mean squared error (RMSE) on validation data


# predict.glinternet.cv(cv_fit,subset_Testing_Data,type="response") # does not work
# glinternet::predict.glinternet.cv(cv_fit,AllTesting_Data,type="response") # does not work
# predict.glinternet(cv_fit,AllTesting_Data,type="response") # does not work
mypred <- predict(cv_fit,x1_test,type="response")
fitted.results_bestglm <- ifelse(mypred > 0.5,1,0)


misClassError(y1_test,fitted.results_bestglm)

# Confusion matrix ####
# a) Confusion matrix manually calculated
obs_pred <- cbind(y1_test,fitted.results_bestglm)
tp <- sum(rowSums(obs_pred)==2)
tn <- sum(rowSums(obs_pred)==0)
fp <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
fn <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab1 <- matrix(c(tp,fn,fp,tn),nrow=2,ncol=2)
con_tab1b <- con_tab1
colnames(con_tab1) <- c('Actual TRUE','Actual FALSE');rownames(con_tab1) <- c('Predicted TRUE','Predicted FALSE')
con_tab1[,1] <- con_tab1[,1]/ sum(y1_test==1)
con_tab1[,2] <-con_tab1[,2]/sum(y1_test==0)
# b) Confusion matrix from package InformationValue
con_tab <- InformationValue::confusionMatrix(y1_test,fitted.results_bestglm)

# Sensitivity and Specificity
# a) Manually calculated
spec <- tn/(tn+fp)
sens <- tp/(tp+fn)
# b) Using package caret
caret::sensitivity(data=as.factor(fitted.results_bestglm),reference=as.factor(y1_test),positive="1",negative="0")
caret::specificity(data=as.factor(fitted.results_bestglm),reference=as.factor(y1_test),positive="1",negative="0")
# c) using package InformationValue
InformationValue::sensitivity(y1_test,fitted.results_bestglm)
InformationValue::specificity(y1_test,fitted.results_bestglm)


# ROC ####
pr <- prediction(mypred, y1_test)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
plotROC(actuals=y1_test,predictedScores=fitted.results_bestglm)
auc2 <- auc(y1_test,fitted.results_bestglm)
auc <- performance(pr, measure = "auc")
auc@y.values[[1]]

# Calculate Sensitivity and Specificity using performance command
spec2 <- performance(pr, measure = "spec")
# spec2b <- performance(pr, measure = "tnr") # equivalent
plot(spec2) # you can see that at cutoff 0.5 it is equal to my result, however I don't know how to extract the value
# sens2 <- performance(pr, measure = "sens") # does now work
sens2 <- performance(pr, measure = "tpr")
sens_spec <- performance(pr, measure="sens", x.measure="spec")
plot(sens_spec) # inverted AUC



# Problems
# 17 months
# ties for some pixels
# no model output

