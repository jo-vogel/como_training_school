library(ncdf4)
library(rgdal)
library(raster)
library(dplyr)
library(glinternet)
library(abind)
library(foreach) 
library(doParallel)


# path_to_NH_files <- "D:/Local/Data/Group project Como"
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"

nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
nh_variables <- lapply(1:length(nh_files),function(x){ncvar_get(nh_data[[x]])})
lati <- ncvar_get(nh_data[[1]],"lat")
long <- ncvar_get(nh_data[[1]],"lon")
lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})  


# ## Adjust growing season #####
gs_start <- nh_variables[[2]][,,1] # Julian date of start of season
gs_length <- matrix(data=nh_variables[[2]],nrow=320*76,ncol=1600)
max_gs_length <- apply(gs_length,1,max) # maximum growing season length
# num_years <- 4
# nh_variables[[1]][180,50,]
# nh_variables[[2]][180,50,1]
# sum(!is.na(max_gs_length))/sum(is.na(max_gs_length)) # only 4% are not NAs and therefore contain wheat yield data


# Calculate threshold
threshold <- 0.05
yields <- matrix(data=nh_variables[[3]],nrow=320*76,ncol=1600) # reshape
# identical(nh_variables[[3]][180,50,],yields[49*320+180,]) # reshaping worked
low_yield <- sapply(1:24320,function(x) {quantile(yields[x,],threshold,na.rm=T)})
# cy <- ifelse(yields<low_yield,0,1)
# myfun <- function(x){ifelse(yields<low_yield,0,1)}
# cy <- apply(yields,1,function(x){ifelse(yields<low_yield,0,1)})
cy <- t(sapply(1:24320,function(x){ifelse(yields[x,]<low_yield[x],0,1)}))
# cy <- sapply(1:24320,function(x){ifelse(yields[x,]<low_yield[x],0,1)})

sum(yields[49*320+180,]<low_yield[49*320+180])
sum(cy[49*320+180,])
a <- rowSums(cy)
sum(!is.na(a)) # all pixels with data: 955
ind <- which(!is.na(a))
table(a) # für die meisten Pixel funktioniert ist, für manche wenige nicht
# problem with ties, for some pixel (e.g. sum(yields[4352,]==0)), there are so many 0s that the treshold is at 0
which(a==1600)
plot(yields[4352,])
# 901/995 pixel are all right

# plot it

# does not work
# coordinates
lon <- rep(long,length(lati))
lat <- rep(lati,each=length(long))
coord <- cbind(lon,lat)
# 
# dat_coord <- cbind(coord,gs_length[,1])
# dat_ras <- rasterFromXYZ(dat_coord)


border <- readOGR('D:/user/vogelj/Data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')	
# download from https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-countries-2/

gs_length_ras <- raster(t(max_gs_length), xmn=min(long), xmx=max(long), ymn=min(lati), ymx=max(lati), crs=CRS(projection(border)))
gs_start_ras <- raster(t(gs_start), xmn=min(long), xmx=max(long), ymn=min(lati), ymx=max(lati), crs=CRS(projection(border)))

projection(border) # projections are equal
projection(gs_length_ras)

# png(file=paste0(getwd(),"/rasterplot.png"),height=8,width=10,unit="cm",
    # pointsize = 6, bg = "white", res = 1000, restoreConsole = TRUE, type = "windows")
x11()
plot(gs_length_ras,asp=1);plot(border,add=T)
plot(gs_start_ras,asp=1);plot(border,add=T)
# dev.off()





# Lasso model building ####

# make it for one pixel first
# then a few
# then parallise it and run globally
# you only need to run 955 pixel, the rest is NA

message('check if dimensions are correct')
# monthly_data <- array(data=c(nh_variables[[5]][,,1:12,],nh_variables[[8]][,,1:12,],nh_variables[[9]][,,1:12,]),dim=c(24320,1600,36))
# Model_data <- abind(cy,monthly_data) # Replace actual crop yield by binary info on fail/success
Data <- array(data=c(cy,nh_variables[[5]][,,1:12,],nh_variables[[8]][,,1:12,],nh_variables[[9]][,,1:12,]),dim=c(24320,1600,37))
Model_data <- Data                
      
# Split data into training and testing data set
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Model_data[,training_indices,]
Testing_Data <- Model_data[,testing_indices,]

all_pixels <- 24320
wheat_pixels <- ind # pixels with wheat yield data
cv_fit_list <- vector(mode='list',length=length(wheat_pixels))
# count <- 1


no_cores <- detectCores() / 2 - 1
cl<-makeCluster(no_cores)
clusterEvalQ(cl, library(glinternet),library(dplyr)) # parallelisation hat eigenes environment, daher m?ssen packages und bestimmte Variablen erneut geladen werden
registerDoParallel(cl)

cv_fit <- foreach (i=wheat_pixels,.combine=list) %dopar% { 
# for (i in 1:wheat_pixels){
  AllData <- as.data.frame(Model_data[i,,])
  AllTraining_Data <- as.data.frame(Training_Data[i,,])
  AllTesting_Data <- as.data.frame(Testing_Data[i,,])
  
  
  X1_train <- AllTraining_Data[,2:length(AllTraining_Data)] # predictors
  y1_train <- AllTraining_Data[,1] # predictand
  # subset_Testing_Data <-  AllTesting_Data[,vec+1]
  x1_test <-  AllTesting_Data[,2:length(AllTesting_Data)] # predictors
  y1_test <- AllTesting_Data[,1] # predictand
  
  # numLevels <- AllTraining_Data[,vec+1] %>% sapply(nlevels)
  numLevels <- AllTraining_Data[,2:length(AllTraining_Data)] %>% sapply(nlevels)
  numLevels[numLevels==0] <- 1 # set to 1 for continuous variables
  
  # Fit model
  cv_fit <- glinternet.cv(X1_train, y1_train, numLevels,family = "binomial")
  glinternet.cv(X1_train, y1_train, numLevels,family = "binomial")
  # cv_fit_list[[i]] <- cv_fit
  # count <- count + 1
}

# plot(cv_fit)
# 
# i_1Std <- which(cv_fit$lambdaHat1Std == cv_fit$lambda) # the preferential lambda (tuning parameter)
# 
# coefs <- coef(cv_fit$glinternetFit)[[i_1Std]] 
# 
# coefs$mainEffects # model part without interactions
# names(numLevels)[coefs$mainEffects$cont] # Main effect variables (without interactions)
# 
# coefs$interactions # model part with interactions pairs
# names(numLevels)[coefs$interactions$contcont] # Main effect variables (with interactions)
# 
# 
# # Assessing performance ####
# sqrt(cv_fit$cvErr[[i_1Std]]) # root mean squared error (RMSE) on validation data
# 
# 
# # predict.glinternet.cv(cv_fit,subset_Testing_Data,type="response") # does not work
# # glinternet::predict.glinternet.cv(cv_fit,AllTesting_Data,type="response") # does not work
# # predict.glinternet(cv_fit,AllTesting_Data,type="response") # does not work
# mypred <- predict(cv_fit,x1_test,type="response")
# fitted.results_bestglm <- ifelse(mypred > 0.5,1,0)
# 
# 
# misClassError(y1_test,fitted.results_bestglm)
# 
# # Confusion matrix ####
# # a) Confusion matrix manually calculated
# obs_pred <- cbind(y1_test,fitted.results_bestglm)
# tp <- sum(rowSums(obs_pred)==2)
# tn <- sum(rowSums(obs_pred)==0)
# fp <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
# fn <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
# con_tab1 <- matrix(c(tp,fn,fp,tn),nrow=2,ncol=2)
# con_tab1b <- con_tab1
# colnames(con_tab1) <- c('Actual TRUE','Actual FALSE');rownames(con_tab1) <- c('Predicted TRUE','Predicted FALSE')
# con_tab1[,1] <- con_tab1[,1]/ sum(y1_test==1)
# con_tab1[,2] <-con_tab1[,2]/sum(y1_test==0)
# # b) Confusion matrix from package InformationValue
# con_tab <- InformationValue::confusionMatrix(y1_test,fitted.results_bestglm)
# 
# # Sensitivity and Specificity
# # a) Manually calculated
# spec <- tn/(tn+fp) 
# sens <- tp/(tp+fn) 
# # b) Using package caret
# caret::sensitivity(data=as.factor(fitted.results_bestglm),reference=as.factor(y1_test),positive="1",negative="0")
# caret::specificity(data=as.factor(fitted.results_bestglm),reference=as.factor(y1_test),positive="1",negative="0")
# # c) using package InformationValue
# InformationValue::sensitivity(y1_test,fitted.results_bestglm)
# InformationValue::specificity(y1_test,fitted.results_bestglm)
# 
# 
# # ROC ####
# pr <- prediction(mypred, y1_test)
# prf <- performance(pr, measure = "tpr", x.measure = "fpr")
# plot(prf)
# plotROC(actuals=y1_test,predictedScores=fitted.results_bestglm)
# auc2 <- auc(y1_test,fitted.results_bestglm)
# auc <- performance(pr, measure = "auc")
# auc@y.values[[1]]
# 
# # Calculate Sensitivity and Specificity using performance command
# spec2 <- performance(pr, measure = "spec") 
# # spec2b <- performance(pr, measure = "tnr") # equivalent
# plot(spec2) # you can see that at cutoff 0.5 it is equal to my result, however I don't know how to extract the value
# # sens2 <- performance(pr, measure = "sens") # does now work
# sens2 <- performance(pr, measure = "tpr")
# sens_spec <- performance(pr, measure="sens", x.measure="spec")
# plot(sens_spec) # inverted AUC



# Problems
# 17 months
# ties for some pixels
# projection, weird spatial plots in general
# order of dimensions
# no model output

