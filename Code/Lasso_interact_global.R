# Lasso with interactions for northern hemisphere

library(ncdf4)
library(glinternet)
library(abind)
library(foreach) 
library(doParallel)
library(InformationValue)
library(ROCR)
library(tictoc)

# Get the data ####
###################
# path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
nh_files <- list.files(path=path_to_NH_files,pattern="NH_yield*") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),
                  FUN = function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
yield <- ncvar_get(nh_data[[1]],"yield")
tasmax <- ncvar_get(nh_data[[1]],"tasmax")
vpd <- ncvar_get(nh_data[[1]],"vpd")
pr <- ncvar_get(nh_data[[1]],"pr")
lat <- ncvar_get(nh_data[[1]],"lat")
lon <- ncvar_get(nh_data[[1]],"lon")
yield_stand <- ncvar_get(nh_data[[2]],"yield")
tasmax_stand <- ncvar_get(nh_data[[2]],"tasmax")
vpd_stand <- ncvar_get(nh_data[[2]],"vpd")
pr_stand <- ncvar_get(nh_data[[2]],"pr")
lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})
coord <- cbind(lat,lon)

# Process data ####
###################

yields_3dim <- array(yield,dim=c(965,1,1600));yields_stand_3dim <- array(yield,dim=c(965,1,1600))
Model_data <- abind(yields_3dim,tasmax,vpd,pr,along=2)
Model_data_stand <- abind(yields_stand_3dim,tasmax_stand,vpd_stand,pr_stand,along=2)


threshold <- 0.05
pix_num <- dim(Model_data)[1]
low_yield <- sapply(1:pix_num,function(x) {quantile(yield[x,],threshold,na.rm=T)})
cy <- t(sapply(1:pix_num,function(x){ifelse(yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield


cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))
Model_data[,1,] <-cy_reshaped
Model_data_stand[,1,] <-cy_reshaped

columnnames <- c("Yield","pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                 "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                 "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2",
                 "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                 "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                 "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
                 "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                 "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                 "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
colnames(Model_data) <- columnnames
colnames(Model_data_stand) <- columnnames



# Exclude NA variable columns
na_col <- vector("logical",52)
na_col <- matrix(data=NA,nrow=pix_num,ncol=52)
for (j in 1:pix_num){
  for (i in 1:52){
    na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)


# Split data into training and testing data set
set.seed(1994)
training_indices <- sort(sample(1:1600, size = floor(1600*0.6)))
testing_indices <- (1:1600)[-training_indices]
Training_Data <- Model_data_stand[,,training_indices]
Testing_Data <- Model_data_stand[,,testing_indices]
pix_in <- 1:pix_num

x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[x,non_na_col[x,],]))}) # predictors
y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[x,1,]}) # predictand

x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[x,non_na_col[x,],]))}) # predictors
y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[x,1,]}) # predictand

var_num <- apply(non_na_col,1,sum)
numLevels_list <- sapply(1:pix_num, function(x){ rep(1,times=var_num[x])})
for (i in 1:pix_num){
  names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
}



# Lasso model ####
##################


pix_with_NA <- which(apply(cy,1,anyNA))
final_pix <- 1:965; 
final_pix <- final_pix[-pix_with_NA]


tic()
no_cores <- detectCores() / 2 - 1
cl<-makeCluster(no_cores)
clusterEvalQ(cl, {
  library(glinternet)
  library(dplyr)
}) # parallelisation has own environment, therefore some packages and variables need be loaded again
registerDoParallel(cl)


# cv_fit <- foreach (i=1:dim(Model_data)[1],.multicombine=TRUE) %dopar% {
cv_fit <- foreach (i=1:5,.multicombine=TRUE) %dopar% {
  # for (i in 1:dim(Model_data)[1][1:3]){
  # Fit model
  # cv_fit <- glinternet.cv(X1_train, y1_train, numLevels,family = "binomial")
  glinternet.cv(x1_train_list[[i]], y1_train_list[[i]], numLevels_list[[i]],family = "binomial")
  # cv_fit_list[[i]] <- cv_fit
  # count <- count + 1
}
stopCluster(cl)
toc()




# Model performance assessment ####
###################################

i_1Std <- lapply(1:5, function(x){ which(cv_fit[[x]]$lambdaHat1Std == cv_fit[[x]]$lambda)}) # the preferential lambda (tuning parameter)

coefs <- lapply(1:5, function(x){coef(cv_fit[[x]]$glinternetFit)[[i_1Std[[x]]]]})

# coefs$mainEffects # model part without interactions
# names(numLevels)[coefs$mainEffects$cont] # Main effect variables (without interactions)
# 
# coefs$interactions # model part with interactions pairs
# names(numLevels)[coefs$interactions$contcont] # Main effect variables (with interactions)


mypred <- lapply(1:5, function(x){predict(cv_fit[[x]],x1_test_list[[x]],type="response")}) 
fitted.results_bestglm <- lapply(1:5, function(x){ifelse(mypred[[x]] > 0.5,1,0)})

mis_clas_err <- lapply(1:5, function(x){misClassError(y1_test_list[[x]],fitted.results_bestglm[[x]])})

con_tab <-  lapply(1:5, function(x){InformationValue::confusionMatrix(y1_test_list[[x]],fitted.results_bestglm[[x]])})
sensi <- sapply(1:5, function(x){InformationValue::sensitivity(y1_test_list[[x]],fitted.results_bestglm[[x]])})
speci <- sapply(1:5, function(x){InformationValue::specificity(y1_test_list[[x]],fitted.results_bestglm[[x]])})



# # Plot specificity and sensitivity on a map

longs <- rep(long,length(lati)) # coordinates rearranged
lats <- rep(lati,each=length(long))
coord_all <- cbind(longs,lats)





spec <- matrix(NA,nrow=320,ncol=320)
sens <- matrix(NA,nrow=320,ncol=320)

pix_in2 <- vector(mode="numeric",length=24320)
pix_in2[pix_in] <- 1 # mark all pixels with values
pix_in2[pix_in[1:50]] <- 1
coord_pix <- cbind(coord,pix_in2)
coord_spec <- coord_pix
# coord_spec[,pix_in[1:50]] <- speci
# coord_spec[,pix_in[x]] <- speci

for (i in 1: 50){
  coord_spec[pix_in[i],3] <- speci[i]
}

coord_spec2 <- matrix(coord_spec[,3],nrow=320,ncol=76)

# dat_ras <- rasterFromXYZ(coord_spec)
spec_ras <- raster(t(coord_spec2[,76:1]), xmn=min(long), xmx=max(long), ymn=min(lati), ymx=max(lati), crs=CRS(projection(border)))

x11()
plot(spec_ras,asp=1);plot(border,add=T)