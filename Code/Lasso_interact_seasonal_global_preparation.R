# Lasso with interactions for northern hemisphere

library(ncdf4)
library(glinternet)
library(abind)
library(foreach) 
library(doParallel)
library(InformationValue)
library(ROCR)
library(tictoc)
library(rgdal)
library(raster)

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
yield_stand <- ncvar_get(nh_data[[2]],"yield")
tasmax_stand <- ncvar_get(nh_data[[2]],"tasmax")
vpd_stand <- ncvar_get(nh_data[[2]],"vpd")
pr_stand <- ncvar_get(nh_data[[2]],"pr")
lat_subset <- ncvar_get(nh_data[[1]],"lat")
lon_subset <- ncvar_get(nh_data[[1]],"lon")
lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})
coord_subset <- cbind(lon_subset,lat_subset)

# load all coordinates of northern hemisphere
nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
lat_all <- ncvar_get(nh_data[[1]],"lat")
lon_all <- ncvar_get(nh_data[[1]],"lon")
lati_all <- rep(lat_all,each=length(lon_all))
long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
coord_all <- cbind(long_all,lati_all)


# Process data ####
###################

yields_3dim <- array(yield,dim=c(965,1,1600));yields_stand_3dim <- array(yield,dim=c(965,1,1600))

seasons <- matrix(c(2:16),nrow=3,ncol=5) # five seasons autumn1, winter, spring, summer, autumn2; first August and last December excluded
tasmax_seas <- array(NA,dim=c(965,5,1600))
pr_seas <- array(NA,dim=c(965,5,1600))
vpd_seas <- array(NA,dim=c(965,5,1600))
tasmax_seas_stand <- array(NA,dim=c(965,5,1600))
pr_seas_stand <- array(NA,dim=c(965,5,1600))
vpd_seas_stand <- array(NA,dim=c(965,5,1600))
for (i in c(1:5)){
  curr_seas <- seasons[,i] # current season (3 months)
  pr_seas[,i,] <- t(sapply(1:965, function(x) {apply(pr[x,curr_seas,],2,mean,na.rm=T)}))
  tasmax_seas[,i,] <- t(sapply(1:965, function(x) {apply(tasmax[x,curr_seas,],2,mean,na.rm=T)}))
  vpd_seas[,i,] <- t(sapply(1:965, function(x) {apply(vpd[x,curr_seas,],2,mean,na.rm=T)}))
  pr_seas_stand[,i,] <- t(sapply(1:965, function(x) {apply(pr_stand[x,curr_seas,],2,mean,na.rm=T)}))
  tasmax_seas_stand[,i,] <- t(sapply(1:965, function(x) {apply(tasmax_stand[x,curr_seas,],2,mean,na.rm=T)}))
  vpd_seas_stand[,i,] <- t(sapply(1:965, function(x) {apply(vpd_stand[x,curr_seas,],2,mean,na.rm=T)}))
}

Model_data <- abind(yields_3dim,tasmax_seas,vpd_seas,pr_seas,along=2)
Model_data_stand <- abind(yields_3dim,tasmax_seas_stand,vpd_seas_stand,pr_seas_stand,along=2)

threshold <- 0.05
pix_num <- dim(Model_data)[1]
low_yield <- sapply(1:pix_num,function(x) {quantile(yield[x,],threshold,na.rm=T)})
cy <- t(sapply(1:pix_num,function(x){ifelse(yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield


cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))
Model_data[,1,] <-cy_reshaped
Model_data_stand[,1,] <-cy_reshaped

columnnames <- c("Yield","tmax_aut_Y1","tmax_win_Y1","tmax_spr_Y1","tmax_sum_Y1","tmax_aut_Y2",
                 "vpd_aut_Y1","vpd_win_Y1","vpd_svpd_Y1","vpd_sum_Y1","vpd_aut_Y2",
                 "pr_aut_Y1","pr_win_Y1","pr_spr_Y1","pr_sum_Y1","pr_aut_Y2")
colnames(Model_data) <- columnnames
colnames(Model_data_stand) <- columnnames



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


# Option 1: Exclude years with NAs (not the whole pixel, just the NA years of it) ####
######################################################################################

# Split data into training and testing data set
vec <- 1:1600

years_with_na <- vector("logical",length=pix_num)
for (i in 1:pix_num){
  years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
}

training_indices <- vector("list",length=pix_num)
testing_indices <- vector("list",length=pix_num)
set.seed(1994)
for (x in 1:pix_num) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((1600-length(na_time[[x]]))*0.6)))
    testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
  } else {
    training_indices[[x]] <- sort(sample(1:1600, size = floor(1600*0.6)))
    testing_indices[[x]] <- (1:1600)[-training_indices[[x]]]    
  }
}
# training_indices2 <- sapply(1:pix_num, function(x) {sort(sample(x=vec[-na_time[[x]]], size = floor((1600-length(na_time[[x]]))*0.6)))})
# testing_indices2 <- sapply(1:pix_num, function(x) {vec[-na_time[[x]]][-training_indices[[x]]]})
# training_indices2 <- sapply(1:pix_num, function(x) {sort(sample(1:(1600-length(na_time[[x]])), size = floor((1600-length(na_time[[x]]))*0.6)))})
# testing_indices2 <- sapply(1:pix_num, function(x) {(1:(1600-length(na_time[[x]])))[-training_indices2[[x]]]})

Training_Data <- lapply(1:pix_num,function(x){Model_data_stand[x,,training_indices[[x]]]})
Testing_Data <- lapply(1:pix_num,function(x){Model_data_stand[x,,testing_indices[[x]]]})

pix_in <- 1:pix_num
# x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[x,non_na_col[x,],]))}) # predictors
# y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[x,1,]}) # predictand
# x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[x,non_na_col[x,],]))}) # predictors
# y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[x,1,]}) # predictand

x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand

var_num <- apply(non_na_col,1,sum)
numLevels_list <- sapply(1:pix_num, function(x){ rep(1,times=var_num[x])})
for (i in 1:pix_num){
  names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
}


# Alternative Option 2: exclude all pixels which have years with NA ####
####################################################################

# # Split data into training and testing data set
# set.seed(1994)
# training_indices <- sort(sample(1:1600, size = floor(1600*0.6)))
# testing_indices <- (1:1600)[-training_indices]
# Training_Data <- Model_data_stand[,,training_indices]
# Testing_Data <- Model_data_stand[,,testing_indices]
# 
# pix_in <- 1:pix_num
# x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[x,non_na_col[x,],]))}) # predictors
# y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[x,1,]}) # predictand
# 
# x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[x,non_na_col[x,],]))}) # predictors
# y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[x,1,]}) # predictand
# 
# var_num <- apply(non_na_col,1,sum)
# numLevels_list <- sapply(1:pix_num, function(x){ rep(1,times=var_num[x])})
# for (i in 1:pix_num){
#   names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
# }

# pix_with_NA <- which(apply(cy,1,anyNA))
# final_pix <- 1:965;
# final_pix <- final_pix[-pix_with_NA]



