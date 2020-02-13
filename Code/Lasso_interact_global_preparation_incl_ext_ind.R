# Lasso with interactions for northern hemisphere including extreme indices: preparation file

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

##### Standardised data
threshold <- 0.05 #threshold for bad yield
segreg_th <- 0.5 #If the fitted proba from the model is below this threshold then it's a bad year


# Get the data ####
###################
# # path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
nh_files <- list.files(path=path_to_NH_files,pattern="NH_yield*") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),
                  FUN = function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
# yield <- ncvar_get(nh_data[[1]],"yield")
# tasmax <- ncvar_get(nh_data[[1]],"tasmax")
# vpd <- ncvar_get(nh_data[[1]],"vpd")
# pr <- ncvar_get(nh_data[[1]],"pr")
# yield_stand <- ncvar_get(nh_data[[2]],"yield")
# tasmax_stand <- ncvar_get(nh_data[[2]],"tasmax")
# vpd_stand <- ncvar_get(nh_data[[2]],"vpd")
# pr_stand <- ncvar_get(nh_data[[2]],"pr")
lat_subset <- ncvar_get(nh_data[[1]],"lat")
lon_subset <- ncvar_get(nh_data[[1]],"lon")
lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})
coord_subset <- cbind(lon_subset,lat_subset)
# 
# # load all coordinates of northern hemisphere
# nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
# nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
# lat_all <- ncvar_get(nh_data[[1]],"lat")
# lon_all <- ncvar_get(nh_data[[1]],"lon")
# lati_all <- rep(lat_all,each=length(lon_all))
# long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
# coord_all <- cbind(long_all,lati_all)


# Process data ####
###################

load("D:/user/vogelj/Data/Group project Como/extremeindices_and_monthlymeteovar_rescaled.Rdata")


yield_3dim <- array(Data_standardized$yield,dim=c(965,1,1600))
dtr_3dim <- array(Data_standardized$dtr,dim=c(965,1,1600))
frs_3dim <- array(Data_standardized$frs,dim=c(965,1,1600))
txx_3dim <- array(Data_standardized$txx,dim=c(965,1,1600))
tnn_3dim <- array(Data_standardized$tnn,dim=c(965,1,1600))
rx5_3dim <- array(Data_standardized$rx5,dim=c(965,1,1600))
tx90p_3dim <- array(Data_standardized$tx90p,dim=c(965,1,1600))
tn10p_3dim <- array(Data_standardized$tn10p,dim=c(965,1,1600))

Model_data_stand <- abind(yield_3dim,dtr_3dim,frs_3dim,txx_3dim,tnn_3dim,rx5_3dim,tx90p_3dim,tn10p_3dim
                          ,Data_standardized$tasmax,Data_standardized$vpd,Data_standardized$pr,along=2)
colnames(Model_data_stand) <- c("Yield", "dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p",
                                "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                                "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                                "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
                                "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                                "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                                "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2",
                                "pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                                "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                                "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")


# Model_data_stand <- abind(yields_stand_3dim,tasmax_stand,vpd_stand,pr_stand,along=2)


# threshold <- 0.05
pix_num <- dim(Model_data_stand)[1]
low_yield <- sapply(1:pix_num,function(x) {quantile(yield_3dim[x,1,],threshold,na.rm=T)})
cy <- t(sapply(1:pix_num,function(x){ifelse(yield_3dim[x,1,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield


cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))
Model_data_stand[,1,] <-cy_reshaped
# Model_data_stand[,1,] <-cy_reshaped

# colnames(Model_data_stand) <- columnnames
# colnames(Model_data_stand) <- columnnames



# Exclude NA variable columns
na_col <- matrix(data=NA,nrow=pix_num,ncol=dim(Model_data_stand)[2])
for (j in 1:pix_num){
  for (i in 1:dim(Model_data_stand)[2]){
    na_col[j,i] <- all(is.na(Model_data_stand[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)


na_time <- vector("list",length=pix_num) # for each pixel, the positions of NAs over time
for (j in 1:pix_num){
  na_time[[j]] <- which(is.na(Model_data_stand[j,1,])) # locations of years with NA values
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
