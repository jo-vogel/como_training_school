# Lasso with interactions for northern hemisphere including extreme indices

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
# lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})
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


threshold <- 0.05
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




# Lasso model ####
##################

# message('This small section if meant for reruns with the preferential lambda')
# load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_complete.RData") # load model output (if you want to recalculate it based on the preferential lambda)
# failed_pixels <- which(sapply(1:965, function(x) {is.character(cv_fit[[x]])})==1)
# work_pix <- pix_in[-failed_pixels] # working pixels
# pref_lam <- rep(NA,965)
# pref_lam[work_pix] <- sapply(work_pix, function(x) cv_fit[[x]]$lambdaHat1Std)


tic()
# cv_fit_list <- vector("list",length=dim(Model_data_stand)[1])

no_cores <- detectCores() / 2 - 1
cl<-makeCluster(no_cores)
clusterEvalQ(cl, {
  library(glinternet)
  library(dplyr)
}) # parallelisation has own environment, therefore some packages and variables need be loaded again
registerDoParallel(cl)


set.seed(100)
cv_fit <- foreach (i=1:dim(Model_data_stand)[1],.multicombine=TRUE) %dopar% {
  # for (i in 1:dim(Model_data)[1]){
  # normal run
  tryCatch(glinternet.cv(x1_train_list[[i]], y1_train_list[[i]], numLevels_list[[i]],family = "binomial",interactionCandidates=""), error=function(e) paste0("Error in iteration ",i))
  # tryCatch(glinternet.cv(x1_train_list[[i]], y1_train_list[[i]], numLevels_list[[i]],family = "binomial"), error=function(e) paste0("Error in iteration ",i))
  
  # run with the preferential lambda from the normal run
  # tryCatch(glinternet.cv(x1_train_list[[i]], y1_train_list[[i]], numLevels_list[[i]],family = "binomial",lambda=pref_lam[i]), error=function(e) paste0("Error in iteration ",i))
  
  # without parallelisation
  # cv_fit <- glinternet.cv(x1_train_list[[i]], y1_train_list[[i]], numLevels_list[[i]],family = "binomial",lambda=pref_lam[i])
  # cv_fit <- tryCatch(glinternet.cv(x1_train_list[[i]], y1_train_list[[i]], numLevels_list[[i]],family = "binomial",lambda=pref_lam[i]), error=function(e) paste0("Error in iteration ",i))
  # cv_fit_list[[i]] <- cv_fit
}
stopCluster(cl)
toc()
# save(cv_fit,file="cv_fit_monthly_without_int_incl_ext.RData")



# Model performance assessment ####
###################################

# Identify pixels with failed runs
failed_pixels <- which(sapply(1:965, function(x) {is.character(cv_fit[[x]])})==1)
# work_pix <- ifelse(length(failed_pixels)==0, pix_in, pix_in[-failed_pixels]) # working pixels
work_pix <- pix_in[-failed_pixels] # working pixels
work_pix <- pix_in


i_1Std <- sapply(work_pix, function(x){ which(cv_fit[[x]]$lambdaHat1Std == cv_fit[[x]]$lambda)}) # the preferential lambda (tuning parameter)
i_1Std_all_pix <- rep(NA,965)
i_1Std_all_pix[work_pix] <- i_1Std # needed as a workaround (to have an object of lenght=965)

coefs <- vector("list",length=965)
coefs[work_pix] <- lapply(work_pix, function(x){coef(cv_fit[[x]]$glinternetFit)[[i_1Std_all_pix[[x]]]]})

# coefs$mainEffects # model part without interactions
# names(numLevels)[coefs$mainEffects$cont] # Main effect variables (without interactions)
# 
# coefs$interactions # model part with interactions pairs
# names(numLevels)[coefs$interactions$contcont] # Main effect variables (with interactions)


#which segregation threshold for the model?
segreg_th <- 0.5
mypred <- lapply(work_pix, function(x){predict(cv_fit[[x]],x1_test_list[[x]],type="response", lambdaType="lambdaHat1Std")}) 
fitted.results_model <- lapply(seq_along(work_pix), function(x){ifelse(mypred[[x]] > segreg_th,1,0)})

y1_test_list_red <- lapply(work_pix,function(work_pix){y1_test_list[[work_pix]]})
mis_clas_err <- rep(NA,965)
# mis_clas_err[work_pix] <- sapply(seq_along(work_pix), function(x){misClassError(y1_test_list_red[[x]],mypred[[x]])})
mis_clas_err[work_pix] <- sapply(seq_along(work_pix), function(x){misClassError(actuals = y1_test_list_red[[x]],
                                                                                predictedScores=mypred[[x]],
                                                                                threshold = segreg_th)})



con_tab <-  lapply(seq_along(work_pix), function(x){InformationValue::confusionMatrix(y1_test_list_red[[x]],fitted.results_model[[x]])})
sensi <- rep(NA,965)
speci <- rep(NA,965)
sensi[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::sensitivity(y1_test_list_red[[x]],fitted.results_model[[x]],
                                                                                         threshold = segreg_th)})
speci[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::specificity(y1_test_list_red[[x]],fitted.results_model[[x]],
                                                                                         threshold = segreg_th)})

# Calculate CSI ####
obs_pred <- lapply(seq_along(work_pix), function(x){cbind(y1_test_list_red[[x]],fitted.results_model[[x]])})
tp <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==2)})
tn <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==0)})
fp <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==0 & obs_pred[[x]][,2]==1)})
fn <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==1 & obs_pred[[x]][,2]==0)})
con_tab1 <- sapply(seq_along(work_pix), function(x){matrix(c(tp[x],fn[x],fp[x],tn[x]),nrow=2,ncol=2)})
# spec <- tn/(tn+fp) 
# sens <- tp/(tp+fn) 
csi <- sapply(seq_along(work_pix), function(x){tn[x]/(tn[x]+fp[x]+fn[x])})



# Plot specificity and sensitivity on a map ####

# Workaround: bind lat and lon to one object, so that you have to look for just one object, and not lat-/lon-pairs
coord_subset_temp <- cbind(coord_subset,paste(coord_subset[,1],coord_subset[,2]))
coord_all_temp <- cbind(coord_all,paste(coord_all[,1],coord_all[,2]))
loc_pix <- which(coord_all_temp[,3] %in% coord_subset_temp [,3]) # locations of our pixels in the whole coordinate set
loc_pix <- loc_pix[-failed_pixels]

coord_all <- cbind(coord_all,rep(NA,24320))
for (i in seq_along(work_pix)){
  coord_all[loc_pix[i],3] <- speci[i]
}

spec_mat <- matrix(as.numeric(coord_all[,3]),nrow=320,ncol=76)

border <- readOGR('D:/user/vogelj/Data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')	
spec_ras <- raster(t(spec_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all), ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))

x11()
plot(spec_ras,asp=1);plot(border,add=T)



