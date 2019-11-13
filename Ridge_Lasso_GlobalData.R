##############################################################################
###########          Ridge and LASSO regression                    ###########
###########   on global standardised crop model data               ###########
##############################################################################


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

# which method? model_name in c("Ridge", "Lasso)
model_name <- "Ridge"
stopifnot(model_name %in% c("Ridge", "Lasso"))

if(model_name=="Lasso"){
  no_model <- 1
}

if(model_name=="Ridge"){
  no_model <- 0
}

#which lambda?
lambda_VALS <- c("lambda.min", "lambda.1se")
lambda_val <- lambda_VALS[1]

#threshold for bad yields in c(0.025,0.05,0.1)
threshold <- 0.05

#which segregation threshold for the model?
segreg_th <- 0.5

##### Initialisation, librairies, data #####

library(ncdf4);library(rgdal);library(raster);library(RColorBrewer);library(viridis)
library(maps);library(mapdata);library(ggplot2)
library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr)
library(foreach);library(doParallel)
library(tictoc)


#############################
##### Standardised data #####
#############################




# Get the data
# path_to_NH_files <- "/scratch3/pauline/Damocles_training_school_Como2019/GroupProject1/Data/NH"
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"

nh_files <- list.files(path=path_to_NH_files,pattern="NH_yield*") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),
                  FUN = function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
yield <- ncvar_get(nh_data[[1]],"yield")
tasmax <- ncvar_get(nh_data[[1]],"tasmax")
vpd <- ncvar_get(nh_data[[1]],"vpd")
pr <- ncvar_get(nh_data[[1]],"pr")
lat_subset <- ncvar_get(nh_data[[1]],"lat")
lon_subset <- ncvar_get(nh_data[[1]],"lon")
yield_stand <- ncvar_get(nh_data[[2]],"yield")
tasmax_stand <- ncvar_get(nh_data[[2]],"tasmax")
vpd_stand <- ncvar_get(nh_data[[2]],"vpd")
pr_stand <- ncvar_get(nh_data[[2]],"pr")
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

lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})


# Process data ####
###################

yields_3dim <- array(yield,dim=c(965,1,1600));yields_stand_3dim <- array(yield,dim=c(965,1,1600))
Model_data <- abind(yields_3dim,tasmax,vpd,pr,along=2)
Model_data_stand <- abind(yields_stand_3dim,tasmax_stand,vpd_stand,pr_stand,along=2)


pix_num <- dim(Model_data)[1]
low_yield <- sapply(1:pix_num,function(x) {quantile(yield[x,],threshold,na.rm=T)})
cy <- t(sapply(1:pix_num,function(x){ifelse(yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield

cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))
Model_data[,1,] <-cy_reshaped
Model_data_stand[,1,] <-cy_reshaped

columnnames <- c("Yield",
                 "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                 "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                 "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
                 "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                 "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                 "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2",
                 "pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                 "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                 "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")
colnames(Model_data) <- columnnames
colnames(Model_data_stand) <- columnnames



# Exclude NA variable columns
na_col <- matrix(data=NA,nrow=pix_num,ncol=52)
for (j in 1:pix_num){
  for (i in 1:52){
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





# Run model ####
################




# no_cores <- detectCores() / 2 - 1
# cl<-makeCluster(no_cores)
# clusterEvalQ(cl, {
#   library(glmnet)
#   library(dplyr)
# }) # parallelisation has own environment, therefore some packages and variables need be loaded again
# registerDoParallel(cl)
# 
# tic()
# # cv_fit <- foreach (i=1:dim(Model_data)[1],.multicombine=TRUE) %dopar% {
# cv_fit <- foreach (i=1:5,.multicombine=TRUE) %dopar% {
#   cv.glmnet(x = as.matrix(x1_train_list[[i]]),
#             y = as.matrix(y1_train_list[[i]]),
#             family = "binomial", alpha = no_model, nfolds = 10)
# }
# stopCluster(cl)
# toc()


#without paralellizing
tic()
model_cv_fitting <- list()
for (pixel in 1:dim(Model_data_stand)[1]) {
# for (pixel in 1:5) {
  model_cv_fitting[[pixel]] <- cv.glmnet(x = as.matrix(x1_train_list[[pixel]]),
                                         y = as.matrix(y1_train_list[[pixel]]),
                                         family = "binomial", alpha = no_model, nfolds = 10)
  print(paste(pixel, "out of", dim(Model_data_stand)[1]))
}#end for pixel

toc()
#15 min for Ridge
save(model_cv_fitting, file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                                    model_name,"_Thresholdbadyield", str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))






# Load the fitted model ####
############################

load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/",
                  model_name,"_thresoldbadyield", str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))







# Model performance assessment ####
###################################
test_length <- length(model_cv_fitting)

coefs <- lapply(1:test_length, function(x){coef(model_cv_fitting[[x]], s=lambda_val)})

# coefs$mainEffects # model part without interactions
# names(numLevels)[coefs$mainEffects$cont] # Main effect variables (without interactions)
# 
# coefs$interactions # model part with interactions pairs
# names(numLevels)[coefs$interactions$contcont] # Main effect variables (with interactions)


mypred <- lapply(1:test_length, function(x){predict(model_cv_fitting[[x]],
                                                    as.matrix(x1_test_list[[x]]),type="response")})

fitted.results_model <- lapply(1:test_length, function(x){ifelse(mypred[[x]] > segreg_th,1,0)})

mis_clas_err <- sapply(1:test_length, function(x){misClassError(actuals = as.matrix(y1_test_list[[x]]),
                                                                predictedScores=mypred[[x]],
                                                                threshold = segreg_th)})

con_tab <-  lapply(1:test_length, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                fitted.results_model[[x]],
                                                                                threshold = segreg_th)})
sensi <- sapply(1:test_length, function(x){InformationValue::sensitivity(as.matrix(y1_test_list[[x]]),
                                                                         fitted.results_model[[x]],
                                                                         threshold = segreg_th)})
speci <- sapply(1:test_length, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                         fitted.results_model[[x]],
                                                                         threshold = segreg_th)})



# Plot specificity and sensitivity on a map ####

# ggplot ####
world <- map_data("world")
DF_miscla <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], miscla = mis_clas_err)

ggplot(data = DF_miscla, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="transparent", color="black", size=0.3) +
  geom_point(shape=15, aes(color=miscla)) +
  scale_color_gradient(limits=c(0,(floor(max(mis_clas_err)*100) + 1)/100),
                        low = "green", high = "red") +
  theme_void()+
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Misclassification\n error",
       title = paste(model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  X11(width = 20, height = 10)


# Workaround: bind lat and lon to one object, so that you have to look for just one object, and not lat-/lon-pairs
coord_subset_temp <- cbind(coord_subset,paste(coord_subset[,1],coord_subset[,2]))
coord_all_temp <- cbind(coord_all,paste(coord_all[,1],coord_all[,2]))
loc_pix <- which(coord_all_temp[,3] %in% coord_subset_temp [,3]) # locations of our pixels in the whole coordinate set

coord_all <- cbind(coord_all,rep(NA,24320))
coord_all_sensi <- cbind(coord_all,rep(NA,24320))
coord_all_classerr <- cbind(coord_all,rep(NA,24320))
for (i in 1:test_length){
  coord_all[loc_pix[i],3] <- speci[i]
  coord_all_sensi[loc_pix[i],3] <- sensi[i]
  coord_all_classerr[loc_pix[i],3] <- mis_clas_err[[i]]
}

spec_mat <- matrix(as.numeric(coord_all[,3]),nrow=320,ncol=76)
sensi_mat <- matrix(as.numeric(coord_all_sensi[,3]),nrow=320,ncol=76)
miscla_mat <- matrix(as.numeric(coord_all_classerr[,3]),nrow=320,ncol=76)

border <- readOGR('C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/ne_50m_admin_0_countries.shp')	
spec_ras_speci <- raster(t(spec_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all),
                         ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))
sensi_ras_speci <- raster(t(sensi_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all),
                          ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))
miscla_ras_speci <- raster(t(miscla_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all),
                          ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))
cuts <- seq(0,1, by = 0.1)
x11(width = 1000, height = 600)
plot(spec_ras_speci,asp=1, breaks=cuts, col=plasma(length(cuts)),
     main=paste("Specificity, ", model_name, " regression\nThreshold bad yield ",
                threshold, " ; threshold segregation", segreg_th));plot(border,add=T)
x11(width = 1000, height = 600)
plot(sensi_ras_speci,asp=1, breaks=cuts, col=plasma(length(cuts)),
     main=paste("Sensitivity, ", model_name, " regression\nThreshold bad yield ",
                threshold, " ; threshold segregation", segreg_th));plot(border,add=T)

tryyyy <- numeric()
for (i in 1:965) {
  tryyyy[i]<-mis_clas_err[[i]]
}
max(tryyyy)
cuts2 <- seq(0,0.2, by = 0.02)
x11(width = 1000, height = 600)
plot(miscla_ras_speci,asp=1, breaks=cuts2, col=viridis(length(cuts2)),
     main=paste("Misclassification Error, ", model_name, " regression\nThreshold bad yield ",
                threshold, " ; threshold segregation", segreg_th));plot(border,add=T)

# Most important coefficients ####
##################################
set.seed(2019)


firstcoeffs <- list()
pixel_i_want <-sample(1:965, size = 1)

rownames(coefs[[pixel_i_want]])[sort(abs(as.numeric(coefs[[pixel_i_want]])),
                          decreasing = T,
                          index.return=T)$ix[-which(sort(abs(as.numeric(coefs[[pixel_i_want]])), 
                                                                        decreasing = T,
                                                                        index.return=T)$ix==1)]]
# Is it the same if I run the model on standardized data (rather than rescaled)?
Stand_meteo_data <- apply(Model_data[pixel_i_want,-1,], FUN = scale, MARGIN = 1)
Stand_data <- t(cbind(Model_data[pixel_i_want,1,], Stand_meteo_data))

training_data_pix <- Stand_data[,training_indices[[pixel_i_want]]]
testing_data_pix <- Stand_data[,testing_indices[[pixel_i_want]]]

pix_in <- 1:pix_num

x1_train_pix <- as.data.frame(t(training_data_pix[non_na_col[pixel_i_want,],])) # predictors
y1_train_pix <- training_data_pix[1,] # predictand
x1_test_pix <- as.data.frame(t(testing_data_pix[non_na_col[pixel_i_want,],])) # predictors
y1_test_pix <- testing_data_pix[1,] # predictand

var_num_pix <- sum(non_na_col[pixel_i_want,])
numLevels_pix <- rep(1,times=var_num[pixel_i_want])
numLevels_pix <-  colnames(x1_test_pix)

model_1pix <- cv.glmnet(x = as.matrix(x1_train_pix),
                        y = as.matrix(y1_train_pix),
                        family = "binomial", alpha = no_model, nfolds = 10)
coef_1pix <- coef(model_1pix)

firstcoeffs[[5]] <- cbind(rownames(coefs[[pixel_i_want]])[sort(abs(as.numeric(coef_1pix)),
                                                       decreasing = T,
                                                       index.return=T)$ix[-which(sort(abs(as.numeric(coef_1pix)), 
                                                                                      decreasing = T,
                                                                                      index.return=T)$ix==1)]][1:10],
                  rownames(coefs[[pixel_i_want]])[sort(abs(as.numeric(coefs[[pixel_i_want]])),
                                                       decreasing = T,
                                                       index.return=T)$ix[-which(sort(abs(as.numeric(coefs[[pixel_i_want]])), 
                                                                                      decreasing = T,
                                                                                      index.return=T)$ix==1)]][1:10])
colnames(firstcoeffs[[5]])=c("Standardization", "Rescale with RANGE")
