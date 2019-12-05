##############################################################################
###########          Ridge and LASSO regression                    ###########
###########   on global standardised crop model data               ###########
##############################################################################


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

# which method? model_name in c("Ridge", "Lasso)
model_name <- "Lasso"
stopifnot(model_name %in% c("Ridge", "Lasso"))

if(model_name=="Lasso"){
  no_model <- 1
}

if(model_name=="Ridge"){
  no_model <- 0
}

#which lambda?
lambda_VALS <- c("lambda.min", "lambda.1se")
lambda_val <- lambda_VALS[2]

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

#which segregation threshold for the model?
segreg_th <- 0.5

fitted.results_model <- lapply(1:test_length, function(x){ifelse(mypred[[x]] > segreg_th,1,0)})

mis_clas_err <- sapply(1:test_length, function(x){misClassError(actuals = as.matrix(y1_test_list[[x]]),
                                                                predictedScores=mypred[[x]],
                                                                threshold = segreg_th)})

con_tab <-  lapply(1:test_length, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                fitted.results_model[[x]],
                                                                                threshold = segreg_th)})
# sensi <- sapply(1:test_length, function(x){InformationValue::sensitivity(as.matrix(y1_test_list[[x]]),
#                                                                          fitted.results_model[[x]],
#                                                                          threshold = segreg_th)})
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

# Plot miscla error ####

world <- map_data("world")
DF_miscla <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], miscla = mis_clas_err)

ggplot(data = DF_miscla, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=miscla)) +
  scale_color_gradient2(limits=c(0,0.2), midpoint = 0.1,
                        low = "yellow", mid = "red3", high = "black") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Misclass.\nerror",
       title = paste("Misclassification error, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


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
       title = paste("Specificity, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


# plot CSI=(hits)/(hits + misses + false alarm) ###
DF_sci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], csi = csi)

ggplot(data = DF_sci, aes(x=lon, y=lat)) +
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
  labs(color="csi",
       title = paste("Critical success index, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)



# Correlation between yield and specificity or miscla error ####

mean_yield <- apply(X=yield, MARGIN = 1, FUN = mean, na.rm=T)
plot(mean_yield, mis_clas_err,
     xlab="Mean Yield (kg/yr)", ylab="Miss-classification error",
     main=paste("Scatterplot mean yield, miss-classification error\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))
plot(mean_yield, csi,
     xlab="Mean Yield (kg/yr)", ylab="CSI",
     main=paste("Scatterplot mean yield, CSI\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))
plot(mean_yield, speci,
     xlab="Mean Yield (kg/yr)", ylab="Specificity",
     main=paste("Scatterplot mean yield, Specificity\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))


pairs(cbind(mean_yield, mis_clas_err, speci, csi),
      main=paste(model_name, " regression, bad yield thr.=", threshold,
                 "\nsegreg. thr.=", segreg_th, sep = ""))

cor(mean_yield, mis_clas_err)
cor(mean_yield, speci)


# Map of the number of coefficients kept (Lasso) #####


source("./Code/get_first_coeff_function.R")

if(model_name=="Lasso"){
  coeff_kep <- numeric()
  
  for (pix in 1:pix_num) {
    # coeff_kep[pix] <- length(coefs[[pix]]$mainEffects$cont)
    coeff_kep[pix] <- sum(coefs[[pix]][row.names(coefs[[pix]])!="(Intercept)"]!=0)
    
  }#end for pix
  
  
  
  # Plot 
  
  DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = coeff_kep)
  
  ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=coeff_kep)) +
    scale_color_gradient(limits=c(1,30),
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
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  
  
  
  
  
  # Plot nb of var kept in every season
  nb_in_aug1 <- numeric()
  nb_in_fall1 <- numeric()
  nb_in_winter1 <- numeric()
  nb_in_spring1 <- numeric()
  nb_in_summer1 <- numeric()
  nb_in_fall2 <- numeric()
  nb_in_dec2 <- numeric()
  
  for (pix in 1:pix_num) {
    all_month_coeff_kept <- get_firstcoeffs(coefs[[pix]], nb_of_coeff = coeff_kep[pix])[,2]
    nb_in_aug1[pix] <- length(which(all_month_coeff_kept=="Aug_Y1"))
    nb_in_fall1[pix] <- length(which(all_month_coeff_kept=="Sep_Y1" | all_month_coeff_kept=="Oct_Y1" 
                                       | all_month_coeff_kept=="Nov_Y1"))
    nb_in_winter1[pix] <- length(which(all_month_coeff_kept=="Dec_Y1" | all_month_coeff_kept=="Jan_Y2" 
                                         | all_month_coeff_kept=="Feb_Y2"))
    nb_in_spring1[pix] <- length(which(all_month_coeff_kept=="Mar_Y2" | all_month_coeff_kept=="Apr_Y2" 
                                         | all_month_coeff_kept=="May_Y2"))
    nb_in_summer1[pix] <- length(which(all_month_coeff_kept=="Jun_Y2" | all_month_coeff_kept=="Jul_Y2" 
                                         | all_month_coeff_kept=="Aug_Y2"))
    nb_in_fall2[pix] <- length(which(all_month_coeff_kept=="Sep_Y2" | all_month_coeff_kept=="Oct_Y2" 
                                       | all_month_coeff_kept=="Nov_Y2"))
    nb_in_dec2[pix] <- length(which(all_month_coeff_kept=="Dec_Y2"))
  }#end for pix
  
  # Plot 
  
  ALL_COL <- scales::seq_gradient_pal("pink", "darkblue", "Lab")(seq(0,1,length.out = 11))
  
  DF_numbcoeff_season <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_var = nb_in_fall2)
  DF_numbcoeff_season$nb_var <- as.factor(DF_numbcoeff_season$nb_var)
  ggplot(data = DF_numbcoeff_season, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=DF_numbcoeff_season$nb_var)) +
    scale_color_manual(values = ALL_COL[(min(as.numeric(DF_numbcoeff_season$nb_var)):max(as.numeric(DF_numbcoeff_season$nb_var)))+1]) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="Nb of var. kept\nin Fall Y2",
         title = paste("Number of variables kept in Fall Y2, simple",model_name,"regression"),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  
  
  
  #plot ratio of nb of variable/nb of variables? because difference in growing season length!
  
  
  
  # #Apparently no correlation between nb coeff kept and mean yield, miscla or specificity
  # plot(coeff_kep, mean_yield)
  # plot(coeff_kep, mis_clas_err)
  # plot(coeff_kep, speci)
  # cor(coeff_kep, mean_yield)
  # cor(coeff_kep, mis_clas_err)
  # cor(coeff_kep, speci)
  # 
  
}#end if Lasso




# Represent the importance of first coefficients Ridge

if(model_name=="Ridge"){
  max_coeff_val <- numeric()
  for (pix in 1:pix_num) {
    max_coeff_val[pix] <- sort(abs(coefs[[pix]][-which(rownames(coefs[[pix]])=="(Intercept)")]),decreasing = T)[1]
  }#end for pix
  
  perc_first_coeffs <- function(coeffs, nb_first_coeff = 1){
    return(round(100*sort(abs(coeffs[-which(rownames(coeffs)=="(Intercept)")]),
                          decreasing = T)[1:nb_first_coeff]/sum(sort(abs(coeffs[-which(rownames(coeffs)=="(Intercept)")]),
                                                                decreasing = T))))
  }#end perc_first_coeffs
  
  coef_1_perc <- sapply(coefs, FUN = perc_first_coeffs)
  coef_1_2_perc <- apply(sapply(coefs, FUN = perc_first_coeffs, nb_first_coeff = 2), FUN = sum, MARGIN = 2)
  coef_1_3_perc <- apply(sapply(coefs, FUN = perc_first_coeffs, nb_first_coeff = 3), FUN = sum, MARGIN = 2)
  coef_1_4_perc <- apply(sapply(coefs, FUN = perc_first_coeffs, nb_first_coeff = 4), FUN = sum, MARGIN = 2)
  coef_1_5_perc <- apply(sapply(coefs, FUN = perc_first_coeffs, nb_first_coeff = 5), FUN = sum, MARGIN = 2)
  

  DF_perc_var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], perc_coef = coef_1_3_perc)
  ggplot(data = DF_perc_var, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=DF_perc_var$perc_coef)) +
    scale_color_gradientn(colours = rainbow(8)[-8], limits = c(13, 84)) +
    # scale_color_gradient2(high = "chartreuse4",low = "darkblue",mid="skyblue",midpoint = 50,limits = c(5, 91))+
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="contribution\n(%)",
         title = paste("Percentage of the 3 first coeff among sum of all coeff, simple",model_name,"regression"),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  
}#end if Ridge





# Map of the most important coefficient ####

source("./Code/get_first_coeff_function.R")

coef_to_plot <- 1 #between 1 and coeff to keep
first_var <- numeric()
month_first_var <- numeric()

for (pix in 1:length(coefs)) {
  first_var[pix] <- get_firstcoeffs(coefs[[pix]], nb_of_coeff = coef_to_plot)[coef_to_plot,1]
  month_first_var[pix] <- get_firstcoeffs(coefs[[pix]], nb_of_coeff = coef_to_plot)[coef_to_plot,2]
}#end for pix


world <- map_data("world")

DF_1var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2],
                       first_var = first_var, month_first_var=month_first_var)

ggplot(data = DF_1var, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=first_var), size=1) +
  scale_color_manual(values = c("pr"="blue", "vpd"="green", "tmax"="red")) +
  # scale_alpha_manual(values = c("Aug_Y1"=(1/14), "Sep_Y1"=(2/14),"Oct_Y1"=(3/14), "Nov_Y1"=4/14,"Dec_Y1"=(5/14),"Jan_Y2"=(6/14),
  #                               "Feb_Y2"=(7/14),"Mar_Y2"=(8/14),"Apr_Y2"=(9/14), "May_Y2"=(10/14),"Jun_Y2"=(11/14),
  #                               "Jul_Y2"=(12/14), "Aug_Y2"=(13/14), "Sep_Y2"=(14/14)),
  #                   breaks = c("Aug_Y1", "Sep_Y1","Oct_Y1", "Nov_Y1","Dec_Y1","Jan_Y2", "Feb_Y2","Mar_Y2",
  #                              "Apr_Y2", "May_Y2","Jun_Y2","Jul_Y2", "Aug_Y2", "Sep_Y2")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="1st var", alpha="Month",
       title = paste("1st predictor, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  guides(color=guide_legend(order = 1), alpha=guide_legend(order = 2, ncol = 2))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 23, height = 7)





# # Plot total % of contribution of most important meteo var ####
# /!\ not correct because variable are not orthogonal, I cannot just some the vaule of coeff to do a %
# 
# contribution_pr <- numeric()
# contribution_vpd <- numeric()
# contribution_tmax <- numeric()
# max_contribution_name <- numeric()
# max_contribution <- numeric()
# for (pix in 1:pix_num) {
#   coeff <- abs(coefs[[pix]])
#   coeff_names <- row.names(coeff)
#   n_coeff <- length(coeff)-1
#   
#   contribution_pr[pix] <- sum(coeff[which(substr(coeff_names, start = 1, stop = 1)=="p")])/sum(coeff[-1])
#   contribution_vpd[pix] <- sum(coeff[which(substr(coeff_names, start = 1, stop = 1)=="v")])/sum(coeff[-1])
#   contribution_tmax[pix] <- sum(coeff[which(substr(coeff_names, start = 1, stop = 1)=="t")])/sum(coeff[-1])
#   
#   max_contribution[pix] <- max(c(contribution_pr[pix],contribution_vpd[pix],contribution_tmax[pix]))
#   
#   max_contribution_name[pix] <- c("pr", "vpd", "tmax")[which.max(c(contribution_pr[pix],
#                                                                    contribution_vpd[pix],contribution_tmax[pix]))]
#   
# }#end for pix
# jet.colors <-
#   colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#                      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# DF_perc_var1 <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], perc_coef = 100*max_contribution)
# ggplot(data = DF_perc_var1, aes(x=lon, y=lat)) +
#   geom_polygon(data = world, aes(long, lat, group=group),
#                fill="white", color="black", size=0.3) +
#   geom_point(shape=15, aes(color=DF_perc_var1$perc_coef)) +
#   scale_color_gradientn(colours = rev(viridis(6))) +
#   theme(panel.ontop = F, panel.grid = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
#   ylab("Lat (°N)") +
#   xlab("Lon (°E)") +
#   coord_fixed(xlim = c(-120, 135),
#               ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
#               ratio = 1.3)+
#   labs(color="contribution\n(%)",
#        title = paste("Most important meteo variable: contribution among all coeff, simple",model_name,"regression"),
#        subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
#   theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
#         legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
#   X11(width = 20, height = 7)




# Plot number of different month in N first variables ####
N_var <- 10
nb_months <- numeric()

for (pix in 1:length(coefs)) {
  if(model_name=="Lasso" & coeff_kep[pix]<N_var){
    nb_months[pix] <- NA
  } else {
    Q <- try(get_firstcoeffs(coefs[[pix]], nb_of_coeff = N_var)[,2], silent = T)
    if(class(Q)=="try-error"){
      nb_months[pix] <- NA
    } else {
    nb_months[pix] <- length(unique(Q))
    }
  }#end if Lasso not enough coeff
}#end for pix


DF_nbdiffmonths <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_months = nb_months)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=nb_months)) +
  scale_color_gradient(low = rgb(1,0.7,0.7), high = rgb(0.1,0,0)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of \ndifferent months",
       title = paste("Number of different months in the ",N_var," first variables, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)




# Plot number of different month in N first variables ####
N_var <- 10
nb_months <- numeric()

for (pix in 1:length(coefs)) {
  if(model_name=="Lasso" & coeff_kep[pix]<N_var){
    nb_months[pix] <- NA
  } else {
    Q <- try(get_firstcoeffs(coefs[[pix]], nb_of_coeff = N_var)[,2], silent = T)
    if(class(Q)=="try-error"){
      nb_months[pix] <- NA
    } else {
      nb_months[pix] <- length(unique(Q))
    }
  }#end if Lasso not enough coeff
}#end for pix


# Plot nb of different meteo variables in 3 first coeff ####
N_var <- 3
nb_meteo <- numeric()

for (pix in 1:length(coefs)) {
  if(model_name=="Lasso" & coeff_kep[pix]<N_var){
    nb_meteo[pix] <- NA
  } else {
    Q <- try(get_firstcoeffs(coefs[[pix]], nb_of_coeff = N_var)[,1], silent = T)
    if(class(Q)=="try-error"){
      nb_meteo[pix] <- NA
    } else {
      nb_meteo[pix] <- length(unique(Q))
    }
  }#end if Lasso not enough coeff
}#end for pix


DF_nbdiffmeteo <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_meteo = nb_meteo)
DF_nbdiffmeteo$nb_meteo <- as.factor(DF_nbdiffmeteo$nb_meteo)

ggplot(data = DF_nbdiffmeteo, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=nb_meteo)) +
  scale_color_manual(values = c("1"=viridis(3)[1], "2"=viridis(3)[2], "3"=viridis(3)[3])) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of different \nmeteo. variables",
       title = paste("Number of different meteorological variables in the ",N_var," first variables, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)



set.seed(2)
sample_try <- sample(1:pix_num, size = 10)
for (pix in sample_try) {
  plot(sort(coefs[[pix]][-1], decreasing = T), main=paste(model_name, "coeff, pixel", pix))
  abline(h=0)
}#end for pix



#Ratio nb pred month compared to growing season length: Take into account that pixels don't haev the same nb of month!!







list_coeff <- list()
nb_season_var <- matrix(data = NA, ncol=2)
for (pix in 1:pix_num) {
  if(length(coefs[[pix]]@i)-1>0){
    list_coeff[[pix]] <- get_firstcoeffs(coefs[[pix]],length(coefs[[pix]]@i)-1)
  } else {
    list_coeff[[pix]] <- NA
  }
}


count_seas_and_var <- function(coeff){
  if (length(coeff@i)-1) {
    coeff_kept <- get_firstcoeffs(coeff,length(coeff@i)-1)
    nb_var <- length(unique(coeff_kept[,1]))
    nb_of_seas <- 0
    
    if(sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Feb", "Dec", "Jan"))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    if(sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("May", "Mar", "Apr"))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    if(sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Jun", "Jul", "Aug"))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    if(sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Sep", "Nov", "Oct"))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    
    return(data.frame("nb_of_seas" = nb_of_seas, "nb_of_var" = nb_var))
  } else {
    return(data.frame("nb_of_seas" = NA, "nb_of_var" = NA))
  }#end if else
  
}

nb_of_seas <- sapply(X=coefs, FUN = count_seas_and_var)[1,]
nb_of_var <- sapply(X=coefs, FUN = count_seas_and_var)[2,]
nb_of_seas <- unlist(nb_of_seas)
nb_of_var <- unlist(nb_of_var)

DF_nbdiffmeteo <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_meteo = nb_of_var)
DF_nbdiffmeteo$nb_meteo <- as.factor(DF_nbdiffmeteo$nb_meteo)

ggplot(data = DF_nbdiffmeteo, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=DF_nbdiffmeteo$nb_meteo)) +
  scale_color_manual(values = c("1"=viridis(3)[1], "2"=viridis(3)[2], "3"=viridis(3)[3])) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of different \nmeteo. variables",
       title = paste("Number of different meteorological variables, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)



DF_nbdiffseason <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_season = nb_of_seas)
DF_nbdiffseason$nb_season <- as.factor(DF_nbdiffseason$nb_season)

ggplot(data = DF_nbdiffseason, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=DF_nbdiffseason$nb_season)) +
  scale_color_manual(values = c("1"=rainbow(4)[1], "2"=rainbow(4)[4], "3"=rainbow(4)[3], "4"=rainbow(4)[2])) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of different \nseasons",
       title = paste("Number of different seasons, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

