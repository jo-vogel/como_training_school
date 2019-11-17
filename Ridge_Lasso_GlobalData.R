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

#which segregation threshold for the model?
segreg_th <- 0.9

fitted.results_model <- lapply(1:test_length, function(x){ifelse(mypred[[x]] > segreg_th,1,0)})

mis_clas_err <- sapply(1:test_length, function(x){misClassError(actuals = as.matrix(y1_test_list[[x]]),
                                                                predictedScores=mypred[[x]],
                                                                threshold = segreg_th)})

# con_tab <-  lapply(1:test_length, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
#                                                                                 fitted.results_model[[x]],
#                                                                                 threshold = segreg_th)})
# sensi <- sapply(1:test_length, function(x){InformationValue::sensitivity(as.matrix(y1_test_list[[x]]),
#                                                                          fitted.results_model[[x]],
#                                                                          threshold = segreg_th)})
speci <- sapply(1:test_length, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                         fitted.results_model[[x]],
                                                                         threshold = segreg_th)})



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

DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], miscla = speci)

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





# Correlation between yield and specificity or miscla error ####

mean_yield <- apply(X=yield, MARGIN = 1, FUN = mean, na.rm=T)
plot(mean_yield, mis_clas_err,
     xlab="Mean Yield (kg/yr)", ylab="Miss-classification error",
     main=paste("Scatterplot mean yield, miss-classification error\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))
plot(mean_yield, speci,
     xlab="Mean Yield (kg/yr)", ylab="Specificity",
     main=paste("Scatterplot mean yield, Specificity\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))


pairs(cbind(mean_yield, mis_clas_err, speci),
      main=paste(model_name, " regression, bad yield thr.=", threshold,
                 "\nsegreg. thr.=", segreg_th, sep = ""))

cor(mean_yield, mis_clas_err)
cor(mean_yield, speci)


# Map of the number of coefficients kept (Lasso) #####
if(model_name=="Lasso"){
  coeff_kep <- numeric()
  
  for (pix in 1:pix_num) {
    coeff_kep[pix] <- sum(coefs[[pix]][row.names(coefs[[pix]])!="(Intercept)"]!=0)
  }#end for pix
  
  # Plot specificity error ####
  
  DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = coeff_kep)
  
  ggplot(data = DF_speci, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=coeff_kep)) +
    scale_color_gradient(limits=c(3,33),
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





# Map of the most important coefficient ####

source("get_first_coeff_function.R")

first_var <- numeric()
month_first_var <- numeric()

for (pix in 1:length(coefs)) {
  first_var[pix] <- get_firstcoeffs(coefs[[pix]])[,1]
  month_first_var[pix] <- get_firstcoeffs(coefs[[pix]])[,2]
}#end for pix

DF_1var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2],
                       first_var = first_var, month_first_var=month_first_var)

ggplot(data = DF_1var, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=first_var, alpha=month_first_var), size=1) +
  scale_color_manual(values = c("pr"="blue", "vpd"="green", "tmax"="red")) +
  scale_alpha_manual(values = c("Aug_Y1"=(1/13), "Sep_Y1"=(2/13),"Oct_Y1"=(3/13), "Nov_Y1"=4/13,"Dec_Y1"=(5/13),"Jan_Y2"=(6/13),
                                "Feb_Y2"=(7/13),"Mar_Y2"=(8/13),"Apr_Y2"=(9/13), "May_Y2"=(10/13),"Jun_Y2"=(11/13),
                                "Jul_Y2"=(12/13), "Aug_Y2"=(13/13)),
                    breaks = c("Aug_Y1", "Sep_Y1","Oct_Y1", "Nov_Y1","Dec_Y1","Jan_Y2", "Feb_Y2","Mar_Y2",
                               "Apr_Y2", "May_Y2","Jun_Y2","Jul_Y2", "Aug_Y2")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="1st var", alpha="Month",
       title = paste("Most important predictor: name and month, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  guides(color=guide_legend(order = 1), alpha=guide_legend(order = 2, ncol = 2))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 23, height = 7)



#Then: Map of the number of important months
#Ratio nb pred