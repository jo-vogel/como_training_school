##############################################################################
###########          Ridge and LASSO regression                    ###########
###########   on global standardised crop model data               ###########
###########                                                        ###########
###########       Author: Pauline Rivoire                          ###########
##############################################################################


# which method? model_name in c("Ridge", "Lasso)
model_name <- "Ridge"

#threshold for bad yields in c(0.025,0.05,0.1)
threshold <- 0.05

#############################
##### Standardised data #####
#############################

##### Initialisation, librairies, data #####

library(glmnet);library(InformationValue)

library(ncdf4)

library(foreach);library(doParallel)


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))




# Get the data
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"

nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),
                  FUN = function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
nh_variables <- lapply(1:length(nh_files),
                       FUN = function(x){ncvar_get(nh_data[[x]])})

lati <- ncvar_get(nh_data[[1]],"lat")
long <- ncvar_get(nh_data[[1]],"lon")
lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})  


#standardisation
T1 <- Sys.time()
yield_stand <- aperm(apply(nh_variables[[which(nh_files=="crop_yield_NH.nc")]],
                           FUN = scale, MARGIN = c(1,2)), perm = c(2,3,1))
T2 <- Sys.time()
difftime(T2,T1)

dim(nh_variables[[which(nh_files=="meteo_pr_NH.nc")]])


#get yield threshold
yield_treshold <- apply(nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]],
                        MARGIN = c(1,2), quantile, probs = threshold, na.rm=T)

# Adjust growing season
sowing_date <- nh_variables[[which(nh_files=="crop_sowing_date_NH.nc")]][,,1] # Julian date of start of season
gs_length <- array(data=nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]],
                   dim = dim(nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]]))
max_gs_length <- apply(gs_length,MARGIN = c(1,2),max) # maximum growing season length
