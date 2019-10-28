##############################################################################
###########          Ridge and LASSO regression                    ###########
###########   on global standardised crop model data               ###########
###########                                                        ###########
###########       Author: Pauline Rivoire                          ###########
##############################################################################



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
                  function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
nh_variables <- lapply(1:length(nh_files),
                       function(x){ncvar_get(nh_data[[x]])})
lati <- ncvar_get(nh_data[[1]],"lat")
long <- ncvar_get(nh_data[[1]],"lon")
lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})  

# ## Adjust growing season #####
gs_start <- nh_variables[[which(nh_files=="crop_sowing_date_NH.nc")]][,,1] # Julian date of start of season
gs_length <- matrix(data=nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]],nrow=320*76,ncol=1600)
max_gs_length <- apply(gs_length,1,max) # maximum growing season length