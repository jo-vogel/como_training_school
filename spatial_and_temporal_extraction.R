###################################################################
###########          Data processing                    ###########
###########   of global standardised crop model data    ###########
###########                                             ###########
###########       Author: Pauline Rivoire               ###########
###################################################################


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

##### Initialisation, librairies, data #####

library(glmnet);library(InformationValue)
library(ncdf4)
library(foreach);library(doParallel)


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

#yield dimension: 320*76*1600
#meteo var dimension: 320*76*17*1600
#Lat and long become only 1 dimension: every slice of constant latitude 
yields <- matrix(data=nh_variables[[which(nh_files=="crop_yield_NH.nc")]],
                 nrow=length(lati)*length(long),
                 ncol=dim(nh_variables[[which(nh_files=="crop_yield_NH.nc")]])[3])

lon <- rep(long,length(lati)) # coordinates rearranged
lat <- rep(lati,each=length(long))

spatial_indices_kept <- !is.na(yields[,1])

yields <- yields[spatial_indices_kept,]

lon_kept <- lon[spatial_indices_kept]
lat_kept <- lat[spatial_indices_kept]

growingseason_length <- matrix(data=nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]],
                               nrow=length(lati)*length(long),
                               ncol=dim(nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]])[3])[spatial_indices_kept,]

sowing_date <- matrix(data=nh_variables[[which(nh_files=="crop_sowing_date_NH.nc")]],
                      nrow=length(lati)*length(long),
                      ncol=dim(nh_variables[[which(nh_files=="crop_sowing_date_NH.nc")]])[3])[spatial_indices_kept,]

precip <- array(data = nh_variables[[which(nh_files=="meteo_pr_NH.nc")]],
                dim = c(length(lati)*length(long),
                        dim(nh_variables[[which(nh_files=="meteo_pr_NH.nc")]][3]),
                        dim(nh_variables[[which(nh_files=="meteo_pr_NH.nc")]][4])))[spatial_indices_kept,,]

tasmax <- array(data = nh_variables[[which(nh_files=="meteo_tasmax_NH.nc")]],
                dim = c(length(lati)*length(long),
                        dim(nh_variables[[which(nh_files=="meteo_tasmax_NH.nc")]][3]),
                        dim(nh_variables[[which(nh_files=="meteo_tasmax_NH.nc")]][4])))[spatial_indices_kept,,]

vpd <- array(data = nh_variables[[which(nh_files=="meteo_vpd_NH.nc")]],
             dim = c(length(lati)*length(long),
                     dim(nh_variables[[which(nh_files=="meteo_vpd_NH.nc")]][3]),
                     dim(nh_variables[[which(nh_files=="meteo_vpd_NH.nc")]][4])))[spatial_indices_kept,,]
