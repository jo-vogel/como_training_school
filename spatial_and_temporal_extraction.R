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
library(ncdf4);library(lubridate)
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

###### Spatial reduction #####

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
nb_pixel_kept <- length(lon_kept)

growingseason_length <- matrix(data=nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]],
                               nrow=length(lati)*length(long),
                               ncol=dim(nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]])[3])[spatial_indices_kept,]

sowing_date <- matrix(data=nh_variables[[which(nh_files=="crop_sowing_date_NH.nc")]],
                      nrow=length(lati)*length(long),
                      ncol=dim(nh_variables[[which(nh_files=="crop_sowing_date_NH.nc")]])[3])[spatial_indices_kept,]

precip <- array(data = nh_variables[[which(nh_files=="meteo_pr_NH.nc")]],
                dim = c(length(lati)*length(long),
                        dim(nh_variables[[which(nh_files=="meteo_pr_NH.nc")]])[3],
                        dim(nh_variables[[which(nh_files=="meteo_pr_NH.nc")]])[4]))[spatial_indices_kept,,]

tasmax <- array(data = nh_variables[[which(nh_files=="meteo_tasmax_NH.nc")]],
                dim = c(length(lati)*length(long),
                        dim(nh_variables[[which(nh_files=="meteo_tasmax_NH.nc")]])[3],
                        dim(nh_variables[[which(nh_files=="meteo_tasmax_NH.nc")]])[4]))[spatial_indices_kept,,]

vpd <- array(data = nh_variables[[which(nh_files=="meteo_vpd_NH.nc")]],
             dim = c(length(lati)*length(long),
                     dim(nh_variables[[which(nh_files=="meteo_vpd_NH.nc")]])[3],
                     dim(nh_variables[[which(nh_files=="meteo_vpd_NH.nc")]])[4]))[spatial_indices_kept,,]

#earn some space
rm(nh_variables)



##### temporal reduction ####
min_sowing_day <- apply(X = sowing_date, MARGIN = 1, FUN = min)
min_sowing_date <- as.Date(min_sowing_day, origin="2019-01-01")
#Let's make sure that all growing season start the same year, before extracting the month
stopifnot(year(min_sowing_date)==2019)
sowing_month <- month(min_sowing_date)

#get the latest harvest day

#for some station, max growing season is ver large and far from the rest of the distribution
plot(ecdf(growingseason_length[406,]))
plot(ecdf(growingseason_length[411,]))
plot(ecdf(growingseason_length[347,]))

# #get rid of too long growing season
# growingseason_length[growingseason_length>365] <- NA

harvest_day <- sowing_date + growingseason_length

max_harvest_day <- apply(X = harvest_day, MARGIN = 1, FUN = max)

max_harvest_date <- as.Date(max_harvest_day, origin="2019-01-01")
stopifnot(year(max_harvest_date)>=2020)
stopifnot(year(max_harvest_date)<=2021)

harvest_month <- numeric(length = nb_pixel_kept)

for (loc in 1:nb_pixel_kept) {
  if(year(max_harvest_date)[loc]==2020){
    harvest_month[loc]<-month(max_harvest_date[loc]) + 12
  } else {
    harvest_month[loc]<-month(max_harvest_date[loc]) + 24
  }#end ifelse
}#end for loc

#remove useless months
months_to_keep <- matrix(data = NA, nrow = nb_pixel_kept,
                         ncol = dim(precip)[2])

for(loc in 1:nb_pixel_kept) {
  months_to_keep[loc,]<-((1:dim(precip)[2])+7>=sowing_month[loc] &
                           (1:dim(precip)[2])+7<=harvest_month[loc])
}#end for loc

#rapid check
for (ind in sample(1:995, size = 4)) {
  print(paste("index", ind))
  print(paste("lat", lat[ind]))
  print(paste("lon", lon[ind]))
  print("min sowing date")
  print(min_sowing_date[ind])
  print("max harvest date")
  print(max_harvest_date[ind])
  print("months to keep")
  print(months_to_keep[ind,])
}#end for ind
plot(ecdf(growingseason_length[517,]))
