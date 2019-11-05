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

max_growingseason_length <- apply(X = growingseason_length, MARGIN = 1, FUN = max)

# #strange years
# plot(max_growingseason_length, xlab="pixel index", ylab="Max gorwing season length")
# strange_gs_pixels <- which(max_growingseason_length>365)
# strange_years_number <- numeric()
# i <- 0
# for (pix in strange_gs_pixels) {
#   i<-i+1
#   strange_years_number[i] <- length(which(growingseason_length[pix,]>365))
# }#end for pix
# 
# plot(strange_years_number/1600*100, ylab="% of years", xlab = "pixels",
#      main="% of years with growing season>365days\n for pixel with some years with growing season>365days")
# plot(growingseason_length[strange_gs_pixel[23],], ylab="growing season length", xlab="year",
#      main="growing season lengths (GSL)\n for 1 pixel with 51% of GSL>365")
# plot(growingseason_length[strange_gs_pixel[25],], ylab="growing season length", xlab="year",
#      main="growing season lengths (GSL)\n for 1 pixel with 5% of GSL>365")
# plot(growingseason_length[815,], ylab="growing season length", xlab="year",
#      main="growing season lengths (GSL)\n for 1 pixel with only 1 year with GSL>365")


#get rid of too long growing season
#I put yield and meteo data as NA for years with growing season>365days
years_to_remove <- growingseason_length>365
growingseason_length_corrected <- growingseason_length
sowing_date_corrected <- sowing_date

for (pix in 1:length(lat_kept)) {
  nb_years_to_rm <- length(which(years_to_remove[pix,]==T))
  
  yields[pix,years_to_remove[pix,]] <- rep(NA, nb_years_to_rm)
  sowing_date_corrected[pix,years_to_remove[pix,]] <- rep(NA, nb_years_to_rm)
  growingseason_length_corrected[pix,years_to_remove[pix,]] <- rep(NA, nb_years_to_rm)
  
  for (mon in 1:dim(vpd)[2]) {#mon=month
    vpd[pix,mon,years_to_remove[pix,]] <- rep(NA, nb_years_to_rm)
    tasmax[pix,mon,years_to_remove[pix,]] <- rep(NA, nb_years_to_rm)
    precip[pix,mon,years_to_remove[pix,]] <- rep(NA, nb_years_to_rm)
  }#end for mon
}#end for pix

min_sowing_day <- apply(X = sowing_date_corrected, MARGIN = 1, FUN = min)
min_sowing_date <- as.Date(min_sowing_day, origin="2019-01-01")
#Let's make sure that all growing season start the same year, before extracting the month
stopifnot(year(min_sowing_date)==2019)
sowing_month <- month(min_sowing_date)

max_growingseason_length <- apply(X = growingseason_length_corrected, MARGIN = 1, FUN = max)

#that's better:
plot(max_growingseason_length, xlab="pixel index", ylab="Max gorwing season length")

harvest_day <- sowing_date_corrected + growingseason_length_corrected

max_harvest_day <- apply(X = harvest_day, MARGIN = 1, FUN = max, na.rm=T)

max_harvest_date <- as.Date(max_harvest_day, origin="2019-01-01")
stopifnot(year(max_harvest_date)==2020)

harvest_month <- numeric(length = nb_pixel_kept)

for (loc in 1:nb_pixel_kept) {
    harvest_month[loc]<-month(max_harvest_date[loc]) + 12
}#end for loc

#remove useless months
months_to_keep <- matrix(data = NA, nrow = nb_pixel_kept,
                         ncol = dim(precip)[2])

for(loc in 1:nb_pixel_kept) {
  months_to_keep[loc,]<-((1:dim(precip)[2])+7>=sowing_month[loc] &
                           (1:dim(precip)[2])+7<=harvest_month[loc])
}#end for loc

#rapid check
for (ind in sample(1:nb_pixel_kept, size = 4)) {
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


for (pix in 1:nb_pixel_kept) {
  for (YY in 1:dim(vpd)[3]) {#YY = year
    month_to_rm <- which(months_to_keep[pix,]==F)
    nb_months_to_rm <- length(month_to_rm)
    vpd[pix,month_to_rm,YY] <- rep(NA, nb_months_to_rm)
    tasmax[pix,month_to_rm,YY] <- rep(NA, nb_months_to_rm)
    precip[pix,month_to_rm,YY] <- rep(NA, nb_months_to_rm)
  }#end for YY
}#end for pix

#Discard pixels with 0.025 percentile yield<1: 31 pixels
pix_to_keep <- which(apply(X=yields, MARGIN = 1, FUN = stats::quantile, probs=0.025, na.rm=T)>1)
pix_to_rm <- setdiff(1:dim(yields)[1], pix_to_keep)
#were these pixels already removed with GSL<=365? Not all of them
#pix_to_rm %in% strange_gs_pixels

lat_kept <- lat_kept[pix_to_keep]
lon_kept <- lon_kept[pix_to_keep]
yields <- yields[pix_to_keep,,]
tasmax <- tasmax[pix_to_keep,,]
vpd <- tasmax[pix_to_keep,,]

# #to plot map
# border <- readOGR('C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/ne_50m_admin_0_countries.shp')	