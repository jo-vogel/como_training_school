###################################################################
###########          Data processing                    ###########
###########   of global standardised crop model data    ###########
###########                                             ###########
###########       Author: Pauline Rivoire               ###########
###################################################################
# It is structured in the following way:
# a) Load the raw Northern hemisphere monthly meteo gridded dataset
# b) Extract gridpoints containing data, and keep in memory mean yield before Data processing
# c) Set to NA the years where growing season>365 days
# d) Get rid of meteo data for months out of the growing season (GS)
# e) Get rid of pixels with meanyield<10th perc and pixels with 0.05 quantile == 0
# f) Scale the data between -1 and 1, thanks to the function FUN = normalize, method = "range", range=c(-1,1)
# g) Build 2 arrays (rescaled and not) with Monthly meteovariables
# h) Store the 2 arrays as Rdata files
# output=stored in drive: /Data/Global_data/monthlymeteovar_V2020-03-20.Rdata and /Data/Global_data/monthlymeteovar_rescaled_V2020-03-20.Rdata



# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

##### Initialisation, librairies, data #####

library(glmnet);library(InformationValue)
library(ncdf4);library(lubridate)
library(foreach);library(doParallel)
library(abind); library(raster); library(BBmisc); library(ggplot2)

# Get the data
#Pauline
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
#Johannes
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
#Cristina
path_to_NH_files <-"C:/Users/39349/Documents/DAMOCLES/Data_global"

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

mean_GS_before_rm_GStoolong <- apply(X = growingseason_length, MARGIN = 1, FUN = mean)

Mean_GSlength_before_rm_GStoolong <- cbind(mean_GS_before_rm_GStoolong, lon_kept, lat_kept)

save(Mean_GSlength_before_rm_GStoolong, file = paste0(path_to_NH_files,"/mean_growingseasonlength_before_process.Rdata"))


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

#save some space
rm(nh_variables)


# Which pixels have with mean yield < 10th perc ####

raw_mean_yield <- apply(yields, MARGIN = 1, FUN=mean, na.rm=T)
raw_sd_yield <- apply(yields, MARGIN = 1, FUN=sd, na.rm=T)

Raw_mean_yield <- cbind(raw_mean_yield, lon_kept, lat_kept)
Raw_sd_yield <- cbind(raw_sd_yield, lon_kept, lat_kept)
colnames(Raw_mean_yield) <- c("mean_yield","longitudes", "latitudes")
colnames(Raw_sd_yield) <- c("sd_yield","longitudes", "latitudes")
save(Raw_mean_yield, file = paste0(path_to_NH_files,"/RawMeanYield_995GP.Rdata"))
save(Raw_sd_yield, file = paste0(path_to_NH_files,"/RawSdYield_995GP.Rdata"))


perc10_yields <- quantile(raw_mean_yield, probs=0.1)

pix_yield_above_10thperc <- which(raw_mean_yield>perc10_yields)


#GP_kept_after_10thperc <- cbind(lon_kept[raw_mean_yield>perc10_yields], lat_kept[raw_mean_yield>perc10_yields])
#colnames(GP_kept_after_10thperc) <- c("longitudes", "latitudes")
#save(GP_kept_after_10thperc, file = paste0(path_to_NH_files,"/895gridpoints_kept_after10thpercyield.Rdata"))

##### temporal reduction ####
min_sowing_day <- apply(X = sowing_date, MARGIN = 1, FUN = min, na.rm=T)
min_sowing_date <- as.Date(min_sowing_day, origin="2019-01-01")
#Let's make sure that all growing season start the same year, before extracting the month
stopifnot(year(min_sowing_date)==2019)
sowing_month <- month(min_sowing_date)
# write.csv2(sowing_month,file="./Code/Workspaces/sowing_month.csv")

max_growingseason_length <- apply(X = growingseason_length, MARGIN = 1, FUN = max)


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

min_sowing_day <- apply(X = sowing_date_corrected, MARGIN = 1, FUN = min, na.rm=T)
min_sowing_date <- as.Date(min_sowing_day, origin="2019-01-01")
#Let's make sure that all growing season start the same year, before extracting the month
stopifnot(year(min_sowing_date[!is.na(min_sowing_date)])==2019)
sowing_month <- month(min_sowing_date)

max_growingseason_length <- apply(X = growingseason_length_corrected, MARGIN = 1, FUN = max, na.rm = T)

#that's better:
plot(max_growingseason_length, xlab="pixel index", ylab="Max growing season length")


mean_GS_after_rm_GStoolong <- apply(X = growingseason_length_corrected, MARGIN = 1, FUN = mean, na.rm=T)
Mean_GSlength_after_rm_GStoolong <- cbind(mean_GS_after_rm_GStoolong, lon_kept, lat_kept)
save(Mean_GSlength_after_rm_GStoolong, file = paste0(path_to_NH_files,"/mean_growingseasonlength_after_rm_GStoolong.Rdata"))

harvest_day <- sowing_date_corrected + growingseason_length_corrected

max_harvest_day <- apply(X = harvest_day, MARGIN = 1, FUN = max, na.rm=T)

max_harvest_date <- as.Date(max_harvest_day, origin="2019-01-01")
# write.csv2(max_harvest_date,file="./Code/Workspaces/max_harvest_date.csv")
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

# Which pixels have 5th percentile yield==0
fifthquantile_yield <- apply(X=yields, MARGIN = 1, FUN = stats::quantile, probs=0.05, na.rm=T)
pix_5th_perc_yield_gt_0 <- which(fifthquantile_yield!=0)


##### Final dataset to store #####


#standardize with range
message('gives warnings due to NAs')
yields_stand_range <- t(apply(yields, FUN = normalize, method = "range", range=c(-1,1), MARGIN = 1))

tasmax_stand_range <- aperm(apply(tasmax, FUN = normalize, method = "range", range=c(-1,1), MARGIN = c(1,2)),
                            perm=c(2,3,1))

vpd_stand_range <- aperm(apply(vpd, FUN = normalize, method = "range", range=c(-1,1), MARGIN = c(1,2)),
                         perm=c(2,3,1))

precip_stand_range <- aperm(apply(precip, FUN = normalize, method = "range", range=c(-1,1), MARGIN = c(1,2)),
                            perm=c(2,3,1))


yields_stand_mean <- t(apply(yields, FUN = base::scale, MARGIN = 1))

tasmax_stand_mean <- aperm(apply(tasmax, FUN = base::scale, MARGIN = c(1,2)), perm=c(2,3,1))

vpd_stand_mean <- aperm(apply(vpd, FUN = base::scale, MARGIN = c(1,2)), perm=c(2,3,1))

precip_stand_mean <- aperm(apply(precip, FUN = base::scale, MARGIN = c(1,2)), perm=c(2,3,1))


colnames(tasmax) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                      "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                      "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")
colnames(tasmax_stand_range) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                            "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                            "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")
colnames(tasmax_stand_mean) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                                  "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                                  "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")

colnames(vpd) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                   "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                   "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
colnames(vpd_stand_range) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                         "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                         "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
colnames(vpd_stand_mean) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                               "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                               "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")

colnames(precip) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                  "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                  "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")
colnames(precip_stand_range) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                        "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                        "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")
colnames(precip_stand_mean) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                                  "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                                  "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")

Data_non_standardized <- list(longitudes=lon_kept, latitudes=lat_kept,
                              yield=yields, precipitation=precip, tasmax = tasmax,
                              vpd=vpd)

Data_standardized <- list(longitudes=lon_kept, latitudes=lat_kept,
                          yield=yields_stand_range, precipitation=precip_stand_range, tasmax = tasmax_stand_range,
                          vpd=vpd_stand_range)

Data_standardized_mean <- list(longitudes=lon_kept, latitudes=lat_kept,
                               yield=yields_stand_mean, precipitation=precip_stand_mean,
                               tasmax = tasmax_stand_mean, vpd=vpd_stand_mean)



save(Data_non_standardized, file = paste0(path_to_NH_files,"/monthlymeteovar_995pixels.Rdata"))
save(Data_standardized, file = paste0(path_to_NH_files,"/monthlymeteovar_rescaled_995pixels.Rdata"))
save(Data_standardized_mean, file = paste0(path_to_NH_files,"/monthlymeteovar_rescaled_with-mean_995pixels.Rdata"))




pix_5th_perc_yield_gt_0_coord  <- data.frame(longitude=lon_kept[pix_5th_perc_yield_gt_0],
                                             latitude=lat_kept[pix_5th_perc_yield_gt_0],
                                             ref_in_995 = pix_5th_perc_yield_gt_0)
pix_yield_above_10thperc_coord  <- data.frame(longitude=lon_kept[pix_yield_above_10thperc],
                                              latitude=lat_kept[pix_yield_above_10thperc],
                                              ref_in_995 = pix_yield_above_10thperc)

final_pix_kept <- (1:nb_pixel_kept)[1:nb_pixel_kept %in% pix_5th_perc_yield_gt_0 & 1:nb_pixel_kept %in% pix_yield_above_10thperc]

final_pixels_coord <- data.frame(longitude=lon_kept[final_pix_kept],
                                 latitude=lat_kept[final_pix_kept],
                                 ref_in_995 = final_pix_kept)

save(final_pixels_coord, file = paste0(path_to_NH_files,"/final_889pixels_coords.Rdata"))
save(pix_5th_perc_yield_gt_0_coord, file = paste0(path_to_NH_files,"/969pixels_5th_perc_yield_gt_0_coords.Rdata"))
save(pix_yield_above_10thperc_coord, file = paste0(path_to_NH_files,"/895pixels_yield_above_10thperc_coords.Rdata"))
