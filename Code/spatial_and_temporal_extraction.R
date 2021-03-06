###################################################################
###########          Data processing                    ###########
###########   of global standardised crop model data    ###########
###########                                             ###########
###########       Author: Pauline Rivoire               ###########
###################################################################
# It is structured in the following way:
# a) Load the raw Northern hemisphere monthly meteo gridded dataset
# b) Extract gridpoints containing data
# c) Set to NA the years where growing season>365 days
# d) Get rid of meteo data for months out of the growing season (GS)
# e) Get rid of pixels where 0.025 quantile == 0
# f) Scale the data between -1 and 1, thanks to the function FUN = normalize, method = "range", range=c(-1,1)
# g) Store the data (rescaled and not) as 2 netcdfiles, structured this way:
#     Variables:
#        yields (dim 965*1600=nb_pixels*nb_of_years)
#        precip (dim 965*17*1600=nb_pixels*max_nb_month*nb_of_years)
#        tasmax  (dim 965*nb_month_GS*1600)
#        vpd      (dim 965*nb_month_GS*1600)
#        latitude   (dim 965)
#     Dimensions:
#        longitude   (length 965)
#        year        (length 1600)
#        month      (length 17)
# output=stored in drive: /Data/Global_data/NH_yield_and_meteovar_range-scaled.nc and /Data/Global_data/NH_yield_and_meteovar.nc


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

##### Initialisation, librairies, data #####

library(glmnet);library(InformationValue)
library(ncdf4);library(lubridate)
library(foreach);library(doParallel)
library(abind); library(raster); library(BBmisc); library(ggplot2)

# Get the data
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"

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
stopifnot(year(min_sowing_date[!is.na(min_sowing_date)])==2019)
sowing_month <- month(min_sowing_date)

max_growingseason_length <- apply(X = growingseason_length_corrected, MARGIN = 1, FUN = max)

#that's better:
plot(max_growingseason_length, xlab="pixel index", ylab="Max growing season length")

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

#Discard pixels with 0.025 percentile yield==0: 31 pixels
# pix_to_keep <- which(apply(X=yields, MARGIN = 1, FUN = stats::quantile, probs=0.025, na.rm=T)>1)
pix_to_keep <- which(apply(X=yields, MARGIN = 1, FUN = stats::quantile, probs=0.025, na.rm=T)!=0)
pix_to_rm <- setdiff(1:dim(yields)[1], pix_to_keep)
#were these pixels already removed with GSL<=365? Not all of them
#pix_to_rm %in% strange_gs_pixels

# #to plot map
# border <- readOGR('C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/ne_50m_admin_0_countries.shp')


##### Final dataset to store #####

lat_kept <- lat_kept[pix_to_keep]
lon_kept <- lon_kept[pix_to_keep]
yields <- yields[pix_to_keep,]
tasmax <- tasmax[pix_to_keep,,]
vpd <- vpd[pix_to_keep,,]
precip <- precip[pix_to_keep,,]


##### Store in netcdf file #####

#standardize with range
message('gives warnings due to NAs')
yields_stand_range <- t(apply(yields, FUN = normalize, method = "range", range=c(-1,1), MARGIN = 1))

tasmax_stand_range <- aperm(apply(tasmax, FUN = normalize, method = "range", range=c(-1,1), MARGIN = c(1,2)),
                            perm=c(2,3,1))

vpd_stand_range <- aperm(apply(vpd, FUN = normalize, method = "range", range=c(-1,1), MARGIN = c(1,2)),
                         perm=c(2,3,1))

precip_stand_range <- aperm(apply(precip, FUN = normalize, method = "range", range=c(-1,1), MARGIN = c(1,2)),
                            perm=c(2,3,1))

setwd("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global")


ncfname1 <- "NH_yield_and_meteovar.nc"
ncfname2 <- "NH_yield_and_meteovar_range-scaled.nc"



# create and write the netCDF file -- ncdf4 version
# define dimensions
londim <- ncdim_def(name = "lon", units = "degrees", vals = lon_kept)
yeardim <- ncdim_def(name = "year", units = "year", vals = 1:1600)
monthdim <- ncdim_def(name = "month", units = "month",
                      vals = (1:17+7))

# define variables
fillvalue <- 1e32

yield_var <- ncvar_def(name = "yield", units = "kg/ha",
                       dim = list(londim, yeardim), missval = fillvalue,
                       prec = "double")

yield_stand_var <- ncvar_def(name = "yield", units = "kg/ha normalized [-1,1]",
                             dim = list(londim, yeardim), missval = fillvalue,
                             prec = "double")

vpd_var <- ncvar_def(name = "vpd", units = "Pa",
                     dim = list(londim, monthdim, yeardim),
                     missval = fillvalue, prec = "double")

vpd_stand_var <- ncvar_def(name = "vpd", units = "Pa normalized [-1,1]",
                           dim = list(londim, monthdim, yeardim),
                           missval = fillvalue, prec = "double")

tasmax_var <- ncvar_def(name = "tasmax", units = "K",
                        dim = list(londim, monthdim, yeardim),
                        missval = fillvalue, prec = "double")

tasmax_stand_var <- ncvar_def(name = "tasmax", units = "K normalized [-1,1]",
                              dim = list(londim, monthdim, yeardim),
                              missval = fillvalue, prec = "double")

precip_var <- ncvar_def(name = "pr", units = "mm/day",
                        dim = list(londim, monthdim, yeardim),
                        missval = fillvalue, prec = "double")

precip_stand_var <- ncvar_def(name = "pr", units = "mm/day normalized [-1,1]",
                              dim = list(londim, monthdim, yeardim),
                              missval = fillvalue, prec = "double")

lat_var <- ncvar_def(name = "lat", units = "degree",
                     dim = list(londim),
                     missval = fillvalue, prec = "double")

# create netCDF file and put arrays
ncout1 <- nc_create(ncfname1,list(yield_var, vpd_var,
                                  tasmax_var, precip_var,
                                  lat_var),force_v4=TRUE)

ncout2 <- nc_create(ncfname2,list(yield_stand_var, vpd_stand_var,
                                  tasmax_stand_var, precip_stand_var,
                                  lat_var),force_v4=TRUE)

# put variables
ncvar_put(nc = ncout1, varid = yield_var, vals = yields)
ncvar_put(ncout1,precip_var,precip)
ncvar_put(ncout1,tasmax_var,tasmax)
ncvar_put(ncout1,vpd_var,vpd)
ncvar_put(ncout1,lat_var,lat_kept)
ncvar_put(ncout2,yield_stand_var,yields_stand_range)
ncvar_put(ncout2,precip_stand_var,precip_stand_range)
ncvar_put(ncout2,tasmax_stand_var,tasmax_stand_range)
ncvar_put(ncout2,vpd_stand_var,vpd_stand_range)
ncvar_put(ncout2,lat_var,lat_kept)

nc_close(ncout1)
nc_close(ncout2)


#Checking these files

TRY <- nc_open(filename = "NH_yield_and_meteovar.nc")
print(TRY$var)
YIELDTRY <- ncvar_get(nc = TRY, varid = "yield")
YIELDTRY[1,]==yields[1,]


TRY2 <- nc_open(filename = "NH_yield_and_meteovar_range-scaled.nc")
PRSTANDTRY <- ncvar_get(nc = TRY2, varid = "pr")
TASMAXSTANDTRY <- ncvar_get(nc = TRY2, varid = "tasmax")

PRSTANDTRY[12,6,]==precip_stand_range[12,6,]

LATTRY2 <- ncvar_get(TRY2, varid = "lat")
LATTRY2 == lat_kept
TASMAXSTANDTRY[1,6,]==tasmax_stand_range[1,6,]

nc_close(TRY)
nc_close(TRY2)


##### Tiff method #####



coords <- cbind(lat_kept,lon_kept)
coords_3d <- array(coords,dim=c(965,2,1600))
yields_3dim <- array(yields,dim=c(965,1,1600))
Model_data <- abind(yields,coords_3d,tasmax,vpd,precip,along=2)

message('gives warnings, maybe due to NAs')
Model_data_stand <- array(data=NA,dim=dim(Model_data))
for (i in 1:dim(Model_data)[1]){
    # Model_data_stand[i,,] <- t(apply(Model_data_wheat[i,,],1,scale)) # reinvert dimensions (inverted by apply)
    Model_data_stand[i,,] <- t(apply(Model_data[i,,],1,normalize, method="range",range=c(-1,1))) # reinvert dimensions (inverted by apply)
}
  


# Create tif raster object
Model_data_reshape <- aperm(Model_data,perm=c(1,3,2))

Model_data_brick <- brick(Model_data_reshape)
Model_data_stack <- stack(Model_data_reshape)


# write it to hard drive
writeRaster(Model_data_brick,filename="Model_data_brick.tif", format="GTiff",overwrite=T)
writeRaster(Model_data_brick,filename="Model_data_stack.tif", format="GTiff",overwrite=T)
save(Model_data,file="Model_data.RData") # works


# load it
load("Model_data.RData") # works
message('tif files do not contain data for unknown reason')
Model_brick <- raster("Model_data_brick.tif")
Model_stack <- raster("Model_data_stack.tif")

































##### Rapid look at the growing season length #####
min_growingseas_length <- apply(X = growingseason_length[pix_to_keep,], MARGIN = 1, FUN = min)
mean_growingseas_length <- apply(X=growingseason_length[pix_to_keep,], MARGIN = 1, FUN = max)
mean_growingseas_length_corrected <- apply(X=growingseason_length_corrected[pix_to_keep,], MARGIN = 1, FUN = max, na.rm=T)
max_growingseas_length <- apply(X = growingseason_length[pix_to_keep,], MARGIN = 1, FUN = max)
max_growingseas_length_corrected <- apply(X = growingseason_length_corrected[pix_to_keep,], MARGIN = 1, FUN = max, na.rm=T)
range_growingseason_length <- max_growingseas_length - min_growingseas_length
range_growingseason_length_corrected <- max_growingseas_length_corrected - min_growingseas_length



# load all coordinates of northern hemisphere
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
lat_all <- ncvar_get(nh_data[[1]],"lat")
lon_all <- ncvar_get(nh_data[[1]],"lon")
lati_all <- rep(lat_all,each=length(lon_all))
long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
coord_all <- cbind(long_all,lati_all)

lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})

coord_subset <- cbind(lon_kept,lat_kept)
world <- map_data("world")


DF_min <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = min_growingseas_length)

ggplot(data = DF_min, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=min_growingseas_length)) +
  scale_color_gradient(low = "black", high = "turquoise") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Min GSL",
       title = paste("Minimum length of the growing season"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


DF_max <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = max_growingseas_length)

ggplot(data = DF_max, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=max_growingseas_length)) +
  scale_color_gradient2(low = "black", mid = "orange", high = "red", midpoint = 365) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Max GSL",
       title = paste("Maximum length of the growing season"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

DF_max_corr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = max_growingseas_length_corrected)

ggplot(data = DF_max_corr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=max_growingseas_length_corrected)) +
  scale_color_gradient(low = "black", high = "orange") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Max GSL",
       title = paste("Maximum length of the growing season (removing years GS>365 days)"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)



DF_mean <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = mean_growingseas_length)

ggplot(data = DF_mean, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=mean_growingseas_length)) +
  scale_color_gradient2(low = "black", mid = "green",
                        high = "white", midpoint = max(mean_growingseas_length_corrected)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Mean GSL",
       title = paste("Mean length of the growing season"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

DF_mean_corr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = mean_growingseas_length_corrected)

ggplot(data = DF_mean_corr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=mean_growingseas_length_corrected)) +
  scale_color_gradient(low = "black", high = "green") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Mean GSL",
       title = paste("Mean length of the growing season (removing years GS>365 days)"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)




DF_range <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = range_growingseason_length)

ggplot(data = DF_range, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=range_growingseason_length)) +
  scale_color_gradientn(colours = c("black", "brown",  "red", "yellow", "purple"),
                        breaks = c(min(range_growingseason_length_corrected),30, 150,
                                   max(range_growingseason_length_corrected),
                                   max(range_growingseason_length))) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Range GSL",
       title = paste("Range of the growing season length"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

DF_range_corr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = range_growingseason_length_corrected)

ggplot(data = DF_range_corr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=range_growingseason_length_corrected)) +
  scale_color_gradientn(colours = c("black", "brown",  "red", "yellow"),
                        breaks = c(min(range_growingseason_length_corrected),30, 150,
                                   max(range_growingseason_length_corrected))) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Range GSL",
       title = paste("Range of the growing season length (removing years GS>365 days)"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
