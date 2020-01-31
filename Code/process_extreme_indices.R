############################################################################
###########                      Process DATA                    ###########
###########         global (standardised) crop model data        ###########
###########         monthly meteo var and extreme indices        ###########
############################################################################


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))


##### Initialisation, librairies, data #####

library(ncdf4);library(BBmisc)


#############################
##### Standardised data #####
#############################




# Get the data
#Pauline's Laptop
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
# path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"


#Monthly Meteovar
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


#extreme indices
nc_extr_ind <- nc_open(paste0(path_to_NH_files,"/extreme_indices_1.nc"))
all_lat_glob <- ncvar_get(nc_extr_ind, varid = "lat")
all_lon_glob <- ncvar_get(nc_extr_ind, varid = "lon")

reshape_longitude <- function(lon){
  if(lon < -180 | lon > 360) return("lon should be in [-180;360]")
  else if(lon>180 ) return(lon-360)
  else return(lon)
}

reshaped_lon <- apply(all_lon_glob, FUN=reshape_longitude, MARGIN = 1)

range_lat <- range(which(round(all_lat_glob,3) %in% round(unique(lat_subset), 3)))
range_lon <- range(which(round(reshaped_lon,3) %in% round(unique(lon_subset), 3)))

dtr <- ncvar_get(nc_extr_ind, varid = "dtr", start = c(range_lon[1], range_lat[1], 1),
                 count=c((range_lon[2]-range_lon[1]+1),( range_lat[2]-range_lat[1]+1), -1))
frs <- ncvar_get(nc_extr_ind, varid = "frs", start = c(range_lon[1], range_lat[1], 1),
                 count=c((range_lon[2]-range_lon[1]+1),( range_lat[2]-range_lat[1]+1), -1))
txx <- ncvar_get(nc_extr_ind, varid = "txx", start = c(range_lon[1], range_lat[1], 1),
                 count=c((range_lon[2]-range_lon[1]+1),( range_lat[2]-range_lat[1]+1), -1))
tnn <- ncvar_get(nc_extr_ind, varid = "tnn", start = c(range_lon[1], range_lat[1], 1),
                 count=c((range_lon[2]-range_lon[1]+1),( range_lat[2]-range_lat[1]+1), -1))
rx5 <- ncvar_get(nc_extr_ind, varid = "rx5", start = c(range_lon[1], range_lat[1], 1),
                 count=c((range_lon[2]-range_lon[1]+1),( range_lat[2]-range_lat[1]+1), -1))
tx90p <- ncvar_get(nc_extr_ind, varid = "tx90p", start = c(range_lon[1], range_lat[1], 1),
                   count=c((range_lon[2]-range_lon[1]+1),( range_lat[2]-range_lat[1]+1), -1))
tn10p <- ncvar_get(nc_extr_ind, varid = "tn10p", start = c(range_lon[1], range_lat[1], 1),
                   count=c((range_lon[2]-range_lon[1]+1),( range_lat[2]-range_lat[1]+1), -1))
nc_close(nc_extr_ind)

lon_partial_glob <- apply(all_lon_glob[range_lon[1]:range_lon[2]], FUN=reshape_longitude, MARGIN = 1)
lat_partial_glob <- apply(all_lat_glob[range_lat[1]:range_lat[2]], FUN=reshape_longitude, MARGIN = 1)

dtr_stand <- matrix(data = NA, nrow = length(lon_subset), ncol = dim(dtr)[3])
frs_stand <- dtr_stand ; txx_stand <- dtr_stand ; tnn_stand <- dtr_stand ; rx5_stand <- dtr_stand ;
tx90p_stand <- dtr_stand ; tn10p_stand <- dtr_stand

dtr_non_stand <- matrix(data = NA, nrow = length(lon_subset), ncol = dim(dtr)[3])
frs_non_stand <- dtr_non_stand ; txx_non_stand <- dtr_non_stand ; tnn_non_stand <- dtr_non_stand ;
rx5_non_stand <- dtr_non_stand ; tx90p_non_stand <- dtr_non_stand ; tn10p_non_stand <- dtr_non_stand

for (pix in 1:length(lon_subset)) {
  print("pixel")
  print(pix)
  num_lon <- which(round(lon_partial_glob,4) == round(lon_subset[pix],4))
  num_lat <- which(round(lat_partial_glob,3) == round(lat_subset[pix],3))
  
  dtr_non_stand[pix,] <- dtr[num_lon,num_lat,];dtr_stand[pix,] <- normalize(dtr_non_stand[pix,], method = "range", range=c(-1,1))
  frs_non_stand[pix,] <- frs[num_lon,num_lat,];frs_stand[pix,] <- normalize(frs_non_stand[pix,], method = "range", range=c(-1,1))
  txx_non_stand[pix,] <- txx[num_lon,num_lat,];txx_stand[pix,] <- normalize(txx_non_stand[pix,], method = "range", range=c(-1,1))
  tnn_non_stand[pix,] <- tnn[num_lon,num_lat,];tnn_stand[pix,] <- normalize(tnn_non_stand[pix,], method = "range", range=c(-1,1))
  rx5_non_stand[pix,] <- rx5[num_lon,num_lat,];rx5_stand[pix,] <- normalize(rx5_non_stand[pix,], method = "range", range=c(-1,1))
  tx90p_non_stand[pix,] <- tx90p[num_lon,num_lat,];tx90p_stand[pix,] <- normalize(tx90p_non_stand[pix,], method = "range", range=c(-1,1))
  tn10p_non_stand[pix,] <- tn10p[num_lon,num_lat,];tn10p_stand[pix,] <- normalize(tn10p_non_stand[pix,], method = "range", range=c(-1,1))
  
  
}#end for pix


colnames(tasmax) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                      "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                      "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")
colnames(tasmax_stand) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                      "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                      "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")

colnames(vpd) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                   "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                   "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
colnames(vpd_stand) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                   "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                   "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")

colnames(pr) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                   "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                   "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")
colnames(pr_stand) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                  "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                  "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")

Data_non_standardized <- list(longitudes=lon_subset, latitudes=lat_subset,
                              yield=yield, precipitation=pr, tasmax = tasmax,
                              vpd=vpd, dtr=dtr_non_stand, frs=frs_non_stand, txx=txx_non_stand,
                              tnn=tnn_non_stand, rx5=rx5_non_stand, tx90p=tx90p_non_stand, tn10p=tn10p_non_stand)

Data_standardized <- list(longitudes=lon_subset, latitudes=lat_subset,
                              yield=yield_stand, precipitation=pr_stand, tasmax = tasmax_stand,
                              vpd=vpd_stand, dtr=dtr_stand, frs=frs_stand, txx=txx_stand,
                              tnn=tnn_stand, rx5=rx5_stand, tx90p=tx90p_stand, tn10p=tn10p_stand)


save(Data_non_standardized, file = paste0(path_to_NH_files,"/extremeindices_and_monthlymeteovar.Rdata"))
save(Data_standardized, file = paste0(path_to_NH_files,"/extremeindices_and_monthlymeteovar_rescaled.Rdata"))
