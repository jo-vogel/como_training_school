############################################################################
###########                      Process DATA                    ###########
###########         global (standardised) crop model data        ###########
###########         monthly meteo var and extreme indices        ###########
###########                                             ###########
###########       Author: Pauline Rivoire               ###########
###################################################################
# It is structured in the following way:
# a) Load the Northern hemisphere monthly meteo ncdf (previously processed with "Data_processing_wo_extreme_indices.R")
# b) Load extreme indices, extract the gripoints corresponding to the one in NH meteovariabeles
# c) Build 2 arrays (rescaled and not) with both Monthly meteovariables and extreme indices
# d) Store the 2 arrays as Rdata files
# output=stored in drive: /Data/Global_data/extremeindices_and_monthlymeteovar_rescaled.Rdata and /Data/Global_data/extremeindices_and_monthlymeteovar.Rdata



# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))


##### Initialisation, librairies, data #####

library(ncdf4);library(BBmisc)


#############################
##### Standardised data #####
#############################




# Get the data: computed from Data_processing_wo_extreme_indices.R
# can be found on the Drive: Data/Global_data/Final_Data
#Pauline's Laptop
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
# path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"


load(paste0(path_to_NH_files,"/monthlymeteovar_995pixels.Rdata"))
load(paste0(path_to_NH_files,"/monthlymeteovar_rescaled_995pixels.Rdata"))
load(paste0(path_to_NH_files,"/monthlymeteovar_rescaled_with-mean_995pixels.Rdata"))

#Monthly Meteovar
yield <- Data_non_standardized$yield
tasmax <- Data_non_standardized$tasmax
vpd <-Data_non_standardized$vpd
pr <- Data_non_standardized$precipitation
lat_subset <- Data_non_standardized$latitudes
lon_subset <- Data_non_standardized$longitudes
yield_stand <- Data_standardized$yield
tasmax_stand <- Data_standardized$tasmax
vpd_stand <- Data_standardized$vpd
pr_stand <- Data_standardized$precipitation
yield_stand_mean <- Data_standardized_mean$yield
tasmax_stand_mean <- Data_standardized_mean$tasmax
vpd_stand_mean <- Data_standardized_mean$vpd
pr_stand_mean <- Data_standardized_mean$precipitation

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

dtr_stand_mean <- matrix(data = NA, nrow = length(lon_subset), ncol = dim(dtr)[3])
frs_stand_mean <- dtr_stand_mean ; txx_stand_mean <- dtr_stand_mean ; tnn_stand_mean <- dtr_stand_mean ;
rx5_stand_mean <- dtr_stand_mean ; tx90p_stand_mean <- dtr_stand_mean ; tn10p_stand_mean <- dtr_stand_mean

dtr_non_stand <- matrix(data = NA, nrow = length(lon_subset), ncol = dim(dtr)[3])
frs_non_stand <- dtr_non_stand ; txx_non_stand <- dtr_non_stand ; tnn_non_stand <- dtr_non_stand ;
rx5_non_stand <- dtr_non_stand ; tx90p_non_stand <- dtr_non_stand ; tn10p_non_stand <- dtr_non_stand

number_NA_years_extremes <- numeric(length = length(lon_subset))

for (pix in 1:length(lon_subset)) {
  num_lon <- which(round(lon_partial_glob,4) == round(lon_subset[pix],4))
  num_lat <- which(round(lat_partial_glob,3) == round(lat_subset[pix],3))
  
  dtr_non_stand[pix,] <- dtr[num_lon,num_lat,];dtr_stand[pix,] <- normalize(dtr_non_stand[pix,], method = "range", range=c(-1,1))
  dtr_stand_mean[pix,] <- base::scale(dtr_non_stand[pix,])
  frs_non_stand[pix,] <- frs[num_lon,num_lat,];frs_stand[pix,] <- normalize(frs_non_stand[pix,], method = "range", range=c(-1,1))
  frs_stand_mean[pix,] <- base::scale(frs_non_stand[pix,])
  txx_non_stand[pix,] <- txx[num_lon,num_lat,];txx_stand[pix,] <- normalize(txx_non_stand[pix,], method = "range", range=c(-1,1))
  txx_stand_mean[pix,] <- base::scale(txx_non_stand[pix,])
  tnn_non_stand[pix,] <- tnn[num_lon,num_lat,];tnn_stand[pix,] <- normalize(tnn_non_stand[pix,], method = "range", range=c(-1,1))
  tnn_stand_mean[pix,] <- base::scale(tnn_non_stand[pix,])
  rx5_non_stand[pix,] <- rx5[num_lon,num_lat,];rx5_stand[pix,] <- normalize(rx5_non_stand[pix,], method = "range", range=c(-1,1))
  rx5_stand_mean[pix,] <- base::scale(rx5_non_stand[pix,])
  tx90p_non_stand[pix,] <- tx90p[num_lon,num_lat,];tx90p_stand[pix,] <- normalize(tx90p_non_stand[pix,], method = "range", range=c(-1,1))
  tx90p_stand_mean[pix,] <- base::scale(tx90p_non_stand[pix,])
  tn10p_non_stand[pix,] <- tn10p[num_lon,num_lat,];tn10p_stand[pix,] <- normalize(tn10p_non_stand[pix,], method = "range", range=c(-1,1))
  tn10p_stand_mean[pix,] <- base::scale(tn10p_non_stand[pix,])
  
  #for some pixels: more NA in extreme indices: set all meteovar to NA
  
  LOC <- which(is.na(dtr_non_stand[pix,]) | is.na(frs_non_stand[pix,]) | is.na(txx_non_stand[pix,])
               | is.na(tnn_non_stand[pix,]) | is.na(rx5_non_stand[pix,]) | is.na(tx90p_non_stand[pix,])
               | is.na(tn10p_non_stand[pix,]))
  
  if(length(LOC) > 0){
    yield[pix,LOC] <- rep(NA, length(LOC))
    tasmax[pix,,LOC] <- rep(NA, length(LOC))
    pr[pix,,LOC] <- rep(NA, length(LOC))
    vpd[pix,,LOC] <- rep(NA, length(LOC))
    yield_stand[pix,LOC] <- rep(NA, length(LOC))
    tasmax_stand[pix,,LOC] <- rep(NA, length(LOC))
    pr_stand[pix,,LOC] <- rep(NA, length(LOC))
    vpd_stand[pix,,LOC] <- rep(NA, length(LOC))
    number_NA_years_extremes[pix] <- length(LOC)
    if(sum(!(LOC %in% which(is.na(dtr_non_stand[pix,]))))){
      print(pix)
      } # if nothing is printed then all info about NA is in dtr
  }
  
}#end for pix


# Store Datasets ####

colnames(tasmax) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                      "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                      "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")
colnames(tasmax_stand) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                      "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                      "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")
colnames(tasmax_stand_mean) <- c("tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                                 "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                                 "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2")

colnames(vpd) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                   "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                   "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
colnames(vpd_stand) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                   "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                   "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
colnames(vpd_stand_mean) <- c("vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                              "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                              "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")

colnames(pr) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                   "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                   "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")
colnames(pr_stand) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                  "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                  "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")
colnames(pr_stand_mean) <- c("pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                             "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                             "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")

Data_xtrm_non_standardized <- list(longitudes=lon_subset, latitudes=lat_subset,
                                   yield=yield, precipitation=pr, tasmax = tasmax,
                                   vpd=vpd, dtr=dtr_non_stand, frs=frs_non_stand, txx=txx_non_stand,
                                   tnn=tnn_non_stand, rx5=rx5_non_stand, tx90p=tx90p_non_stand, tn10p=tn10p_non_stand)

Data_xtrm_standardized <- list(longitudes=lon_subset, latitudes=lat_subset,
                               yield=yield_stand, precipitation=pr_stand, tasmax = tasmax_stand,
                               vpd=vpd_stand, dtr=dtr_stand, frs=frs_stand, txx=txx_stand,
                               tnn=tnn_stand, rx5=rx5_stand, tx90p=tx90p_stand, tn10p=tn10p_stand)

Data_xtrm_standardized_mean <- list(longitudes=lon_subset, latitudes=lat_subset,
                                    yield=yield_stand_mean, precipitation=pr_stand_mean, tasmax = tasmax_stand_mean,
                                    vpd=vpd_stand_mean, dtr=dtr_stand_mean, frs=frs_stand_mean, txx=txx_stand_mean,
                                    tnn=tnn_stand_mean, rx5=rx5_stand_mean, tx90p=tx90p_stand_mean, tn10p=tn10p_stand_mean)


save(Data_xtrm_non_standardized, file = paste0(path_to_NH_files,"/extremeindices_and_monthlymeteovar_995pix.Rdata"))
save(Data_xtrm_standardized, file = paste0(path_to_NH_files,"/extremeindices_and_monthlymeteovar_rescaled_995pix.Rdata"))
save(Data_xtrm_standardized_mean, file = paste0(path_to_NH_files,"/extremeindices_and_monthlymeteovar_rescaled_mean_995pix.Rdata"))






# Look at problematic years for pixels we will actully keep ####
pixels_with_na_in_xtrm <- cbind(number_na = number_NA_years_extremes[number_NA_years_extremes>0],
                                longitude = lon_subset[number_NA_years_extremes>0],
                                latitude = lat_subset[number_NA_years_extremes>0])

save(pixels_with_na_in_xtrm, file = paste0(path_to_NH_files,"/gridpoints_with_NA_in_extremes.Rdata"))
saveRDS(pixels_with_na_in_xtrm, file = paste0(path_to_NH_files,"/gridpoints_with_NA_in_extremes.rds"))
write.csv(pixels_with_na_in_xtrm, file = paste0(path_to_NH_files,"/gridpoints_with_NA_in_extremes.csv"))

read.csv(paste0(path_to_NH_files,"/gridpoints_with_NA_in_extremes.csv"))

load(paste0(path_to_NH_files,"/final_889pixels_coords.Rdata"))

both_criteria <- (1:995 %in% which(number_NA_years_extremes>0) & 1:995 %in% final_pixels_coord$ref_in_995)

pixels_kept_with_na_in_xtrm <- data.frame(number_NA = number_NA_years_extremes[both_criteria],
                                     longitude = lon_subset[both_criteria],
                                     latitude = lat_subset[both_criteria])

number_NA_years_extremes[final_pixels_coord$ref_in_995]
sum(number_NA_years_extremes[final_pixels_coord$ref_in_995])
sum(number_NA_years_extremes[final_pixels_coord$ref_in_995]>0)

NA_in_yield <- numeric()
ind <- 0
for (pix in 1:889) {
  if(sum(is.na(Data_non_standardized$yield[final_pixels_coord$ref_in_995[pix],]))>0) {
    ind <- ind + 1
    NA_in_yield[ind] <- final_pixels_coord$ref_in_995[pix]
  }
}

NA_in_yield %in% which(number_NA_years_extremes>0)[which(number_NA_years_extremes>0) %in% final_pixels_coord$ref_in_995]
#79 pixels with NA due to GS too long, and these pixels are contained in the 82 pixels with NA years in extreme indices
length(NA_in_yield)

years_with_na_yield <- list()
for (piiix in 1:length(NA_in_yield)) {
  years_with_na_yield[[piiix]] <- which(is.na(Data_non_standardized$yield[NA_in_yield[piiix],]))
}
print(paste("total nb years with na because GS too long:", length(unlist(years_with_na_yield))))


#Check if years rm because of GS too long are included in years with NA in extreme
for (piiix in 1:length(NA_in_yield)) {
  if(sum(1-(years_with_na_yield[[piiix]] %in% which(is.na(dtr_non_stand[NA_in_yield[piiix],]) | is.na(frs_non_stand[NA_in_yield[piiix],]) |
                                                    is.na(txx_non_stand[NA_in_yield[piiix],]) | is.na(tnn_non_stand[NA_in_yield[piiix],]) |
                                                    is.na(rx5_non_stand[NA_in_yield[piiix],]) | is.na(tx90p_non_stand[NA_in_yield[piiix],]) |
                                                    is.na(tn10p_non_stand[NA_in_yield[piiix],]))
  )
  )
  )
    print(NA_in_yield[piiix])
}
#No output=>each of the 79 pixels, the years removed because of GS too long are included in the years removed because of NA in extreme indices


#Check in the 79 pixels if years with NA in extreme are included in years rm because of GS too long
for (piiix in 1:length(NA_in_yield)) {
  print(sum(1-(which(is.na(dtr_non_stand[NA_in_yield[piiix],]) | is.na(frs_non_stand[NA_in_yield[piiix],]) |
                       is.na(txx_non_stand[NA_in_yield[piiix],]) | is.na(tnn_non_stand[NA_in_yield[piiix],]) |
                       is.na(rx5_non_stand[NA_in_yield[piiix],]) | is.na(tx90p_non_stand[NA_in_yield[piiix],]) |
                       is.na(tn10p_non_stand[NA_in_yield[piiix],])) %in% years_with_na_yield[[piiix]])))
}
#No, there are 7+6+4+8+1+31 years for which the GS season too long doesn't explain everything. (same for the 82-79=3 remaining pixels)

#Let's check the 3 pixels remaining

pix_NA_extr_final <- which(number_NA_years_extremes>0)[which(number_NA_years_extremes>0) %in% final_pixels_coord$ref_in_995]

pix_NA_extr_final[!(pix_NA_extr_final %in% NA_in_yield)]
number_NA_years_extremes[pix_NA_extr_final[!(pix_NA_extr_final %in% NA_in_yield)]]
# 1 and 2 and 1 years

print(paste("total nb years with na because of extreme indices:", sum(number_NA_years_extremes[pix_NA_extr_final])))
