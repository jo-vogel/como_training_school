##############################################################################
###########          Ridge and LASSO regression                    ###########
###########   on global standardised crop model data               ###########
###########                                                        ###########
###########       Author: Pauline Rivoire                          ###########
##############################################################################


# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

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

image(nh_variables[[which(nh_files=="crop_yield_NH.nc")]][,,1])

##### Reduce the size (spatial and growing season length) #####
#index of lat with at least one non NA value
index_lati_kept <- (1:length(lati))[apply(nh_variables[[which(nh_files=="crop_yield_NH.nc")]],
                                          MARGIN = 2, FUN = max, na.rm=T)>0]
#index of lon with at least one non NA value
index_long_kept <- (1:length(long))[apply(nh_variables[[which(nh_files=="crop_yield_NH.nc")]],
                                          MARGIN = 1, FUN = max, na.rm=T)>0]


#index of growing season starting month (1=August)
sowing_month <- (floor(apply(nh_variables[[which(nh_files=="crop_sowing_date_NH.nc")]],
                             MARGIN = c(1,2), FUN = min)/30.5)-7)[index_long_kept,index_lati_kept]
month_length_growing <- (floor(apply(nh_variables[[which(nh_files=="crop_growingseason_length_NH.nc")]],
                                    MARGIN = c(1,2), FUN = max)/30.5)+1)[index_long_kept,index_lati_kept]


#reduce the variables
yield <-nh_variables[[which(nh_files=="crop_yield_NH.nc")]][index_long_kept,index_lati_kept,]
vpd <- nh_variables[[which(nh_files=="meteo_vpd_NH.nc")]][index_long_kept,index_lati_kept,,]
prec <- nh_variables[[which(nh_files=="meteo_pr_NH.nc")]][index_long_kept,index_lati_kept,,]
tasmax <- nh_variables[[which(nh_files=="meteo_tasmax_NH.nc")]][index_long_kept,index_lati_kept,,]
nb_years <- dim(prec)[4]
for (LAT in 1:length(index_lati_kept)) {
  for (LON in 1:length(index_long_kept)) {
    month_index <- 1:dim(vpd)[3]<sowing_month[LON,LAT] | 1:dim(vpd)[3]>(sowing_month[LON,LAT]+month_length_growing[LON,LAT]-1)
    if (!is.na(sum(month_index))){
      for (year in 1:nb_years) {
        vpd[LON, LAT,month_index,year] <- rep(NA, sum(month_index))
        prec[LON, LAT,month_index,year] <- rep(NA, sum(month_index))
        tasmax[LON, LAT,month_index,year] <- rep(NA, sum(month_index))
      }#enf for year
    }#END IF
  }#end for LON
}#end for LAT

#Explore the data
#nb points with 
sum(month_length_growing<=3, na.rm=T)
sum(month_length_growing==3, na.rm=T)

# Look at the erroneous pixel (see Erroneous pixel.R)
month_length_growing[which(index_long_kept==180),which(index_lati_kept==12)]
prec[which(index_long_kept==180),which(index_lati_kept==12),,1]





##### Standardisation #####
T1 <- Sys.time()
yield_stand <- aperm(apply(yield, FUN = scale, MARGIN = c(1,2)), perm = c(2,3,1))
vpd_stand <- aperm(apply(vpd, FUN = scale, MARGIN = c(1,2,3)), perm=c(2,3,4,1))
pr_stand <- aperm(apply(prec, FUN = scale, MARGIN = c(1,2,3)), perm=c(2,3,4,1))
tasmax_stand <- aperm(apply(tasmax, FUN = scale, MARGIN = c(1,2,3)), perm=c(2,3,4,1))
T2 <- Sys.time()
difftime(T2,T1)



#get yield threshold
yield_treshold <- apply(yield_stand,MARGIN = c(1,2), quantile, probs = threshold, na.rm=T)

#quick check: compare
yield_th_not_stand <- apply(yield,MARGIN = c(1,2), quantile, probs = threshold, na.rm=T)
image(yield_th_not_stand)
image((yield_treshold*apply(yield,MARGIN = c(1,2),FUN=sd, na.rm=T))+apply(yield,MARGIN = c(1,2),FUN=mean, na.rm=T))
