# Load data ####
library(ncdf4)
setwd("C:/Promotion/Training school") # path should be your local path to the file


nc_cy <- nc_open("C:/Promotion/Training school/Data/crop_yield_US.nc")
nc_gsl <- nc_open("C:/Promotion/Training school/Data/crop_growingseason_length_US.nc")
nc_dps<- nc_open("C:/Promotion/Training school/Data/meteo_dps_US.nc")
nc_pr <- nc_open("C:/Promotion/Training school/Data/meteo_pr_US.nc")
nc_sfcwind <- nc_open("C:/Promotion/Training school/Data/meteo_sfcwind_US.nc")
nc_rsds <- nc_open("C:/Promotion/Training school/Data/meteo_rsds_US.nc")
nc_tmax <- nc_open("C:/Promotion/Training school/Data/meteo_tasmax_US.nc")
nc_tmin <- nc_open("C:/Promotion/Training school/Data/meteo_tasmin_US.nc") 

cy_US <- ncvar_get(nc_cy)
gsl_US <- ncvar_get(nc_gsl)
dps_US <- ncvar_get(nc_dps)
pr_US <- ncvar_get(nc_pr)
sfcwind_US <- ncvar_get(nc_sfcwind)
rsds_US <- ncvar_get(nc_rsds)
tasmax_US <- ncvar_get(nc_tmax)
tasmin_US <- ncvar_get(nc_tmin)

nc_close(nc_cy)
nc_close(nc_gsl)
nc_close(nc_dps)
nc_close(nc_pr)
nc_close(nc_sfcwind)
nc_close(nc_rsds)
nc_close(nc_tmax)
nc_close(nc_tmin)

list_US <- list(dps_US,pr_US,sfcwind_US,rsds_US,tasmax_US,tasmin_US)

# plot first five years of the data (daily values)
plot(tasmax_US[,1,1], type ="l")



## Adjust growing season #####
gs_start <- 276 # Julian date of start of season
gs_length <- 315 # growing season length
# gs_end <- gs_start+gs_length-365-1 # Julian date of end of season
num_years <- 4

cube_US_2d <- array(NA,c(1826,400,6)) # 365*4+366: 5 years
cube_US_gs_4years <- array(NA,c(gs_length*num_years,400,6))
cube_US <- array(NA,c(1600,gs_length,6))

for (m in 1:length(list_US)){ # all 6 daily variables
  
  
  # reshape from 3 dimensions to 2 dimensions
  # tasmax_US_2d <- matrix(tasmax_US,dim(tasmax_US)[1],dim(tasmax_US)[2]*dim(tasmax_US)[3]) # reshape from 3d to 2d
  cube_US_2d[,,m] <- matrix(list_US[[m]],dim(list_US[[m]])[1],dim(list_US[[m]])[2]*dim(list_US[[m]])[3]) # reshape from 3d to 2d
  
  days <-c(365,366,365,365,365) # 2011-2015 (leap year 2012)
  days <- cumsum(days) # cumulative sum
  startvec <- c(gs_start,days[1:3]+gs_start) # start of growing season for each year
  endvec <- startvec+gs_length-1 # end of growing season for each year
  plot(startvec,ylim=c(1,1826))
  points(endvec,col='green')
  vec_temp <- sapply(1:4, function(x) {c((startvec[x]):(endvec[x]))}) # indices of the days of each growing season of the years 2011-2015
  vec_temp <- as.vector(vec_temp)
  cube_US_gs_4years[,,m] <-  cube_US_2d[vec_temp,,m] # cut to growing season
  
  
  # reshape: time series of 4 years per column are separated into 1 year per row
  M2 <- NULL
  for (j in 1:400){ # simulations
    
    M1<-matrix(NA,4,gs_length)
    k<-1
    s<-gs_length
    for (i in 1:4) # years
    {
      M1[i,]<-cube_US_gs_4years[k:s,j,m]
      k<-s+1 # start of next year
      s<-s+gs_length # end of next year
    }
    
    M2<-rbind(M2,M1)
  }
  cube_US[,,m] <- M2
  
}

# crop yield
cy_US_2d <- matrix(cy_US,dim(cy_US)[1],dim(cy_US)[2]*dim(cy_US)[3]) # reshape from 3d to 2d
cy_US_final <- as.vector(cy_US_2d)

# growing season length
gsl_US_2d <- matrix(gsl_US,dim(gsl_US)[1],dim(gsl_US)[2]*dim(gsl_US)[3]) # reshape from 3d to 2d
gsl_US_final <- as.vector(gsl_US_2d)

cy_gsl_US <- cbind(cy_US_final,gsl_US_final)
colnames(cy_gsl_US) <- c("Crop yield","Growing season length")

# write data into files
# All data in an RData-file
save(cube_US,cy_gsl_US,file='cube_US.RData')
# Daily data
write.table(cube_US[,,1],"Dps_US.txt")
write.table(cube_US[,,2],"Pr_US.txt")
write.table(cube_US[,,3],"Wind_US.txt")
write.table(cube_US[,,4],"Rsds_US.txt")
write.table(cube_US[,,5],"T_max_US.txt")
write.table(cube_US[,,6],"T_min_US.txt")
# Yearly data
write.table(cy_gsl_US,file="Crop_yield_Growing_season_length_US.txt")

# Plotting

# Plot one exemplary year
plot(cube_US[1,,1]) # dps
plot(cube_US[1,,2]) # pr
plot(cube_US[1,,3]) # wind
plot(cube_US[1,,4]) # rsds
plot(cube_US[1,,5]) # tmax
plot(cube_US[1,,6]) # tmin

# Plot one exemplary time step
plot(cube_US[,1,1])