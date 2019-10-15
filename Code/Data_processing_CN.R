# Load data ####
library(ncdf4)
setwd("C:/Promotion/Training school") # path should be your local path to the file


nc_cy <- nc_open("C:/Promotion/Training school/Data/crop_yield_CN.nc")
nc_gsl <- nc_open("C:/Promotion/Training school/Data/crop_growingseason_length_CN.nc")
nc_dps<- nc_open("C:/Promotion/Training school/Data/meteo_dps_CN.nc")
nc_pr <- nc_open("C:/Promotion/Training school/Data/meteo_pr_CN.nc")
nc_sfcwind <- nc_open("C:/Promotion/Training school/Data/meteo_sfcwind_CN.nc")
nc_rsds <- nc_open("C:/Promotion/Training school/Data/meteo_rsds_CN.nc")
nc_tmax <- nc_open("C:/Promotion/Training school/Data/meteo_tasmax_CN.nc")
nc_tmin <- nc_open("C:/Promotion/Training school/Data/meteo_tasmin_CN.nc") 

cy_CN <- ncvar_get(nc_cy)
gsl_CN <- ncvar_get(nc_gsl)
dps_CN <- ncvar_get(nc_dps)
pr_CN <- ncvar_get(nc_pr)
sfcwind_CN <- ncvar_get(nc_sfcwind)
rsds_CN <- ncvar_get(nc_rsds)
tasmax_CN <- ncvar_get(nc_tmax)
tasmin_CN <- ncvar_get(nc_tmin)

nc_close(nc_cy)
nc_close(nc_gsl)
nc_close(nc_dps)
nc_close(nc_pr)
nc_close(nc_sfcwind)
nc_close(nc_rsds)
nc_close(nc_tmax)
nc_close(nc_tmin)

list_CN <- list(dps_CN,pr_CN,sfcwind_CN,rsds_CN,tasmax_CN,tasmin_CN)

# plot first five years of the data (daily values)
plot(tasmax_CN[,1,1], type ="l")



## Adjust growing season #####
gs_start <- 272 # Julian date of start of season
gs_length <- 294 # growing season length
# gs_end <- gs_length-(365-gs_start) # Julian date of end of season
num_years <- 4

cube_CN_2d <- array(NA,c(1826,400,6)) # 365*4+366: 5 years
cube_CN_gs_4years <- array(NA,c(gs_length*num_years,400,6))
cube_CN <- array(NA,c(1600,gs_length,6))

for (m in 1:length(list_CN)){ # all 6 daily variables
  
  
  # reshape from 3 dimensions to 2 dimensions
  # tasmax_CN_2d <- matrix(tasmax_CN,dim(tasmax_CN)[1],dim(tasmax_CN)[2]*dim(tasmax_CN)[3]) # reshape from 3d to 2d
  cube_CN_2d[,,m] <- matrix(list_CN[[m]],dim(list_CN[[m]])[1],dim(list_CN[[m]])[2]*dim(list_CN[[m]])[3]) # reshape from 3d to 2d
  
  days <-c(365,366,365,365,365) # 2011-2015 (leap year 2012)
  days <- cumsum(days) # cumulative sum
  startvec <- c(gs_start,days[1:3]+gs_start) # start of growing season for each year
  endvec <- startvec+gs_length-1 # end of growing season for each year
  plot(startvec,ylim=c(1,1826))
  points(endvec,col='green')
  vec_temp <- sapply(1:4, function(x) {c((startvec[x]):(endvec[x]))}) # indices of the days of each growing season of the years 2011-2015
  vec_temp <- as.vector(vec_temp)
  cube_CN_gs_4years[,,m] <-  cube_CN_2d[vec_temp,,m] # cut to growing season
  
  
  # reshape: time series of 4 years per column are separated into 1 year per row
  M2 <- NULL
  for (j in 1:400){ # simulations
    
    M1<-matrix(NA,4,gs_length)
    k<-1
    s<-gs_length
    for (i in 1:4) # years
    {
      M1[i,]<-cube_CN_gs_4years[k:s,j,m]
      k<-s+1 # start of next year
      s<-s+gs_length # end of next year
    }
    
    M2<-rbind(M2,M1)
  }
  cube_CN[,,m] <- M2
  
}

# crop yield
cy_CN_2d <- matrix(cy_CN,dim(cy_CN)[1],dim(cy_CN)[2]*dim(cy_CN)[3]) # reshape from 3d to 2d
cy_CN_final <- as.vector(cy_CN_2d)

# growing season length
gsl_CN_2d <- matrix(gsl_CN,dim(gsl_CN)[1],dim(gsl_CN)[2]*dim(gsl_CN)[3]) # reshape from 3d to 2d
gsl_CN_final <- as.vector(gsl_CN_2d)

cy_gsl_CN <- cbind(cy_CN_final,gsl_CN_final)
colnames(cy_gsl_CN) <- c("Crop yield","Growing season length")

# write data into files
# All data in an RData-file
save(cube_CN,cy_gsl_CN,file='cube_CN.RData')
# Daily data
write.table(cube_CN[,,1],"Dps_CN.txt")
write.table(cube_CN[,,2],"Pr_CN.txt")
write.table(cube_CN[,,3],"Wind_CN.txt")
write.table(cube_CN[,,4],"Rsds_CN.txt")
write.table(cube_CN[,,5],"T_max_CN.txt")
write.table(cube_CN[,,6],"T_min_CN.txt")
# Yearly data
write.table(cy_gsl_CN,file="Crop_yield_Growing_season_length_CN.txt")

# Plotting

# Plot one exemplary year
plot(cube_CN[1,,1]) # dps
plot(cube_CN[1,,2]) # pr
plot(cube_CN[1,,3]) # wind
plot(cube_CN[1,,4]) # rsds
plot(cube_CN[1,,5]) # tmax
plot(cube_CN[1,,6]) # tmin

# Plot one exemplary time step
plot(cube_CN[,1,1])