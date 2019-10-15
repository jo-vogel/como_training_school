# Load data ####
warning('Upload changed data (small changes')
library(ncdf4)
setwd("C:/Promotion/Training school") # path should be your local path to the file


nc_cy <- nc_open("C:/Promotion/Training school/Data/crop_yield_FR.nc")
nc_gsl <- nc_open("C:/Promotion/Training school/Data/crop_growingseason_length_FR.nc")
nc_dps<- nc_open("C:/Promotion/Training school/Data/meteo_dps_FR.nc")
nc_pr <- nc_open("C:/Promotion/Training school/Data/meteo_pr_FR.nc")
nc_sfcwind <- nc_open("C:/Promotion/Training school/Data/meteo_sfcwind_FR.nc")
nc_rsds <- nc_open("C:/Promotion/Training school/Data/meteo_rsds_FR.nc")
nc_tmax <- nc_open("C:/Promotion/Training school/Data/meteo_tasmax_FR.nc")
nc_tmin <- nc_open("C:/Promotion/Training school/Data/meteo_tasmin_FR.nc")  

cy_FR <- ncvar_get(nc_cy)
gsl_FR <- ncvar_get(nc_gsl)
dps_FR <- ncvar_get(nc_dps)
pr_FR <- ncvar_get(nc_pr)
sfcwind_FR <- ncvar_get(nc_sfcwind)
rsds_FR <- ncvar_get(nc_rsds)
tasmax_FR <- ncvar_get(nc_tmax)
tasmin_FR <- ncvar_get(nc_tmin)

nc_close(nc_cy)
nc_close(nc_gsl)
nc_close(nc_dps)
nc_close(nc_pr)
nc_close(nc_sfcwind)
nc_close(nc_rsds)
nc_close(nc_tmax)
nc_close(nc_tmin)

list_FR <- list(dps_FR,pr_FR,sfcwind_FR,rsds_FR,tasmax_FR,tasmin_FR)

# plot first five years of the data (daily values)
plot(tasmax_FR[,1,1], type ="l")



## Adjust growing season #####
gs_start <- 304 # Julian date of start of season
gs_length <- 310 # growing season length
# gs_end <- gs_length-(365-gs_start) # Julian date of end of season
num_years <- 4

cube_FR_2d <- array(NA,c(1826,400,6))  # 365*4+366: 5 years
cube_FR_gs_4years <- array(NA,c(gs_length*num_years,400,6)) # Number of obs., sim., var.
cube_FR <- array(NA,c(1600,gs_length,6))

for (m in 1:length(list_FR)){ # all 6 daily variables
  
  
  # reshape from 3 dimensions to 2 dimensions
  cube_FR_2d[,,m] <- matrix(list_FR[[m]],dim(list_FR[[m]])[1],dim(list_FR[[m]])[2]*dim(list_FR[[m]])[3]) # reshape from 3d to 2d
  
  days <-c(365,366,365,365,365) # 2011-2015 (leap year 2012)
  days <- cumsum(days) # cumulative sum
  startvec <- c(gs_start,days[1:3]+gs_start) # start of growing season for each year
  endvec <- startvec+gs_length-1 # end of growing season for each year
  plot(startvec,ylim=c(1,1826))
  points(endvec,col='green')
  vec_temp <- sapply(1:4, function(x) {c((startvec[x]):(endvec[x]))}) # indices of the days of each growing season of the years 2011-2015
  vec_temp <- as.vector(vec_temp)
  cube_FR_gs_4years[,,m] <-  cube_FR_2d[vec_temp,,m] # cut to growing season
  
  
  # reshape: time series of 4 years per column are separated into 1 year per row
  M2 <- NULL
  for (j in 1:400){ # simulations
    
    M1<-matrix(NA,4,gs_length)
    k<-1
    s<-gs_length
    for (i in 1:4) # years
    {
      M1[i,]<-cube_FR_gs_4years[k:s,j,m]
      k<-s+1 # start of next year
      s<-s+gs_length # end of next year
    }
    
    M2<-rbind(M2,M1)
  }
  cube_FR[,,m] <- M2
  
}

# crop yield
cy_FR_2d <- matrix(cy_FR,dim(cy_FR)[1],dim(cy_FR)[2]*dim(cy_FR)[3]) # reshape from 3d to 2d
cy_FR_final <- as.vector(cy_FR_2d)

# growing season length
gsl_FR_2d <- matrix(gsl_FR,dim(gsl_FR)[1],dim(gsl_FR)[2]*dim(gsl_FR)[3]) # reshape from 3d to 2d
gsl_FR_final <- as.vector(gsl_FR_2d)

cy_gsl_FR <- cbind(cy_FR_final,gsl_FR_final)
colnames(cy_gsl_FR) <- c("Crop yield","Growing season length")


# write data into files
# All data in an RData-file
save(cube_FR,cy_gsl_FR,file='cube_FR.RData')
# Daily data
write.table(cube_FR[,,1],"Dps_FR.txt")
write.table(cube_FR[,,2],"Pr_FR.txt")
write.table(cube_FR[,,3],"Wind_FR.txt")
write.table(cube_FR[,,4],"Rsds_FR.txt")
write.table(cube_FR[,,5],"T_max_FR.txt")
write.table(cube_FR[,,6],"T_min_FR.txt")
# Yearly data
write.table(cy_gsl_FR,file="Crop_yield_Growing_season_length_FR.txt")


# Plotting

# Plot one exemplary year
plot(cube_FR[1,,1]) # dps
plot(cube_FR[1,,2]) # pr
plot(cube_FR[1,,3]) # wind
plot(cube_FR[1,,4]) # rsds
plot(cube_FR[1,,5]) # tmax
plot(cube_FR[1,,6]) # tmin

# Plot one exemplary time step
plot(cube_FR[,1,1])