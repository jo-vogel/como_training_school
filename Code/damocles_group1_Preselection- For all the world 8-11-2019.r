  ### LEILA ###

# Lasso with interactions for northern hemisphere

library(ncdf4)
library(glinternet)
library(abind)
library(foreach)
library(doParallel)
library(InformationValue)
library(ROCR)
library(tictoc)

# Get the data ####
###################
# path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"

#path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"

path_to_NH_files <- "E:/POLIMI-RESEARCH/18. Training School/project-1/Big_MeteoData"

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

# load all coordinates of northern hemisphere
nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})

lat_all <- ncvar_get(nh_data[[1]],"lat")
lon_all <- ncvar_get(nh_data[[1]],"lon")
lati_all <- rep(lat_all,each=length(lon_all))
long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
coord_all <- cbind(long_all,lati_all)


# Process data ####
###################

yields_3dim <- array(yield,dim=c(965,1,1600));yields_stand_3dim <- array(yield,dim=c(965,1,1600))
Model_data <- abind(yields_3dim,tasmax,vpd,pr,along=2)
Model_data_stand <- abind(yields_stand_3dim,tasmax_stand,vpd_stand,pr_stand,along=2)


threshold <- 0.05
pix_num <- dim(Model_data)[1]
low_yield <- sapply(1:pix_num,function(x) {quantile(yield[x,],threshold,na.rm=T)})
cy <- t(sapply(1:pix_num,function(x){ifelse(yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield


cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))
Model_data[,1,] <-cy_reshaped
Model_data_stand[,1,] <-cy_reshaped

columnnames <- c("Yield","pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                 "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                 "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2",
                 "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                 "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                 "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
                 "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                 "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                 "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
colnames(Model_data) <- columnnames
colnames(Model_data_stand) <- columnnames



# Exclude NA variable columns
na_col <- matrix(data=NA,nrow=pix_num,ncol=52)
for (j in 1:pix_num){
  for (i in 1:52){
    na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)


na_time <- vector("list",length=pix_num) # for each pixel, the positions of NAs over time
for (j in 1:pix_num){
  na_time[[j]] <- which(is.na(Model_data[j,1,])) # locations of years with NA values
}




###################################################################################################
##################################    Preselections    ############################################
threshold_low <- 0.05                                                                              
threshold_high <- 0.95                                                                             
pix_num <- dim(Model_data)[1]                                                                     
                                                                                                   
low_yield <- sapply(1:pix_num,function(x) {quantile(Model_data[x,1,],threshold_low,na.rm=T)})      
high_yield <- sapply(1:pix_num,function(x) {quantile(Model_data[x,1,],threshold_high,na.rm=T)})    
                                                                                                   
Nrow_low_yield <-sapply(1:pix_num,function(x) {which(Model_data[x,1,]<=low_yield[x])})             
Nrow_high_yield <-sapply(1:pix_num,function(x) {which(Model_data[x,1,]>=high_yield[x])})           
                                                                                                                                                                                               #  
###################################################################################################


Botomm_Up_Low <-  array(NA,c(pix_num,52,2))              # Botomm_Up for low yeild boxplot
Botomm_Up_High <- array(NA,c(pix_num,52,2))              # Botomm_Up for high yeild boxplot

Up <- 0.75
Bot <- 0.25

   Not_NA <- sapply(1:pix_num,function(x) {which(na_col[x,]==F)[2:length(which(na_col[x,]==F))]})   # False means entire column isnot NA
   
   columnnames <- c("Yield","pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                 "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                 "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2",
                 "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                 "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                 "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
                 "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                 "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                 "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2")
                 
   Separate_variable <- matrix(NA,ncol =52, nrow=pix_num) 
   Dif_Driver <- matrix(NA,ncol =52, nrow=pix_num)
   colnames(Separate_variable) <- columnnames
   colnames(Dif_Driver) <- columnnames
          
   for (x in 1:pix_num) {
   
    for(mm in 1:length(unlist(Not_NA[x]))){
   
   
    Botomm_Up_Low [x,Not_NA[[x]][mm],1]<-quantile(Model_data[x,Not_NA[[x]][mm],Nrow_low_yield[[x]]],Up)          # Up_low_yield 
     
    Botomm_Up_Low [x,Not_NA[[x]][mm],2]<-quantile(Model_data[x,Not_NA[[x]][mm],Nrow_low_yield[[x]]],Bot)         # Botomm_Low_yield
 
    Botomm_Up_High [x,Not_NA[[x]][mm],1]<-quantile(Model_data[x,Not_NA[[x]][mm],Nrow_high_yield[[x]]],Up)        # Up_high_yield 
    
    Botomm_Up_High [x,Not_NA[[x]][mm],2]<-quantile(Model_data[x,Not_NA[[x]][mm],Nrow_high_yield[[x]]],Bot)       # Botomm_high_yield
    
    
    if(Botomm_Up_High [x,Not_NA[[x]][mm],2]> Botomm_Up_Low[x,Not_NA[[x]][mm],1] | Botomm_Up_Low[x,Not_NA[[x]][mm],2]> Botomm_Up_High [x,Not_NA[[x]][mm],1]){Separate_variable[x,Not_NA[[x]][mm]]<-1
                                                                                                         Dif_Driver[x,Not_NA[[x]][mm]]<- max((Botomm_Up_High [x,Not_NA[[x]][mm],2]- Botomm_Up_Low[x,Not_NA[[x]][mm],1]),(Botomm_Up_Low[x,Not_NA[[x]][mm],2]- Botomm_Up_High [x,Not_NA[[x]][mm],1]))}else{
                                                                                                          Separate_variable[x,Not_NA[[x]][mm]]<-0
                                                                                                          Dif_Driver[x,Not_NA[[x]][mm]]<- 0}
                     }
                            
                     }



Driver_selected <-sapply(1:pix_num,function(x) {
                    
    if (length(which(Separate_variable[x,]==1))>0){Model_data[x,which(Separate_variable[x,]==1),]}})       # Preselected driver
    
   
    
                  