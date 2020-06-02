##############################################################################
###########                 Final plots (glmnet)                   ###########
###########   NH crop model, monthly meteovar, XtrM indices        ###########
###########               Author: Pauline Rivoire                  ###########
##############################################################################
# It WAS structured in the following way: (structure to be updated)
# a) Load the desired model
# b) plot map CSI
# c) plot map specificity
# d) plot map nb of variable
# e) plot histo nb of variable
# f) plots map nb of extremes kept
# g) plot map nb of seasons
# h) plot barplot nb pix for each var
# i) plot difference CSI w/o xtrm indices





# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))


model_name <- "Lasso"

#which lambda?
lambda_VALS <- c("lambda.min", "lambda.1se")
lambda_NAMES <- c("lambdamin", "lambda1se")

#threshold for bad yields
threshold <- 0.05


##### Initialisation, librairies, data #####

library(ncdf4);library(glmnet);library(InformationValue);library(ROCR);library(ggpubr)
library(abind);library(stringr);library(tictoc);library(ggplot2);library(viridis)






##### Load standardized Data #####

# in the drive folder Data/Global_Data/Final_Data
# Pauline
path <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
# Johannes
path <- "D:/user/vogelj/Group_project/Data"
# Cristina
path <- "C:/Users/39349/Documents/DAMOCLES/Data_global"

load(paste0(path,"/extremeindices_and_monthlymeteovar_rescaled_995pix.Rdata"))
load(paste0(path,"/extremeindices_and_monthlymeteovar_995pix.Rdata"))


##### Process data #####
total_nb_pix <- length(Data_xtrm_non_standardized$longitudes)
yield_3dim <- array(Data_xtrm_standardized$yield,dim=c(total_nb_pix,1,1600))
dtr_3dim <- array(Data_xtrm_standardized$dtr,dim=c(total_nb_pix,1,1600))
frs_3dim <- array(Data_xtrm_standardized$frs,dim=c(total_nb_pix,1,1600))
txx_3dim <- array(Data_xtrm_standardized$txx,dim=c(total_nb_pix,1,1600))
tnn_3dim <- array(Data_xtrm_standardized$tnn,dim=c(total_nb_pix,1,1600))
rx5_3dim <- array(Data_xtrm_standardized$rx5,dim=c(total_nb_pix,1,1600))
tx90p_3dim <- array(Data_xtrm_standardized$tx90p,dim=c(total_nb_pix,1,1600))
tn10p_3dim <- array(Data_xtrm_standardized$tn10p,dim=c(total_nb_pix,1,1600))

Model_data <- abind(yield_3dim,dtr_3dim,frs_3dim,txx_3dim,tnn_3dim,rx5_3dim,tx90p_3dim,tn10p_3dim
                    ,Data_xtrm_standardized$tasmax,Data_xtrm_standardized$vpd,Data_xtrm_standardized$pr,along=2)
colnames(Model_data) <- c("Yield", "dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p",
                          "tmax_Aug_Y1","tmax_Sep_Y1","tmax_Oct_Y1","tmax_Nov_Y1","tmax_Dec_Y1","tmax_Jan_Y2",
                          "tmax_Feb_Y2","tmax_Mar_Y2","tmax_Apr_Y2","tmax_May_Y2","tmax_Jun_Y2","tmax_Jul_Y2",
                          "tmax_Aug_Y2","tmax_Sep_Y2","tmax_Oct_Y2","tmax_Nov_Y2","tmax_Dec_Y2",
                          "vpd_Aug_Y1","vpd_Sep_Y1","vpd_Oct_Y1","vpd_Nov_Y1","vpd_Dec_Y1","vpd_Jan_Y2",
                          "vpd_Feb_Y2","vpd_Mar_Y2","vpd_Apr_Y2","vpd_May_Y2","vpd_Jun_Y2","vpd_Jul_Y2",
                          "vpd_Aug_Y2","vpd_Sep_Y2","vpd_Oct_Y2","vpd_Nov_Y2","vpd_Dec_Y2",
                          "pr_Aug_Y1","pr_Sep_Y1","pr_Oct_Y1","pr_Nov_Y1","pr_Dec_Y1","pr_Jan_Y2",
                          "pr_Feb_Y2","pr_Mar_Y2","pr_Apr_Y2","pr_May_Y2","pr_Jun_Y2","pr_Jul_Y2",
                          "pr_Aug_Y2","pr_Sep_Y2","pr_Oct_Y2","pr_Nov_Y2","pr_Dec_Y2")


Yield <- Data_xtrm_standardized$yield
low_yield <- apply(Yield, MARGIN = 1, FUN=quantile, probs=threshold, na.rm=T)
cy <- t(sapply(1:total_nb_pix,function(x){ifelse(Yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield

cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))

Model_data[,1,] <-cy_reshaped


# Exclude NA variable columns
na_col <- matrix(data=NA,nrow=total_nb_pix,ncol=dim(Model_data)[2])
for (j in 1:total_nb_pix){
  for (i in 1:dim(Model_data)[2]){
    na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)

na_time <- vector("list",length=total_nb_pix) # for each pixel, the positions of NAs over time
for (j in 1:total_nb_pix){
  na_time[[j]] <- which(is.na(Model_data[j,1,])) # locations of years with NA values
}


##### Split data into training and testing data set #####
vec <- 1:1600

years_with_na <- vector("logical",length=total_nb_pix)
for (i in 1:total_nb_pix){
  years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
}

training_indices <- vector("list",length=total_nb_pix)
testing_indices <- vector("list",length=total_nb_pix)



seed=1994
train_size <- 70

set.seed(seed)
for (x in 1:total_nb_pix) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((1600-length(na_time[[x]]))*(train_size/100))))
    testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
  } else {
    training_indices[[x]] <- sort(sample(1:1600, size = floor(1600*(train_size/100))))
    testing_indices[[x]] <- (1:1600)[-training_indices[[x]]]    
  }
}


Training_Data <- lapply(1:total_nb_pix,function(x){Model_data[x,,training_indices[[x]]]})
Testing_Data <- lapply(1:total_nb_pix,function(x){Model_data[x,,testing_indices[[x]]]})

pix_in <- 1:total_nb_pix


x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand


var_num <- apply(non_na_col,1,sum)
numLevels_list <- sapply(1:total_nb_pix, function(x){ rep(1,times=var_num[x])})
for (i in 1:total_nb_pix){
  names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
}


##### Load the model #####

# On the Drive you can find my data in:
# Models/LASSO-Ridge regression/Lasso_glmnet_final_results/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_995pixels.RData
# Models/LASSO-Ridge regression/Lasso_glmnet_final_results/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_995pixels.RData


# Pauline
load(file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
seed, "_train", train_size,"_995pixels.Rdata"))
# load(file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
#                    seed, "_train", train_size,"_995pixels.Rdata"))


# Johannes
load(paste0("D:/user/vogelj/Group_project/Code/Workspaces/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
             seed, "_train", train_size,"_995pixels.Rdata"))
# load(paste0("D:/user/vogelj/Group_project/Code/Workspaces/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
#              seed, "_train", train_size,"_995pixels.Rdata"))



#Cristina

load(paste0("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
            seed, "_train", train_size,"_995pixels.Rdata"))
load(paste0("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
            seed, "_train", train_size,"_995pixels.Rdata"))


Model_chosen <- lasso_model_lambda1se
lambda_val <- lambda_VALS[2]
lambda_name <- lambda_NAMES[2]


##### Plot Raw mean yield, and pixels that were removed #####

# in the drive folder Data/Global Data/Final_Data

#Pixels kept
load(paste0(path,"/final_889pixels_coords.Rdata"))
#Raw mean yield
load(paste0(path,"/RawMeanYield_995GP.Rdata"))

world <- map_data("world")

DF_meanY <- data.frame(lon=Raw_mean_yield[,"longitudes"], lat = Raw_mean_yield[,"latitudes"],
                       meany = Raw_mean_yield[,"mean_yield"]/1000)

pixels_excluded <- as.logical(1-(1:total_nb_pix %in% final_pixels_coord$ref_in_995))
DF_excluded_pix <- data.frame(lon = Raw_mean_yield[pixels_excluded,"longitudes"],
                              lat = Raw_mean_yield[pixels_excluded,"latitudes"])

ggplot(data = DF_meanY, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3)+  geom_tile(aes(fill=DF_meanY$meany)) +
  scale_fill_gradient2(midpoint = max(DF_meanY$meany, na.rm = T)/2,
                       limits=c(0,max(DF_meanY$meany)),
                       low = "#f7fcb9", mid = "#addd8e", high = "#31a354") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-115, 130),
              ylim = c(min(DF_meanY$lat), max(DF_meanY$lat)),
              ratio = 1.3)+
  labs(fill="Mean yield\n(t/ha)"  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))+
  geom_point(data = DF_excluded_pix, aes(x = DF_excluded_pix$lon, y = DF_excluded_pix$lat),
             color = "black", size = 0.89, pch=4) +
  X11(width = 20, height = 6)



# ##### Map mean growing season length #####
# # load mean growing season, in the drive folder Data/Global Data
# load(paste0(path,"/mean_growingseasonlength_before_process.Rdata"))
# load(paste0(path,"/mean_growingseasonlength_after_rm_GStoolong.Rdata"))
# 
# world <- map_data("world")
# 
# #Choose before or after rm GS too long
# GSlength <- Mean_GSlength_before_rm_GStoolong
# GSlength <- Mean_GSlength_after_rm_GStoolong
# 
# world <- map_data("world")
# 
# DF_GSl <- data.frame(lon=GSlength[,"lon_kept"], lat = GSlength[,"lat_kept"],
#                      GSl = GSlength[,1])
# 
# ggplot(data = DF_GSl, aes(x=DF_GSl$lon, y=DF_GSl$lat)) +
#   geom_polygon(data = world, aes(long, lat, group=group),
#                fill="white", color="black", size=0.3) +
#   geom_tile(aes(fill=DF_GSl$GSl)) +
#   scale_fill_gradient2(midpoint = max(DF_GSl$GSl, na.rm = T)/2,
#                        #limits=c(0,1),
#                        low = "#bfd3e6", mid = "#8c96c6", high = "#810f7c") +
#   theme(panel.ontop = F, panel.grid = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
#   ylab("Lat (°N)") +
#   xlab("Lon (°E)") +
#   coord_fixed(xlim = c(-120, 135),
#               ylim = c(min(DF_GSl$lat)-1, max(DF_GSl$lat+1)),
#               ratio = 1.3)+
#   labs(fill="Mean growing\nseason length\n(days)"
#        #,title = paste("CSI, simple",model_name,"regression, "),
#        #subtitle = paste("Monthly meteo var + extreme indices, cutoff level=", round(segreg_th,3),", ",lambda_val, sep = "")
#   )+
#   theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
#         legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
#   X11(width = 20, height = 6)


##### Map number of month used in the analysis #####
world <- map_data("world")

coord_subset <- cbind(final_pixels_coord$longitude, final_pixels_coord$latitude)

final_pix_num <- length(final_pixels_coord$latitude)
nb_month_GS <- integer(length = final_pix_num)
for (pix in 1:final_pix_num) {
  pixel <- final_pixels_coord$ref_in_995[pix]
  nb_month_GS[pix] <- sum(substr(colnames(x1_train_list[[pixel]]), start = 1, stop = 3)=="pr_")
}

levels_nb_month <- cut(nb_month_GS, breaks = c(5,8,11,14), right = F)
levels_nb_month <- gsub(","," - ",levels_nb_month,fixed=TRUE)
DF_GSmonth <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], levels_nb_month = levels_nb_month)
DF_GSmonth$levels_nb_month <- gsub("\\[|\\)","",levels_nb_month)


# DF_GSmonth <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2],
#                          levels_nb_month = nb_month_GS)
# DF_GSmonth$levels_nb_month <- as.integer(DF_GSmonth$levels_nb_month)

ggplot(data = DF_GSmonth, aes(x=DF_GSmonth$lon, y=DF_GSmonth$lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF_GSmonth$levels_nb_month)) +
  scale_fill_manual(values=c("5 - 8" = "#edf8b1", "8 - 11" = "#41b6c4", "11 - 14" = "#2c7fb8" ),
                    breaks=c("5 - 8", "8 - 11", "11 - 14"),
                    label=c("5 - 8"="5 - 7", "8 - 11"="8 - 10", "11 - 14"="11 - 13")) +
  # scale_fill_gradient2(midpoint = 12,
  #                      #limits=c(0,1),
  #                      low = "#edf8b1", mid = "#41b6c4", high = "#081d58") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(DF_GSmonth$lat)-1, max(DF_GSmonth$lat+1)),
              ratio = 1.3)+
  labs(fill="Growing season\nlength (months)"
       #,title = paste("CSI, simple",model_name,"regression, "),
       #subtitle = paste("Monthly meteo var + extreme indices, cutoff level=", round(segreg_th,3),", ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 5)


##### Adjust cutoff level #####

# source("./Code/Simple_Lasso_Ridge_ElasticNet/cutoff_adj_glmnet_lambda1se.R")
# y1_train_list_simple_lasso <- y1_train_list
# x1_train_list_simple_lasso <- x1_train_list
# Model_chosen_889 <- list()
# 
# y1_train_list_simple_lasso <- list()
# x1_train_list_simple_lasso <- list()
# work_pix_tmp <- numeric()
# for (pixel in 1:final_pix_num) {
#   pix_in_995 <- final_pixels_coord$ref_in_995[pixel]
#   y1_train_list_simple_lasso[[pixel]] <- y1_train_list[[pix_in_995]]
#   x1_train_list_simple_lasso[[pixel]] <- x1_train_list[[pix_in_995]]
#   Model_chosen_889[[pixel]] <- Model_chosen[[pix_in_995]]
#   if(is.character(Model_chosen[[pix_in_995]])){work_pix_tmp[pixel]<-0} else {work_pix_tmp[pixel]<-1}
# }#end for pixel
# 
# 
# 
# cost_fp_simple_lasso <- 100 # Misses: this should be associated with a higher cost, as it is more detrimental
# cost_fn_simple_lasso <- 100 # False alarms
# 
# work_pix <- which(work_pix_tmp==1)
# library(pbapply)
# #return the mean value, over all pixels, of the adjusted cutoff
# cutoff_simple_lasso <- adjust_cutoff(model_vector = Model_chosen_889,x1_train_list = x1_train_list_simple_lasso, y1_train_list = y1_train_list_simple_lasso,
#                                      work_pix = work_pix, cost_fp = cost_fp_simple_lasso, cost_fn= cost_fn_simple_lasso)


segreg_th <- 0.6582418

##### Model performance assessment #####
world <- map_data("world")

coord_subset <- cbind(final_pixels_coord$longitude, final_pixels_coord$latitude)

final_pix_num <- length(final_pixels_coord$latitude)



coeff  <-list()
speci <- rep(NA, final_pix_num)
sensi <- rep(NA, final_pix_num)
csi <- rep(NA, final_pix_num)
mypred <- vector("list", final_pix_num)
fitted_bad_yield <- vector("list", final_pix_num)
nb_training_years <- rep(NA, final_pix_num)
nb_testing_years <- rep(NA, final_pix_num)

for (pixel in 1:final_pix_num) {
  pix <- final_pixels_coord$ref_in_995[pixel]
  coeff[[pixel]] <- coefficients(Model_chosen[[pix]])
  
  mypred[[pixel]] <- predict(Model_chosen[[pix]], as.matrix(x1_test_list[[pix]]),type="response")
  
  fitted_bad_yield[[pixel]] <- ifelse(mypred[[pixel]] > segreg_th,1,0)
  
  speci[pixel] <- InformationValue::specificity(actuals = as.matrix(y1_test_list[[pix]]),
                                                predictedScores = fitted_bad_yield[[pixel]],
                                                threshold = segreg_th)
  sensi[pixel] <- InformationValue::sensitivity(actuals = as.matrix(y1_test_list[[pix]]),
                                                predictedScores = fitted_bad_yield[[pixel]],
                                                threshold = segreg_th)
  
  con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(y1_test_list[[pix]]),
                                               predictedScores = fitted_bad_yield[[pixel]],
                                               threshold = segreg_th)
  csi[pixel] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
  if(is.na(con_tab["0","0"])){
    csi[pixel] <- 0
  }
  
  nb_training_years[pixel] <- length(training_indices[[pix]])
  nb_testing_years[pixel] <- length(testing_indices[[pix]])
  
}#end pixel
pixel_with_pb <- which(is.na(csi))


# plot(sort((1600 - (nb_training_years + nb_testing_years))[which(nb_training_years + nb_testing_years != 1600)]),
#      ylab="number of years with GS too long",
#      xlab="pixels with some GS too long\n(sorted by increasing nb)")

extreme_in_coeff <- function(coeff_list){ #function to check how many extreme indeices are kept as predictors
  extreme_indices <- c("dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p")
  if(max(abs(coeff_list[extreme_indices,]))==0){
    return(0)
  } else {
    return(length(which(abs(coeff_list[extreme_indices,])>0)))
  }
}#end func extreme_in_coeff


number_coeff_kept <- function(coeff_list){#give number of coeff !=0
  return(length(which(abs(coeff_list[-1,])>0)))
}

nb_extr_kept <- numeric()
nb_coeff_kept <- numeric()

for (pixel in 1:final_pix_num) {
  nb_extr_kept[pixel] <- extreme_in_coeff(coeff[[pixel]])
  nb_coeff_kept[pixel] <- number_coeff_kept(coeff[[pixel]])
}#end pixel



var_yield <- apply(Data_xtrm_non_standardized$yield,MARGIN = 1, FUN = var, na.rm=T)


##### Plot results ######
# # load all coordinates of northern hemisphere
# path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global" # Pauline
# # path_to_NH_files <- "D:/user/vogelj/Data/Group project Como" # Johannes
# nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
# nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
# lat_all <- ncvar_get(nh_data[[1]],"lat")
# lon_all <- ncvar_get(nh_data[[1]],"lon")
# lati_all <- rep(lat_all,each=length(lon_all))
# long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
# coord_all <- cbind(long_all,lati_all)
# 
# lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})

coord_subset <- cbind(final_pixels_coord$longitude, final_pixels_coord$latitude)

# # save final coordinates on Pauline's Laptop
# save(coord_subset, file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Final_889GP_coordinates.Rdata")
# saveRDS(coord_subset, file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Final_889GP_coordinates.rds")
# write.csv(coord_subset, file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Final_889GP_coordinates.csv")

world <- map_data("world")



# plot CSI=(hits)/(hits + misses + false alarm) ####
DF_csi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], csi = csi)

ggplot(data = DF_csi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=csi)) +
  scale_fill_gradient2(midpoint = max(csi, na.rm = T)/2,
                       #limits=c(0,1),
                       low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="CSI"
       #,title = paste("CSI, simple",model_name,"regression, "),
       #subtitle = paste("Monthly meteo var + extreme indices, cutoff level=", round(segreg_th,3),", ",lambda_val, sep = "")
       )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 5)

cor.test(Raw_mean_yield[final_pixels_coord$ref_in_995,"mean_yield"], csi)
cor.test(Raw_mean_yield[final_pixels_coord$ref_in_995,"mean_yield"], csi, method = "kendall")
par(mar=c(4.1,4.1,1,1))

plot(Raw_mean_yield[final_pixels_coord$ref_in_995,"mean_yield"]/1000, csi,
     col="black",pch=19, xlab="Mean yield [t/ha]", ylab="CSI", cex.lab=1.2, cex.axis=1.2,  cex.sub=1.2,
     font.lab=2) 
#plot(mean_yield[number_pix_to_keep_in_969], csi,col="black",pch=19, xlab="Mean yield [kg/ha]", ylab="CSI", cex.lab=1.2, cex.axis=1.2,  cex.sub=1.2) 




# number of season kept


source("./Code/get_first_coeff_function.R")

count_seas_and_var <- function(coefff){
  if (length(which((coefff)!=0))-1>0) {
    coeff_kept <- get_firstcoeffs(coefff, nb_of_coeff = length(which((coefff)!=0))-1)
    nb_of_seas <- 0
    
    if(sum(!is.na(coeff_kept[,2])) & (sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Feb", "Dec", "Jan")))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    if(sum(!is.na(coeff_kept[,2])) & sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("May", "Mar", "Apr"))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    if(sum(!is.na(coeff_kept[,2])) & sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Jun", "Jul", "Aug"))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    if(sum(!is.na(coeff_kept[,2])) & sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Sep", "Nov", "Oct"))){
      nb_of_seas <- nb_of_seas + 1
    }
    
    if (nb_of_seas>0){
      return(nb_of_seas)
    } else {
      return("No met. var")
    }
    
  } else {
    return("No met. var")
  }#end if else
  
}

nb_of_seas <- numeric()
for (pix in 1:final_pix_num) {
  nb_of_seas[pix] <- count_seas_and_var(coeff[[pix]])
}



nb_meteo <- numeric(length=length(coeff))
met_type <- numeric(length=length(coeff))
meteo_type <- as.character(met_type)
met_strings <- vector("list",length=length(coeff))
met_vpd <- numeric(length=length(coeff))
met_pr <- numeric(length=length(coeff))
met_temp <- numeric(length=length(coeff))

source("./Code/get_first_coeff_function.R")
for (pix in 1:length(coeff)) {
  if((length(which(coeff[[pix]]!=0))-1)<1){
    nb_meteo[pix] <- 0
    met_type[pix] <- 8
    meteo_type[pix] <- "None"
  } else {
    coeff_kept <- get_firstcoeffs(coeff = coeff[[pix]],
                                  nb_of_coeff = (length(which(coeff[[pix]]!=0))-1))
    nb_meteo[pix] <- length(unique(coeff_kept[,1]))
    for (extr in c("dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p")) {
      if (extr %in% unique(get_firstcoeffs(coeff = coeff[[pix]],
                                           nb_of_coeff = length(coeff_kept[,1]))[,1])){
        nb_meteo[pix] <- nb_meteo[pix] - 1
        
      }#end fi
    }#end for xtrm
    met_string <- (c("vpd","pr","tmax") %in% unique(coeff_kept[,1])) 
    met_strings[[pix]] <- met_string
    if (identical(met_string, c(T,F,F))) {
      met_type[pix] <- 1
      meteo_type[pix] <- "VPD"
    } else if (identical(met_string, c(F,T,F))) {
      met_type[pix] <- 2
      meteo_type[pix] <- "Pr"
    } else if (identical(met_string, c(F,F,T))) {
      met_type[pix] <- 3
      meteo_type[pix] <- "Tmax"
    } else if (identical(met_string, c(T,T,F))) {
      met_type[pix] <- 4
      meteo_type[pix] <- "VPD & Pr"
    } else if (identical(met_string, c(T,F,T))) {
      met_type[pix] <- 5
      meteo_type[pix] <- "VPD & Tmax"
    } else if (identical(met_string, c(F,T,T))) {
      met_type[pix] <- 6
      meteo_type[pix] <- "Pr & Tmax"
    } else if (identical(met_string, c(T,T,T))) {
      met_type[pix] <- 7
      meteo_type[pix] <- "All"
    } else if (identical(met_string, c(F,F,F)) | is.null(met_string)) {
      met_type[pix] <- 8
      meteo_type[pix] <- "None"
    }
    if ("vpd" %in% unique(coeff_kept[,1])) {met_vpd[pix] <- 1} 
    if ("pr" %in% unique(coeff_kept[,1])) {met_temp[pix] <- 2} 
    if ("tmax" %in% unique(coeff_kept[,1])) {met_pr[pix] <- 3}
  }#end ifelse
}#end for pix


# 4 plots combined: number of var kept, number of extreme indices, combination of meteovar and nb seasons ####

# Plot number of variables kept
levels_nb_var <- cut(nb_coeff_kept, breaks = c(0,5,10,15,20,25,30), right = F)
levels_nb_var <- gsub(","," - ",levels_nb_var,fixed=TRUE)
DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = levels_nb_var)
DF_numbcoeff$levels_nb_var <- gsub("\\[|\\)","",levels_nb_var)

P1 <- ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF_numbcoeff$levels_nb_var)) +
  scale_fill_manual(values=c("0 - 5"="#f1eef6", "5 - 10"="#d4b9da", "10 - 15"="#c994c7",
                             "15 - 20"="#df65b0", "20 - 25"="#dd1c77", "25 - 30"="#980043"),
                    breaks=c("0 - 5", "5 - 10", "10 - 15", "15 - 20", "20 - 25", "25 - 30"),
                    label=c("0 - 4", "5 - 9", "10 - 14", "15 - 10", "20 - 24", "25 - 29"))+
  # scale_fill_gradient(low = "#e7e1ef", high = "#dd1c77") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Nb of var."
       #,title = paste("Number of variables kept, simple",model_name,"regression, "),
       #subtitle = paste("Monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))


# Plot number of extreme indices kept
DF_numbextr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_extr_kept)
DF_numbextr$coeff_kep <- as.factor(DF_numbextr$coeff_kep)

# mycolors <- c("orange", rgb(0.3,0.3,0.5), rgb(0,0,0.7),
# rgb(0.2,0.2,0.8), rgb(0.4,0.4,0.9), rgb(0.7,0.7,1), rgb(0.8,0.8,1), rgb(0.9,0.9,1))
mycolors <- rev(hcl.colors(n=8,palette="Viridis"))
# mycolors <- rev(hcl.colors(n=8,palette="Plasma"))

P2 <- ggplot(data = DF_numbextr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=coeff_kep)) +
  scale_fill_manual(values = mycolors, ) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Nb extr.\n ind."
       #,title = paste("Number of exteme indices kept, simple",model_name,"regression,"),
       #subtitle = paste("Monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))
#   + X11(width = 20, height = 5)


# Combinations of meteorological variables
DF_meteo_type <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], met_type = meteo_type)
#DF_meteo_type$met_type <- as.factor(DF_meteo_type$met_type)
#cols <- c("1" = "#7FC97F", "2" = "cadetblue2", "3" = "#386CB0", "4" = "#824D99", "5" = "#F0027F", "6" = "darkred" , "7" = "#FDC086", "8" = "#FFFF99")

cols <- c("VPD" = "#7FC97F", "Pr" = "cadetblue2", "Tmax" = "#386CB0", "VPD & Pr" = "#824D99",
          "VPD & Tmax" = "#F0027F", "Pr & Tmax" = "darkred" , "All" = "#FDC086", "None" = "#FFFF99")

P3 <- ggplot(data = DF_meteo_type, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=met_type)) +
  scale_fill_manual(values = cols,breaks=c("VPD","Pr","Tmax","VPD & Pr","VPD & Tmax","Pr & Tmax", "All","None")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Combination\nof met.\nvariables"
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))
   # +X11(width = 20, height = 5)




DF_nbseason <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_season = nb_of_seas)
DF_nbseason$nb_season <- as.factor(DF_nbseason$nb_season)


#Plot nb of seasons
P4 <- ggplot(data = DF_nbseason, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF_nbseason$nb_season)) +
  scale_fill_manual(values = c("0"=viridis(6)[1], "1"=viridis(6)[3], "2"=viridis(6)[4],
                               "3"=viridis(6)[5], "4"=viridis(6)[6], "No met. var"="pink")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Nb of seas."
       #,title = paste("Number of different seasons, simple",model_name,"regression"),
       #subtitle = paste("monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))
# + X11(width = 20, height = 5)


L1 <- get_legend(P1+ theme(legend.title = element_text(size=11),
                           legend.text = element_text(size=10),
                           legend.key.size = unit(0.6,"line")))
L2 <- get_legend(P2+ theme(legend.title = element_text(size=11),
                           legend.text = element_text(size=10),
                           legend.key.size = unit(0.6,"line")))
L3 <- get_legend(P3+ theme(legend.title = element_text(size=11),
                           legend.text = element_text(size=10),
                           legend.key.size = unit(0.6,"line")))
L4 <- get_legend(P4+ theme(legend.title = element_text(size=11),
                           legend.text = element_text(size=10),
                           legend.key.size = unit(0.6,"line")))

ggarrange(P1 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8)),
          L1,
          P2 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8)),
          L2,
          P3 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8)),
          L3,
          P4 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8)),
          L4,
          nrow = 4, ncol=2,labels = c("(a)", "", "(b)", "",
                                      "(c)", "", "(d)", ""),
          label.x = -0.015
          ,widths=c(8,1.5), heights=c(1,1,1), font.label = list(size = 14, face = "plain", color ="black")
          )+
  X11(width = 30, height = 25)


# Three most kept variables: dtr, frs and vpd_Jun_Y2 ####
coeff_dtr <- numeric()
coeff_frs <- numeric()
coeff_vpd_Jun_Y2 <- numeric()

TEST_coeff_vpd_Jun_Y2 <- logical()

coeff_test <- sapply(1:length(coeff), function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])

for (pix in 1:length(coeff)) {
  coeff_dtr[pix] <- coeff[[pix]]["dtr",]
  coeff_frs[pix] <- coeff[[pix]]["frs",]
  if("vpd_Jun_Y2" %in% rownames(coeff[[pix]])){
  coeff_vpd_Jun_Y2[pix] <- coeff[[pix]]["vpd_Jun_Y2",]
  } else {
    coeff_vpd_Jun_Y2[pix] <- NA
  }
  
  TEST_coeff_vpd_Jun_Y2[pix] <- "vpd_Jun_Y2" %in% coeff_test[[pix]]
}#end for pix

# I don't obtain at all the same numbers as the barplot of Figure 8, for meteo variables:
sum(coeff_vpd_Jun_Y2!=0, na.rm = T)
sum(TEST_coeff_vpd_Jun_Y2)

cbind((coeff_vpd_Jun_Y2!=0),TEST_coeff_vpd_Jun_Y2)

#ex pix 472 doesn't give the same result
coeff[[472]]
names(numLevels_list[[472]])[coeff[[472]][-1]!=0]#code barplot: vpd_Jun_Y2 not kept but...
coeff[[pix]]["vpd_Jun_Y2",]

#Origine of the pb: shift in coeff... there are more months kept in  numLevels_list?! Problem in the definition of numLevels_list
cbind(names(numLevels_list[[472]]), rownames(coeff[[472]])[-1]) #warnong because not same vector length
Data_xtrm_standardized$tasmax[final_pixels_coord$ref_in_995[472],,1] #Aug_Y2" and "Sept_Y2" should not be included



#plot anyway, to be changed:

#dtr
DF1_dtr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sign_coeff = sign(coeff_dtr))
DF1_dtr$sign_coeff <- as.factor(DF1_dtr$sign_coeff)


#Plot sign of coeff
P1 <- ggplot(data = DF1_dtr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF1_dtr$sign_coeff)) +
  scale_fill_manual(values = c("0"="gray", "1"="#ef8a62", "-1"="#67a9cf"),
                    labels=c("0"="not kept", "1"="positive", "-1"="negative"),
                    breaks=c("-1", "1", "0")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Sign of\ndtr coeff."
       #,title = paste("Number of different seasons, simple",model_name,"regression"),
       #subtitle = paste("monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))


#frs
DF2_frs <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sign_coeff = sign(coeff_frs))
DF2_frs$sign_coeff <- as.factor(DF2_frs$sign_coeff)


#Plot sign of coeff
P2 <- ggplot(data = DF2_frs, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF2_frs$sign_coeff)) +
  scale_fill_manual(values = c("0"="gray", "1"="#ef8a62", "-1"="#67a9cf"),
                    labels=c("0"="not kept", "1"="positive", "-1"="negative"),
                    breaks=c("-1", "1", "0")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Sign of\nfrs coeff."
       #,title = paste("Number of different seasons, simple",model_name,"regression"),
       #subtitle = paste("monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))



#frs
DF3_vpdJun <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sign_coeff = sign(coeff_vpd_Jun_Y2))
DF3_vpdJun$sign_coeff[is.na(DF3_vpdJun$sign_coeff)] <- rep(0, sum(is.na(DF3_vpdJun$sign_coeff)))
DF3_vpdJun$sign_coeff <- as.factor(DF3_vpdJun$sign_coeff)


#Plot sign of coeff
P3 <- ggplot(data = DF3_vpdJun, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF3_vpdJun$sign_coeff)) +
  scale_fill_manual(na.value="#d9d9d9", values = c("0"="#969696", "1"="#ef8a62", "-1"="#67a9cf"),
                    labels=c("0"="not kept", "1"="positive", "-1"="negative")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Sign of\nVPD_Jun_Y2 coeff."
       #,title = paste("Number of different seasons, simple",model_name,"regression"),
       #subtitle = paste("monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))
# +  X11(width = 20, height = 5)


L1 <- get_legend(P1+ theme(legend.title = element_text(size=12),
                           legend.text = element_text(size=12),
                           legend.key.size = unit(0.6,"line")))
L2 <- get_legend(P2+ theme(legend.title = element_text(size=12),
                           legend.text = element_text(size=12),
                           legend.key.size = unit(0.6,"line")))
L3 <- get_legend(P3+ theme(legend.title = element_text(size=12),
                           legend.text = element_text(size=12),
                           legend.key.size = unit(0.6,"line")))

ggarrange(P1 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10)),
          L1,
          P2 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10)),
          L2,
          P3 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10)),
          L3,
          nrow = 3, ncol=2,labels = c("(a)", "", "(b)", "",
                                      "(c)", ""),
          label.x = -0.015
          ,widths=c(8,1.5), heights=c(1,1,1), font.label = list(size = 14, face = "plain", color ="black")
)+
  X11(width = 25, height = 18)



# Histogram nb variables ####
hist(nb_coeff_kept, main = paste0("Nb of variables kept, ",lambda_val), xlab = "number of variables", breaks=14)






# ##### Differences wo extreme indices #####
# 
# plot(csi - csi_wo_extremes, nb_coeff_kept - nb_coeff_kept_wo_extremes, xlab="Diff in CSI", ylab="Diff in nb var")
# abline(v=0, col="red")
# abline(h=0, col="chartreuse4")
# hist(csi - csi_wo_extremes)
# abline(v=mean(csi - csi_wo_extremes, na.rm=T), col="blue")
# abline(v=median(csi - csi_wo_extremes, na.rm=T), col="blue4", lty=2)
# hist(nb_coeff_kept - nb_coeff_kept_wo_extremes)
# abline(v=mean(nb_coeff_kept - nb_coeff_kept_wo_extremes, na.rm=T), col="blue")
# abline(v=median(nb_coeff_kept - nb_coeff_kept_wo_extremes, na.rm=T), col="blue4", lty=2)
# 
# which((csi - csi_wo_extremes)<0)
# plot(csi - csi_wo_extremes, nb_extr_kept, xlab="Diff in CSI", ylab="nb of extremes kept")
# abline(v=0, col="red")
# mean(nb_extr_kept[which((csi - csi_wo_extremes)<0)])
# mean(nb_extr_kept)


# # Difference nb var kept w/o extreme indices ####
# nb_coeff_kept_extremes <- nb_coeff_kept
# # load(file=paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/nb_coeff_wo_xtrms_",
# #                  lambda_name,"_V2020-03-20.RData")) # Pauline
# load(file=paste0("D:/user/vogelj/Group_project/Code/Workspaces/nb_coeff_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData")) # Johannes
# world <- map_data("world")
# 
# nbcoeff_sub <- nb_coeff_kept_extremes - nb_coeff_kept_wo_extremes
# 
# DF_sub <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sub_score = nbcoeff_sub)
# 
# ggplot(data = DF_sub, aes(x=DF_sub$lon, y=DF_sub$lat)) +
#   geom_polygon(data = world, aes(long, lat, group=group),
#                fill="white", color="black", size=0.3) +
#   geom_tile(aes(fill=DF_sub$sub_score)) +
#   scale_fill_gradient2(low = "blue", high = "red") +
#   theme(panel.ontop = F, panel.grid = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
#   ylab("Lat (°N)") +
#   xlab("Lon (°E)") +
#   coord_fixed(xlim = c(-120, 135),
#               ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
#               ratio = 1.3)+
#   labs(fill="Diff in nb var",
#        title = "Difference number of variable kept with xtrm - w/o xtrm",
#        subtitle = paste("Adjusted cut-off", sep = ""))+
#   theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
#         legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
#   X11(width = 20, height = 7)
# 
# 
# mean_yield <- apply(X = Data_xtrm_non_standardized$yield, MARGIN = 1, FUN = mean, na.rm=T)
# plot(mean_yield, csi_sub, ylab="csi with xtrms - w/o xtrms (lambdamin)")






# # Difference CSI w/o extreme indices ####
# csi_extremes <- csi
# csi_sub <- csi_extremes - csi_wo_extremes
# 
# 
# # load(file=paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/csi_wo_xtrms_",
# #                  lambda_name,"_V2020-03-20.RData")) # Pauline
# load(file=paste0("D:/user/vogelj/Group_project/Code/Workspaces/csi_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData")) # Johannes
# 
# world <- map_data("world")
# 
# DF_sub <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sub_score = csi_sub)
# 
# ggplot(data = DF_sub, aes(x=DF_sub$lon, y=DF_sub$lat)) +
#   geom_polygon(data = world, aes(long, lat, group=group),
#                fill="white", color="black", size=0.3) +
#   geom_tile(aes(fill=DF_sub$sub_score)) +
#   scale_fill_gradient2(low = "blue", high = "red") +
#   theme(panel.ontop = F, panel.grid = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
#   ylab("Lat (°N)") +
#   xlab("Lon (°E)") +
#   coord_fixed(xlim = c(-120, 135),
#               ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
#               ratio = 1.3)+
#   labs(fill="Diff in CSI",
#        title = "Difference CSI with xtrm - w/o xtrm",
#        subtitle = paste("Adjusted cut-off", sep = ""))+
#   theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
#         legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
#   X11(width = 20, height = 7)
# hist(csi_sub, main="csi with xtrms - w/o xtrms (lambdamin)")
# 
# mean_yield <- apply(X = Data_xtrm_non_standardized$yield, MARGIN = 1, FUN = mean, na.rm=T)
# plot(mean_yield, csi_sub, ylab="csi with xtrms - w/o xtrms (lambdamin)")

#### Barplot: combinations of months and variables ####
coefs_seas <- sapply(1:length(coeff), function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])
coefs_seas_vec <- unlist(coefs_seas)

pdf(file="C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Images/Final_plots/lambda1se/barplot_variables_lambda1se.pdf",
    width=6,height=8)

par(mar=c(5,7,1,1))
barplot(sort(table(coefs_seas_vec)),horiz=T,las=1,col="grey",
        xlab="Number of grid points, where variable is included in the model",cex.names=0.6)
dev.off()





##### Plot of 1 variable (in the 10 most present variables) #####
top10variables <- names(sort(table(coefs_seas_vec), decreasing = T)[1:10])
allvariables <- colnames(Model_data)[-1]
allvariables_adj <- allvariables
allvariables_adj <- gsub(x=allvariables_adj, pattern="vpd", replacement = "VPD")
allvariables_adj <- gsub(x=allvariables_adj, pattern="tmax", replacement = "Tmax")
allvariables_adj <- gsub(x=allvariables_adj, pattern="pr", replacement = "Pr")
allvariables_adj <- gsub(x=allvariables_adj, pattern="APr", replacement = "Apr") # recorrect April
varia_names_extr <- c("dtr","frs","TXx","TNn","Rx5day","TX90p","TN10p")

# par(mfrow=c(4,2))
# for (varia in 1:10) {
for (varia in 1:length(allvariables)) {
# for (varia in 1:7) {
    
    # varia_name <- top10variables[varia]
    varia_name <- allvariables[varia]
    varia_in_pix <- numeric()
    plots <- vector("list",length=(length(allvariables)))
    for (pix in 1:final_pix_num) {
      varia_in_pix[pix] <- (varia_name %in% row.names(coeff[[pix]])[which(coeff[[pix]]!=0)])
    }#end for pix
    # varia_name <- varia_names_extr[varia] # only for ext. ind. naming
    varia_name <- allvariables_adj[varia]
    
    DF_var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], var_in = varia_in_pix)
    DF_var$var_in <- as.factor(DF_var$var_in)
    
    ggplot(data = DF_var, aes(x=lon, y=lat)) +
      geom_polygon(data = world, aes(long, lat, group=group),
                   fill="white", color="black", size=0.3) +
      geom_tile(aes(fill=DF_var$var_in)) +
      scale_fill_manual(drop=F,values = c("1"="#fc8d62", "0"="#8da0cb"),
                        label= c("1"="Yes", "0"="No"),
                        breaks=c("1","0"), limits=c(0,1)) +
      theme(panel.ontop = F, panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
      ylab("Lat (°N)") +
      xlab("Lon (°E)") +
      coord_fixed(xlim = c(-120, 135),
                  ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                  ratio = 1.3)+
      labs(fill="",title = varia_name)+
      theme(plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15),
            legend.text = element_text(size = 14)) +
      X11(width = 20, height = 6)
      ggsave(filename=paste0("D:/user/vogelj/Group_project/Output/Plots/All_variables/",varia_name,".jpg"))
}







# Indiviual monthly variable predictors
DF_meteo_cat <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], met_cat = met_vpd)
# DF_meteo_cat <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], met_cat = met_temp)
# DF_meteo_cat <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], met_cat = met_pr)
DF_meteo_cat$met_cat <- as.factor(DF_meteo_cat$met_cat)
cols <- c("1" = colors2[1], "2" = colors2[3], "3" = colors2[5],  "4" = colors2[7])

ggplot(data = DF_meteo_cat, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=met_cat)) +
  # scale_fill_manual(values = cols,labels=c("VPD","Pr","T","VPD & Pr","VPD & T","Pr & T", "All","None")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Occurrence of variable"
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 6)


source("./Code/regions_barplot.r") # create barplot figure with variable number of all grid points categorized by continent



# Numbers mentioned in the manuscript text
mean(csi,na.rm=T)
max(csi[loc_eur_pixels_num],na.rm=T)
max(csi[loc_no_am_pixels_num],na.rm=T)
failed_pix <- sapply(1:final_pix_num, function(x) !any(fitted_bad_yield[[x]] == 0) )
sum(failed_pix)
summary(nb_coeff_kept)
summary(nb_extr_kept)
sum(meteo_type=="All")/final_pix_num
sum(meteo_type=="VPD" | meteo_type=="VPD & Tmax" | meteo_type=="VPD & Pr" |meteo_type=="All")/final_pix_num
