##############################################################################
###########                 Final plots (glmnet)                   ###########
###########   NH crop model, monthly meteovar, XtrM indices        ###########
###########               Author: Pauline Rivoire                  ###########
##############################################################################
# It is structured in the following way:
# a) Load the desired model
# b) Load data for initial 995 pixels
# c) plot map CSI
# d) Scatterplot CSIvs mean/variance/Sd Yield with excluded pixels
# e) plot map mean yield with 2 criteria of excluded pixels
# f) plot map length of growing season with 2 criteria of excluded pixels
# g) plot map specificity
# h) plot map nb of variable
# i) plot histo nb of variable
# l) plots map nb of extremes kept
# m) plot map nb of seasons
# n) plot barplot nb pix for each var
# o) plot difference CSI w/o xtrm indices





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

library(ncdf4);library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr);library(tictoc);library(ggplot2);library(viridis);library(maps)




##### Load standardized Data #####

# in the drive folder Data/Global Data
#Pauline's Laptop
#load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/extremeindices_and_monthlymeteovar_rescaled_V2020-03-20.Rdata")
#load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/extremeindices_and_monthlymeteovar_V2020-03-20.Rdata")
# Johannes
# load("D:/user/vogelj/Group_project/Data/extremeindices_and_monthlymeteovar_rescaled_V2020-03-20.Rdata")
# load("D:/user/vogelj/Group_project/Data/extremeindices_and_monthlymeteovar_V2020-03-20.Rdata")
#Cristina
load("C:/Users/39349/Documents/DAMOCLES/Data_global/extremeindices_and_monthlymeteovar_rescaled_V2020-03-20.Rdata")
load("C:/Users/39349/Documents/DAMOCLES/Data_global/extremeindices_and_monthlymeteovar_V2020-03-20.Rdata")


##### Process data #####
yield_3dim <- array(Data_xtrm_standardized$yield,dim=c(969,1,1600))
dtr_3dim <- array(Data_xtrm_standardized$dtr,dim=c(969,1,1600))
frs_3dim <- array(Data_xtrm_standardized$frs,dim=c(969,1,1600))
txx_3dim <- array(Data_xtrm_standardized$txx,dim=c(969,1,1600))
tnn_3dim <- array(Data_xtrm_standardized$tnn,dim=c(969,1,1600))
rx5_3dim <- array(Data_xtrm_standardized$rx5,dim=c(969,1,1600))
tx90p_3dim <- array(Data_xtrm_standardized$tx90p,dim=c(969,1,1600))
tn10p_3dim <- array(Data_xtrm_standardized$tn10p,dim=c(969,1,1600))

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

pix_num <- dim(Model_data)[1]
Yield <- Data_xtrm_standardized$yield
low_yield <- apply(Yield, MARGIN = 1, FUN=quantile, probs=threshold, na.rm=T)
cy <- t(sapply(1:pix_num,function(x){ifelse(Yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield

cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))

Model_data[,1,] <-cy_reshaped


# Exclude NA variable columns
na_col <- matrix(data=NA,nrow=pix_num,ncol=dim(Model_data)[2])
for (j in 1:pix_num){
  for (i in 1:dim(Model_data)[2]){
    na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)

na_time <- vector("list",length=pix_num) # for each pixel, the positions of NAs over time
for (j in 1:pix_num){
  na_time[[j]] <- which(is.na(Model_data[j,1,])) # locations of years with NA values
}


##### Split data into training and testing data set #####
vec <- 1:1600

years_with_na <- vector("logical",length=pix_num)
for (i in 1:pix_num){
  years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
}

training_indices <- vector("list",length=pix_num)
testing_indices <- vector("list",length=pix_num)



seed=1994
train_size <- 70

set.seed(seed)
for (x in 1:pix_num) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((1600-length(na_time[[x]]))*(train_size/100))))
    testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
  } else {
    training_indices[[x]] <- sort(sample(1:1600, size = floor(1600*(train_size/100))))
    testing_indices[[x]] <- (1:1600)[-training_indices[[x]]]    
  }
}


Training_Data <- lapply(1:pix_num,function(x){Model_data[x,,training_indices[[x]]]})
Testing_Data <- lapply(1:pix_num,function(x){Model_data[x,,testing_indices[[x]]]})

pix_in <- 1:pix_num


x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand


var_num <- apply(non_na_col,1,sum)
numLevels_list <- sapply(1:pix_num, function(x){ rep(1,times=var_num[x])})
for (i in 1:pix_num){
  names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
}


##### Load the model #####

# On the Drive you can find my data in:
# Models/LASSO-Ridge regression/regression_output_Global_data/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_969GP.RData
# Models/LASSO-Ridge regression/regression_output_Global_data/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed1994_train70_969GP.RData

load(file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/SensitivityAnalysis/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
                   seed, "_train", train_size,"_969GP.RData"))
# load(file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/SensitivityAnalysis/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
#                    seed, "_train", train_size,"_969GP.RData"))



# Johannes
# load(paste0("D:/user/vogelj/Group_project/Code/Workspaces/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
#              seed, "_train", train_size,"_969GP.RData"))
# load(paste0("D:/user/vogelj/Group_project/Code/Workspaces/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
#              seed, "_train", train_size,"_969GP.RData"))


#Cristina

load(paste0("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
             seed, "_train", train_size,"_969GP.RData"))
load(paste0("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
             seed, "_train", train_size,"_969GP.RData"))


##### Adjust cutoff level #####

Model_chosen <- lasso_model_lambda1se
lambda_val <- lambda_VALS[2]
lambda_name <- lambda_NAMES[2]


#with 969 pixels, seed 1994, train size 70%
segreg_th_adj_1se <- 0.6500294

segreg_th <- segreg_th_adj_1se


######### Load Data for the initials 995 Pixels #########
# You can find it in Drive: Figures/Data-Processing #####
load("C:/Users/39349/Documents/DAMOCLES/Data_global/Data_995_pixR.RData")

Mean_yields_995<- apply(yields_995,MARGIN = 1, FUN = mean, na.rm=T)  #Median for the row

Mean_grow_season_995<- apply(length_growseas_995,1,MARGIN = 1, FUN = mean, na.rm=T) #Median for the row


#Remove pixels with low mean yield, lower than the10th percentile on the 995 VERY initial pixels
mean_yield <- apply(Data_xtrm_non_standardized$yield,MARGIN = 1, FUN = mean, na.rm=T)

pix_to_keep <- which(mean_yield > quantile_10_raw_data)
excluded_pixel<-which(mean_yield <= quantile_10_raw_data)

final_pix_num <- length(pix_to_keep)


##### Model performance assessment #####


coeff  <-list()
speci <- rep(NA, final_pix_num)
sensi <- rep(NA, final_pix_num)
csi <- rep(NA, final_pix_num)

for (pixel in 1:final_pix_num) {
  pix <- pix_to_keep[pixel]
  coeff[[pixel]] <- coefficients(Model_chosen[[pix]])
  
  mypred <- predict(Model_chosen[[pix]], as.matrix(x1_test_list[[pix]]),type="response")
  
  fitted_bad_yield <- ifelse(mypred > segreg_th,1,0)
  
  speci[pixel] <- InformationValue::specificity(actuals = as.matrix(y1_test_list[[pix]]),
                                                predictedScores = fitted_bad_yield,
                                                threshold = segreg_th)
  sensi[pixel] <- InformationValue::sensitivity(actuals = as.matrix(y1_test_list[[pix]]),
                                                predictedScores = fitted_bad_yield,
                                                threshold = segreg_th)
  
  con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(y1_test_list[[pix]]),
                                               predictedScores = fitted_bad_yield,
                                               threshold = segreg_th)
  csi[pixel] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
  if(is.na(con_tab["0","0"])){
    csi[pixel] <- 0
  }
  
}#end pixel

# CSI for all the pixels #

csi_all <- rep(NA, final_pix_num)

for (pixel in 1:pix_num) {
  pix <- pixel
  coeff[[pixel]] <- coefficients(Model_chosen[[pix]])
  
  mypred <- predict(Model_chosen[[pix]], as.matrix(x1_test_list[[pix]]),type="response")
  
  fitted_bad_yield <- ifelse(mypred > segreg_th,1,0)
  
  con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(y1_test_list[[pix]]),
                                               predictedScores = fitted_bad_yield,
                                               threshold = segreg_th)
  
  csi_all[pixel] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
  if(is.na(con_tab["0","0"])){
    csi_all[pixel] <- 0
  }  }
  
#######

# load(file=paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/coeff_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData")) # Pauline
# load(file=paste0("D:/user/vogelj/Group_project/Code/Workspaces/coeff_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData")) # Johannes

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

# load(file=paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/nb_coeff_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData")) # Pauline
# load(file=paste0("D:/user/vogelj/Group_project/Code/Workspaces/nb_coeff_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData"))# Johannes

# load(file=paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/csi_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData")) # Pauline
# load(file=paste0("D:/user/vogelj/Group_project/Code/Workspaces/csi_wo_xtrms_",
#                  lambda_name,"_V2020-03-20.RData")) # Johannes

mean_yield <- apply(Data_xtrm_non_standardized$yield,MARGIN = 1, FUN = mean, na.rm=T)
var_yield <- apply(Data_xtrm_non_standardized$yield,MARGIN = 1, FUN = var, na.rm=T)


##### Plot results ######
# load all coordinates of northern hemisphere
#path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global" # Pauline
# path_to_NH_files <- "D:/user/vogelj/Data/Group project Como" # Johannes
path_to_NH_files <- "C:/Users/39349/Documents/DAMOCLES/Data_global" # Cristina

nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
lat_all <- ncvar_get(nh_data[[1]],"lat")
lon_all <- ncvar_get(nh_data[[1]],"lon")
lati_all <- rep(lat_all,each=length(lon_all))
long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
coord_all <- cbind(long_all,lati_all)

lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})

coord_subset <- cbind(Data_xtrm_standardized$longitudes[pix_to_keep],
                      Data_xtrm_standardized$latitudes[pix_to_keep])
world <- map_data("world")


# plot CSI=(hits)/(hits + misses + false alarm) ####
DF_sci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], csi = csi)

ggplot(data = DF_sci, aes(x=lon, y=lat)) +
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

plot(mean_yield[pix_to_keep],csi)
cor.test(mean_yield[pix_to_keep], csi)
cor.test(mean_yield[pix_to_keep], csi, method = "kendall")
plot(var_yield[pix_to_keep],csi)
cor.test(var_yield[pix_to_keep], csi)
plot(mean_yield[pix_to_keep],var_yield[pix_to_keep])
cor.test(mean_yield[pix_to_keep], var_yield[pix_to_keep])

##### ScatterPlot CSI all pixel and mean yield/ yield variance #########
          ## Plot with 969 pixels: Rejected the pixel with Yield<..### 


#Criterion for excluded pixels that are in the plot
#Mean Yield< quantile 0.1 calculated on the raw data Yield (before excluding year with GrowS>365)#
  
  
#Mean Yield and CSI
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(mean_yield, csi_all,col="black",pch=19, xlab="Mean Yield", ylab="CSI", cex.lab=1.2, cex.axis=1.2,  cex.sub=1.2)
points(mean_yield[excluded_pixel], csi_all[excluded_pixel],col="red",pch=19)
legend("topleft",c("Excluded pixels","Included pixels"),col=c("red","black"), pch=19)

# Yield Variance and CSI
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(var_yield,csi_all,col="black",pch=19, xlab="Yield variance", ylab="CSI", cex.lab=1.2, cex.axis=1.2,  cex.sub=1.2)
points( var_yield[excluded_pixel],csi_all[excluded_pixel],col="red",pch=19)
legend("topleft",c("Excluded pixels","Included pixels"),col=c("red","black"), pch=19)

# Yield Standard Dev and CSI
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(sqrt(var_yield),csi_all,col="black",pch=19, xlab="Yield Standard Deviation", ylab="CSI", cex.lab=1.2, cex.axis=1.2,  cex.sub=1.2)
points( sqrt(var_yield[excluded_pixel]),csi_all[excluded_pixel],col="red",pch=19)
legend("topleft",c("Excluded pixels","Included pixels"),col=c("red","black"), pch=19)


#Just excluded pixels

plot(var_yield[excluded_pixel], csi_all[excluded_pixel],col="red",pch=19, main="Only excluded pixel")



# MAP Plot Mean Yield ####
# Plot with the initial 995 pixels #


excluded_pixel_995<-which(Mean_yields_995<= quantile_10_raw_data)
pix_to_keep_995<- which(Mean_yields_995> quantile_10_raw_data)

#excluded_pixel_995<-which(Mean_yields_995<= quantile(Mean_yields_995,0.10))
#pix_to_keep_995<- which(Mean_yields_995>quantile(Mean_yields_995,0.10))



pixel<-rep(1:995,1)

#Plot data Processing: scatterplot pixels and mean yields

plot(pixel,Mean_yields_995, xlab="Pixels", ylab="Mean Yield")
points(pixel[excluded_pixel_995],Mean_yields_995[excluded_pixel_995], col="red",pch=16)
points(pixel[pix_to_rm],Mean_yields_995[pix_to_rm],col="green",pch=16)
legend("topleft",c("5th percentile yield==0","Mean Yield< 10th percentile"),col=c("green","red"), pch=19)

# MAP Plot Mean Yield 

DF_meanY <- data.frame(lon=coord_all_995[,1], lat = coord_all_995[,2], meany = Mean_yields_995)

DF_meanY_kept<-DF_meanY [pix_to_keep_995,]
DF_meanY_EX2<-DF_meanY [excluded_pixel_995,]
DF_meanY_EX1<- DF_meanY [pix_to_rm,]

p1<-ggplot(data = DF_meanY, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3)+  geom_tile(aes(fill=meany)) +
  scale_fill_gradient2(midpoint = max(Mean_yields_995, na.rm = T)/2,
                       limits=c(0,max(Mean_yields_995)),
                       low = "black", mid = "red3", high = "yellow") +
  #scale_shape_manual(values=c("Excluded pixels"=4, "Exc2=15"))+
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-115, 130),
              ylim = c(min(coord_all_995[,2]), max(coord_all_995[,2])),
              ratio = 1.3)+
  labs(fill="Mean Yield"  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))+
  geom_point(data = DF_meanY_EX2, aes(x = lon, y = lat), color = "yellow", size = 0.8, pch=4)+
  geom_point(data = DF_meanY_EX1, aes(x = lon, y = lat), color = "yellow", size = 0.4, pch=15)


#ggsave("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/MeanYield_2crit_1.png", units="in", dpi=400)
ggsave("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/MeanYield_2crit_1.png", units="in", width=16.7,height=7.52, dpi=400)

 #+X11(width = 20, height = 6)



### MAP Plot Growing Season ####
# Plot with the initial 995 pixels #


DF_meanGrSeas<-data.frame(lon=coord_all_995[,1], lat = coord_all_995[,2], meany = Mean_grow_season_995)

DF_meanGrSeas_kept<-DF_meanY [pix_to_keep_995,]
DF_meanGrSeas_EX2<-DF_meanY [excluded_pixel_995,]
DF_meanGrSeas_EX1<- DF_meanY [pix_to_rm,]


p2<-ggplot(data = DF_meanGrSeas, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3)+  geom_tile(aes(fill=meany)) +
  scale_fill_gradient2(midpoint = max(Mean_grow_season_995, na.rm = T)/2,
                       limits=c(0,max(Mean_grow_season_995)),
                       low = "black", mid = "red3", high = "blue") +
  #scale_shape_manual(values=c("Excluded pixels"=4))+
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-115, 130),
              ylim = c(min(coord_all_995[,2]), max(coord_all_995[,2])),
              ratio = 1.3)+
  labs(fill="Mean growing season"  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))+
  geom_point(data = DF_meanGrSeas_EX2, aes(x = lon, y = lat), color = "yellow", size = 0.8, pch=4, show.legend=TRUE)+
  geom_point(data = DF_meanGrSeas_EX1, aes(x = lon, y = lat), color = "yellow", size = 0.4, pch=15, show.legend=TRUE)

#ggsave("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/Grow_season_2crit_col2.png", units="in", dpi=400)
ggsave("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/Grow_season_2crit_col2.png", units="in", width=16.7,height=7.52, dpi=400)

ggarrange(p1, p2, nrow = 2,ncol=1,labels = c("A)", "B)"))
ggsave("C:/Users/39349/Documents/DAMOCLES/Final Workspace LASSO/MATCH.png", units="in", dpi=400)


# # Plot specificity error ####
# 
# DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = speci)
# 
# ggplot(data = DF_speci, aes(x=lon, y=lat)) +
#   geom_polygon(data = world, aes(long, lat, group=group),
#                fill="white", color="black", size=0.3) +
#   geom_tile(aes(fill=speci)) +
#   scale_fill_gradient2(limits=c(0,1), midpoint = 0.5,
#                        low = "black", mid = "red3", high = "yellow") +
#   theme(panel.ontop = F, panel.grid = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
#   ylab("Lat (°N)") +
#   xlab("Lon (°E)") +
#   coord_fixed(xlim = c(-120, 135),
#               ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
#               ratio = 1.3)+
#   labs(fill="Specif."
#        #,title = paste("Specificity")
#        #,subtitle = paste("Monthly meteo var + extreme indices, cutoff level=", round(segreg_th,3),", ",lambda_val, sep = "")
#   )+
#   theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
#         legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
#   X11(width = 20, height = 6)
# 
# plot(mean_yield,speci)
# cor.test(mean_yield, speci)
# cor.test(mean_yield, speci, method = "kendall")



# Plot number of variables kept ####
DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_coeff_kept)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=nb_coeff_kept)) +
  scale_fill_gradient(low = "pink", high = "darkblue") +
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
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 5)

# Histogram nb variables ####
hist(nb_coeff_kept, main = paste0("Nb of variables kept, ",lambda_val), xlab = "number of variables", breaks=14)


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




# Plot number of extreme indices kept ####
DF_numbextr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_extr_kept)
DF_numbextr$coeff_kep <- as.factor(DF_numbextr$coeff_kep)

# mycolors <- c("orange", rgb(0.3,0.3,0.5), rgb(0,0,0.7),
# rgb(0.2,0.2,0.8), rgb(0.4,0.4,0.9), rgb(0.7,0.7,1), rgb(0.8,0.8,1), rgb(0.9,0.9,1))
mycolors <- rev(hcl.colors(n=8,palette="Viridis"))
# mycolors <- rev(hcl.colors(n=8,palette="Plasma"))

ggplot(data = DF_numbextr, aes(x=lon, y=lat)) +
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
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 5)
# ggsave(filename="D:/user/vogelj/Group_project/Output/Plots/num_extr_ind_viridis.pdf")

mean_yield <- apply(X=Data_xtrm_non_standardized$yield, MARGIN = 1, FUN = mean, na.rm=T)

pairs(cbind(mean_yield[pix_to_keep], csi, speci, nb_coeff_kept, nb_extr_kept), lower.panel = NULL)




# Map of the number of season kept #####


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
      return("No meteo var")
    }
    
  } else {
    return("No meteo var")
  }#end if else
  
}

nb_of_seas <- numeric()
for (pix in 1:final_pix_num) {
  nb_of_seas[pix] <- count_seas_and_var(coeff[[pix]])
}

DF_nbseason <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_season = nb_of_seas)
DF_nbseason$nb_season <- as.factor(DF_nbseason$nb_season)

ggplot(data = DF_nbseason, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF_nbseason$nb_season)) +
  scale_fill_manual(values = c("0"=viridis(6)[1], "1"=viridis(6)[3], "2"=viridis(6)[4],
                               "3"=viridis(6)[5], "4"=viridis(6)[6], "No meteo var"="pink")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Nb of seasons"
       #,title = paste("Number of different seasons, simple",model_name,"regression"),
       #subtitle = paste("monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 5)


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





# Plot nb of different meteo variables ####

nb_meteo <- numeric()
source("./Code/get_first_coeff_function.R")

for (pix in 1:length(coeff)) {
  if((length(which(coeff[[pix]]!=0))-1)<1){
    nb_meteo[pix] <- 0
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
    
  }#end ifelse
}#end for pix

DF_nbmeteo <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_meteo = nb_meteo)
DF_nbmeteo$nb_meteo <- as.factor(DF_nbmeteo$nb_meteo)

ggplot(data = DF_nbmeteo, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=nb_meteo)) +
  scale_fill_manual(values = c("0"="#fbb4b9", "1"="#f768a1", "2"="#c51b8a", "3"="#7a0177")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(fill="Nb of\nmeteo var."
       #,title = paste("Number of variables kept, simple",model_name,"regression, "),
       #subtitle = paste("Monthly meteo var + extreme indices, ",lambda_val, sep = "")
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 5)





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
varia_names_extr <- c("dtr","frs","TXx","TNn","Rx5day","TX90p","TN10p")

# par(mfrow=c(4,2))
# for (varia in 1:10) {
for (varia in 1:length(allvariables)) {
  # for (varia in 1:7) {
  
  # varia_name <- top10variables[varia]
  varia_name <- allvariables[varia]
  varia_in_pix <- numeric()
  plots <- vector("list",length=(length(allvariables)))
  for (pix in 1:pix_num) {
    varia_in_pix[pix] <- (varia_name %in% row.names(coeff[[pix]])[which(coeff[[pix]]!=0)])
  }#end for pix
  varia_name <- varia_names_extr[varia] # only for ext. ind. naming
  
  DF_var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], var_in = varia_in_pix)
  DF_var$var_in <- as.factor(DF_var$var_in)
  
  # plots[[varia]] <- 
  ggplot(data = DF_var, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_tile(aes(fill=DF_var$var_in)) +
    scale_fill_manual(values = c("1"="#fc8d62", "0"="#8da0cb"),
                      label= c("1"="Yes", "0"="No"),
                      breaks=c("1","0")) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(fill=paste0(varia_name,"\nselected")
         ,title = varia_name
         #,title = paste("Number of variables kept, simple",model_name,"regression, "),
         #subtitle = paste("Monthly meteo var + extreme indices, ",lambda_val, sep = "")
    )+
    theme(plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 6)
  ggsave(filename=paste0("D:/user/vogelj/Group_project/Output/Plots/",varia_name,".jpg"))
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
    
  }#end ifelse
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
    meteo_type[pix] <- "T"
  } else if (identical(met_string, c(T,T,F))) {
    met_type[pix] <- 4
    meteo_type[pix] <- "VPD & Pr"
  } else if (identical(met_string, c(T,F,T))) {
    met_type[pix] <- 5
    meteo_type[pix] <- "VPD & T"
  } else if (identical(met_string, c(F,T,T))) {
    met_type[pix] <- 6
    meteo_type[pix] <- "Pr & T"
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
  
}#end for pix


# Combinations of meteorological variables
DF_meteo_type <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], met_type = meteo_type)
#DF_meteo_type$met_type <- as.factor(DF_meteo_type$met_type)
#cols <- c("1" = "#7FC97F", "2" = "cadetblue2", "3" = "#386CB0", "4" = "#824D99", "5" = "#F0027F", "6" = "darkred" , "7" = "#FDC086", "8" = "#FFFF99")

cols <- c("VPD" = "#7FC97F", "Pr" = "cadetblue2", "T" = "#386CB0", "VPD & Pr" = "#824D99",
          "VPD & T" = "#F0027F", "Pr & T" = "darkred" , "All" = "#FDC086", "None" = "#FFFF99")

ggplot(data = DF_meteo_type, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=met_type)) +
  scale_fill_manual(values = cols,breaks=c("VPD","Pr","T","VPD & Pr","VPD & T","Pr & T", "All","None")) +
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
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 5)
#ggsave("D:/user/vogelj/Group_project/Output/Plots/meteotype.pdf")


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
ggsave("D:/user/vogelj/Group_project/Output/Plots/met_vpd.pdf")
# ggsave("D:/user/vogelj/Group_project/Output/Plots/met_temp.pdf")
# ggsave("D:/user/vogelj/Group_project/Output/Plots/met_pr.pdf")


