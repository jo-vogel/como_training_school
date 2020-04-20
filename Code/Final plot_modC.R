##############################################################################
###########                 Final plots (glmnet)                   ###########
###########   NH crop model, monthly meteovar, XtrM indices        ###########
###########               Author: Pauline Rivoire                  ###########
##############################################################################
# It is structured in the following way:
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



#Remove pixels with low mean yield, lower than the10th percentile on the 995 VERY initial pixels
mean_yield <- apply(Data_xtrm_non_standardized$yield,MARGIN = 1, FUN = mean, na.rm=T)
pix_to_keep <- which(mean_yield > 434.24)
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

# CSI for all the pixel #

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

##### Plot csi all pixel and mean yield/ yield variance #########

excluded_pixel<-which(mean_yield <= 434.24)


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



