##############################################################################
###########          Comparison of seasonal and monthly results    ###########
###########          Ridge, LASSO, elasticnet regression           ###########
###########      global standardised crop model data               ###########
##############################################################################

# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

# which method? model_name in c("Ridge", "Lasso", "ElasticNet")
model_name <- "ElasticNet"
stopifnot(model_name %in% c("Ridge", "Lasso", "ElasticNet"))


if(model_name=="Lasso"){
  no_model <- 1
}

if(model_name=="Ridge"){
  no_model <- 0
}

if(model_name=="ElasticNet"){
  no_model <- 0.5
}


#############################
#### Load seasonal model ####
#############################

#threshold for bad yields in c(0.025,0.05,0.1)
threshold <- 0.05

#which segregation threshold for the model?
segreg_th <- 0.5

##### Initialisation, librairies, data #####

library(ncdf4);library(rgdal);library(raster);library(RColorBrewer);library(viridis)
library(maps);library(mapdata);library(ggplot2)
library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr)
library(foreach);library(doParallel)
library(tictoc)


# Get the data ####
###################
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
# path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
nh_files <- list.files(path=path_to_NH_files,pattern="NH_yield*") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),
                  FUN = function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
yield <- ncvar_get(nh_data[[1]],"yield")
tasmax <- ncvar_get(nh_data[[1]],"tasmax")
vpd <- ncvar_get(nh_data[[1]],"vpd")
pr <- ncvar_get(nh_data[[1]],"pr")
yield_stand <- ncvar_get(nh_data[[2]],"yield")
tasmax_stand <- ncvar_get(nh_data[[2]],"tasmax")
vpd_stand <- ncvar_get(nh_data[[2]],"vpd")
pr_stand <- ncvar_get(nh_data[[2]],"pr")
lat_subset <- ncvar_get(nh_data[[1]],"lat")
lon_subset <- ncvar_get(nh_data[[1]],"lon")
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

seasons <- matrix(c(2:16),nrow=3,ncol=5) # five seasons autumn1, winter, spring, summer, autumn2; first August and last December excluded
tasmax_seas <- array(NA,dim=c(965,5,1600))
pr_seas <- array(NA,dim=c(965,5,1600))
vpd_seas <- array(NA,dim=c(965,5,1600))
tasmax_seas_stand <- array(NA,dim=c(965,5,1600))
pr_seas_stand <- array(NA,dim=c(965,5,1600))
vpd_seas_stand <- array(NA,dim=c(965,5,1600))
for (i in c(1:5)){
  curr_seas <- seasons[,i] # current season (3 months)
  pr_seas[,i,] <- t(sapply(1:965, function(x) {apply(pr[x,curr_seas,],2,mean,na.rm=T)}))
  tasmax_seas[,i,] <- t(sapply(1:965, function(x) {apply(tasmax[x,curr_seas,],2,mean,na.rm=T)}))
  vpd_seas[,i,] <- t(sapply(1:965, function(x) {apply(vpd[x,curr_seas,],2,mean,na.rm=T)}))
  pr_seas_stand[,i,] <- t(sapply(1:965, function(x) {apply(pr_stand[x,curr_seas,],2,mean,na.rm=T)}))
  tasmax_seas_stand[,i,] <- t(sapply(1:965, function(x) {apply(tasmax_stand[x,curr_seas,],2,mean,na.rm=T)}))
  vpd_seas_stand[,i,] <- t(sapply(1:965, function(x) {apply(vpd_stand[x,curr_seas,],2,mean,na.rm=T)}))
}

Model_data <- abind(yields_3dim,tasmax_seas,vpd_seas,pr_seas,along=2)
Model_data_stand <- abind(yields_3dim,tasmax_seas_stand,vpd_seas_stand,pr_seas_stand,along=2)

threshold <- 0.05
pix_num <- dim(Model_data)[1]
low_yield <- sapply(1:pix_num,function(x) {quantile(yield[x,],threshold,na.rm=T)})
cy <- t(sapply(1:pix_num,function(x){ifelse(yield[x,]<low_yield[x],0,1)})) # identical for standardised and non-standardised yield


cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,1600))
Model_data[,1,] <-cy_reshaped
Model_data_stand[,1,] <-cy_reshaped

columnnames <- c("Yield","tmax_aut_Y1","tmax_win_Y1","tmax_spr_Y1","tmax_sum_Y1","tmax_aut_Y2",
                 "vpd_aut_Y1","vpd_win_Y1","vpd_svpd_Y1","vpd_sum_Y1","vpd_aut_Y2",
                 "pr_aut_Y1","pr_win_Y1","pr_spr_Y1","pr_sum_Y1","pr_aut_Y2")
colnames(Model_data) <- columnnames
colnames(Model_data_stand) <- columnnames



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


# Option 1: Exclude years with NAs (not the whole pixel, just the NA years of it) ####
######################################################################################

# Split data into training and testing data set
vec <- 1:1600

years_with_na <- vector("logical",length=pix_num)
for (i in 1:pix_num){
  years_with_na[i] <- ifelse(length(na_time[[i]] ) ==0,F,T)
}

training_indices <- vector("list",length=pix_num)
testing_indices <- vector("list",length=pix_num)
set.seed(1994)
for (x in 1:pix_num) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x=vec[-na_time[[x]]], size = floor((1600-length(na_time[[x]]))*0.6)))
    testing_indices[[x]] <- vec[-c(na_time[[x]], training_indices[[x]])]
  } else {
    training_indices[[x]] <- sort(sample(1:1600, size = floor(1600*0.6)))
    testing_indices[[x]] <- (1:1600)[-training_indices[[x]]]    
  }
}

Training_Data <- lapply(1:pix_num,function(x){Model_data_stand[x,,training_indices[[x]]]})
Testing_Data <- lapply(1:pix_num,function(x){Model_data_stand[x,,testing_indices[[x]]]})

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




load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/seas_",
                  model_name,"_lambda1se_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))
load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/seas_",
                  model_name,"_lambdamin_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))
load(file = paste("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Seasonal",
                  model_name,"_lambdamin_alpha",str_pad(alpha_elastic*10, 2, pad = "0"),"_threshbadyield",
                  str_pad(threshold*100, 3, pad = "0"),".RData", sep = ""))



#### Model performance ####

test_length <- length(model_fit)

coefs <- lapply(1:test_length, function(x){coef(model_fit[[x]])})



mypred <- lapply(1:test_length, function(x){predict(model_fit[[x]],
                                                    as.matrix(x1_test_list[[x]]),type="response")})

#which segregation threshold for the model?
segreg_th <- 0.5

fitted.results_model <- lapply(1:test_length, function(x){ifelse(mypred[[x]] > segreg_th,1,0)})

mis_clas_err <- sapply(1:test_length, function(x){misClassError(actuals = as.matrix(y1_test_list[[x]]),
                                                                predictedScores=mypred[[x]],
                                                                threshold = segreg_th)})

con_tab <-  lapply(1:test_length, function(x){InformationValue::confusionMatrix(as.matrix(y1_test_list[[x]]),
                                                                                fitted.results_model[[x]],
                                                                                threshold = segreg_th)})

speci <- sapply(1:test_length, function(x){InformationValue::specificity(as.matrix(y1_test_list[[x]]),
                                                                         fitted.results_model[[x]],
                                                                         threshold = segreg_th)})

# CSI #
csi <- numeric()
for(pix in 1:length(coefs)){
  csi[pix] <- con_tab[[pix]]["0","0"]/(con_tab[[pix]]["0","0"] + con_tab[[pix]]["1","0"] + con_tab[[pix]]["0","1"])
  if(is.na(con_tab[[pix]]["0","0"])){
    csi[pix] <- 0
  }
}#end for pix







# Create comparison maps ####
#############################
world <- map_data("world")
substract_score_plot <- function(score_name, score_1, model1_name, score_2, model2_name){#score name in c("speci","csi", EDI")
  sub_score <- score_1 - score_2
  DF_sub <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sub_score = sub_score)
  
  ggplot(data = DF_sub, aes(x=DF_sub$lon, y=DF_sub$lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=DF_sub$sub_score)) +
    scale_color_gradient2(low = "blue", high = "red") +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color=paste("Diff in", score_name),
         title = paste("Difference",score_name,model1_name, "-", model2_name),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  
}


# Specificity map #####
#######################
pairs(cbind(speci_bestglm, speci_lwi, speci_ridge, speci_simplelasso, speci_elastic))



possible_pairs <- c("BestGLM-Ridge",
                    "Lasso_with_interact-Ridge",
                    "LassoSimple-Ridge",
                    "Best GLM-Lasso_with_interact",
                    "Best GLM-LassoSimple",
                    "Lasso_with_interact-LassoSimple",
                    "Elastic-BestGLM",
                    "Elastic-Ridge",
                    "Elastic-LassoSimple",
                    "Elastic-Lasso_with_interact")

score <- "Specificity"
model_1 <- "Elastic"
score_1 <- speci_elastic
model_2 <- "Lasso_with_interact"
score_2 <- speci_lwi

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)

# CSI map ####
##############
pairs(cbind(csi_bestglm, csi_lwi, csi_ridge, csi_simplelasso, csi_elastic))

score <- "CSI"
model_1 <- "Elastic"
score_1 <- csi_elastic
model_2 <- "BestGLM"
score_2 <- csi_bestglm

substract_score_plot(score_name = score,
                     score_1 = score_1, model1_name = model_1,
                     score_2 = score_2, model2_name = model_2)

ggsave(paste(output_path, score,"_",model_1,"VS",model_2,".png", sep = ""), width = 18, height = 6)