##############################################################################
###########          Simple Lasso regression (glmnet)              ###########
###########   NH crop model, monthly meteovar, XtrM indices        ###########
###########               Author: Pauline Rivoire                  ###########
##############################################################################
# It is structured in the following way:
# a) Load the data, build training and testing dataset (monthly meteo var and extreme indices)
# b) run the cross validation on training data for Lasso glmnet, to find lambdamin and lambda 1se
# c) run Lasso glmnet on training data withlambdamin and lambda 1se to get the coefficients of the logistic regression
# d) Model assessment with testing data: plot of CSI, Specificity, nb of coefficients and exteme indices kept
# e) extract Lambda min and lambda 1se




# Clean everything
print("Are you sure you want to run the next line? ;) everything in the environment will be removed")
rm(list=ls(all=TRUE))

# which method? model_name in c("Ridge", "Lasso)
model_name <- "Lasso"
stopifnot(model_name %in% c("Ridge", "Lasso"));if(model_name=="Lasso"){no_model <- 1};if(model_name=="Ridge"){no_model <- 0}

#which lambda?
lambda_VALS <- c("lambda.min", "lambda.1se")
lambda_val <- lambda_VALS[2]

#threshold for bad yields in c(0.025,0.05,0.1)
threshold <- 0.05

#which cutoff level?
segreg_th <- 0.5

##### Initialisation, librairies, data #####

library(ncdf4);library(glmnet);library(InformationValue);library(ROCR)
library(abind);library(stringr);library(tictoc);library(ggplot2)




##### Load standardized Data #####

#Pauline's Laptop
load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/extremeindices_and_monthlymeteovar_rescaled.Rdata")
load("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global/extremeindices_and_monthlymeteovar.Rdata")



##### Process data #####
yield_3dim <- array(Data_standardized$yield,dim=c(965,1,1600))
dtr_3dim <- array(Data_standardized$dtr,dim=c(965,1,1600))
frs_3dim <- array(Data_standardized$frs,dim=c(965,1,1600))
txx_3dim <- array(Data_standardized$txx,dim=c(965,1,1600))
tnn_3dim <- array(Data_standardized$tnn,dim=c(965,1,1600))
rx5_3dim <- array(Data_standardized$rx5,dim=c(965,1,1600))
tx90p_3dim <- array(Data_standardized$tx90p,dim=c(965,1,1600))
tn10p_3dim <- array(Data_standardized$tn10p,dim=c(965,1,1600))

Model_data <- abind(yield_3dim,dtr_3dim,frs_3dim,txx_3dim,tnn_3dim,rx5_3dim,tx90p_3dim,tn10p_3dim
                    ,Data_standardized$tasmax,Data_standardized$vpd,Data_standardized$pr,along=2)
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
Yield <- Data_standardized$yield
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





# #### Run the CrossValidation #####
# tic()
# model_cv_fitting <- list()
# nbyears_final_training_data <- numeric()
# for (pixel in 1:pix_num) {
#   # for (pixel in 1:5) {
#   
#   var_pix <- as.matrix(x1_train_list[[pixel]])
#   yield_pix <- as.matrix(y1_train_list[[pixel]])
#   which_na_xtrm <- which(is.na(var_pix[,1]))
#   nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))
#   
#     if (sum(yield_pix[-which_na_xtrm,])<=nbyears_final_training_data[pixel] & 
#         sum(yield_pix[-which_na_xtrm,])>=nbyears_final_training_data[pixel]-8){
#       #if all the year kept are actually good years, not prediction possible, or there are only 1 or 2 bad years : not possible
#       model_cv_fitting[[pixel]] <- "Training years w/o na in extremes have less than 8 bad years"
#     } else {
#       if(length(which_na_xtrm)>0){
#       model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix[-which_na_xtrm,],
#                                              y = yield_pix[-which_na_xtrm,],
#                                              family = "binomial", alpha = no_model, nfolds = 10)
#       } else {
#         model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix,
#                                                y = yield_pix,
#                                                family = "binomial", alpha = no_model, nfolds = 10)
#         
#       }#end if exists na else
# 
#     }#end if all years kept are good else
# 
#   print(paste(pixel, "out of", pix_num))
# }#end for pixel
# 
# toc() #3.6h for Lasso
# 
# save(model_cv_fitting, file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/cv_month_xtrm_",
#                                    model_name,"_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData"))


# ##### Run the model with lambda1se and lambda min #####
# 
# load(file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/cv_month_xtrm_",
#                    model_name,"_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData"))
# 
# tic()
# 
# if (model_name == "Lasso"){
#   lasso_model_lambdamin <- list()
#   lasso_model_lambda1se <- list()
# 
#   for (pixel in 1:pix_num) {
#     if(is.character(model_cv_fitting[[pixel]])){
#       lasso_model_lambdamin[[pixel]] <- model_cv_fitting[[pixel]]
#       lasso_model_lambda1se[[pixel]] <- model_cv_fitting[[pixel]]
#     } else {
#       var_pix <- as.matrix(x1_train_list[[pixel]])
#       yield_pix <- as.matrix(y1_train_list[[pixel]])
#       which_na_xtrm <- which(is.na(var_pix[,1]))
#       nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))
#       
#       training_years_wo_na <- which(!is.na(x1_train_list[[pixel]]$dtr))
#       
#       if(length(which_na_xtrm)>0){
#         
#         lasso_model_lambdamin[[pixel]] <- glmnet(x = var_pix[-which_na_xtrm,],
#                                                  y = yield_pix[-which_na_xtrm,],
#                                                  family = "binomial", alpha = no_model,
#                                                  lambda = model_cv_fitting[[pixel]]$lambda.min)
#         
#         lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix[-which_na_xtrm,],
#                                                  y = yield_pix[-which_na_xtrm,],
#                                                  family = "binomial", alpha = no_model,
#                                                  lambda = model_cv_fitting[[pixel]]$lambda.1se)
#       } else {
#         
#         lasso_model_lambdamin[[pixel]] <- glmnet(x = var_pix, y = yield_pix,
#                                                  family = "binomial", alpha = no_model,
#                                                  lambda = model_cv_fitting[[pixel]]$lambda.min)
#         
#         lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix, y = yield_pix,
#                                                  family = "binomial", alpha = no_model,
#                                                  lambda = model_cv_fitting[[pixel]]$lambda.1se)
#       }#end if exists na else
# 
#     } #end ifelse
#     print(paste(pixel, "out of", pix_num))
# 
#   }#end for pixel
# 
#   save(lasso_model_lambdamin, file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambdamin_month_xtrm_",
#                                             model_name,"_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData"))
#   save(lasso_model_lambda1se, file = paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambda1se_month_xtrm_",
#                                             model_name,"_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData"))
# }
# toc()
# # 2min for Lasso




##### Load the model #####

# On the Drive you can find my data in:
# Models/LASSO-Ridge regression/regression_results_Global_wo_interactions/Lasso_lambda1se_month_xtrm_Lasso_threshbadyield005.RData
# Models/LASSO-Ridge regression/regression_results_Global_wo_interactions/Lasso_lambdamin_month_xtrm_Lasso_threshbadyield005.RData

load(paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambda1se_month_xtrm_",
            model_name,"_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData"))
load(paste0("C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Lasso_lambdamin_month_xtrm_",
            model_name,"_threshbadyield", str_pad(threshold*100, 3, pad = "0"),".RData"))







##### Extract lambda value #####
lambda1se <- numeric()
lambdamin <- numeric()
for (pix in 1:pix_num) {
  if(is.character(lasso_model_lambda1se[[pix]])) {
    lambda1se[pix] <- NA
    lambdamin[pix] <- NA
  } else {
    lambda1se[pix] <- lasso_model_lambda1se[[pix]]$lambda
    lambdamin[pix] <- lasso_model_lambdamin[[pix]]$lambda
  }#end ifelse
}#end for pix
colnames(coord_subset) <- c("Longitude", "Latitude")
LAMBDAS <- list(lambda_min = lambdamin, lambda_1se = lambda1se, coordinates = coord_subset )
save(LAMBDAS,
     file = "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/lambda_values_Lasso_glmnet_monthly_xtrm_indices.Rdata")





##### Model performance assessment #####


Model_chosen <- lasso_model_lambda1se

pix_model_failed <- numeric(length = pix_num)
coeff  <-list()
speci <- rep(NA, pix_num)
sensi <- rep(NA, pix_num)
csi <- rep(NA, pix_num)

for (pixel in 1:pix_num) {
  if(is.character(Model_chosen[[pixel]])){
    pix_model_failed[pixel] <- 1
    coeff[[pixel]] <- "No model"
  } else {
    
    coeff[[pixel]] <- coefficients(Model_chosen[[pixel]])
    
    mypred <- predict(Model_chosen[[pixel]], as.matrix(x1_test_list[[pixel]]),type="response")
    
    fitted_bad_yield <- ifelse(mypred > segreg_th,1,0)
    
    speci[pixel] <- specificity(actuals = as.matrix(y1_test_list[[pixel]]),
                                predictedScores = fitted_bad_yield,
                                threshold = segreg_th)
    sensi[pixel] <- sensitivity(actuals = as.matrix(y1_test_list[[pixel]]),
                                predictedScores = fitted_bad_yield,
                                threshold = segreg_th)
    
    con_tab <- confusionMatrix(actuals = as.matrix(y1_test_list[[pixel]]),
                                  predictedScores = fitted_bad_yield,
                                  threshold = segreg_th)
    csi[pixel] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
    if(is.na(con_tab["0","0"])){
      csi[pixel] <- 0
    }
      
    }#end if else pb model
}#end pixel



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

for (pixel in 1:pix_num) {
  if(is.character(Model_chosen[[pixel]])){
    nb_extr_kept[pixel] <- NA
    nb_coeff_kept[pixel] <- NA
  } else {
    
    nb_extr_kept[pixel] <- extreme_in_coeff(coeff[[pixel]])
    nb_coeff_kept[pixel] <- number_coeff_kept(coeff[[pixel]])
    
  }#end if else pb model
}#end pixel








##### Plot results ######
# load all coordinates of northern hemisphere
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
lat_all <- ncvar_get(nh_data[[1]],"lat")
lon_all <- ncvar_get(nh_data[[1]],"lon")
lati_all <- rep(lat_all,each=length(lon_all))
long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
coord_all <- cbind(long_all,lati_all)

lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})

coord_subset <- cbind(Data_standardized$longitudes,Data_standardized$latitudes)
world <- map_data("world")

# Plot specificity error ####

DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = speci)

ggplot(data = DF_speci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=speci)) +
  scale_color_gradient2(limits=c(0,1), midpoint = 0.5,
                        low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Specif.",
       title = paste("Specificity, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", cutoff level=", segreg_th,", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


# plot CSI=(hits)/(hits + misses + false alarm) ###
DF_sci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], csi = csi)

ggplot(data = DF_sci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=csi)) +
  scale_color_gradient2(limits=c(0,1), midpoint = 0.5,
                        low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="CSI",
       title = paste("CSI, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", cutoff level=", segreg_th,", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


#Plot number of variables kept
DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_coeff_kept)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=nb_coeff_kept)) +
  scale_color_gradient(low = "pink", high = "darkblue") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of var.",
       title = paste("Number of variables kept, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold, ", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

#Plot number of extreme indices kept
DF_numbextr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_extr_kept)
DF_numbextr$coeff_kep <- as.factor(DF_numbextr$coeff_kep)

mycolors <- c("orange", rgb(0.3,0.3,0.5), rgb(0,0,0.7),
              rgb(0.2,0.2,0.8), rgb(0.4,0.4,0.9), rgb(0.7,0.7,1), rgb(0.8,0.8,1), rgb(0.9,0.9,1))

ggplot(data = DF_numbextr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=coeff_kep)) +
  scale_color_manual(values = mycolors, ) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb ind.",
       title = paste("Number of exteme indices kept, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold, ", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

mean_yield <- apply(X=Data_non_standardized$yield, MARGIN = 1, FUN = mean, na.rm=T)

pairs(cbind(mean_yield, csi, speci, nb_coeff_kept, nb_extr_kept), lower.panel = NULL)


#Plot percentage of extreme indices in var kept
DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], perc_kept = nb_extr_kept/nb_coeff_kept)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=DF_numbcoeff$perc_kept)) +
  scale_color_gradient(low = "yellow", high = "blue") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of var.",
       title = paste("Perc. of extreme indices in variables kept, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold, ", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)






#### Barplot: combinations of months and variables ####
coefs_seas <- sapply(1:length(coeff), function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])
coefs_seas_vec <- unlist(coefs_seas)

png(filename="C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/RidgeRegression/Global_results/Images/Monthly_and_extreme_indices/barplot_variables_lambda1se.png",
    res=2000,units="cm",width=15,height=20)

par(mar=c(5,7,1,1))
barplot(sort(table(coefs_seas_vec)),horiz=T,las=1,col="ForestGreen",
        xlab="Number of pixels, where variable is included in the model\nLasso glmnet lambda 1se",cex.names=0.6)
dev.off()
 










##### Adjust cutoff level #####

# Model_chosen <- lasso_model_lambda1se
# lambda_val <- lambda_VALS[2]

Model_chosen <- lasso_model_lambdamin
lambda_val <- lambda_VALS[1]

source("./Code/Simple_Lasso_Ridge_ElasticNet/cutoff_adj_glmnet_lambda1se.R")
y1_train_list_simple_lasso <- y1_train_list
x1_train_list_simple_lasso <- x1_train_list
cost_fp_simple_lasso <- 100 # Misses: this should be associated with a higher cost, as it is more detrimental
cost_fn_simple_lasso <- 100 # False alarms


work_pix_tmp <- numeric()
for (pix in 1:pix_num) {
  if(is.character(Model_chosen[[pix]])){work_pix_tmp[pix]<-0} else {work_pix_tmp[pix]<-1}
}#end for pix

work_pix <- which(work_pix_tmp==1)
library(pbapply)
cutoff_simple_lasso <- adjust_cutoff(model_vector = Model_chosen,x1_train_list = x1_train_list_simple_lasso, y1_train_list = y1_train_list_simple_lasso,
                                     work_pix = work_pix, cost_fp = cost_fp_simple_lasso, cost_fn= cost_fn_simple_lasso)
# segreg_th_adj_1se <- cutoff_simple_lasso # replace the default threshold = 0.5, by the calculated optimal cutoff
# segreg_th_adj_min <- cutoff_simple_lasso # replace the default threshold = 0.5, by the calculated optimal cutoff

segreg_th_adj_1se <- 0.6599446
segreg_th_adj_min <- 0.5975667

##### Model performance assessment #####

segreg_th <- segreg_th_adj_1se

pix_model_failed <- numeric(length = pix_num)
coeff  <-list()
speci <- rep(NA, pix_num)
sensi <- rep(NA, pix_num)
csi <- rep(NA, pix_num)

for (pixel in 1:pix_num) {
  if(is.character(Model_chosen[[pixel]])){
    pix_model_failed[pixel] <- 1
    coeff[[pixel]] <- "No model"
  } else {
    
    coeff[[pixel]] <- coefficients(Model_chosen[[pixel]])
    
    mypred <- predict(Model_chosen[[pixel]], as.matrix(x1_test_list[[pixel]]),type="response")
    
    fitted_bad_yield <- ifelse(mypred > segreg_th,1,0)
    
    speci[pixel] <- InformationValue::specificity(actuals = as.matrix(y1_test_list[[pixel]]),
                                predictedScores = fitted_bad_yield,
                                threshold = segreg_th)
    sensi[pixel] <- InformationValue::sensitivity(actuals = as.matrix(y1_test_list[[pixel]]),
                                predictedScores = fitted_bad_yield,
                                threshold = segreg_th)
    
    con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(y1_test_list[[pixel]]),
                               predictedScores = fitted_bad_yield,
                               threshold = segreg_th)
    csi[pixel] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
    if(is.na(con_tab["0","0"])){
      csi[pixel] <- 0
    }
    
  }#end if else pb model
}#end pixel



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

for (pixel in 1:pix_num) {
  if(is.character(Model_chosen[[pixel]])){
    nb_extr_kept[pixel] <- NA
    nb_coeff_kept[pixel] <- NA
  } else {
    
    nb_extr_kept[pixel] <- extreme_in_coeff(coeff[[pixel]])
    nb_coeff_kept[pixel] <- number_coeff_kept(coeff[[pixel]])
    
  }#end if else pb model
}#end pixel








##### Plot results ######
# load all coordinates of northern hemisphere
path_to_NH_files <- "C:/Users/admin/Documents/Damocles_training_school_Como/GroupProject1/Data/Global"
nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
lat_all <- ncvar_get(nh_data[[1]],"lat")
lon_all <- ncvar_get(nh_data[[1]],"lon")
lati_all <- rep(lat_all,each=length(lon_all))
long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
coord_all <- cbind(long_all,lati_all)

lapply(1:length(nh_files),function(x){nc_close(nh_data[[x]])})

coord_subset <- cbind(Data_standardized$longitudes,Data_standardized$latitudes)
world <- map_data("world")

# Plot specificity error ####

DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], speci = speci)

ggplot(data = DF_speci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=speci)) +
  scale_color_gradient2(limits=c(0,1), midpoint = 0.5,
                        low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Specif.",
       title = paste("Specificity, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", cutoff level=", round(segreg_th,3),", ",lambda_val, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


# plot CSI=(hits)/(hits + misses + false alarm) ###
DF_sci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], csi = csi)

ggplot(data = DF_sci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=csi)) +
  scale_color_gradient2(limits=c(0,1), midpoint = 0.5,
                        low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="CSI",
       title = paste("CSI, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", cutoff level=", round(segreg_th,3),", ",lambda_val, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


#Plot number of variables kept
DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_coeff_kept)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=nb_coeff_kept)) +
  scale_color_gradient(low = "pink", high = "darkblue") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of var.",
       title = paste("Number of variables kept, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold, ", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

#Plot number of extreme indices kept
DF_numbextr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_extr_kept)
DF_numbextr$coeff_kep <- as.factor(DF_numbextr$coeff_kep)

mycolors <- c("orange", rgb(0.3,0.3,0.5), rgb(0,0,0.7),
              rgb(0.2,0.2,0.8), rgb(0.4,0.4,0.9), rgb(0.7,0.7,1), rgb(0.8,0.8,1), rgb(0.9,0.9,1))

ggplot(data = DF_numbextr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=coeff_kep)) +
  scale_color_manual(values = mycolors, ) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb ind.",
       title = paste("Number of exteme indices kept, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold, ", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)

mean_yield <- apply(X=Data_non_standardized$yield, MARGIN = 1, FUN = mean, na.rm=T)

pairs(cbind(mean_yield, csi, speci, nb_coeff_kept, nb_extr_kept), lower.panel = NULL)


#Plot percentage of extreme indices in var kept
DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], perc_kept = nb_extr_kept/nb_coeff_kept)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=DF_numbcoeff$perc_kept)) +
  scale_color_gradient(low = "yellow", high = "blue") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of var.",
       title = paste("Perc. of extreme indices in variables kept, simple",model_name,"regression, monthly meteo var + extreme indices"),
       subtitle = paste("Bad yield threshold=", threshold, ", lambda 1se", sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)



















# miscellaneous plots to check NA

nb_na_dtr <- numeric()
nb_na_frs <- numeric()
nb_na_txx <- numeric()
nb_na_tnn <- numeric()
nb_na_rx5 <- numeric()
nb_na_tx90p <- numeric()
nb_na_tn10p <- numeric()

for (pixel in 1:pix_num) {
  nb_na_dtr[pixel] <- length(which(is.na(Model_data[pixel,"dtr",])))
  nb_na_frs[pixel] <- length(which(is.na(Model_data[pixel,"frs",])))
  nb_na_txx[pixel] <- length(which(is.na(Model_data[pixel,"txx",])))
  nb_na_tnn[pixel] <- length(which(is.na(Model_data[pixel,"tnn",])))
  nb_na_rx5[pixel] <- length(which(is.na(Model_data[pixel,"rx5",])))
  nb_na_tx90p[pixel] <- length(which(is.na(Model_data[pixel,"tx90p",])))
  nb_na_tn10p[pixel] <- length(which(is.na(Model_data[pixel,"tn10p",])))
}



DF_numbna <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nbna = nb_na_dtr/16)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=DF_numbna$nbna)) +
  scale_color_gradient(low = "yellow", high = "blue") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="% of years",
       title = paste("% years (out of 1600) with NA in extreme indices"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)


nb_GS_toolong <- numeric()
for (pixel in 1:pix_num) {
  nb_GS_toolong[pixel] <- length(which(is.na(Yield[pixel,])))
}

DF_numbna <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nbna = nb_GS_toolong/16)
ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=DF_numbna$nbna)) +
  scale_color_gradient(low = "yellow", high = "blue", limits=c(0,100)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="% of years",
       title = paste("% years (out of 1600) with GS>365 days"))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
