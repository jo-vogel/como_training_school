##########################################
# Visualisation of results 
# by Pauline, partially edited by Johannes
##########################################

library(ggplot2)
library(scales)
library(viridis)
library(maps)
library(glmnet)

# Load model output ####
########################

# source('./Code/Lasso_interact_global.R') # takes ca. 3 hours
source('./Code/Lasso_interact_global_preparation.R')
load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_complete.RData")



# Model performance assessment ####
###################################

# Identify pixels with failed runs
failed_pixels <- which(sapply(1:965, function(x) {is.character(cv_fit[[x]])})==1)
work_pix <- pix_in[-failed_pixels] # working pixels


i_1Std <- sapply(work_pix, function(x){ which(cv_fit[[x]]$lambdaHat1Std == cv_fit[[x]]$lambda)}) # the preferential lambda (tuning parameter)
i_1Std_all_pix <- rep(NA,965)
i_1Std_all_pix[work_pix] <- i_1Std # needed as a workaround (to have an object of lenght=965)

coefs <- vector("list",length=965)
coefs[work_pix] <- lapply(work_pix, function(x){coef(cv_fit[[x]]$glinternetFit)[[i_1Std_all_pix[[x]]]]})

# coefs$mainEffects # model part without interactions
# names(numLevels)[coefs$mainEffects$cont] # Main effect variables (without interactions)
# 
# coefs$interactions # model part with interactions pairs
# names(numLevels)[coefs$interactions$contcont] # Main effect variables (with interactions)


#which segregation threshold for the model?
segreg_th <- 0.5
mypred <- lapply(work_pix, function(x){predict(cv_fit[[x]],x1_test_list[[x]],type="response")}) 
fitted.results_model <- lapply(seq_along(work_pix), function(x){ifelse(mypred[[x]] > segreg_th,1,0)})

y1_test_list_red <- lapply(work_pix,function(work_pix){y1_test_list[[work_pix]]})
mis_clas_err <- rep(NA,965)
# mis_clas_err[work_pix] <- sapply(seq_along(work_pix), function(x){misClassError(y1_test_list_red[[x]],mypred[[x]])})
mis_clas_err[work_pix] <- sapply(seq_along(work_pix), function(x){misClassError(actuals = y1_test_list_red[[x]],
                                                                                predictedScores=mypred[[x]],
                                                                                threshold = segreg_th)})



con_tab <-  lapply(seq_along(work_pix), function(x){InformationValue::confusionMatrix(y1_test_list_red[[x]],fitted.results_model[[x]])})
sensi <- rep(NA,965)
speci <- rep(NA,965)
sensi[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::sensitivity(y1_test_list_red[[x]],fitted.results_model[[x]],
                                                                                         threshold = segreg_th)})
speci[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::specificity(y1_test_list_red[[x]],fitted.results_model[[x]],
                                                                                         threshold = segreg_th)})

obs_pred <- lapply(seq_along(work_pix), function(x){cbind(y1_test_list_red[[x]],fitted.results_model[[x]])})
tp <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==2)})
tn <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==0)})
fp <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==0 & obs_pred[[x]][,2]==1)})
fn <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==1 & obs_pred[[x]][,2]==0)})
con_tab1 <- sapply(seq_along(work_pix), function(x){matrix(c(tp[x],fn[x],fp[x],tn[x]),nrow=2,ncol=2)})
# spec <- tn/(tn+fp) 
# sens <- tp/(tp+fn) 
csi <- rep(NA,965)
csi[work_pix] <- sapply(seq_along(work_pix), function(x){tn[x]/(tn[x]+fp[x]+fn[x])})



# Visualisation ####
####################

model_name <- "Lasso with interactions"
world <- map_data("world")


# Plot miscla error ####

DF_miscla <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], miscla = mis_clas_err)

ggplot(data = DF_miscla, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.1) +
  geom_point(shape=15, aes(color=miscla),size=0.7) +
  scale_color_gradient2(limits=c(min(mis_clas_err,na.rm=T),max(mis_clas_err,na.rm=T)),midpoint=min(mis_clas_err,na.rm=T)+(max(mis_clas_err,na.rm=T)-min(mis_clas_err,na.rm=T))/2,
                        low = "yellow", mid = "red3", high = "black") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Misclass.\nerror",
       title = paste("Misclassification error, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group project/Output/Plots/Mis_class_error_lasso_interact_map.png")


# Plot specificity ####

DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], specificity = speci)

ggplot(data = DF_speci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=speci),size=0.7) +
  scale_color_gradient2(limits=c(min(speci,na.rm=T),max(speci,na.rm=T)),midpoint=min(speci,na.rm=T)+(max(speci,na.rm=T)-min(speci,na.rm=T))/2,
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
       title = paste("Specificity, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group project/Output/Plots/Specificity_lasso_interact_map.png")


# Plot sensitivity ####

DF_sensi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sensitivity = sensi)

ggplot(data = DF_sensi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=sensi),size=0.7) +
  scale_color_gradient2(limits=c(min(sensi,na.rm=T),max(sensi,na.rm=T)),midpoint=min(sensi,na.rm=T)+(max(sensi,na.rm=T)-min(sensi,na.rm=T))/2,
                        low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Sensi.",
       title = paste("Sensitivity, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))  +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group project/Output/Plots/Sensitivity_lasso_interact_map.png")


# Plot critical succes index (CSI) ####

DF_csi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], Critical_success_index = csi)

ggplot(data = DF_csi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=csi),size=0.7) +
  scale_color_gradient2(limits=c(min(csi,na.rm=T),max(csi,na.rm=T)),midpoint=min(csi,na.rm=T)+(max(csi,na.rm=T)-min(csi,na.rm=T))/2,
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
       title = paste("Critical succes index, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))  +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group project/Output/Plots/CSI_lasso_interact_map.png")


# Correlation between yield and specificity or miscla error ####

mean_yield <- apply(X=yield, MARGIN = 1, FUN = mean, na.rm=T)
plot(mean_yield, mis_clas_err,
     xlab="Mean Yield (kg/yr)", ylab="Miss-classification error",
     main=paste("Scatterplot mean yield, miss-classification error\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))
plot(mean_yield, speci,
     xlab="Mean Yield (kg/yr)", ylab="Specificity",
     main=paste("Scatterplot mean yield, Specificity\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))


pairs(cbind(mean_yield, mis_clas_err, speci),
      main=paste(model_name, " regression, bad yield thr.=", threshold,
                 "\nsegreg. thr.=", segreg_th, sep = ""))

cor(mean_yield, mis_clas_err)
cor(mean_yield, speci)


# Map of the number of coefficients kept (Lasso) #####
# if(model_name=="Lasso"){
  coeff_kep <- numeric()
  
  for (pix in 1:pix_num) {
    coeff_kep[pix] <- length(coefs[[pix]]$mainEffects$cont)
    # coeff_kep[pix] <- sum(coefs[[pix]][row.names(coefs[[pix]])!="(Intercept)"]!=0)
  }#end for pix
  

  DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = coeff_kep)
  
  ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=coeff_kep),size=0.7) +
    # scale_color_gradientn(limits=c(0,max(DF_numbcoeff[,3])), 
    # colours=c(gray.colors(1),topo.colors(23)[-c(1,3,5,16:23)],rev(heat.colors(10))) ,values=rescale(0:22,c(0,1))) + # mixed visualisation
    # colours=c(gray.colors(1),topo.colors(20)[-c(13:20)],rev(heat.colors(10))) ,values=rescale(0:22,c(0,1))) + # mixed visualisation v1
    scale_color_gradientn(limits=c(0,15),
    colours=c(gray.colors(1),topo.colors(9)[-c(8,9)],rev(heat.colors(7))) ,values=rescale(0:15,c(0,1)),
    breaks=c(0,3,6,9,12,15),labels=c("0","3","6","9","12",">=15")) + # mixed visualisation with cutoff
    # scale_color_gradientn(limits=c(0,15), # cut off high values for better visualisation of the rest
                          # colours=c(gray.colors(1),topo.colors(14)) ,values=rescale(0:15,c(0,1))) + # topo color scheme
    #                      colours=c(gray.colors(1),rainbow(14)) ,values=rescale(0:15,c(0,1))) + # rainbow color scheme
    # scale_color_gradient2(limits=c(0,15),midpoint=15/2, # cut off high values for better visualisation of the rest
    #                      low = "blue", mid = "yellow", high = "red3") +
    # scale_color_gradient2(limits=c(0,max(DF_numbcoeff[,3])),midpoint=max(DF_numbcoeff[,3])/2,
                          # low = "blue", mid = "yellow", high = "red3") +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="Nb of coefficients",
         title = paste("Number of coefficients kept, simple",model_name,"regression"),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
    # d + scale_fill_discrete(breaks=c(2,4,6,8,10), labels = c("A", "B", "C", "D", "E"))

  ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_coefficients_lasso_interact_map.png")
# }
  plot(table(coeff_kep))  # overview of distribution of the number of coefficients
  

# Map of number of interactions ####
  
  num_interact <- numeric()
  
  for (pix in 1:pix_num) {
    if (is.null(dim(coefs[[pix]]$interactions$contcont)[1])){
      num_interact[pix] <- 0
    } else {
      
    num_interact[pix] <-  dim(coefs[[pix]]$interactions$contcont)[1]
    }
  }
  

  DF_numb_interact <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], num_interact = num_interact)
  
  ggplot(data = DF_numb_interact, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=num_interact),size=0.7) +
    # scale_color_gradient2(limits=c(0,max(DF_numb_interact[,3])),midpoint=max(DF_numb_interact[,3]/2),
                          # low = "blue", mid = "yellow", high = "red3") +
    scale_color_gradientn(limits=c(0,16), breaks=c(0,4,8,12,16),labels=c("0","4","8","12",">=16"),
                          colours=c(gray.colors(1),topo.colors(10)[-c(9,10)],rev(heat.colors(7))) ,values=rescale(0:16,c(0,1))) + # mixed visualisation with cutoff
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="Nb of int.",
         title = paste("Number of coefficient interactions from",model_name,"regression"),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_coefficient_interactions_lasso_interact_map.png")
  plot(table(num_interact))  # overview of distribution of the number of interactions of coefficients

  
  ####################
  # Number of variables and seasons per pixel
  ####################
  
  # number of variables: which are included
  nb_of_var <- vector("numeric",length=length(coefs))
  for (j in  1:length(coefs)){
    coefs_str <- names(numLevels_list[[j]])[coefs[[j]]$mainEffects$cont]
    str_var <- vector(mode="character",length=length(coefs_str))
    for (i in seq_along(coefs_str)){
      str_var[i] <- strsplit(coefs_str[i],split="_")[[1]][1]
      nb_of_var[j] <- length(unique(str_var))
    }
  }
  plot(table(nb_of_var))
  
  
  # number of months: list of variables, select months
  num_months <- vector("numeric",length=length(coefs))
  months_all_pix <- vector("list",length=length(coefs))
  for (j in  1:length(coefs)){
    coefs_str <- names(numLevels_list[[j]])[coefs[[j]]$mainEffects$cont]
    str_months <- vector(mode="character",length=length(coefs_str))
    str_months_short <- vector("character",length(coefs_str))
    for (i in seq_along(coefs_str)){
      str_months[i] <- paste(strsplit(coefs_str[i],split="_")[[1]][2],strsplit(coefs_str[i],split="_")[[1]][3])
      num_months[j] <- length(unique(str_months))
      str_months_short[i] <- strsplit(coefs_str[i],split="_")[[1]][2]
    }
    months_all_pix[[j]] <- str_months_short
  }
  plot(table(num_months))
  
  
  nb_of_seas <- vector("numeric",length=length(coefs))
  for (i in 1:length(coefs)){
    if(sum(months_all_pix[[i]] %in% c("Feb", "Dec", "Jan"))){
      nb_of_seas[i] <- nb_of_seas[i] + 1
    }
    
    if(sum(months_all_pix[[i]] %in% c("May", "Mar", "Apr"))){
      nb_of_seas[i] <- nb_of_seas[i] + 1
    }
    
    if(sum(months_all_pix[[i]] %in% c("Jun", "Jul", "Aug"))){
      nb_of_seas[i] <- nb_of_seas[i] + 1
    }
    
    if(sum(months_all_pix[[i]] %in% c("Sep", "Nov", "Oct"))){
      nb_of_seas[i] <- nb_of_seas[i] + 1
    }
  }  
  count_seans_and_var <- data.frame("nb_of_seas" = nb_of_seas, "nb_of_var" = nb_of_var)
  
  
  # Plot number of variables
  DF_nbdiffmeteo <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_meteo = nb_of_var)
  DF_nbdiffmeteo$nb_meteo <- as.factor(DF_nbdiffmeteo$nb_meteo)
  
  ggplot(data = DF_nbdiffmeteo, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=DF_nbdiffmeteo$nb_meteo),size=0.7) +
    scale_color_manual(values = c("0"=rainbow(4)[1],"1"=rainbow(4)[2], "2"=rainbow(4)[3], "3"=rainbow(4)[4])) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="Nb of different \nmeteo. variables",
         title = paste("Number of different meteorological variables, simple",model_name,"regression"),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  ggsave(file="D:/user/vogelj/Group project/Output/Plots/Number_of_variables_vpd_tmax_prec.png")
  
  # Plot number of months  
  DF_nbdiffseason <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_season = nb_of_seas)
  DF_nbdiffseason$nb_season <- as.factor(DF_nbdiffseason$nb_season)
  
  ggplot(data = DF_nbdiffseason, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=DF_nbdiffseason$nb_season),size=0.7) +
    scale_color_manual(values = c("0"=rainbow(5)[1],"1"=rainbow(5)[2], "2"=rainbow(5)[3], "3"=rainbow(5)[4], "4"=rainbow(5)[5])) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="Nb of different \nseasons",
         title = paste("Number of different seasons, simple",model_name,"regression"),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  ggsave(file="D:/user/vogelj/Group project/Output/Plots/Number_of_seasons.png")
  
  
  
  
  # ROC ####
  ########## 
  
  pred <- sapply(seq_along(work_pix), function(x) prediction(mypred[[x]], y1_test_list_red[[x]]))
  prf <- sapply(seq_along(work_pix), function(x) performance(pred[[x]], measure = "tpr", x.measure = "fpr"))
  par(mfrow=c(5,5))
  # plot(prf[[3]])
  for (i in sample(1:963,25)) (plot(prf[[3]]))
  ROCs <- sapply(seq_along(work_pix), function(x) plotROC(actuals=y1_test_list_red[[x]],predictedScores=fitted.results_model[[x]]))
  for (i in sample(1:963,25)) (plot(ROCs[[3]]))
  auc2 <- sapply(seq_along(work_pix), function(x) auc(y1_test_list_red[[x]],fitted.results_model[[x]]))
  auc <- sapply(seq_along(work_pix), function(x) performance(pred[[x]], measure = "auc"))
  # auc@y.values[[1]]