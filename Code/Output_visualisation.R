##########################################
# Visualisation of results 
# by Pauline, partially edited by Johannes
##########################################


# Load model output ####
########################

# source('./Code/Lasso_interact_global.R') # takes ca. 3 hours
source('./Code/Lasso_interact_global_preparation.R')
load('D:/user/vogelj/Group project/Code/Workspaces/cv_fit_complete.RData ')



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

# Plot miscla error ####

world <- map_data("world")
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
  
  # Plot specificity error ####
  
  DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = coeff_kep)
  
  ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_point(shape=15, aes(color=coeff_kep)) +
    scale_color_gradient(limits=c(3,33),
                         low = "pink", high = "darkblue") +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="Nb of var.",
         title = paste("Number of variables kept, simple",model_name,"regression"),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
# }