
# Visualisation ####
####################



world <- map_data("world")


# Rename to fit variable names (lwi)
csi <- csi_lwi
csi_adj <- csi_adj_lwi
speci <- speci_lwi
speci_adj <- speci_adj_lwi
sensi <- sensi_lwi
sensi_adj <- sensi_adj_lwi
coeff_kep <- coeff_kep_lwi
coefs <- coefs_lwi
cutoff_avg <- cutoff_lwi
work_pix <- work_pix_lwi


# # Plot miscla error ####
# 
# DF_miscla <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], miscla = mis_clas_err)
# 
# ggplot(data = DF_miscla, aes(x=lon, y=lat)) +
#   geom_polygon(data = world, aes(long, lat, group=group),
#                fill="white", color="black", size=0.1) +
#   geom_point(shape=15, aes(color=miscla),size=0.7) +
#   scale_color_gradient2(limits=c(min(mis_clas_err,na.rm=T),max(mis_clas_err,na.rm=T)),midpoint=min(mis_clas_err,na.rm=T)+(max(mis_clas_err,na.rm=T)-min(mis_clas_err,na.rm=T))/2,
#                         low = "yellow", mid = "red3", high = "black") +
#   theme(panel.ontop = F, panel.grid = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
#   ylab("Lat (°N)") +
#   xlab("Lon (°E)") +
#   coord_fixed(xlim = c(-120, 135),
#               ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
#               ratio = 1.3)+
#   labs(color="Misclass.\nerror",
#        title = paste("Misclassification error, simple",model_name,"regression"),
#        subtitle = paste("Bad yield threshold=", threshold,
#                         ", segregation threshold=", segreg_th, sep = ""))+
#   theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
#         legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
#   X11(width = 20, height = 7)
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Mis_class_error_lasso_interact_map.png")
# # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Mis_class_error_lasso_interact_seasonal_map.png")


# Plot specificity ####

DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], specificity = speci)

ggplot(data = DF_speci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=speci),size=0.7) +
  # scale_color_gradient2(limits=c(min(speci,na.rm=T),max(speci,na.rm=T)),midpoint=min(speci,na.rm=T)+(max(speci,na.rm=T)-min(speci,na.rm=T))/2,
  #                       low = "black", mid = "red3", high = "yellow") +
  geom_tile(aes(fill=speci)) +
  scale_fill_viridis(na.value="grey50")+
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Specificity_lasso_interact_map.png")
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Specificity_lasso_interact_seasonal_map.png")


DF_speci_adj <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], specificity = speci_adj)
ggplot(data = DF_speci_adj, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=speci_adj),size=0.7) +
  # scale_color_gradient2(limits=c(min(speci_adj,na.rm=T),max(speci_adj,na.rm=T)),midpoint=min(speci_adj,na.rm=T)+(max(speci_adj,na.rm=T)-min(speci_adj,na.rm=T))/2,
  #                       low = "black", mid = "red3", high = "yellow") +
  geom_tile(aes(fill=speci_adj)) +
  scale_fill_viridis(na.value="grey50")+
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Specif.",
       title = paste("Specificity (adj.), simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Specificity_adj_lasso_interact_map.png")


# Plot sensitivity ####

DF_sensi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], sensitivity = sensi)

ggplot(data = DF_sensi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=sensi),size=0.7) +
  # scale_color_gradient2(limits=c(min(sensi,na.rm=T),max(sensi,na.rm=T)),midpoint=min(sensi,na.rm=T)+(max(sensi,na.rm=T)-min(sensi,na.rm=T))/2,
  #                       low = "black", mid = "red3", high = "yellow") +
  geom_tile(aes(fill=sensi)) +
  scale_fill_viridis(na.value="grey50")+
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Sensitivity_lasso_interact_map.png")
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Sensitivity_lasso_interact_seasonal_map.png")


# Plot critical succes index (CSI) ####

DF_csi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], Critical_success_index = csi)

ggplot(data = DF_csi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=csi),size=0.7) +
  # scale_color_gradient2(limits=c(min(csi,na.rm=T),max(csi,na.rm=T)),midpoint=min(csi,na.rm=T)+(max(csi,na.rm=T)-min(csi,na.rm=T))/2,
  #                       low = "black", mid = "red3", high = "yellow") +
  geom_tile(aes(fill=csi)) +
  scale_fill_viridis(na.value="grey50")+
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/CSI_lasso_interact_map.png")
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/CSI_lasso_interact_seasonal_map.png")


DF_csi_adj <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], Critical_success_index = csi_adj)
ggplot(data = DF_csi_adj, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=csi_adj)) +
  scale_fill_viridis(na.value="grey50")+
  # scale_fill_gradient2(limits=c(min(csi_adj,na.rm=T),max(csi_adj,na.rm=T)),midpoint=min(csi_adj,na.rm=T)+(max(csi_adj,na.rm=T)-min(csi_adj,na.rm=T))/2,
  # low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="CSI",
       title = paste("Critical succes index (adj.), simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))  +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/CSI_adj_lasso_interact_map.png")


# # Correlation between yield and specificity or miscla error ####

if (model_name == "cv_fit_monthly_without_int_incl_ext") {yield <- matrix(Data_standardized$yield,965,1600)}

mean_yield <- apply(X=yield, MARGIN = 1, FUN = mean, na.rm=T)
# plot(mean_yield, mis_clas_err,
#      xlab="Mean Yield (kg/yr)", ylab="Miss-classification error",
#      main=paste("Scatterplot mean yield, miss-classification error\nbad yield threshold=", threshold,
#                 "\nsegregation threshold=", segreg_th, sep = ""))
plot(mean_yield, speci,
     xlab="Mean Yield (kg/yr)", ylab="Specificity",
     main=paste("Scatterplot mean yield, Specificity\nbad yield threshold=", threshold,
                "\nsegregation threshold=", segreg_th, sep = ""))


# pairs(cbind(mean_yield, mis_clas_err, speci),
#       main=paste(model_name, " regression, bad yield thr.=", threshold,
#                  "\nsegreg. thr.=", segreg_th, sep = ""))

# cor(mean_yield, mis_clas_err)
cor(mean_yield, speci)


# Map of the number of coefficients kept (Lasso) #####


DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = coeff_kep)

ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=coeff_kep),size=0.7) +
  # scale_color_gradientn(limits=c(0,max(DF_numbcoeff[,3])), 
  # colours=c(gray.colors(1),topo.colors(23)[-c(1,3,5,16:23)],rev(heat.colors(10))) ,values=rescale(0:22,c(0,1))) + # mixed visualisation
  # colours=c(gray.colors(1),topo.colors(20)[-c(13:20)],rev(heat.colors(10))) ,values=rescale(0:22,c(0,1))) + # mixed visualisation v1
  
  # lambdaHat1Std
  # scale_color_gradientn(limits=c(0,15),
  #                       colours=c(gray.colors(1),topo.colors(9)[-c(8,9)],rev(heat.colors(7))) ,values=rescale(0:15,c(0,1)),
  #                       breaks=c(0,3,6,9,12,15),labels=c("0","3","6","9","12",">=15")) + # mixed visualisation with cutoff
  
  # lambdaHat
  # scale_color_gradientn(limits=c(0,40),
  #                       colours=c(gray.colors(1),topo.colors(9)[-c(8,9)],rev(heat.colors(7))) ,values=rescale(0:40,c(0,1)),
  #                       breaks=c(0,10,20,30,40),labels=c("0","10","20","30","40")) + # mixed visualisation with cutoff
  
  # scale_color_gradientn(limits=c(0,15), # cut off high values for better visualisation of the rest
  # colours=c(gray.colors(1),topo.colors(14)) ,values=rescale(0:15,c(0,1))) + # topo color scheme
  #                      colours=c(gray.colors(1),rainbow(14)) ,values=rescale(0:15,c(0,1))) + # rainbow color scheme
  # scale_color_gradient2(limits=c(0,15),midpoint=15/2, # cut off high values for better visualisation of the rest
  #                      low = "blue", mid = "yellow", high = "red3") +
# scale_color_gradient2(limits=c(0,max(DF_numbcoeff[,3])),midpoint=max(DF_numbcoeff[,3])/2,
# low = "blue", mid = "yellow", high = "red3") +
  geom_tile(aes(fill=coeff_kep)) +
  scale_fill_viridis(na.value="grey50")+
theme(panel.ontop = F, panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Nb of coefficients",
       title = paste("Number of coefficients kept",model_name),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
# d + scale_fill_discrete(breaks=c(2,4,6,8,10), labels = c("A", "B", "C", "D", "E"))

ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_coefficients_lasso_interact_map.png")
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_coefficients_lasso_interact_seasonal_map.png")
# }
plot(table(coeff_kep))  # overview of distribution of the number of coefficients


# Map of number of interactions ####

DF_numb_interact <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], num_interact = num_interact)

ggplot(data = DF_numb_interact, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=num_interact),size=0.7) +
  # scale_color_gradient2(limits=c(0,max(DF_numb_interact[,3])),midpoint=max(DF_numb_interact[,3]/2),
  # low = "blue", mid = "yellow", high = "red3") +
  # scale_color_gradientn(limits=c(0,16), breaks=c(0,4,8,12,16),labels=c("0","4","8","12",">=16"), # lamdaHat1Std
  #                       colours=c(gray.colors(1),topo.colors(10)[-c(9,10)],rev(heat.colors(7))) ,values=rescale(0:16,c(0,1))) + # mixed visualisation with cutoff
  # scale_color_gradientn(limits=c(0,64), breaks=c(0,12,24,36,48),labels=c("0","12","24","36",">=48"), # lambdaHat
  #                       colours=c(gray.colors(1),topo.colors(10)[-c(9,10)],rev(heat.colors(7))) ,values=rescale(0:64,c(0,1))) + # mixed visualisation with cutoff
  geom_tile(aes(fill=num_interact)) +
  scale_fill_viridis(na.value="grey50")+
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
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_coefficient_interactions_lasso_interact_seasonal_map.png")
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
  }
  
  # account for extreme indicators  
  str_var[str_var %in% c("dtr","frs","tn10p","tnn","tx90p","txx")] <- "tmax"     # Group all extreme temperature indicators into same group as tmax
  str_var[str_var %in% "rx5"] <- "pr" # Group all extreme precipitation indicators into same group as precipitation
  
  nb_of_var[j] <- length(unique(str_var))
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

# combinations of months and variables
coefs_seas <- sapply(1:length(coefs), function(x) names(numLevels_list[[x]])[coefs[[x]]$mainEffects$cont])
coefs_seas_vec <- unlist(coefs_seas)
png(filename="D:/user/vogelj/Group_project/Output/Plots/Combination_variable_season_barplot.png",res=2000,units="cm",width=15,height=22)
par(mar=c(5,7,4,2))
barplot(sort(table(coefs_seas_vec)),horiz=T,las=1,col="ForestGreen",xlab="Number of pixels, where variable is included in the model",cex.names=0.8)
dev.off()

nb_of_seas <- vector("numeric",length=length(coefs))
for (i in 1:length(coefs)){
  if(sum(months_all_pix[[i]] %in% c("Feb", "Dec", "Jan","win"))){
    nb_of_seas[i] <- nb_of_seas[i] + 1
  }
  
  if(sum(months_all_pix[[i]] %in% c("May", "Mar", "Apr","spr"))){
    nb_of_seas[i] <- nb_of_seas[i] + 1
  }
  
  if(sum(months_all_pix[[i]] %in% c("Jun", "Jul", "Aug","sum"))){
    nb_of_seas[i] <- nb_of_seas[i] + 1
  }
  
  if(sum(months_all_pix[[i]] %in% c("Sep", "Nov", "Oct","aut"))){
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
  # geom_point(shape=15, aes(color=DF_nbdiffmeteo$nb_meteo),size=0.7) +
  # scale_color_manual(values = c("0"=rainbow(4)[1],"1"=rainbow(4)[2], "2"=rainbow(4)[3], "3"=rainbow(4)[4])) +
  geom_tile(aes(fill=DF_nbdiffmeteo$nb_meteo)) +
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_variables_vpd_tmax_prec.png")
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_variables_vpd_tmax_prec_seasonal.png")

# Plot number of months  
DF_nbdiffseason <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_season = nb_of_seas)
DF_nbdiffseason$nb_season <- as.factor(DF_nbdiffseason$nb_season)

ggplot(data = DF_nbdiffseason, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=DF_nbdiffseason$nb_season),size=0.7) +
  # scale_color_manual(values = c("0"=rainbow(5)[1],"1"=rainbow(5)[2], "2"=rainbow(5)[3], "3"=rainbow(5)[4], "4"=rainbow(5)[5])) +
  geom_tile(aes(fill=DF_nbdiffseason$nb_season)) +
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_seasons.png")
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_seasons_seasonal.png")



# Cutoff level plot
cutoffs <- rep(NA,965)
cutoffs[work_pix] <- sapply(work_pix, function(x){cutoff_avg[x]})
DF_cutoff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], cutoff = cutoffs)

ggplot(data = DF_cutoff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_point(shape=15, aes(color=cutoffs),size=0.7) +
  # scale_color_gradient2(limits=c(min(cutoffs,na.rm=T),max(cutoffs,na.rm=T)),midpoint=min(cutoffs,na.rm=T)+(max(cutoffs,na.rm=T)-min(cutoffs,na.rm=T))/2,
  #                       low = "black", mid = "red3", high = "yellow") +
  geom_tile(aes(fill=cutoffs)) +
  scale_fill_viridis(na.value="grey50")+
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Cutoff",
       title = paste("Cutoff",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Cutoff_lasso_interact_map.png")


# Plots of singular extreme indicators ####
if (model_name == "cv_fit_monthly_without_int_incl_ext"){
  

  # number of variables: which are included
  dtr_all <- rep(NA,965);  frs_all <- rep(NA,965); rx5_all <- rep(NA,965)
  tn10p_all <- rep(NA,965); tnn_all <- rep(NA,965);  tx90p_all <- rep(NA,965);  
  txx_all <- rep(NA,965)
  for (j in  1:length(coefs)){
    coefs_str <- names(numLevels_list[[j]])[coefs[[j]]$mainEffects$cont]
    
    dtr_all[j] <- ifelse(sum(coefs_str %in% "dtr")==1,T,F)
    frs_all[j] <- ifelse(sum(coefs_str %in% "frs")==1,T,F)
    rx5_all[j] <- ifelse(sum(coefs_str %in% "rx5")==1,T,F)
    tn10p_all[j] <- ifelse(sum(coefs_str %in% "tn10p")==1,T,F)
    tnn_all[j] <- ifelse(sum(coefs_str %in% "tnn")==1,T,F)
    tx90p_all[j] <- ifelse(sum(coefs_str %in% "tx90p")==1,T,F)
    txx_all[j] <- ifelse(sum(coefs_str %in% "txx")==1,T,F)
  }
  
  DF_ext_ind<- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], dtr = dtr_all, frs=frs_all,rx5=rx5_all,tn10p=tn10p_all,tnn=tnn_all,tx90p=tx90p_all,txx=txx_all)
  apply(DF_ext_ind,2,sum)
  
  ggplot(data = DF_ext_ind, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    # geom_point(shape=15, size=0.7, 
               # aes(color=DF_ext_ind$dtr)) +
               # aes(color=DF_ext_ind$frs)) +
               # aes(color=DF_ext_ind$rx5)) +
               # aes(color=DF_ext_ind$tn10p)) +
               # aes(color=DF_ext_ind$tnn)) +
               # aes(color=DF_ext_ind$tx90p)) +
               # aes(color=DF_ext_ind$txx)) +
    geom_tile(aes(fill=DF_ext_ind$dtr)) +
    # geom_tile(aes(fill=DF_ext_ind$frs)) +
    # geom_tile(aes(fill=DF_ext_ind$rx5)) +
    # geom_tile(aes(fill=DF_ext_ind$tn10p)) +
    # geom_tile(aes(fill=DF_ext_ind$tnn)) +
    # geom_tile(aes(fill=DF_ext_ind$tx90p)) +
    # geom_tile(aes(fill=DF_ext_ind$txx)) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(color="Nb of different \nmeteo. variables",
         title = paste("Pixels with dtr",model_name),
         # title = paste("Pixels with frs",model_name),
         # title = paste("Pixels with rx5",model_name),
         # title = paste("Pixels with tn10p",model_name),
         # title = paste("Pixels with tnn",model_name),
         # title = paste("Pixels with tx90p",model_name),
         # title = paste("Pixels with txx",model_name),
         subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 7)
  ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_dtr.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_frs.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_rx5.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_tn10p.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_tnn.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_tx90p.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_txx.png")

}


# source("./Code/Problematic_pixels.r")


# ROC ####
########## 

# pred <- sapply(seq_along(work_pix), function(x) prediction(mypred[[x]], y1_test_list_red[[x]]))
# prf <- sapply(seq_along(work_pix), function(x) performance(pred[[x]], measure = "tpr", x.measure = "fpr"))
# par(mfrow=c(5,5))
# # plot(prf[[3]])
# for (i in sample(1:963,25)) (plot(prf[[3]]))
# # takes long
# # ROCs <- sapply(seq_along(work_pix), function(x) plotROC(actuals=y1_test_list_red[[x]],predictedScores=fitted.results_model[[x]]))
# # for (i in sample(1:963,25)) (plot(ROCs[[3]]))
# auc2 <- sapply(seq_along(work_pix), function(x) auc(y1_test_list_red[[x]],fitted.results_model[[x]]))
# auc <- sapply(seq_along(work_pix), function(x) performance(pred[[x]], measure = "auc"))
# # auc@y.values[[1]]



