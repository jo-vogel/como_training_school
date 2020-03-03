
# Investigate the link between pixels with low yield, only positive predictions and only 1 or 2 variable-models

model_name <- "Lasso (glmnet)"

world <- map_data("world")
if (model_name == "cv_fit_monthly_without_int_incl_ext" | model_name == "Lasso (glmnet)") {
  yield <- matrix(Data_standardized$yield,965,1600)
  }
mean_yield <- apply(X=yield, MARGIN = 1, FUN = mean, na.rm=T)


# Rename to fit variable names (lwi) for glinternet
csi <- csi_lwi
csi_adj <- csi_adj_lwi
speci <- speci_lwi
speci_adj <- speci_adj_lwi
sensi <- sensi_lwi
sensi_adj <- sensi_adj_lwi
coeff_kep <- coeff_kep_lwi
# coefs <- coefs_lwi
# cutoff_avg <- cutoff_lwi
work_pix <- work_pix_lwi
fn_adj <- fn_adj_lwi
tn_adj <- tn_adj_lwi
mypred <- mppred_lwi


# Rename to fit variable names for glmnet
csi <- csi_simplelasso
csi_adj <- csi_simplelasso_adj
speci <- speci_simplelasso
speci_adj <- speci_simplelasso_adj
sensi <- sensi_simplelasso
sensi_adj <- sensi_simplelasso_adj
mypred <- pred_simplelasso
work_pix <- 965

Model_chosen <- lasso_model_lambda1se
pix_model_failed <- numeric(length = pix_num)
coeff  <-list()
for (pixel in 1:pix_num) {
  if(is.character(Model_chosen[[pixel]])){
    pix_model_failed[pixel] <- 1
    coeff[[pixel]] <- "No model"
  } else {
    
    coeff[[pixel]] <- coefficients(Model_chosen[[pixel]])
  }#end if else pb model
}#end pixel

number_coeff_kept <- function(coeff_list){#give number of coeff !=0
  return(length(which(abs(coeff_list[-1,])>0)))
}
nb_coeff_kept <- numeric()
for (pixel in 1:pix_num) {
  if(is.character(Model_chosen[[pixel]])){
    nb_coeff_kept[pixel] <- NA
  } else {
    nb_coeff_kept[pixel] <- number_coeff_kept(coeff[[pixel]])
  }#end if else pb model
}#end pixel
coeff_kep <- nb_coeff_kept 

tn_adj_simpleLasso <- numeric(965)
fn_adj_simpleLasso <- numeric(965)
for(pix in 1:nb_pix_simplelasso){
  tn_adj_simpleLasso[pix] <- con_tab_simplelasso_adj[[pix]]["0","0"] 
  fn_adj_simpleLasso[pix] <- con_tab_simplelasso_adj[[pix]]["0","1"]
  message('note that rows are predicted and columns observed here (reverse in case of glinternet)')
  if(is.na(con_tab_simplelasso_adj[[pix]]["0","0"])){ #no extreme event forecasted => no first line in contengency table
    tn_adj_simpleLasso[pix] <- 0
    fn_adj_simpleLasso[pix] <- 0
  }
}
tn_adj <- tn_adj_simpleLasso
fn_adj <- fn_adj_simpleLasso




# Pixel with only positive predictions
# all_n_adj <- rep(NA,965)
# all_n_adj[work_pix] <- fn_adj + tn_adj # all (true and false) predicted negatives
all_n_adj <- fn_adj + tn_adj # all (true and false) predicted negatives
no_neg <- rep(NA,965)
# no_neg[work_pix] <- all_n_adj == 0  # pixels with no negative predictions
no_neg <- all_n_adj == 0  # pixels with no negative predictions
sum(no_neg) # there are 112 pixels with no negative predictions

DF_no_neg <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], none_neg = no_neg)
ggplot(data = DF_no_neg, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=none_neg)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Pixel with no predicted negatives",
       title = paste("Pixel with no predicted negatives",model_name),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_pred_glinternet.png")
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_pred_glmnet.png")

pred_min <- sapply(mypred,min)
sum(pred_min>cutoff_simple_lasso) # should be identical to no_neg

# png(filename = "D:/user/vogelj/Group_project/Output/Plots/hist_min_pred_val_glinternet.png")
png(filename = "D:/user/vogelj/Group_project/Output/Plots/hist_min_pred_val_glmnet.png")
# hist(pred_min,main="Minimum predicted value for glinternet")
hist(pred_min,main="Minimum predicted value for glmnet")
dev.off()
# png(filename = "D:/user/vogelj/Group_project/Output/Plots/boxplot_min_pred_val_glinternet.png")
png(filename = "D:/user/vogelj/Group_project/Output/Plots/boxplot_min_pred_val_glmnet.png")
# boxplot(pred_min,main="Minimum predicted value for glinternet")
boxplot(pred_min,main="Minimum predicted value for glmnet")
dev.off()


# Pixels with only 1 or 2 variable-models

plot(table(coeff_kep))

coeff_kep_low <- coeff_kep <= 2
# coeff_kep_low <- coeff_kep == 0
sum(coeff_kep_low) # there are 289 pixels with only 1 or 2 coefficients in the model

DF_coeff_kep_low <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_low = coeff_kep_low)
ggplot(data = DF_coeff_kep_low, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill= coeff_low)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Pixel with only 1 or 2 variable-models",
       title = paste("Pixel with only 1 or 2 variable-models",model_name),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
# ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_coeff_kep_low_glinternet.png")
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_coeff_kep_0_glmnet.png")



# Pixels with low yield

# pick the lowest 112 pixels (to have a number that is in the range of the number of pixel with only positive predictions)

hist(mean_yield)

yield_sorted <- sort(mean_yield)
pix_yield_low <- mean_yield <= yield_sorted[[112]]

DF_yield_low <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], yield_low = pix_yield_low)
ggplot(data = DF_yield_low, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill= pix_yield_low)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Pixel with low yield",
       title = paste("Pixel with low yield",model_name),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_low yield.png")




# compare how the 3 criteria are connected


all_3 <- no_neg + coeff_kep_low + pix_yield_low
all_bad <- rep(F,965)
all_bad[which(all_3 == 3)] <- T
sum(all_bad)

all_con <- cbind(no_neg , coeff_kep_low , pix_yield_low)
no_neg_yield_low <- (all_con[,1]==T & all_con[,3]==T)
sum(no_neg_yield_low)
no_neg_coeff_low <- (all_con[,1]==T & all_con[,2]==T)
sum(no_neg_coeff_low)
coeff_low_yield_low <- (all_con[,2]==T & all_con[,3]==T)
sum(coeff_low_yield_low)

# Map

DF_bad_pix <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], all_bad = all_bad, no_neg_coeff_low= no_neg_coeff_low, no_neg_yield_low=no_neg_yield_low)
ggplot(data = DF_yield_low, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_tile(aes(fill= all_bad)) +
  geom_tile(aes(fill= no_neg_coeff_low)) +
  # geom_tile(aes(fill= no_neg_yield_low)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  # labs(title = paste("Pixel with only positive predictions, low yield and only 1 or 2 variable-models for",model_name))+
  labs(title = paste("Pixel with only positive predictions and only 1 or 2 variable-models for",model_name))+
  # labs(title = paste("Pixel with only positive predictions and low yield for",model_name))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_low yield_low_coeff_glinternet.png")
  ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_low_coeff_glinternet.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_low yield_glinternet.png")

  
  