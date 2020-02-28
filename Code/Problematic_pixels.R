
# Investigate the link between pixels with low yield, only positive predictions and only 1 or 2 variable-models


# Rename to fit variable names (lwi)
world <- map_data("world")
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



# Pixel with only positive predictions
all_n_adj_lwi <- rep(NA,965)
all_n_adj_lwi[work_pix_lwi] <- fn_adj_lwi + tn_adj_lwi # all (true and false) predicted negatives
no_neg <- rep(NA,965)
no_neg[work_pix_lwi] <- all_n_adj_lwi == 0  # pixels with no negative predictions
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_pred.png")



# Pixels with only 1 or 2 variable-models

plot(table(coeff_kep))

coeff_kep_low <- coeff_kep <= 2
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_coeff_kep_low.png")



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


# Map

DF_bad_pix <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], all_bad = all_bad, no_neg_coeff_low= no_neg_coeff_low, no_neg_yield_low=no_neg_yield_low)
ggplot(data = DF_yield_low, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  # geom_tile(aes(fill= all_bad)) +
  # geom_tile(aes(fill= no_neg_coeff_low)) +
  geom_tile(aes(fill= no_neg_yield_low)) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  # labs(title = paste("Pixel with only positive predictions, low yield and only 1 or 2 variable-models for",model_name))+
  # labs(title = paste("Pixel with only positive predictions and only 1 or 2 variable-models for",model_name))+
  labs(title = paste("Pixel with only positive predictions and low yield for",model_name))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_low yield_low_coeff.png")
  # ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_low_coeff.png")
  ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Pixels_no_neg_low yield.png")
