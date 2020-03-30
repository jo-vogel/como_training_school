# Author: Johannes Vogel
# Investigate differences between glinternet model with and without extreme indices

# Run with glinternet section from Model_comparison.r

hist(csi_adj_lwi_monthly_no_int); mean(csi_adj_lwi_monthly_no_int,na.rm=T);
hist(csi_adj_lwi_monthly_without_int_incl_ext); mean(csi_adj_lwi_monthly_without_int_incl_ext)
hist(coeff_kep_lwi_monthly_no_int); mean(coeff_kep_lwi_monthly_no_int)
hist(coeff_kep_lwi_monthly_without_int_incl_ext); mean(coeff_kep_lwi_monthly_without_int_incl_ext)

csi_diff <- csi_adj_lwi_monthly_without_int_incl_ext - csi_adj_lwi_monthly_no_int; mean(csi_diff,na.rm=T)
coefs_diff <- coeff_kep_lwi_monthly_without_int_incl_ext - coeff_kep_lwi_monthly_no_int; mean(coefs_diff)


DF_csi_adj <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], Critical_success_index = csi_diff)
ggplot(data = DF_csi_adj, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=csi_diff)) +
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
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/CSI_adj_glinternet_diff_with_without_int_map.png")




DF_numbcoeff_diff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = coefs_diff)

ggplot(data = DF_numbcoeff_diff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
geom_tile(aes(fill=coefs_diff)) +
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

ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Number_of_coefficients_glinternet_diff_with_without_int_map.png")