source("get_first_coeff_function.R")

nb_coeff_to_keep <- 3 #if you want to plot the 1st, second, third variable...

coef_to_plot <- 1 #between 1 and coeff to keep
first_var <- numeric()
month_first_var <- numeric()

for (pix in 1:length(coefs)) {
  first_var[pix] <- get_firstcoeffs(coefs[[pix]], nb_of_coeff = nb_coeff_to_keep)[coef_to_plot,1]
  month_first_var[pix] <- get_firstcoeffs(coefs[[pix]], nb_of_coeff = nb_coeff_to_keep)[coef_to_plot,2]
}#end for pix

DF_1var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2],
                      first_var = first_var, month_first_var=month_first_var)

ggplot(data = DF_1var, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_point(shape=15, aes(color=first_var, alpha=month_first_var), size=1) +
  scale_color_manual(values = c("pr"="blue", "vpd"="green", "tmax"="red")) +
  scale_alpha_manual(values = c("Aug_Y1"=(1/13), "Sep_Y1"=(2/13),"Oct_Y1"=(3/13), "Nov_Y1"=4/13,"Dec_Y1"=(5/13),"Jan_Y2"=(6/13),
                                "Feb_Y2"=(7/13),"Mar_Y2"=(8/13),"Apr_Y2"=(9/13), "May_Y2"=(10/13),"Jun_Y2"=(11/13),
                                "Jul_Y2"=(12/13), "Aug_Y2"=(13/13)),
                     breaks = c("Aug_Y1", "Sep_Y1","Oct_Y1", "Nov_Y1","Dec_Y1","Jan_Y2", "Feb_Y2","Mar_Y2",
                                "Apr_Y2", "May_Y2","Jun_Y2","Jul_Y2", "Aug_Y2")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="1st var", alpha="Month",
       title = paste("Most important predictor: name and month, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold, sep = ""))+
  guides(color=guide_legend(order = 1), alpha=guide_legend(order = 2, ncol = 2))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 23, height = 7)