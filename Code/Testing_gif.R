# Version 1: with ggplot2 function

# https://cran.r-project.org/web/packages/magick/vignettes/intro.html

library(gapminder)
library(ggplot2)
library(magick)

top10variables <- names(sort(table(coefs_seas_vec), decreasing = T)[1:10])
allvariables <- colnames(Model_data)[-1]

DF_list <- vector("list",length=10)

for (varia in 1:10) {
# for (varia in 1:length(allvariables)) {
  
  varia_name <- top10variables[varia]
  # varia_name <- allvariables[varia]
  varia_in_pix <- numeric()
  # plots <- vector("list",length=(length(allvariables)))
  for (pix in 1:pix_num) {
    varia_in_pix[pix] <- (varia_name %in% row.names(coeff[[pix]])[which(coeff[[pix]]!=0)])
  }#end for pix
  
  DF_var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], var_in = varia_in_pix)
  DF_var$var_in <- as.factor(DF_var$var_in)
  DF_list[[varia]] <- DF_var
} 

img <- image_graph(600, 340, res = 96)
out <- lapply(1:10, function(i){
  p <- ggplot(data = DF_list[[i]], aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_tile(aes(fill=DF_var$var_in)) +
    scale_fill_manual(values = c("1"="#fc8d62", "0"="#8da0cb"),
                      label= c("1"="Yes", "0"="No"),
                      breaks=c("1","0")) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(fill=paste0(top10variables[varia],"\nselected")
         #,title = paste("Number of variables kept, simple",model_name,"regression, "),
         #subtitle = paste("Monthly meteo var + extreme indices, ",lambda_val, sep = "")
    )+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) #+
    # X11(width = 20, height = 6)
  # ggsave(filename=paste0("D:/user/vogelj/Group_project/Output/Plots/",top10variables[varia],".jpg"))
    print(p)
})
dev.off()
animation <- image_animate(img, fps = 2, optimize = TRUE)
print(animation)
magick::image_write(animation, path="D:/user/vogelj/Group_project/Output/Plots/myanimation.gif")




# Version 2: with saveGIF from animation librabry


top10variables <- names(sort(table(coefs_seas_vec), decreasing = T)[1:10])
allvariables <- colnames(Model_data)[-1]

library(animation)
saveGIF(expr= {for (varia in 1:10) {
  # for (varia in 1:length(allvariables)) {
  
  varia_name <- top10variables[varia]
  # varia_name <- allvariables[varia]
  varia_in_pix <- numeric()
  plots <- vector("list",length=(length(allvariables)))
  for (pix in 1:pix_num) {
    varia_in_pix[pix] <- (varia_name %in% row.names(coeff[[pix]])[which(coeff[[pix]]!=0)])
  }#end for pix
  
  DF_var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], var_in = varia_in_pix)
  DF_var$var_in <- as.factor(DF_var$var_in)
  
  # plots[[varia]] <-
  ggplot(data = DF_var, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_tile(aes(fill=DF_var$var_in)) +
    scale_fill_manual(values = c("1"="#fc8d62", "0"="#8da0cb"),
                      label= c("1"="Yes", "0"="No"),
                      breaks=c("1","0")) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    ylab("Lat (°N)") +
    xlab("Lon (°E)") +
    coord_fixed(xlim = c(-120, 135),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1.3)+
    labs(fill=paste0(varia_name,"\nselected")
         #,title = paste("Number of variables kept, simple",model_name,"regression, "),
         #subtitle = paste("Monthly meteo var + extreme indices, ",lambda_val, sep = "")
    )+
    theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    X11(width = 20, height = 6)
  # ggsave(filename=paste0("D:/user/vogelj/Group_project/Output/Plots/",varia_name,".pdf"))
}},
movie.name="D:/user/vogelj/Group_project/Output/Plots/animation.gif"
)

