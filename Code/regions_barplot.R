# Author Johannes Vogel

# Run Final_Plots.r

library(raster)
library(rgdal)

# Create matrix with all coordinates
path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
nh_files <- list.files(path=path_to_NH_files,pattern="*NH.nc") # all files from northern hemisphere
nh_data <- lapply(1:length(nh_files),function(x){nc_open(paste0(path_to_NH_files,"/",nh_files[x]))})
lat_all <- ncvar_get(nh_data[[1]],"lat")
lon_all <- ncvar_get(nh_data[[1]],"lon")
lati_all <- rep(lat_all,each=length(lon_all))
long_all <- rep(lon_all,length(lat_all)) # coordinates rearranged
coord_all <- cbind(long_all,lati_all)

# Connect variables to spatial locations
coord_subset_temp <- cbind(coord_subset,paste(coord_subset[,1],coord_subset[,2]))
coord_all_temp <- cbind(coord_all,paste(coord_all[,1],coord_all[,2]))
loc_pix <- which(coord_all_temp[,3] %in% coord_subset_temp [,3]) # locations of our pixels in the whole coordinate set

# Create continent polygons
border <- readOGR('D:/user/vogelj/Data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')	
continents <- readOGR("D:/user/vogelj/Data/continent_shapefile/continent.shp") # from https://www.arcgis.com/home/item.html?id=5cf4f223c4a642eb9aa7ae1216a04372
africa <- subset(continents,subset=continents@data[["CONTINENT"]]=="Africa")
europe <- subset(continents,subset=continents@data[["CONTINENT"]]=="Europe")
no_am <- subset(continents,subset=continents@data[["CONTINENT"]]=="North America")
asia <- subset(continents,subset=continents@data[["CONTINENT"]]=="Asia")

# connect 969 pixels to their coordinates
work_pix <- 1:969
coord_assigned <- cbind(coord_all,rep(NA,320*76))
for (i in seq_along(work_pix)){
  coord_assigned[loc_pix[i],3] <- i
}

# Extract pixels by continent
loc_mat <- matrix(as.numeric(coord_assigned[,3]),nrow=320,ncol=76)
loc_ras <- raster(t(loc_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all), ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))
loc_afr_pixels <- extract(loc_ras,africa)
loc_eur_pixels <- extract(loc_ras,europe)
loc_no_am_pixels <- extract(loc_ras,no_am)
loc_asia_pixels <- extract(loc_ras,asia)

sum(!is.na(loc_eur_pixels[[1]])) # 233 pixels
sum(!is.na(loc_no_am_pixels[[1]])) # 421 pixels
sum(!is.na(loc_afr_pixels[[1]])) # 45 pixels
sum(!is.na(loc_asia_pixels[[1]])) # 265 pixels
# 45+233+421+265=964; 969-964: 5 pixel are missing

loc_afr_pixels_num <- loc_afr_pixels[[1]][!is.na(loc_afr_pixels[[1]])] 
loc_eur_pixels_num <- loc_eur_pixels[[1]][!is.na(loc_eur_pixels[[1]])]
loc_no_am_pixels_num <- loc_no_am_pixels[[1]][!is.na(loc_no_am_pixels[[1]])]
loc_asia_pixels_num <- loc_asia_pixels[[1]][!is.na(loc_asia_pixels[[1]])]

# add missing points: see extra file Missing_points.r
loc_eur_pixels_num <- c(loc_eur_pixels_num,529)
loc_no_am_pixels_num <- c(loc_no_am_pixels_num,470)
loc_asia_pixels_num <- c(loc_asia_pixels_num,33, 372, 642)

# Get data into right format for the plot (there should be a simpler approach than this)
coefs_seas <- sapply(1:length(coeff), function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])
coefs_seas_vec <- unlist(coefs_seas)

coefs_seas_afr <- sapply(loc_afr_pixels_num, function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])
coefs_seas_afr_vec <- unlist(coefs_seas_afr )
coefs_seas_eur <- sapply(loc_eur_pixels_num, function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])
coefs_seas_eur_vec <- unlist(coefs_seas_eur )
coefs_seas_no_am <- sapply(loc_no_am_pixels_num, function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])
coefs_seas_no_am_vec <- unlist(coefs_seas_no_am )
coefs_seas_asia <- sapply(loc_asia_pixels_num, function(x) names(numLevels_list[[x]])[coeff[[x]][-1]!=0])
coefs_seas_asia_vec <- unlist(coefs_seas_asia )
# coefs_seas_mat <- cbind(coefs_seas_afr_vec,coefs_seas_eur_vec,coefs_seas_no_am_vec,coefs_seas_asia_vec)

coefs_seas_afr_tab <- as.data.frame(table(coefs_seas_afr_vec))
coefs_seas_eur_tab <-  as.data.frame(table(coefs_seas_eur_vec))
coefs_seas_no_am_tab <- as.data.frame(table(coefs_seas_no_am_vec))
coefs_seas_asia_tab <- as.data.frame(table(coefs_seas_asia_vec))
coefs_seas_vec_tab <- as.data.frame(table(coefs_seas_vec))

# correct variable names
list_coefs <- list(coefs_seas_afr_tab,coefs_seas_eur_tab,coefs_seas_no_am_tab,coefs_seas_asia_tab,coefs_seas_vec_tab)
for (i in 1:length(list_coefs)){
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="vpd", replacement = "VPD")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="pr", replacement = "Pr")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="APr", replacement = "Apr") # recorrect April
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tmax", replacement = "Tmax")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="txx", replacement = "TXx")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tnn", replacement = "TNn")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="rx5", replacement = "Rx5day")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tx90p", replacement = "TX90p")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tn10p", replacement = "TN10p")
}
coefs_seas_asia_tab[,1] <- coefs_seas_vec_tab[,1]

colnames(coefs_seas_afr_tab) <- c("Variables","Freq_Afr")
colnames(coefs_seas_eur_tab) <- c("Variables","Freq_Eur")
colnames(coefs_seas_no_am_tab) <- c("Variables","Freq_No_Am")
colnames(coefs_seas_asia_tab) <- c("Variables","Freq_Asia")
colnames(coefs_seas_vec_tab) <- c("Variables","Freq_All")


coefs_all_cont <- merge(coefs_seas_afr_tab,coefs_seas_eur_tab,all=T,by="Variables")
coefs_all_cont <- merge(coefs_all_cont,coefs_seas_no_am_tab,all=T,by="Variables")
coefs_all_cont <- merge(coefs_all_cont,coefs_seas_asia_tab,all=T,by="Variables")
coefs_all_cont <- merge(coefs_all_cont,coefs_seas_vec_tab,all=T,by="Variables")
coefs_all_cont <- coefs_all_cont[order(coefs_all_cont$Freq_All,decreasing=F),]
coefs_all_cont_mat <- as.matrix(coefs_all_cont)
row.names(coefs_all_cont_mat) <- coefs_all_cont_mat[,1]
coefs_all_cont_mat[which(is.na(coefs_all_cont_mat))] <- 0
# coefs_all_cont_mat <- coefs_all_cont_mat[,-1]
coefs_all_cont_mat <- coefs_all_cont_mat[,-c(1,6)]




# Plot coefficients by continent
pdf(file="D:/user/vogelj/Group_project/Output/Plots/barplot_variables_lambda1se_all_continents.pdf",
    width=6,height=8)

par(mar=c(5,7,1,1))
barplot(sort(table(coefs_seas_afr_vec)),horiz=T,las=1,col="grey",main="Africa",
        xlab="Number of grid points, where variable is included in the model",cex.names=0.6)
barplot(sort(table(coefs_seas_eur_vec)),horiz=T,las=1,col="grey",main="Europe",
        xlab="Number of grid points, where variable is included in the model",cex.names=0.6)
barplot(sort(table(coefs_seas_no_am_vec)),horiz=T,las=1,col="grey",main="North America",
        xlab="Number of grid points, where variable is included in the model",cex.names=0.6)
barplot(sort(table(coefs_seas_asia_vec)),horiz=T,las=1,col="grey",main="Asia",
        xlab="Number of grid points, where variable is included in the model",cex.names=0.6)

barplot(t(coefs_all_cont_mat),horiz=T,las=1,col=c("brown3","DarkOrange2","goldenrod3","burlywood1"),
        xlab="Number of grid points, where variable is included in the model",cex.names=0.6,
        legend.text=c("Africa","Europe","North America","Asia"),args.legend =list(x= "bottomright"))
dev.off()

# only the last plot
pdf(file="D:/user/vogelj/Group_project/Output/Plots/barplot_variables_lambda1se_all_continents_summary.pdf",
    width=6,height=8)
par(mar=c(5,7,1,1))
barplot(t(coefs_all_cont_mat),horiz=T,las=1,col=c("brown3","DarkOrange2","goldenrod3","burlywood1"),
        xlab="Number of grid points, where variable is included in the model",cex.names=0.6,
        legend.text=c("Africa","Europe","North America","Asia"),args.legend =list(x= "bottomright"))
dev.off()



# 4 Plots in one line (all, N.America, Europe, Asia) ####

line2user <- function(line, side) { # from https://stackoverflow.com/questions/14660372/common-main-title-of-a-figure-panel-compiled-with-parmfrow
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}


# sort them in the same way: sort by one column; they are already sorted in coefs_all_cont
# coefs_all_cont2 <- coefs_all_cont[order(coefs_all_cont$Freq_All),]
pdf(file="D:/user/vogelj/Group_project/Output/Plots/barplot_variables_lambda1se_separate_continents.pdf",
    width=16,height=8)
# x11(width=16,height=8)
par(mfrow=c(1,4),mar=c(c(5, 4, 1, 0.5)),oma=c(0,4,2,0))
# title("Number of grid points, where variable is included in the model")
# text(x=0,y=0,"Number of grid points, where variable is included in the model")
barplot(t(coefs_all_cont_mat),horiz=T,las=1,col=c("brown3","DarkOrange2","goldenrod3","burlywood1"),
        xlab="",cex.names=1,font=2,
        legend.text=c("Africa","Europe","North America","Asia"),args.legend =list(x= "bottomright",cex=1.5,text.font=2),main="")
mtext("All continents",side=3,font=2,line=-0.8, cex=1.2)
mtext("a)",side=3,font=2,line=-0.2, adj=0.1,cex=1.6)
barplot(coefs_all_cont$Freq_No_Am,horiz=T,las=1,col="LightCyan3",main="",
        xlab="",cex.names=0.6)
mtext("North America",side=3,font=2,line=-0.8, cex=1.2)
mtext("b)",side=3,font=2,line=-0.2, adj=0.1,cex=1.6)
barplot(coefs_all_cont$Freq_Eur,horiz=T,las=1,col="LightCyan3",main="",
        xlab="",cex.names=0.6)
mtext("Europe",side=3,font=2,line=-0.8, cex=1.3)
mtext("c)",side=3,font=2,line=-0.2, adj=0.1,cex=1.6)
text(line2user(line=mean(par('mar')[c(2, 3)]), side=2), 
     line2user(line=3.5, side=1), 'Number of grid points, where variable is included in the model', xpd=NA, cex=1.5, font=2)
barplot(coefs_all_cont$Freq_Asia,horiz=T,las=1,col="LightCyan3",main="",
        xlab="",cex.names=0.6)
mtext("Asia",side=3,font=2,line=-0.8, cex=1.2)
mtext("d)",side=3,font=2,line=-0.2, adj=0.1,cex=1.6)
dev.off()

