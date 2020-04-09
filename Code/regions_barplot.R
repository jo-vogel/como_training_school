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