# Regional comparisons

# Author: Johannes Vogel

message("run with Model_comparison.r")


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

coord_all_var <- cbind(coord_all,rep(NA,320*76),rep(NA,320*76),rep(NA,320*76))
for (i in seq_along(work_pix)){
  coord_all_var[loc_pix[i],3] <- speci_adj[i]
  coord_all_var[loc_pix[i],4] <- csi_adj[i]
  coord_all_var[loc_pix[i],5] <- coeff_kep[i]
}

# Create continent polygons
border <- readOGR('D:/user/vogelj/Data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')	
continents <- readOGR("D:/user/vogelj/Data/continent_shapefile/continent.shp") # from https://www.arcgis.com/home/item.html?id=5cf4f223c4a642eb9aa7ae1216a04372
africa <- subset(continents,subset=continents@data[["CONTINENT"]]=="Africa")
europe <- subset(continents,subset=continents@data[["CONTINENT"]]=="Europe")
no_am <- subset(continents,subset=continents@data[["CONTINENT"]]=="North America")
so_am <- subset(continents,subset=continents@data[["CONTINENT"]]=="South America")
oceania <- subset(continents,subset=continents@data[["CONTINENT"]]=="Oceania")
asia <- subset(continents,subset=continents@data[["CONTINENT"]]=="Asia")
australia <- subset(continents,subset=continents@data[["CONTINENT"]]=="Australia")
antarctica <- subset(continents,subset=continents@data[["CONTINENT"]]=="Antarctica")


# Specificity
spec_mat <- matrix(as.numeric(coord_all_var[,3]),nrow=320,ncol=76)
spec_ras <- raster(t(spec_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all), ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))
spec_afr_pixels <- extract(spec_ras,africa)
spec_eur_pixels <- extract(spec_ras,europe)
spec_no_am_pixels <- extract(spec_ras,no_am)
spec_asia_pixels <- extract(spec_ras,asia)
# so_am_pixels <- extract(spec_ras,so_am)
# oce_pixels <- extract(spec_ras,oceania)
# aus_pixels <- extract(spec_ras,australia)
# ant_pixels <- extract(spec_ras,antarctica)

sum(!is.na(spec_eur_pixels[[1]])) # 233 pixels
sum(!is.na(spec_no_am_pixels[[1]])) # 419 pixels
sum(!is.na(spec_afr_pixels[[1]])) # 45 pixels
sum(!is.na(spec_asia_pixels[[1]])) # 263 pixels
# 45+233+419+263=960; 965-960: 5 pixel are missing

boxplot(spec_eur_pixels[[1]],spec_no_am_pixels[[1]],spec_afr_pixels[[1]], spec_asia_pixels[[1]],speci_adj,names=c("Europe (233)","No. America (419)","Africa (45)","Asia (263)","World (965)"), main="Adjusted specificity for glinternet wit extr. indicators",sub="Number of pixels in brackets")


# CSI
csi_mat <- matrix(as.numeric(coord_all_var[,4]),nrow=320,ncol=76)
csi_ras <- raster(t(csi_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all), ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))

csi_afr_pixels <- extract(csi_ras,africa)
csi_eur_pixels <- extract(csi_ras,europe)
csi_no_am_pixels <- extract(csi_ras,no_am)
csi_asia_pixels <- extract(csi_ras,asia)
boxplot(csi_eur_pixels[[1]],csi_no_am_pixels[[1]],csi_afr_pixels[[1]], csi_asia_pixels[[1]],csi_adj,names=c("Europe (233)","No. America (419)","Africa (45)","Asia (263)","World (965)"), main="Adjusted CSI for glinternet wit extr. indicators",sub="Number of pixels in brackets")


# Number of variables coeff_kep
coeff_kep_mat <- matrix(as.numeric(coord_all_var[,5]),nrow=320,ncol=76)
coeff_kep_ras <- raster(t(coeff_kep_mat[,76:1]), xmn=min(lon_all), xmx=max(lon_all), ymn=min(lat_all), ymx=max(lat_all), crs=CRS(projection(border)))

coeff_kep_afr_pixels <- extract(coeff_kep_ras,africa)
coeff_kep_eur_pixels <- extract(coeff_kep_ras,europe)
coeff_kep_no_am_pixels <- extract(coeff_kep_ras,no_am)
coeff_kep_asia_pixels <- extract(coeff_kep_ras,asia)
boxplot(coeff_kep_eur_pixels[[1]],coeff_kep_no_am_pixels[[1]],coeff_kep_afr_pixels[[1]], coeff_kep_asia_pixels[[1]],coeff_kep,names=c("Europe (233)","No. America (419)","Africa (45)","Asia (263)","World (965)"), main="Number of coefficients for glinternet wit extr. indicators",sub="Number of pixels in brackets")
