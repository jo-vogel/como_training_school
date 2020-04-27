# Author Johannes Vogel

# subfile for regions_barplot.r
# Missing_points is a small script to identify, which points are not extracted correctly using the continental shapefiles. 

final_pix_num <- 889

all_num <- c(loc_afr_pixels_num,loc_eur_pixels_num,loc_no_am_pixels_num,loc_asia_pixels_num)
miss_points <- !(c(1:final_pix_num) %in% all_num) # five points are missing from the overlay with the shapefiles
which(miss_points)
coord_assigned[loc_pix[miss_points],] 
# 33, 372, 642 are in asia, 470 is in north america, 529 is in Europe (see map below)
coord_subset[miss_points,]

afr_loc_points <- (c(1:final_pix_num) %in% loc_afr_pixels_num)

# does not work: they are inverted
# missing_points <- xyFromCell(loc_ras,loc_pix[miss_points])
# afr_points <- xyFromCell(loc_ras2,loc_pix[afr_loc_points])
# all_points <- xyFromCell(loc_ras,loc_pix)

# any_points <- xyFromCell(loc_ras,sample(1:24000,50))


plot(loc_ras)
# points(missing_points,pch=3,col='blue') # does not work: they are inverted
# points(afr_points,pch=3,col='green')
# points(all_points,pch=3,col='brown')
points(coord_subset[miss_points,],col="red")
# points(coord_subset[afr_loc_points,],col="green")
# points(any_points,pch=3,col='red')
plot(africa,add=T)
plot(no_am,add=T)
plot(europe,add=T)
plot(asia,add=T)