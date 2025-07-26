library(raster)
library(marmap)

add_bathy = function(data){
### Bathymetry
# read combined data
NOAA.bathy = getNOAA.bathy(-180,180,-90,-30)
data$LON[which(data$LON < -179.96)] = -179.96
data$LON[which(data$LON > 179.96)] = 179.96
data$BATHY_DEPTH_m[is.finite(data$LON)] = get.depth(NOAA.bathy,data$LON[is.finite(data$LON) ] ,data$LAT[is.finite(data$LON)] ,locator = FALSE)$depth
bathy_rast = as.raster(NOAA.bathy)
cont_500 = rasterToContour(bathy_rast,levels = -500)
# data(wrld_simpl, package = "maptools")
# wrld_subset <- crop(wrld_simpl, extent(-180, 180, -90, -30)) 
points = SpatialPoints(data[is.finite(data$LON),c("LON","LAT")])
dist.mat <- geosphere::dist2Line(p = points, line = cont_500)
data$coast_dist = NA
data$coast_dist[is.finite(data_CTD$LON)] = dist.mat[,1]/1000
data$coast_dist[data$BATHY_DEPTH_m > - 500] = 0
data
}
