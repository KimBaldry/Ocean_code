library(dplyr)
library(raadfiles)
library(raadtools)

add_ice = function(data){
  coords = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d"))
#  DIE = extract(distance_to_ice, coords[is.finite(coords$Lon),])
  IC = raadtools::extract(readice_daily, coords[is.finite(coords$Lon),])
  data$ice_con[is.finite(data$LON)] = IC
#  data$die[is.finite(data$LON)] = DIE/1000
  #data$die[which(data$ice_con>15)] = -1*data$die[which(data$ice_con>15)]
  data
}
