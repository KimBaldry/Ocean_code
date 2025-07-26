library(dplyr)
library(raadfiles)
library(raadtools)

add_wind = function(data){
  coords1 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d"), id = 1:nrow(data))
  coords2 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d") - 18*3600, id = 1:nrow(data))
  coords3 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600, id = 1:nrow(data))
  coords4 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600, id = 1:nrow(data))
  coords5 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600, id = 1:nrow(data))
  coords = rbind(coords1, coords2, coords3, coords4, coords5)
  ws = raadtools::extract(readwind, coords[is.finite(coords$Lon),],magonly = TRUE)
  coords$ws = ws
  ws = stats::aggregate(coords[,4:5],by = list(coords[,4]), FUN = function(x){mean(x, na.rm = T)})
  data$ws[is.finite(data$LON)] = ws[,3]
  data
}


add_wind_storm = function(data){
  coords1 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d"), id = 1:nrow(data))
  coords2 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d") - 6*3600, id = 1:nrow(data))
  coords3 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600, id = 1:nrow(data))
  coords4 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d") - 18*3600, id = 1:nrow(data))
  coords5 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600, id = 1:nrow(data))
  coords6 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600, id = 1:nrow(data))
  coords7 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600, id = 1:nrow(data))
  coords8 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600, id = 1:nrow(data))
  coords9 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600, id = 1:nrow(data))

  coords = rbind(coords1, coords2, coords3, coords4, coords5,coords6, coords7, coords8, coords9)
  ws = raadtools::extract(readwind, coords[is.finite(coords$Lon),],magonly = TRUE)
  coords$ws = ws
  ws = stats::aggregate(coords[,4:5],by = list(coords[,4]), FUN = function(x){any(x > 10, na.rm = T)})
  data$ws_storm[is.finite(data$LON)] = ws[,3]
  data
}

add_wind_stormT = function(data){
  coords1 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d"), id = 1:nrow(data))
  coords2 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600, id = 1:nrow(data))
  coords3 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600, id = 1:nrow(data))
  coords4 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-")  , format = "%Y-%m-%d")- 18*3600, id = 1:nrow(data))
  coords5 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600, id = 1:nrow(data))
  coords6 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600, id = 1:nrow(data))
  coords7 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d") - 12*3600- 24*3600, id = 1:nrow(data))
  coords8 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600, id = 1:nrow(data))
  coords9 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600, id = 1:nrow(data))

  coords10 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d") - 6*3600- 24*3600*2, id = 1:nrow(data))
  coords11 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*2, id = 1:nrow(data))
  coords12 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*2, id = 1:nrow(data))
  coords13 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*2, id = 1:nrow(data))
  coords14 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*3, id = 1:nrow(data))
  coords15 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*3, id = 1:nrow(data))
  coords16 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*3, id = 1:nrow(data))
  coords17 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*3, id = 1:nrow(data))
  coords18 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*4, id = 1:nrow(data))
  coords19 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*4, id = 1:nrow(data))
  coords20 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"),  format = "%Y-%m-%d")- 18*3600 - 24*3600*4, id = 1:nrow(data))
  coords21 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*4, id = 1:nrow(data))
  coords22 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*5, id = 1:nrow(data))
  coords23 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*5, id = 1:nrow(data))
  coords24 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*5, id = 1:nrow(data))
  coords25 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d") - 24*3600- 24*3600*5, id = 1:nrow(data))
  coords26 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*6, id = 1:nrow(data))
  coords27 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*6, id = 1:nrow(data))
  coords28 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*6, id = 1:nrow(data))
  coords29 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*6, id = 1:nrow(data))
  coords30 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*7, id = 1:nrow(data))
  coords31 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*7, id = 1:nrow(data))
  coords32 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*7, id = 1:nrow(data))
  coords33 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*7, id = 1:nrow(data))
  coords30 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*8, id = 1:nrow(data))
  coords31 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*8, id = 1:nrow(data))
  coords32 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*8, id = 1:nrow(data))
  coords33 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*8, id = 1:nrow(data))
  coords30 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*9, id = 1:nrow(data))
  coords31 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*9, id = 1:nrow(data))
  coords32 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*9, id = 1:nrow(data))
  coords33 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*9, id = 1:nrow(data))
  coords30 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 6*3600- 24*3600*10, id = 1:nrow(data))
  coords31 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 12*3600- 24*3600*10, id = 1:nrow(data))
  coords32 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 18*3600 - 24*3600*10, id = 1:nrow(data))
  coords33 = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-") , format = "%Y-%m-%d")- 24*3600- 24*3600*10, id = 1:nrow(data))
  
  coords = rbind(coords1, coords2, coords3, coords4, coords5,coords6, coords7, coords8, coords9, coords10,coords11, coords12, coords13, coords14, coords15,coords16, coords17, coords18, coords19, coords20,coords21, coords22, coords23, coords24, coords25, coords26, coords27, coords28, coords29, coords30, coords31, coords32, coords33)
  ws = raadtools::extract(readwind, coords[is.finite(coords$Lon),],magonly = TRUE)
  coords$ws = ws
  ws = stats::aggregate(coords[,4:5],by = list(coords[,4]), FUN = function(x){((which(x > 10)[1]-1)*6)/24})
  data$ws_tstorm[is.finite(data$LON)] = ws[,3]
  data
}


