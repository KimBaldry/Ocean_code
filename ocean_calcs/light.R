
add_light = function(data){
  coords = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d"))
  pr = raadtools::extract(read_par, coords[is.finite(coords$Lon),])
  data$par[is.finite(data$LON)] = pr
  data
}
