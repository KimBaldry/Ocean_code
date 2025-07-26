library(dplyr)
library(raadfiles)
library(raadtools)

add_chl = function(data,res = "monthly"){
  coords = data.frame("Lon" = data$LON, "Lat" = data$LAT, "date"=as.POSIXct(paste(data$YYYY,data$MM,data$DD,sep = "-"), format = "%Y-%m-%d"))
  if(res == "daily"){cl = raadtools::extract(readCHL_daily, coords[is.finite(coords$Lon),])
  data$sat_chl_day[is.finite(data$LON)] = cl}
  if(res == "monthly"){  cl = raadtools::extract(readCHL_monthly, coords[is.finite(coords$Lon),])
  data$sat_chl_month[is.finite(data$LON)] = cl}
  if(res == "weekly"){  cl = raadtools::extract(readCHL_weekly, coords[is.finite(coords$Lon),])
  data$sat_chl_week[is.finite(data$LON)] = cl}
  data
}
