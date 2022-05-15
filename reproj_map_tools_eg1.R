### L. Bach - plot k, wind speed, sea-ice climatology (south 60S) 
### and SOCAT OF and K ###
# Author: K. Baldry
# Date of creation: 30012020
#
#
#
# 
# colour scheme: #8c510a #d8b365 #f6e8c3 #c7eae5 #5ab4ac #01665e
library(ggplot2)
library(mapproj)
library(raster)
library(ggpubr)
library(sf)
library(rgdal)
library(dplyr)
library(SOmap)
library(ggnewscale)
library(ggmap)
library(scales)
library(viridis)
library(marmap)
library(sp)
library(maptools)
library(rgeos)
library(polyclip)
library(GISTools)
library(animation)
library(gridExtra)
library(gam)
library(data.table)

mainDir = "C:/Users/kabaldry/Documents/PhD/Collaborations/L. Bach"
# read data
path = "C:/Users/kabaldry/Documents/FluxEngine-master/SOCATv2"
data_CTD = fread(file.path(path,"FluxEngine_out_30012020.tsv"), header = T, stringsAsFactors = F)
data_SOCAT = data_CTD %>% filter(OF != -999)
nc_path = "C:/Users/kabaldry/OneDrive - University of Tasmania/PhD/Collaborations/L. Bach/Data/K_calc"
SO_proj = "+proj=laea +lat_0=-90 +lon_0=147 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
# set clip latitude here
clip_lat = -60

# create extent references for two projections (cartesian and stereographic)
extent_p = data.frame("x"=c(0,90,180,-90), "y" = rep(clip_lat,4))
coordinates(extent_p) <- ~x+y
projection(extent_p) <- "+proj=longlat +datum=WGS84"



# plot data tailored to data
plot_data = map_plot_proj(extent_p,proj_str = SO_proj,resolution = 0.5,circle = T)


names = c("Dec","Jan","Feb","March")
i = 1
plotlist = list()
of_socat_l = list()
wind_l = list()
k_sd_l = list()
k_av_l = list()
ok_socat_l = list()
density_socat_l = list()
si_l = list()
sst_l = list()

for(m in c(12,1:3))
{
  nc_file = file.path(nc_path,paste("All_data_for_k_month_",m,".nc",sep = ""))
  ## SOCAT  density
  points = data_SOCAT[which(data_SOCAT$mon == m),c("lon","lat","OF","OK3")]
  colnames(points) = c("long","lat","OF","OK3")
  coordinates(points) <- ~long+lat #check right way round
  projection(points) <- "+proj=longlat +datum=WGS84"
  points_ster = spTransform(points,SOmap()$projection)
  rast <- raster(x = extent_p2, nrows = 60, ncols = 60)
  gridded_counts_profile = rasterize(points_ster, rast, field = 1, fun = "count", background = 0)
  gridded_counts_profile <- as(gridded_counts_profile, "SpatialPointsDataFrame")
  gridded_counts_profile <-as.data.frame(gridded_counts_profile)
  # background to NA
  gridded_counts_profile$layer[which(gridded_counts_profile$layer == 0)] = NA
  
  ## BGCArgo density
  
  ## Chl + SE
  
  ## SI
  
  ## OF_SOCAT
  socat_of = rasterize(points_ster, rast, field = "OF", fun = mean, background = 0)
  socat_of <- as(socat_of, "SpatialPointsDataFrame")
  socat_of <-as.data.frame( socat_of)
  # background to NA
  socat_of$layer[which(socat_of$layer == 0)] = NA
  
  ## OF_BGCArgo + SE
  
  ## K
  socat_ok = rasterize(points_ster, rast, field = "OK3", fun = mean, background = 0)
  socat_ok <- as(socat_ok, "SpatialPointsDataFrame")
  socat_ok <-as.data.frame( socat_ok)
  # background to NA
  socat_ok$layer[which(socat_ok$layer == 0)] = NA
  
  
  ## WS
  ws = stack(nc_file,varname = "av_ws")
  av_ws = calc(ws, fun = function(x){mean(x,na.rm = T)})
  av_ws = change_lon(av_ws)
  ws_ster = projectRaster(av_ws,crs = SOmap()$projection)
  ws_ster <- as(ws_ster, "SpatialPointsDataFrame")
  ws_ster <-as.data.frame(ws_ster)
  ## SI
  si = stack(nc_file,varname = "sea_ice")
  av_si = calc(si, fun = function(x){mean(x,na.rm = T)})
  av_si = change_lon(av_si)
  si_ster = projectRaster(av_si,crs = SOmap()$projection)
  si_ster <- as(si_ster, "SpatialPointsDataFrame")
  si_ster <-as.data.frame(si_ster)
  ## SST - convert to K
  sstC = calc(stack(nc_file,varname ="SST"), fun = function(x){x - 0.17})
  av_sst = calc(sstC, fun = function(x){mean(x,na.rm = T)})
  av_sst = change_lon(av_sst)
  sst_ster = projectRaster(av_sst,crs = SOmap()$projection)
  sst_ster <- as(sst_ster, "SpatialPointsDataFrame")
  sst_ster <-as.data.frame(sst_ster)
  
  
  
  sst = stack(nc_file,varname ="SST") - 0.17
  
  sc = calc(sst, fun = function(x){2116.8 + (-136.25*x) + (4.7353 * (x^2)) + (-0.092307  *(x^3)) + (0.0007555*(x^4))}, forceapply = T, forcefun = T)
  ## K
  k1 =  calc(ws, fun = function(x){0.251 * (x^2)}, forceapply = T, forcefun = T)
  sc2 = calc(sc,fun = function(x){sqrt(660.0 * (x^-1))}, forceapply = T, forcefun = T)
  
  for(sp in 1:nlayers(k1)){
    a =  raster(k1, layer=sp)
    b = raster(sc2, layer = sp)
    c = calc(raster(si, layer = sp), fun = function(x){(100-x) * 0.01}, forceapply = T, forcefun = T)
    
    k_sub = overlay(a,b,fun = function(x,y){x * y}, forceapply = T, forcefun = T)
    k_sub2 = overlay(k_sub, c,fun = function(x,y){x * y }, forceapply = T, forcefun = T)
    if(sp>1){k = stack(k,k_sub2)}else{k = k_sub2}
  }
  rm(a,b,c,k_sub)
  
  
  av_k = overlay(k, fun = mean, na.rm = T)
  sd_k = overlay(k, fun = sd, na.rm =T) 
  av_k = change_lon(av_k)
  sd_k = change_lon(sd_k)
  k_sd = projectRaster(sd_k,crs = SOmap()$projection)
  k_sd <- as(k_sd, "SpatialPointsDataFrame")
  k_sd <-as.data.frame(k_sd)
  
  k_av = projectRaster(av_k,crs = SOmap()$projection)
  k_av <- as(k_av, "SpatialPointsDataFrame")
  k_av <-as.data.frame(k_av)
  #
  if(m == 12){
    lim_dens = quantile(gridded_counts_profile$layer, c(0,0.9), na.rm = T)
    lim_of_socat = quantile(socat_of$layer, c(0,0.9), na.rm = T)
    lim_ok_socat = quantile(socat_ok$layer, c(0,0.9), na.rm = T)
    lim_k = quantile(k_av$layer, c(0,0.9), na.rm = T)
    lim_sd =  quantile(k_sd$layer, c(0,0.9), na.rm = T)
    lim_w = quantile(ws_ster$layer, c(0,0.99), na.rm = T)
    lim_si = quantile(si_ster$layer, c(0,1), na.rm = T)
    lim_sst = quantile(sst_ster$layer, c(0,0.9), na.rm = T)
  }
  density_socat =  base_map_gg +  
    geom_raster(data = gridded_counts_profile, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_dens, na.value = "#FFFFFF00",oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(density_socat )
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  density_socat_l[[i]] = plot.lab   
  
  of_socat = base_map_gg +  
    geom_raster(data = socat_of, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_of_socat, na.value = "#FFFFFF00", oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(of_socat )
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  of_socat_l[[i]] = plot.lab                        
  
  ok_socat = base_map_gg +  
    geom_raster(data = socat_ok, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_ok_socat, na.value = "#FFFFFF00", oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(ok_socat )
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  ok_socat_l[[i]] = plot.lab  
  
  
  k_average = base_map_gg +  
    geom_raster(data = k_av, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_k, na.value = "#FFFFFF00", oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(k_average )
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  k_av_l[[i]] = plot.lab  
  
  k_stdev = base_map_gg +  
    geom_raster(data = k_sd, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_sd, na.value = "#FFFFFF00", oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(k_stdev )
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  k_sd_l[[i]] = plot.lab  
  
  wind = base_map_gg +  
    geom_raster(data = ws_ster, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_w, na.value = "#FFFFFF00", oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(wind )
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  wind_l[[i]] = plot.lab  
  
  SI = base_map_gg +  
    geom_raster(data = si_ster, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_si, na.value = "#FFFFFF00", oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(SI )
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  si_l[[i]] = plot.lab  
  
  SST = base_map_gg +  
    geom_raster(data = sst_ster, aes(x=x, y=y, fill = layer)) + 
    scale_fill_viridis(name = "No. measurements", limits = lim_sst, na.value = "#FFFFFF00", oob = squish) +
    theme(legend.position = "none", plot.title = element_text(size =12)) +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(-1000,0), col = "lightgrey") +
    stat_contour(data = bathySOdf,aes(x=x,y=y,z=layer),breaks = c(0), col = "black", fill = "white", geom = "polygon") +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
    geom_polygon(data = fortify(mask),aes(x=long,y=lat), fill = "white", col = "NA")+
    geom_path(data = trim_points, aes(x=x,y=y))
  
  plot.lab = ggplotGrob(SST)
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  sst_l[[i]] = plot.lab  
  
  print(i)
  
  i = i+1
}

# plot obs density
dummy = ggplot(data = gridded_counts_profile, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "No. measurements", limits = lim_dens, oob = squish) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
density_socat_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(NA,i,i,NA))

jpeg(file.path(mainDir, "obs_density.jpg"),width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs =density_socat_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()

# plot OF
dummy_of = ggplot(data = socat_of, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "Flux [g C m-2 day-1]", limits = lim_of_socat, oob = squish) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy_of)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
of_socat_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(i,i,i,i))


jpeg(file.path(mainDir, "of_plots.jpg"),width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs = of_socat_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()

## plot SOCAT K
dummy = ggplot(data = socat_ok, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "k [m/s]", limits = lim_ok_socat, oob = squish) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
ok_socat_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(NA,i,i,NA))

jpeg(file.path(mainDir, "ok_socat_plots.jpg"), width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs =ok_socat_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()

## plot K climatology
dummy = ggplot(data = k_av, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "k [cm/h]", limits = lim_k, oob = squish) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
k_av_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(NA,i,i,NA))

jpeg(file.path(mainDir, "k_clim_plots.jpg"), width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs = k_av_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()

dummy = ggplot(data = k_sd, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "stdev k [cm/h]", limits = lim_sd, oob = squish) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
k_sd_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(NA,i,i,NA))

jpeg(file.path(mainDir, "k_clim_sd_plots.jpg"), width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs = k_sd_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()

# ws
dummy = ggplot(data = ws_ster, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "wind speed [m/s]", limits = lim_w, oob = squish ) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
wind_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(NA,i,i,NA))

jpeg(file.path(mainDir, "wind_plots.jpg"),width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs = wind_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()

# sst
dummy = ggplot(data = sst_ster, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "wind speed [m/s]", limits = lim_sst, oob = squish ) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
sst_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(NA,i,i,NA))

jpeg(file.path(mainDir, "sst_plots.jpg"), width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs = sst_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()


# sea-ice
dummy = ggplot(data = si_ster, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "sea-ice [%]", limits = lim_si, oob = squish ) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
si_l[[i]] = legend
lay = rbind(matrix(c(1:4),nrow = 1, ncol = 4, byrow = T),c(NA,i,i,NA))

jpeg(file.path(mainDir, "si_plots.jpg"), width = 1200, height = 500, units = "px", quality = 200,res = 150)
grid.arrange(grobs = si_l,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()
