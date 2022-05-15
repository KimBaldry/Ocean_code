#Read in SOM output nc file

mainDir <-"C:/Users/kabaldry/OneDrive - University of Tasmania/Documents/PhD/Active Collaborations/Segment stereographic plots"
in_dir  <- mainDir
f_name  <- "SOM_daily_MSLP_3_4.nc"

library(fpc)
library(tools)
library(ncdf4)
library(clValid)
library(maptools)
library(sf)
library(ggplot2)
#library(tidyverse)  #R package for data science include ggplot
#library(maps)
library(ggcorrplot)   #for correlation between SOM nodes and daily MSLP field
library(viridis)   #color palettes
library(reshape)
library(raster)
library(maps)
library(rgdal)
library(rgeos)
library(polyclip)
library(gridExtra)
pkg_dir = "C:/Users/kabaldry/OneDrive - University of Tasmania/Documents/PhD/Argo_code_development/Hard_plotting"
source(file.path(path,"transform_rast.R"))
# install SOmap
# remotes::install_github("AustralianAntarcticDivision/SOmap")
#read in netcdf

nc  <- nc_open(file.path(mainDir,f_name), readunlim=FALSE)
# what is the resolution of the lat/lon grid?
# the middle longitude is?
midlon = 105
#proj string 
DanisProj = paste("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=",midlon," +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",sep = "")
#DanisProj = paste("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=",midlon," +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",sep = "")

#define variables
SOM_msl <- ncvar_get(nc,varid="node")
lon    <- ncvar_get(nc, "Lon")
lat    <- ncvar_get(nc, "Lat")
SOMextent = extent(c(range(lon), range(lat)))
# plot data tailored to data
plot_lims = range(SOM_msl)
plot_data = map_plot_proj(SOMextent,proj_str = DanisProj,resolution = 0.75)

# create an empty raster and set values in loop
SOMdata = raster(SOMextent, ncol = length(lon), nrow = length(lat))
projection(SOMdata) = "+proj=longlat +datum=WGS84" # add CRS string

plotlist = list()
for(sp in ncvar_get(nc, "nodeN")){
  # put data from SOM1
  SOMdata = setValues(SOMdata,t(SOM_msl[,,sp]))
  # reproject raster (function produces an interpolated grid)
  SOMdata_reproj = transform_rast(rast = SOMdata,proj = DanisProj, plot_extent = extent(plot_data$border))
  # plot
  SOMplot = ggplot(data = SOMdata_reproj, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) +
  geom_polygon(data = plot_data$border, aes(x=x, y = y), col = "black",fill = NA,lwd = 1) + 
  scale_fill_viridis(name = "MSLP", limits = plot_lims)+
  geom_contour(aes(z= layer)) + 
  geom_path(data = plot_data$coastlines, aes(x=x,y=y, group = group),  col = "black") +
  theme(legend.position = "none", legend.key.width = unit(2,"lines"),
        plot.title = element_text(size =12),panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title  = element_blank(), axis.ticks = element_blank()) + 
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))
# list of plots
  plot.lab = ggplotGrob(SOMplot + ggtitle(paste("SOM",sp)))
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  plotlist[[sp]] = plot.lab                        
}
# guide
dummy = ggplot(data = SOMdata_reproj, aes(x=x,y=y)) +geom_raster(aes(fill= layer)) + 
  scale_fill_viridis(name = "MSLP", limits = plot_lims) +
  guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "bottom"))+
  theme(legend.position = "bottom",legend.key.width = unit(3,"line"))

g = ggplotGrob(dummy)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
plotlist[[sp+1]] = legend
lay = rbind(matrix(c(1:sp),nrow = 3, ncol = 4, byrow = T),c(NA,sp+1,sp+1,NA))

jpeg(file.path(mainDir, "SOM_plots.jpg"), width = 1000, height = 500, units = "px", quality = 90)
grid.arrange(grobs = plotlist,layout_matrix = lay, heights = c(rep(2,nrow(lay)-1),1))
dev.off()
