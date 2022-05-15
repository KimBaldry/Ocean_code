# Title: create coastline data
# Author: K. Baldry - IMAS/UTAS
# Created: 11/05/2020
#
# Convert shapefiles into a r data object and save this object in the package data (reproj_map_tools)


# import coastlines
coast_path = file.path("E:/PhD_Data","Coastlines_USGS_gen01")
coastlines = readOGR(dsn = coast_path,layer = "Coastlines_USGS_gen01")
pkg_dir = "C:/Users/kabaldry/OneDrive - University of Tasmania/Documents/PhD/Argo_code_development/reproj_map_tools"

save(coastlines, file=file.path(pkg_dir,"data","Coastlines_USGS_gen01.rda"))
