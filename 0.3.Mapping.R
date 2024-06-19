### Mapping of Antarctic sites 

library(ggOceanMapsData)
library(ggOceanMaps)
library(sf)
library(ggspatial)


pt <- read.csv("sitemetadata.csv")

pdf("mapping/OverviewV2.pdf",height = 8,width = 6)
basemap(limits = c(-80, -46, -80, -58), bathy.style = "rcb",rotate = TRUE,glaciers = TRUE)+
  #basemap(limits = c(-25, -12, 63, 67),rotate = TRUE)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt[1:3,], color ="black",cex=5,crs = 4326)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt[1:3,], color = c("#009e73","#D55E00","#0072B2"),cex=4,crs = 4326)+
  #ggspatial::geom_spatial_text(aes(x = Lon, y = Lat+0.04,label=CoreID),data = metadata, color = "darkred",cex=2)+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_scale(location = "bl")+ 
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()

pdf("mapping/OverviewV2.zoom.pdf",height = 8,width = 6)
basemap(limits = c(-70, -56, -70, -60), bathy.style = "rcb",rotate = TRUE,glaciers = TRUE)+
  #basemap(limits = c(-25, -12, 63, 67),rotate = TRUE)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt[1:3,], color ="black",cex=3,crs = 4326)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt[1:3,], color = c("#009e73","#D55E00","#0072B2"),cex=2,crs = 4326)+
  #ggspatial::geom_spatial_text(aes(x = Lon, y = Lat+0.04,label=CoreID),data = metadata, color = "darkred",cex=2)+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_scale(location = "bl")+ 
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()


pdf("mapping/OverviewV2.2.zoom.pdf",height = 8,width = 6)
basemap(limits = c(-70, -56, -70, -60), bathy.style = "rcb",rotate = TRUE,glaciers = TRUE)+
  #basemap(limits = c(-25, -12, 63, 67),rotate = TRUE)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt[1:3,], color ="black",cex=5,crs = 4326)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt[1:3,], color = c("#009e73","#D55E00","#0072B2"),cex=4,crs = 4326)+
  #ggspatial::geom_spatial_text(aes(x = Lon, y = Lat+0.04,label=CoreID),data = metadata, color = "darkred",cex=2)+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.9, "cm"), width = unit(0.9, "cm"))+
                           annotation_scale(location = "bl")
dev.off()

pdf("mapping/OverviewV2.zoom.newcores.pdf",height = 8,width = 6)
basemap(limits = c(-70, -56, -70, -60), bathy.style = "rcb",rotate = TRUE,glaciers = TRUE)+
  #basemap(limits = c(-25, -12, 63, 67),rotate = TRUE)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt, color ="black",cex=7,crs = 4326)+
  ggspatial::geom_spatial_point(aes(x = Lon, y = Lat),data = pt, color = c("#009e73","#D55E00","#0072B2","darkorchid4","darkorchid4"),cex=6,crs = 4326)+
  #ggspatial::geom_spatial_text(aes(x = Lon, y = Lat+0.04,label=CoreID),data = metadata, color = "darkred",cex=2)+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_scale(location = "bl")+ 
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.9, "cm"), width = unit(0.9, "cm"))
dev.off()





pdf("mapping/Antarctica.pdf",height = 4,width = 4)
dt <- data.frame(lon = c(-80, -46, -46, -80), lat = c(-80, -80, -58, -58))
basemap(limits =-56,gla.col = "black",shapefiles = "Antarctic",glaciers = T,grid.size = 0) +
  ggspatial::geom_spatial_polygon(
    data = dt, 
    aes(x = lon, y = lat), fill = NA, color = "darkred",lwd=2)+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
dev.off()


pdf("mapping/Antarctica2.pdf",height = 4,width = 4)
dt <- data.frame(lon = c(-70, -56, -56, -70), lat = c(-70, -70, -60, -60))
basemap(limits =-56,gla.col = "black",shapefiles = "Antarctic",glaciers = T,grid.size = 0) +
  ggspatial::geom_spatial_polygon(
    data = dt, 
    aes(x = lon, y = lat), fill = NA, color = "darkred",lwd=2)+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

dev.off()




pdf("mapping/AntarcticaOLD.pdf",height = 4,width = 4)
basemap(limits =-60,gla.col = "black",shapefiles = "Antarctic",glaciers = T,grid.size = 0)
dev.off()



  
  
  
#archive basement 



# libraries
library(ncdf4)
library(raster)
set.seed("123456")



# Data in 
metadata <- read.csv("metadata.csv")

## A little map 

#load GEBCO_2022 netcdf downloaded on 110822 - https://www.gebco.net/data_and_products/gridded_bathymetry_data/
gebco <- raster("mapping/gebco_2022_sub_ice_topo/GEBCO_2022_sub_ice_topo.nc")
# Create extent (our map area)
ant.ex <- extent(-80, -46, -75, -56)
# Create a crop of the bathymetric data
ant.crop <- crop(gebco, ant.ex)
#Save file 
#saveRDS(gebco.crop,file="mapBuilding/GalapBathy.rds")

#here is a crs centred around the antartic pole  
newproj1 <- crs("+init=epsg:3031")
# the same projection but now the prime meridian is roughly in line with the Ant. peninsula 
newproj2 <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-64 +x_0=0 +y_0=0 +datum=WGS84 +units=m"

#let's project this bad boi
ant.crop.p <- projectRaster(ant.crop,crs=newproj2)


#function to calculate break points for colours from https://www.benjaminbell.co.uk/2019/08/bathymetric-maps-in-r-colour-palettes.html
## x = raster, b1 & b2 = number of divisions for each sequence, r1 & r2 = rounding value
colbr <- function(x, b1=50, b2=50, r1=-2, r2=-2) {
  # Min/max values of the raster (x)
  mi <- cellStats(x, stat="min")-100
  ma <- cellStats(x, stat="max")+100
  # Create sequences, but only use unique numbers
  s1 <- unique(round(seq(mi, 0, 0-mi/b1),r1))
  s2 <- unique(round(seq(0, ma, ma/b2),r2))
  # Combine sequence for our break points, removing duplicate 0
  s3 <- c(s1, s2[-1])
  # Create a list with the outputs
  # [[1]] = length of the first sequence minus 1 (water)
  # [[2]] = length of the second sequence minus 1 (land)
  # [[3]] = The break points
  x <- list(length(s1)-1, length(s2)-1, s3)
}

ant.br <- colbr(ant.crop.p,b2=1)
# Colour palette
blue.col <- colorRampPalette(c("darkblue", "lightblue"))


#lets project the sites, 
lat <- -62.623554 
lon <- -60.291926


pt <- read.csv("sitemetadata.csv")
coordinates(pt) <- ~Lat+Lon
crs(pt) <- crs(gebco)
pt_projected <- spTransform(pt,newproj2)

plot(ant.crop.p, col=c(blue.col(ant.br[[1]]), grey.colors(ant.br[[2]])), breaks=ant.br[[3]],axes = FALSE,box=F,legend = FALSE)
points(pt_projected,col="red",pch=16)
points(pt,col="red",pch=16)

points(lat,lon,col="red",pch=16)
#First plot

pdf(width = 8,height=6.5,file="mapping/test.pdf")

plot(ant.crop.p, col=c(blue.col(ant.br[[1]]), grey.colors(ant.br[[2]])), breaks=ant.br[[3]],axes = FALSE,box=F,legend = FALSE)
points(pt_projected,pch=16, col='black',cex=1)
points(pt_projected,pch=16, col='white',cex=0.5)

dev.off()







