library(sp)
library(rgeos)
library(rgdal)
library(raster)

Path.GIS <- "/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/MeshModeling/Data/GISData/"
Path.city <- paste0(Path.GIS,"City_Boundary")
city.shp <- readOGR(Path.city,"City_Boundary") 
cellsize.x <- 660
cellsize.y <- 660
X_range <- city.shp@bbox[1,]
Y_range <- city.shp@bbox[2,]
grd.full <- expand.grid(list(X_COORD=seq(X_range[1],X_range[2],by=cellsize.x),
                             Y_COORD=seq(Y_range[1],Y_range[2],by=cellsize.y)))
coordinates(grd.full) = ~X_COORD+Y_COORD # convert to SpatialPoints

prj <- paste("+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 +k=0.9999749999999999",
             "+x_0=300000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
proj4string(grd.full) <- prj

# rasterize the city spatial polygon to get a grid template
grd.full <- SpatialPixels(grd.full)
r <- raster(ncol=grd.full@grid@cells.dim[1],nrow=grd.full@grid@cells.dim[2],
            xmn=grd.full@bbox[1,1],xmx=grd.full@bbox[1,2],ymn=grd.full@bbox[2,1],ymx=grd.full@bbox[2,2],crs=CRS(prj))

city.raster <- rasterize(city.shp,r,0)
city.df_full <- as.data.frame(city.raster,xy=TRUE)
city.df_full <- city.df_full[,1:2]
names(city.df_full) <- c("X_COORD","Y_COORD")
RegGrd.full <- city.df_full

coordinates(city.df_full) <- c("X_COORD", "Y_COORD") 
proj4string(city.df_full) <- prj
BoundedOverFullGrd <- over(city.df_full, city.shp)
isInCity <- !is.na(BoundedOverFullGrd$OBJECTID)
RegGrd <- RegGrd.full[isInCity,]

# KDEgrd.full is identical to RegGrd.full except the Y_COORD is in reverse order 
# KDEgrd.full is used for KDE and mesh generation (to be compatiable with the function output coordinate order)
KDEgrd.full <- expand.grid(list(X_COORD=seq(X_range[1],X_range[2],by=cellsize.x),
                                Y_COORD=seq(Y_range[1],Y_range[2],by=cellsize.y)))
KDEgrd.sp <- KDEgrd.full 
coordinates(KDEgrd.sp) <- c("X_COORD", "Y_COORD") 
proj4string(KDEgrd.sp) <- prj
BoundedOverFullGrd2 <- over(KDEgrd.sp, city.shp)
isInCity2 <- !is.na(BoundedOverFullGrd2$OBJECTID)
KDEgrd <- KDEgrd.full[isInCity2,] 

# the order relationship between RegGrd and KDEgrd
ord <- rep(NA,nrow(KDEgrd))
for (i in 1:nrow(KDEgrd)){
  ord[i] <- which(RegGrd[,1]==KDEgrd[i,1]&RegGrd[,2]==KDEgrd[i,2])
}

ord.full <- rep(NA,nrow(KDEgrd.full))
for (i in 1:nrow(KDEgrd.full)){
  ord.full[i] <- which(RegGrd.full[,1]==KDEgrd.full[i,1]&RegGrd.full[,2]==KDEgrd.full[i,2])
}