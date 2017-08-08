require(sp)
require(raster)
require(spatstat)
require(KernSmooth)

generateDenIm <- function(CrimeData,r,grd,bw=NA){
  # Generate density image of the raw crime data via KDE
  
  # subset data
  CrimeData <- subset(CrimeData,select=c("X_COORD","Y_COORD","INC_CNT"))
  
  # KDE bandwidth
  X_range <- range(grd$X_COORD)
  Y_range <- range(grd$Y_COORD)
  # find the optimal kernel
  if (any(is.na(bw))){
    # rasterize the data to an image
    CrimeData.Raster <- rasterize(CrimeData[,c("X_COORD","Y_COORD")], r, CrimeData$INC_CNT, fun=sum)
    rasterMat <- as.matrix(CrimeData.Raster@data)
    rasterMat[is.na(rasterMat)] <- 0 
    
    CrimeData.pp <- ppp(grd$X_COORD,grd$Y_COORD,window=owin(xrange=X_range,yrange=Y_range),marks=as.vector(rasterMat))
    bw <- bw.diggle(CrimeData.pp)
    bw <- c(round(as.numeric(bw)),round(as.numeric(bw)))
  }
  
  # Kernel Smoothing
  M <- r@nrows
  N <- r@ncols
  KDE.df <- data.frame(X_COORD=grd$X_COORD,Y_COORD=grd$Y_COORD,VALUE=rep(NA,N*M))
  kernSm <- bkde2D(data.matrix(CrimeData[,c("X_COORD","Y_COORD")]), bandwidth=bw, 
                   gridsize=c(N, M), range.x=list(X_range,Y_range))
  pixelIm <- matrix(as.vector(kernSm$fhat),N,M) 
  
  return(pixelIm)
}


generatePixGrd <- function(NumPix,RegGrd,r){
  # generate pixel-based grid given the expected number of pixels and the reference grid
  
  # The ratio of the coarse grid (with number of pixels for evaluation) over the reference (fine) grid (raster r)
  scaleRatio <- sqrt(NumPix/nrow(RegGrd))
  
  Grd <- expand.grid(list(X_COORD=seq(r@extent@xmin,r@extent@xmax,length.out=round(scaleRatio*r@ncols)),
                          Y_COORD=seq(r@extent@ymin,r@extent@ymax,length.out=round(scaleRatio*r@nrows))))
  Raster <- raster(ncol=round(scaleRatio*r@ncols),nrow=round(scaleRatio*r@nrows),
                   xmn=r@extent@xmin,xmx=r@extent@xmax,ymn=r@extent@ymin,ymx=r@extent@ymax,crs=r@crs) 
  return(list(grd=Grd,raster=Raster))
}