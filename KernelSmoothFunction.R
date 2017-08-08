filterConv <- function(x,y){
  # Note that in R, the usual definition of convolution of two sequences x and y is given by convolve(x, rev(y), type = "open")
  # "filter" returns the middle sub-vector of "open". By zero-padding x, it can have the same effect of conv(x,y,'same') as in MATLAB
  return(zapsmall(convolve(c(rep(0,floor(length(y)/2)),x,rep(0,floor(length(y)/2))),rev(y),type="filter")))
}

truncConv <- function(x,y){
  # x is of length M, y is of length N, linear convolution returns a sequence z of length M+N-1
  # truncConv returns the first M samples of z
  z <- zapsmall(convolve(x,rev(y),type="open"))
  return(z[1:length(x)])
}

SpatialKernSmCrime <- function(CrimeData,Grid,raster,kernel.x,kernel.y,isInPoly,prj){
  Crime.NeighborArray <- matrix(NA,nrow=raster@nrows,ncol=raster@ncols)
  NeighborArray_KS.xy <- matrix(NA,nrow=raster@nrows,ncol=raster@ncols)
  KernSm.df_full <- Grid[,c("X_COORD","Y_COORD")]
  
  CrimePts <- subset(CrimeData,select=c("X_COORD","Y_COORD","INC_CNT"))
  CrimePts.raster <- rasterize(CrimePts[,c("X_COORD","Y_COORD")], raster, CrimePts$INC_CNT, fun=sum)
  Crime.NeighborArray <- as.matrix(CrimePts.raster)
  Crime.NeighborArray[is.na(Crime.NeighborArray)] <- 0
  
  NeighborArray_KS.x <- apply(Crime.NeighborArray,1,filterConv,kernel.x)
  NeighborArray_KS.xy <- apply(NeighborArray_KS.x,1,filterConv,kernel.y)
  
  KernSm.df_full$KS_VAL <- as.vector(t(NeighborArray_KS.xy))
  
  KernSm.sp <- KernSm.df_full
  coordinates(KernSm.sp) <- c("X_COORD", "Y_COORD") 
  proj4string(KernSm.sp) <- prj
  KernSm.sp <- as(KernSm.sp,"SpatialPointsDataFrame")
  
  KernSm.sp_inPoly <- KernSm.sp
  KernSm.sp_inPoly@data <- data.frame(KS_VAL=KernSm.sp@data$KS_VAL[isInPoly])
  KernSm.sp_inPoly@coords <- KernSm.sp@coords[isInPoly,]
  
  KernSm.df_inPoly <- as.data.frame(KernSm.sp_inPoly)
  
  return(list(KernSm.df_inPoly=KernSm.df_inPoly,KernSm.sp_inPoly=KernSm.sp_inPoly))
} 

minmaxScale <- function(x,center=min(x),scale=max(x)-min(x)){
  x <- (x-center)/scale
  return(x)
}