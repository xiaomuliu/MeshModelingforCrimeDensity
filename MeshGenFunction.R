require(sp)
require(raster)
require(spatstat)
require(KernSmooth)
require(geometry)
source("KernelSmoothFunction.R")
source("Dither.R")

generateMesh <- function(CrimeData,r,grd,ENode,bw=NA,plot=TRUE){

  # tryCatch(r@nrows*r@ncols!=nrow(grd), error = function(e) print(e), finally=print("Aborted")) 
  
  # subset data
  CrimeData <- subset(CrimeData,select=c("X_COORD","Y_COORD","INC_CNT"))
  
  ## KDE
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
  KDE.df$VALUE <- as.vector(kernSm$fhat)      
#   # scale KDE values to [0,1]
#   KDE.df$val_scale <- minmaxScale(KDE.df$VALUE)
  KDE.df$val_scale <- KDE.df$VALUE
  
  # error diffusion 
  threshold <- sum(KDE.df$val_scale)/(2*ENode)
  pixelIm <- matrix(KDE.df$val_scale,N,M) 
  
  b <- errorDiffusion(pixelIm, 1, threshold)
  
  # node placement
  Reg_x <- matrix(grd$X_COORD,N,M)
  Reg_y <- matrix(grd$Y_COORD,N,M)
  
  temp <- t(Reg_x)
  Node_x <- temp[t(b)==1]
  temp <- t(Reg_y)
  Node_y <- temp[t(b)==1]
  
  Node_x[Node_x>X_range[2]] <- X_range[2]
  Node_y[Node_y>Y_range[2]] <- Y_range[2]
  Node_x[max(Node_x)==Node_x] <- X_range[2]
  Node_y[max(Node_y)==Node_y] <- Y_range[2]
  
  # Delaunay triangulation
  pts <- cbind(Node_x,Node_y)
  tri <- delaunayn(pts)
  triSides <- rbind(tri[, -1], tri[, -2], tri[, -3])
  
  if (plot==TRUE){
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    par(mfrow=c(1,2),oma=c(0,0,2,0))
    image(x=seq(X_range[1],X_range[2],length.out=N),y=seq(Y_range[1],Y_range[2],length.out=M),
          z=pixelIm,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
    plot.new()
    plot.window(xlim=c(X_range[1],X_range[2]),ylim=c(Y_range[1],Y_range[2]))
    axis(1)
    axis(2)
    box()
    segments(pts[triSides[, 1], 1], pts[triSides[, 1], 2], pts[triSides[, 2], 1], pts[triSides[, 2], 2], col="blue")
  }
  
  return(list(Tri=tri,meshNode=pts,triSides=triSides,pixelIm=pixelIm))
}


generateConstraintMesh <- function(CrimeData,r,grd,ENode,BndyCoord,BndyNode=round(0.1*ENode),bw=NA,plot=TRUE){
  # Generate nodes uniformly along the boundary to form constrained mesh grids
  
  # subset data
  CrimeData <- subset(CrimeData,select=c("X_COORD","Y_COORD","INC_CNT"))
  
  ## KDE
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
  KDE.df$VALUE <- as.vector(kernSm$fhat)      
  # #scale KDE values to [0,1]
  # KDE.df$val_scale <- minmaxScale(KDE.df$VALUE)
  KDE.df$val_scale <- KDE.df$VALUE

  # error diffusion 
  threshold <- sum(KDE.df$val_scale)/(2*ENode)
  pixelIm <- matrix(KDE.df$val_scale,N,M) 
  
  b <- errorDiffusion(pixelIm, 1, threshold)
  
  # node placement
  Reg_x <- matrix(grd$X_COORD,N,M)
  Reg_y <- matrix(grd$Y_COORD,N,M)
  
  temp <- t(Reg_x)
  Node_x <- temp[t(b)==1]
  temp <- t(Reg_y)
  Node_y <- temp[t(b)==1]
  
  Node_x[Node_x>X_range[2]] <- X_range[2]
  Node_y[Node_y>Y_range[2]] <- Y_range[2]
  Node_x[max(Node_x)==Node_x] <- X_range[2]
  Node_y[max(Node_y)==Node_y] <- Y_range[2]
  
  # Uniformly sample nodes on the boundary
  BndyCoord.Sample <- BndyCoord[seq(1,nrow(BndyCoord),length.out=BndyNode),]
  # Delaunay triangulation
  pts <- cbind(c(Node_x,BndyCoord.Sample[,1]),c(Node_y,BndyCoord.Sample[,2]))
  tri.bndy <- delaunayn(pts) 
  triSides <- rbind(tri.bndy[, -1], tri.bndy[, -2], tri.bndy[, -3])
  
  if (plot==TRUE){
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    par(mfrow=c(1,2),oma=c(0,0,0,0))
    image(x=seq(X_range[1],X_range[2],length.out=N),y=seq(Y_range[1],Y_range[2],length.out=M),
          z=pixelIm,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
    plot.new()
    plot.window(xlim=c(X_range[1],X_range[2]),ylim=c(Y_range[1],Y_range[2]))
    axis(1)
    axis(2)
    box()
    segments(pts[triSides[, 1], 1], pts[triSides[, 1], 2], pts[triSides[, 2], 1], pts[triSides[, 2], 2], col="blue")
    
    par(mfrow=c(1,1),oma=c(0,0,0,0))
    image(x=seq(X_range[1],X_range[2],length.out=N),y=seq(Y_range[1],Y_range[2],length.out=M),
          z=pixelIm,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
    trimesh(tri.bndy,pts,add=TRUE,axis=FALSE,boxed=TRUE,col="red")
  }
  
  return(list(Tri=tri.bndy,meshNode=pts,triSides=triSides,pixelIm=pixelIm))
}



generateConstraintMesh2 <- function(CrimeData,r,grd,ENode,isOnBndy,BndyNode=round(0.1*ENode),bw=NA,plot=TRUE){
  # The KDE of raw crime data is digitized as a image. Therefore when sampling boundary points, some of them may fall out of pixels.
  # This function generates constrained mesh grids and forces sampled boundary points in pixels.
  
  # subset data
  CrimeData <- subset(CrimeData,select=c("X_COORD","Y_COORD","INC_CNT"))
  
  ## KDE
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
  KDE.df$VALUE <- as.vector(kernSm$fhat)      
  #   #scale KDE values to [0,1]
  #   KDE.df$val_scale <- minmaxScale(KDE.df$VALUE)
  KDE.df$val_scale <- KDE.df$VALUE
  
  # error diffusion 
  threshold <- sum(KDE.df$val_scale)/(2*ENode)
  pixelIm <- matrix(KDE.df$val_scale,N,M) 
  
  b <- errorDiffusion(pixelIm, 1, threshold)
  
  # node placement
  Reg_x <- matrix(grd$X_COORD,N,M)
  Reg_y <- matrix(grd$Y_COORD,N,M)
  
  temp <- t(Reg_x)
  Node_x <- temp[t(b)==1]
  temp <- t(Reg_y)
  Node_y <- temp[t(b)==1]
  
  Node_x[Node_x>X_range[2]] <- X_range[2]
  Node_y[Node_y>Y_range[2]] <- Y_range[2]
  Node_x[max(Node_x)==Node_x] <- X_range[2]
  Node_y[max(Node_y)==Node_y] <- Y_range[2]
  
  # Uniformly sample nodes on the boundary pixels
  BndyCoord <- grd[isOnBndy,]
  BndyCoord.Sample <- BndyCoord[seq(1,nrow(BndyCoord),length.out=BndyNode),]
  
  # Delaunay triangulation
  pts <- cbind(c(Node_x,BndyCoord.Sample[,1]),c(Node_y,BndyCoord.Sample[,2]))
  tri.bndy <- delaunayn(pts) 
  triSides <- rbind(tri.bndy[, -1], tri.bndy[, -2], tri.bndy[, -3])
  
  if (plot==TRUE){
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    par(mfrow=c(1,2),oma=c(0,0,0,0))
    image(x=seq(X_range[1],X_range[2],length.out=N),y=seq(Y_range[1],Y_range[2],length.out=M),
          z=pixelIm,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
    plot.new()
    plot.window(xlim=c(X_range[1],X_range[2]),ylim=c(Y_range[1],Y_range[2]))
    axis(1)
    axis(2)
    box()
    segments(pts[triSides[, 1], 1], pts[triSides[, 1], 2], pts[triSides[, 2], 1], pts[triSides[, 2], 2], col="blue")
    
    par(mfrow=c(1,1),oma=c(0,0,0,0))
    image(x=seq(X_range[1],X_range[2],length.out=N),y=seq(Y_range[1],Y_range[2],length.out=M),
          z=pixelIm,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
    trimesh(tri.bndy,pts,add=TRUE,axis=FALSE,boxed=TRUE,col="red")
  }
  
  return(list(Tri=tri.bndy,meshNode=pts,triSides=triSides,pixelIm=pixelIm))
}