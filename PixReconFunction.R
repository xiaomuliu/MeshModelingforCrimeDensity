# REMARK: For pixel-based (regular grid) reconstruction, all work is done over the full grids instead of the grids inside the city
# The reconstruction pixels can be constrained inside the city area in the end

EucDist <- function(x1, x2) {sqrt(sum((x1 - x2)^2))} # vector-wise Euclidean distance
VecL2Norm <- function(x){return(sqrt(t(x)%*%x))} # vector L2 norm

PixInterpBasis <- function(Grd1,Grd2,method=c('linear','nearest','constant')){
  # Generate interpolation basis matrix 
  # Grd1: coarse grid 
  # Grd2: fine grid (full)
  # method: 
  # linear: bilinear interpolation; 
  # nearest: interpolation using the nearest neighor's value
  # constant: piecewise constant using the average of four vertices values
  
  InterpBasis <- matrix(0,nrow=nrow(Grd2),ncol=nrow(Grd1))
  if (method=='linear'){
    for (i in 1:nrow(InterpBasis)){
      # Find four nearest vertices (only works if assuming grid sizes are equal in both x and y directions)
      # VertexIdx <- sort.int(apply(Grd1,MARGIN=1,FUN=EucDist,Grd2[i,]),index.return=TRUE)$ix[1:4]
      # x1 <- min(Grd1$X_COORD[VertexIdx])
      # y1 <- min(Grd1$X_COORD[VertexIdx])
      # x2 <- max(Grd1$X_COORD[VertexIdx])
      # y2 <- max(Grd1$X_COORD[VertexIdx])
      
      # Rigorously, when coarse grid cells are not square, 
      # one should first find two nearest vertices in x direction then find two nearest vertices in y direction
      V.x <- unique(Grd1$X_COORD[order(abs(Grd2$X_COORD[i]-Grd1$X_COORD))])[1:2]
      V.y <- unique(Grd1$Y_COORD[order(abs(Grd2$Y_COORD[i]-Grd1$Y_COORD))])[1:2]
      x1 <- V.x[1]
      x2 <- V.x[2]
      y1 <- V.y[1]
      y2 <- V.y[2]
      
      # reorder VertexIdx in the order of [ (x1,y1), (x1,y2), (x2,y1), (x2,y2) ]
      VertexIdx <- rep(0,4)
      VertexIdx[1] <- which(Grd1$X_COORD==x1&Grd1$Y_COORD==y1)
      VertexIdx[2] <- which(Grd1$X_COORD==x1&Grd1$Y_COORD==y2)
      VertexIdx[3] <- which(Grd1$X_COORD==x2&Grd1$Y_COORD==y1)
      VertexIdx[4] <- which(Grd1$X_COORD==x2&Grd1$Y_COORD==y2)
      
      x <- Grd2$X_COORD[i]
      y <- Grd2$Y_COORD[i]
      
      # The below method leads to an ill-conditioned system matrix as x,y coordinates are large integers (not preferred)
      # sysMat <- matrix(c(1,1,1,1,x1,x1,x2,x2,y1,y2,y1,y2,x1*y1,x1*y2,x2*y1,x2*y2),ncol=4,byrow=FALSE)
      # InterpBasis[i,VertexIdx] <- t(solve(sysMat)) %*% c(1,x,y,xy)
      
      # Explicitly express the interpolation instead
      a1 <- (x2-x)*(y2-y)/((x2-x1)*(y2-y1))
      a2 <- (x-x1)*(y2-y)/((x2-x1)*(y2-y1))
      a3 <- (x2-x)*(y-y1)/((x2-x1)*(y2-y1))
      a4 <- (x-x1)*(y-y1)/((x2-x1)*(y2-y1))
      InterpBasis[i,VertexIdx] <- c(a1,a2,a3,a4)
    }
  }else if (method=='nearest'){
    for (i in 1:nrow(InterpBasis)){
      # This method below is slow
      VertexIdx <- which.min(apply(Grd1,MARGIN=1,FUN=EucDist,Grd2[i,]))
      InterpBasis[i,VertexIdx] <- 1
    }
  }else if (method=='constant'){
    for (i in 1:nrow(InterpBasis)){
      # This method below is slow
      # VertexIdx <- sort.int(apply(Grd1,MARGIN=1,FUN=EucDist,Grd2[i,]),index.return=TRUE)$ix[1:4]
      
      # Use this method instead
      # Rigorously, when coarse grid cells are not square, 
      # one should first find two nearest vertices in x direction then find two nearest vertices in y direction
      V.x <- unique(Grd1$X_COORD[order(abs(Grd2$X_COORD[i]-Grd1$X_COORD))])[1:2]
      V.y <- unique(Grd1$Y_COORD[order(abs(Grd2$Y_COORD[i]-Grd1$Y_COORD))])[1:2]

      # reorder VertexIdx in the order of [ (x1,y1), (x1,y2), (x2,y1), (x2,y2) ]
      VertexIdx <- rep(0,4)
      VertexIdx[1] <- which(Grd1$X_COORD==V.x[1]&Grd1$Y_COORD==V.y[1])
      VertexIdx[2] <- which(Grd1$X_COORD==V.x[1]&Grd1$Y_COORD==V.y[2])
      VertexIdx[3] <- which(Grd1$X_COORD==V.x[2]&Grd1$Y_COORD==V.y[1])
      VertexIdx[4] <- which(Grd1$X_COORD==V.x[2]&Grd1$Y_COORD==V.y[2])
      
      InterpBasis[i,VertexIdx] <- 0.25
    }
  }

  # For points in Grd1 that are outside Grd2:
  # Put extraplotation code here (Leave for future work)
  return(InterpBasis)
}


PixVal_MLE_EM <- function(Grd,InterpBasis,ImPix,reltol=1e-8,iterMax=20,eps=1e-10){
  converge.flag <- FALSE
  PixValEst_old <- runif(nrow(Grd),min=0,max=max(ImPix))
  PixValEst_new <- PixValEst_old
  
  cnt <- 0  
  while (!converge.flag & cnt<iterMax){
    pixArray <- InterpBasis %*% PixValEst_old
    pixArray[pixArray==0] <- eps
    for (n in 1:nrow(Grd)){
      PixValEst_new[n] <- PixValEst_old[n]/sum(InterpBasis[,n])*
        sum((InterpBasis[,n]*ImPix)/pixArray) 
    }
    cnt <- cnt + 1
    converge.flag <- VecL2Norm(PixValEst_new-PixValEst_old)<reltol*VecL2Norm(PixValEst_old)
    PixValEst_old <- PixValEst_new
  }
  return(PixValEst_new)
}


source("AdjList.R")

PixDiffAdjMat <- function(Grd,r,type=c("queen","rook")){
  # Calculate the derivative of potential function (Euclidean-distance based) (up to a scale of 1/2)
  # If ith and jth have an edge, then C(i,j)=-1. Diagnal C(i,i)= (total # of neighbors) so that sum(C(i,))=0
  
  NumGrd <- nrow(Grd)
  DiffAdjMat <- matrix(0,NumGrd,NumGrd)
  AdjList <- RasterAdjList(r,type=type)
  
  for (i in 1:NumGrd){
    DiffAdjMat[i,AdjList[[i]]] <- -1
  }
  diag(DiffAdjMat) <- -rowSums(DiffAdjMat)
  return(DiffAdjMat)
}


PixVal_MAP_EM <- function(Grd,InterpBasis,ImPix,Raster,type=c("queen","rook"),SmoothParam=1,reltol=1e-8,iterMax=20,eps=1e-10){
  converge.flag <- FALSE
  PixValEst_old <- runif(nrow(Grd),min=0,max=max(ImPix))
  PixValEst_new <- PixValEst_old
  DiffAdjMat <- PixDiffAdjMat(Grd,Raster,type)
  
  cnt <- 0
  while (!converge.flag & cnt<iterMax){
    pixArray <- InterpBasis %*% PixValEst_old
    pixArray[pixArray==0] <- eps
    for (n in 1:nrow(Grd)){
      PixValEst_new[n] <- PixValEst_old[n]/(sum(InterpBasis[,n])+2*SmoothParam*(DiffAdjMat[n,] %*% PixValEst_old))*
        sum((InterpBasis[,n]*ImPix)/pixArray)
    }
    cnt <- cnt + 1
    if (any(PixValEst_new<0)){
      # MAP-EM alogrithm may yield negative estimates which violates Poisson assumption
      PixValEst_new <- NA 
      break
    }else{
      converge.flag <- VecL2Norm(PixValEst_new-PixValEst_old)<reltol*VecL2Norm(PixValEst_old)
      PixValEst_old <- PixValEst_new
    }
  }
  return(PixValEst_new)
}


PixRecon <- function(CrimeObsPts,PixGrd,KDEgrd,iterMax=20,reltol=1e-8,eps=1e-10,Estimation=c("ML","MAP"),
                     interpMethod=c('linear','nearest','constant'),Raster=NULL,type=c("queen","rook"),SmoothParam=1){  
  ## Reconstruction
  # Interpolation basis
  InterpBasis <- PixInterpBasis(PixGrd,KDEgrd,method=interpMethod)
  
  reltol <- 1e-8
  if (Estimation=="ML"){
    PixValEst <- PixVal_MLE_EM(PixGrd,InterpBasis,ImPix=CrimeObsPts$INC_CNT,reltol=reltol,iterMax=iterMax,eps=eps)
  }else if (Estimation=="MAP"){
    PixValEst <- PixVal_MAP_EM(PixGrd,InterpBasis,ImPix=CrimeObsPts$INC_CNT,Raster=Raster,type=type,
                               SmoothParam=SmoothParam,reltol=reltol,iterMax=iterMax,eps=eps)
  }
  
  ## Interpolation
  Recon.city <- KDEgrd
  if (any(is.na(PixValEst))){
    Recon.city$VALUE<-NA
  }else{
    Recon.city$VALUE <- rep(0,nrow(Recon.city))
    Recon.city$VALUE <- InterpBasis %*% PixValEst
  }
  
  return(Recon.city)
}


################################################################################
# Theoratically, ML estimates of pixel values are the observed values themselves
################################################################################
# PixVal_ML <- function(ImPix,NumPix,RegGrd){
#   
#   
#   # The ratio of the coarse grid (with number of pixels for evaluation) over the reference (fine) grid
#   scaleRatio <- sqrt(NumPix[i]/nrow(RegGrd))
#   
#   PixRaster <- raster(ncol=round(scaleRatio*grd.full@grid@cells.dim[1]),nrow=round(scaleRatio*grd.full@grid@cells.dim[2]),
#                       xmn=grd.full@bbox[1,1],xmx=grd.full@bbox[1,2],ymn=grd.full@bbox[2,1],ymx=grd.full@bbox[2,2],crs=CRS(prj)) 
#   CrimeObsPixRaster <- rasterize(CrimeObsPts.df_inCity[,c("X_COORD","Y_COORD")],PixRaster,CrimeObsPts.df_inCity$INC_CNT,fun=sum)
#   
#   CrimeObsPixProj <- resample(CrimeObsPixRaster,r,method='ngb')  # disaggregate to resolution of raster r
#   # CrimeObsPixProj <- disaggregate(CrimeObsPixRaster, round(1/scaleRatio), method='bilinear')  # disaggregate to resolution of raster r
#   
#   CrimeObsPix.df_inCity <- as.data.frame(CrimeObsPixProj,xy=TRUE)[isInCity,]
#   names(CrimeObsPix.df_inCity) <- c("X_COORD","Y_COORD","VALUE")
#   CrimeObsPix.df_inCity$VALUE[is.na(CrimeObsPix.df_inCity$VALUE)] <- 0
#   # # The summation of values in disaggregated cells should be equal to the value in the bigger cell
#   # # This step is unnecessary as the density is calculated in the end
#   # CrimeObsPix.df_inCity$VALUE <- CrimeObsPix.df_inCity$VALUE * scaleRatio
#   
#   # re-arrange the order according to KDE grid
#   CrimeObsPix.df_inCity <- CrimeObsPix.df_inCity[ord,]
# }

# Interp2RefGrd <- function(Grd1,Grd2,data,method=c('linear','nearest'),isInCity=rep(TRUE,nrow(Grd2))){
#   # Interpolate data frome Grd1 to Grd2
#   
#   # Grd1: coarse grid 
#   # Grd2: fine grid
#   interpX <- unique(Grd2$X_COORD) 
#   interpY <- unique(Grd2$Y_COORD)
#   
#   # method 1 (supports linear interpolation only)
#   library(akima)
#   # data must be in the form of a vector with length corresponding to Grd1
#   interpMat <- interp(x=Grd1$X_COORD, y=Grd2$Y_COORD, z=data,xo=interpX, yo=interpY,
#                       linear=TRUE, extrap=TRUE, duplicate="mean")$z
#   
#   # method 2 (supports linear interpolation only)
#   library(fields)
#   # data must be in the form of matrix with dimension corresponding to Grd1
#   interpMat <- interp.surface.grid(obj=list(x=Grd1$X_COORD,y=Grd1$Y_COORD,z=data), grid.list=list(x=interpX,y=interpY))$z
#   
#   # method 3
#   library(pracma)
#   # x, y must be vectors with monotonically increasing elements 
#   # z must be a length(x)-by-length(y) matrix
#   xi <- sort(unique(Grd1$X_COORD))
#   yi <- sort(unique(Grd1$Y_COORD))
#   interpMat <- interp2(x=xi, y=yi, Z=data, xp=interpX, yp=interpY, method=method)
#   
#   interpDf <- Grd2[,c("X_COORD","Y_COORD")]
#   interpDf$Val <- as.vector(interpMat)
#   interpDf <- interpDf[isInCity,]
#   return(interpDf)
# }