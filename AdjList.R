# NOTE: Rasterization scans an image (matrix) row-wise (from left to right) while converting a matrix to a vector 
#       is done through concatinating elements column-wise.
require(spdep)

# Given a raster index, the function returns the corresponding matrix index
RasterIdx2MatVecIdx <- function(RasterIdx,ncol,nrow){
  MatVecIdx <- as.vector(matrix(1:(ncol*nrow),ncol=ncol,nrow=nrow,byrow=T))
  return(which(MatVecIdx==RasterIdx))
}

# Given a matrix index, the function returns the corresponding raster index
MatVecIdx2RasterIdx <- function(MatVecIdx,ncol,nrow){
  RasterIdx <- as.vector(matrix(1:(ncol*nrow),ncol=nrow,nrow=ncol,byrow=T))
  return(which(RasterIdx==MatVecIdx))
}

RasterAdjList <- function(r,type=c("queen","rook")){
  # Given a raster, returns the neighbor index list based on vectorized grids
  # r: a raster object
  # type: connectivity type
  
  AdjList.mat <- cell2nb(ncol=r@ncols,nrow=r@nrows,type=type)
  raster_idx <- 1:(r@ncols*r@nrows)
  AdjList <- vector("list",r@ncols*r@nrows)
  
  for (i in 1:length(raster_idx)){
    mat_idx <- RasterIdx2MatVecIdx(raster_idx[i],ncol=r@ncols,nrow=r@nrows)
    AdjList[[i]] <- sapply(AdjList.mat[[mat_idx]],FUN=MatVecIdx2RasterIdx,ncol=r@ncols,nrow=r@nrows,USE.NAMES=F)
  }

  return(AdjList)
}
