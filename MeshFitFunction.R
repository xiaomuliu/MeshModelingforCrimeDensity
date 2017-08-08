InterpolationBasis <- function(grd,MeshNode,MeshTri,BaryCoordList){
  # Generate interpolation basis matrix using conversion between barycentric and Cartesian coordinates on triangles
  InterpBasis <- matrix(0,nrow=nrow(grd),ncol=nrow(MeshNode))
  for (i in 1:nrow(InterpBasis)){
    if (is.na(BaryCoordList$idx[i])){
      # The point to be interpolated is not in any elements (triangles)
      next
    }else{
      VertexIdx <- MeshTri[BaryCoordList$idx[i],]
      InterpBasis[i,VertexIdx] <- BaryCoordList$p[i,]
    }  
  }
  return(InterpBasis)
}


FixBasisErr <- function(InterpBasis,grd,MeshNode){
  # Deal with the issue that a point to be interpolated is not in any elements (triangles) so that the bases are zeros
  # They all lie on the boundary corners. We ignore these pixels.
  if (any(rowSums(InterpBasis)==0)){
    OutPix <- rowSums(InterpBasis)==0
  }else{
    OutPix <- rep(FALSE,nrow(InterpBasis))
  }
  grd2 <- grd[!OutPix,]
  
  InterpBasis2 <- InterpBasis[!OutPix,]
  
  # Deal with the issue that a vetex is not in any pixels.
  # They all lie on the boundary. We shift them to the nearest pixels inside boundary
  if (any(colSums(InterpBasis2)==0)){
    OutVertex <- which(colSums(InterpBasis2)==0)
    for (i in OutVertex){
      NodeCoord <- MeshNode[i,]
      NodeCoord_mat <- matrix(rep(NodeCoord,nrow(grd2)),nrow=nrow(grd2),byrow=TRUE)
      distVec <- apply(as.matrix(grd2)-NodeCoord_mat, 1, FUN=function(x){return(sqrt(x[1]^2+x[2]^2))})
      MeshNode[i,] <- as.vector(grd2[which.min(distVec),],mode="numeric")
      
      # Some mesh nodes are single points that are not associated with any triangles 
      # or all the pixels inside those triangle are at places out of boundary (convex hull problem)
      InterpBasis2[which(grd2$X_COORD==MeshNode[i,1]&grd2$Y_COORD==MeshNode[i,2]),i] <- 1
    } 
  }else{
    OutVertex <- 0
  }
  
  return(list(InterpBasis2=InterpBasis2,MeshNode=MeshNode,OutPix=OutPix,OutVertex=OutVertex,grd2=grd2))
}


MeshLSfit <- function(KDE,MeshNode,MeshTri,KDEgrd){  
  ## Mesh image representation (Least-square fit)
  # Interpolation basis
  # linear interpolation through Barycentric coordinate system
  BaryCoordList <- tsearch(MeshNode[,1],MeshNode[,2],MeshTri,KDEgrd$X_COORD,KDEgrd$Y_COORD,bary=TRUE)
  
  InterpBasis <- InterpolationBasis(KDEgrd,MeshNode,MeshTri,BaryCoordList)
  
  # Deal with these issues:
  # (i) A point to be interpolated is not in any elements (triangles) so that the bases are zeros;
  # (ii) A vetex is not in any pixels.
  if (any(rowSums(InterpBasis)==0) | any(colSums(InterpBasis)==0)){
    FixBasis <- FixBasisErr(InterpBasis,KDEgrd,MeshNode)
    InterpBasis2 <- FixBasis$InterpBasis2
    MeshNode <- FixBasis$MeshNode
    OutPix <- FixBasis$OutPix
  }else{
    InterpBasis2 <- InterpBasis
    OutPix <- rep(FALSE,nrow(KDEgrd))
  }
  KDE2 <- KDE[!OutPix,]
  
  # LS fit
  # KDE_hat <- solve(t(InterpBasis2) %*% InterpBasis2) %*% t(InterpBasis2) %*% KDE2$DENVAL
  KDE_fit <- lsfit(InterpBasis2, KDE2$DENVAL,intercept=FALSE)
  
  ## Interpolation 
  KDE_hat.city <- KDEgrd
  KDE_hat.city$VALUE <- rep(0,nrow(KDE_hat.city))
  KDE_hat.city$VALUE[!OutPix] <- InterpBasis2 %*% KDE_fit$coef
  
  return(list(KDE_hat.city=KDE_hat.city,MeshNode=MeshNode))
}