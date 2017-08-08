VecL2Norm <- function(x){return(sqrt(t(x)%*%x))} # vector L2 norm

MeshInterpBasis <- function(grd,MeshNode,MeshTri,BaryCoordList){
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


NodalVal_MLE_EM <- function(MeshNode,InterpBasis,ImPix,reltol=1e-8,iterMax=20,eps=1e-10){
  converge.flag <- FALSE
  NodalValEst_old <- runif(nrow(MeshNode),min=0,max=max(ImPix))
  NodalValEst_new <- NodalValEst_old
  
  cnt <- 0  
  while (!converge.flag & cnt<iterMax){
    pixArray <- InterpBasis %*% NodalValEst_old
    pixArray[pixArray==0] <- eps
    for (n in 1:nrow(MeshNode)){
      NodalValEst_new[n] <- NodalValEst_old[n]/sum(InterpBasis[,n])*
        sum((InterpBasis[,n]*ImPix)/pixArray) 
    }
    cnt <- cnt + 1
    converge.flag <- VecL2Norm(NodalValEst_new-NodalValEst_old)<reltol*VecL2Norm(NodalValEst_old)
    NodalValEst_old <- NodalValEst_new
  }
  return(NodalValEst_new)
}

MeshDiffAdjMat <- function(MeshNode,MeshTri){
  # Calculate the derivative of potential function (Euclidean-distance based) (up to a scale of 1/2)
  # If ith and jth have an edge, then C(i,j)=-1. Diagnal C(i,i)= (total # of neighbors) so that sum(C(i,))=0 
  NumMeshNode <- nrow(MeshNode)
  DiffAdjMat <- matrix(0,NumMeshNode,NumMeshNode)
  
  for (i in 1:NumMeshNode){  
    # find node's neighbors
    simplexInd <- which(MeshTri==i,arr.ind=TRUE)
    neighborNodeIdx <- as.vector(MeshTri[simplexInd[,1],])
    neighborNodeIdx <- unique(neighborNodeIdx[which(neighborNodeIdx!=i)])
    
    DiffAdjMat[i,neighborNodeIdx] <- -1
  }
  diag(DiffAdjMat) <- -rowSums(DiffAdjMat)
  return(DiffAdjMat)
}


NodalVal_MAP_EM <- function(MeshNode,InterpBasis,ImPix,MeshTri,DiffAdjMat=MeshDiffAdjMat(MeshNode,MeshTri),SmoothParam=1,
                            reltol=1e-8,iterMax=20,eps=1e-10){
  converge.flag <- FALSE
  NodalValEst_old <- runif(nrow(MeshNode),min=0,max=max(ImPix))
  NodalValEst_new <- NodalValEst_old
  
  cnt <- 0  
  while (!converge.flag & cnt<iterMax){
    pixArray <- InterpBasis %*% NodalValEst_old
    pixArray[pixArray==0] <- eps
    for (n in 1:nrow(MeshNode)){
      NodalValEst_new[n] <- NodalValEst_old[n]/(sum(InterpBasis[,n])+2*SmoothParam*(DiffAdjMat[n,] %*% NodalValEst_old))*
        sum((InterpBasis[,n]*ImPix)/pixArray) 
    }
    cnt <- cnt + 1
    if (any(NodalValEst_new<0)){
      # MAP-EM alogrithm may yield negative estimates which violates Poisson assumption
      NodalValEst_new <- NA 
      break
    }else{
      converge.flag <- VecL2Norm(NodalValEst_new-NodalValEst_old)<reltol*VecL2Norm(NodalValEst_old)
      NodalValEst_old <- NodalValEst_new
    }
  }
  return(NodalValEst_new)
}


MeshRecon <- function(CrimeObsPts,MeshNode,MeshTri,KDEgrd,iterMax=20,reltol=1e-8,eps=1e-10,
                      Estimation=c("ML","MAP"),DiffAdjMat=NA,SmoothParam=1){  
  ## Reconstruction
  # Interpolation basis
  # linear interpolation through Barycentric coordinate system
  BaryCoordList <- tsearch(MeshNode[,1],MeshNode[,2],MeshTri,KDEgrd$X_COORD,KDEgrd$Y_COORD,bary=TRUE)
  # Due to machine precision, some barycentric coordinates returns tiny negative values
  # Force them to be zero (to avoid divergence to negative values later in the EM algorithm)
  BaryCoordList$p[BaryCoordList$p<0] <- 0
  
  InterpBasis <- MeshInterpBasis(KDEgrd,MeshNode,MeshTri,BaryCoordList)
  
  # Deal with the issues (a) a point to be interpolated is not in any elements (triangles) so that the bases are zeros;
  # (b) a vetex is not in any pixels.
  if (any(rowSums(InterpBasis)==0) | any(colSums(InterpBasis)==0)){
    FixBasis <- FixBasisErr(InterpBasis,KDEgrd,MeshNode)
    InterpBasis2 <- FixBasis$InterpBasis2
    MeshNode <- FixBasis$MeshNode
    OutPix <- FixBasis$OutPix
  }else{
    InterpBasis2 <- InterpBasis
    OutPix <- rep(FALSE,nrow(KDEgrd))
  }
  CrimeObsPts2 <- CrimeObsPts[!OutPix,]
  
  if (Estimation=="ML"){
    NodalValEst <- NodalVal_MLE_EM(MeshNode,InterpBasis2,ImPix=CrimeObsPts2$INC_CNT,
                                   reltol=reltol,iterMax=iterMax,eps=eps)
  }else if (Estimation=="MAP"){
    NodalValEst <- NodalVal_MAP_EM(MeshNode,InterpBasis2,ImPix=CrimeObsPts2$INC_CNT,MeshTri,
                                   DiffAdjMat=MeshDiffAdjMat(MeshNode,MeshTri),SmoothParam=SmoothParam,
                                   reltol=reltol,iterMax=iterMax,eps=eps)
  }

  ## Interpolation 
  Recon.city <- KDEgrd
  if (any(is.na(NodalValEst))){
    Recon.city$VALUE<-NA
  }else{
    Recon.city$VALUE <- rep(0,nrow(Recon.city))
    Recon.city$VALUE[!OutPix] <- InterpBasis2 %*% NodalValEst
  }
  return(list(Recon.city=Recon.city,MeshNode=MeshNode,MeshTri=MeshTri))
}

