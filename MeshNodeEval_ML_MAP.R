# Evaluate Mesh reconstruction of crime density map
# Nodal values are estimated by ML and MAP estimation 

setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/MeshModeling/")

PredTarget <- "ViolentCrime"
source("SetupGrid.R")
source("SetupCrimeData.R")

## Crime density and its corresponding mesh
source("MeshGenFunction.R")
source("IMfunction.R")
source("MeshReconFunction.R")
source("MeshEvalFunction.R")

startDate.mesh <- as.Date("2001-01-01")
endDate.mesh <- as.Date("2010-12-31")
CrimeHistPts <- subset(CrimeData,DATEOCC>=startDate.mesh & DATEOCC<=endDate.mesh,select=c("X_COORD","Y_COORD","INC_CNT"))

startDate.obs <- as.Date("2011-01-01")
endDate.obs <- as.Date("2011-12-31")
CrimeObsPts <- subset(CrimeData,DATEOCC>=startDate.obs & DATEOCC<=endDate.obs,select=c("X_COORD","Y_COORD","INC_CNT"))

CrimeObsPts.raster <- rasterize(CrimeObsPts[,c("X_COORD","Y_COORD")], r, CrimeObsPts$INC_CNT, fun=sum)
CrimeObsPts.df_inCity <- as.data.frame(CrimeObsPts.raster,xy=TRUE)[isInCity,]
names(CrimeObsPts.df_inCity) <- c("X_COORD","Y_COORD","INC_CNT")
CrimeObsPts.df_inCity$INC_CNT[is.na(CrimeObsPts.df_inCity$INC_CNT)] <- 0
CrimeObsPts.df_inCity <- CrimeObsPts.df_inCity[ord,]

# Extract City Boundary 
isInCity2_mat <- matrix(isInCity2,nrow=r@ncols,ncol=r@nrows)
Bndy_mat <- DetectEdge(isInCity2_mat)
isOnBndy <- as.vector(Bndy_mat)  

h <- grd.full@grid@cellsize  # kernel bandwidth for KDE
# EnodeSet <- round(2^seq(8.5,11,by=0.5)) # expected number of nodes
EnodeSet <- seq(200,2000,by=200)
bndyPercent <- 0.1 # the percent of nodes on boundary
# SmParamSet <- 10^seq(-4,-3,by=0.5) # Gibbs prior parameter
SmParamSet <- seq(0.0005,0.0015,by=0.0005)
iterMax <- 100 # maximal EM iterations
reltol <- 1e-8 # relevance tolerance for convergence
eps <- 1e-20
NumNode <- rep(0,length(EnodeSet))
Mesh.ML_PSNR <- rep(0,length(EnodeSet))
Mesh.ML_repErr <- rep(0,length(EnodeSet))
Mesh.MAP_PSNR <- matrix(0,nrow=length(EnodeSet),ncol=length(SmParamSet))
Mesh.MAP_repErr <- matrix(0,nrow=length(EnodeSet),ncol=length(SmParamSet))

ptm <- proc.time()
for (i in 1:length(EnodeSet)){
  ENode <- EnodeSet[i] 
  BndyNode <- round(bndyPercent*ENode)
  
  ## generate constrained Mesh
  meshList <- generateConstraintMesh2(CrimeHistPts,r,KDEgrd.full,ENode,isOnBndy,BndyNode,bw=h,plot=FALSE)
  MeshNode_raw <- meshList$meshNode
  MeshTri <- meshList$Tri
  DenIm <- meshList$pixelIm
  
  KDE.Hist_full <- KDEgrd.full[,c("X_COORD","Y_COORD")]
  KDE.Hist_full$VALUE <- as.vector(DenIm)
  KDE.Hist_inCity <- KDE.Hist_full[isInCity2,]  
  KDE.Hist_inCity$DENVAL <- KDE.Hist_inCity$VALUE/sum(KDE.Hist_inCity$VALUE) 
  
  #   # assumed "true" density of observed data
  #   KDE.Obs_full <- KDEgrd.full
  #   kernSm <- bkde2D(data.matrix(CrimeObsPts[,c("X_COORD","Y_COORD")]), bandwidth=h, 
  #                    gridsize=c(r@ncols, r@nrows), range.x=list(X_range,Y_range))
  #   KDE.Obs_full$VALUE <- as.vector(kernSm$fhat)      
  #   KDE.Obs_inCity <- KDE.Obs_full[isInCity2,] 
  #   KDE.Obs_inCity$DENVAL <- KDE.Obs_inCity$VALUE/sum(KDE.Obs_inCity$VALUE)
  
  # ML
  ReconList <- MeshRecon(CrimeObsPts.df_inCity,MeshNode_raw,MeshTri,KDEgrd,iterMax=iterMax,Estimation="ML",reltol=reltol,eps=eps)
  Recon.city <- ReconList$Recon.city
  MeshNode <- ReconList$MeshNode
  
  NumNode[i] <- nrow(MeshNode)
  Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE) 
  
  # Evaluation
  err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL)) 
  Mesh.ML_PSNR[i] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2) # peak SNR
  Mesh.ML_repErr[i] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0)) # representation error (percent RMSE)
  
  # MAP
  for (j in 1:length(SmParamSet)){
    ReconList <- MeshRecon(CrimeObsPts.df_inCity,MeshNode_raw,MeshTri,KDEgrd,iterMax=iterMax,Estimation="MAP",SmoothParam=SmParamSet[j])
    Recon.city <- ReconList$Recon.city
    MeshNode <- ReconList$MeshNode
    
    Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE)  
    
    err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL)) 
    Mesh.MAP_PSNR[i,j] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2) # peak SNR
    Mesh.MAP_repErr[i,j] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0)) # representation error (percent RMSE)
  }
  
}
proc.time()-ptm

# For some bigger smoothing parameters,EM does not converge (produces negative estimates that violents Poisson assumption)
# Suppress these results by NA
# Future work: use gradient-descent-based methods for constrained optimizations

source("PlotRecon_mesh.R")
