setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/MeshModeling/")

PredTarget <- "ViolentCrime"
source("SetupGrid.R")
source("SetupCrimeData.R")

source("IMfunction.R")
source("PixGenFunction.R")
source("PixReconFunction.R")
source("MeshGenFunction.R")
source("MeshReconFunction.R")
source("MeshEvalFunction.R")

startDate.pix <- as.Date("2001-01-01")
endDate.pix <- as.Date("2010-12-31")
CrimeHistPts <- subset(CrimeData,DATEOCC>=startDate.pix & DATEOCC<=endDate.pix,select=c("X_COORD","Y_COORD","INC_CNT"))

startDate.obs <- as.Date("2011-01-01")
endDate.obs <- as.Date("2012-12-31")
CrimeObsPts <- subset(CrimeData,DATEOCC>=startDate.obs & DATEOCC<=endDate.obs,select=c("X_COORD","Y_COORD","INC_CNT"))

CrimeObsPts.raster <- rasterize(CrimeObsPts[,c("X_COORD","Y_COORD")], r, CrimeObsPts$INC_CNT, fun=sum)
# Used for mesh reconstruction (only needs grids inside the city area)
CrimeObsPts.df_inCity <- as.data.frame(CrimeObsPts.raster,xy=TRUE)[isInCity,]
names(CrimeObsPts.df_inCity) <- c("X_COORD","Y_COORD","INC_CNT")
CrimeObsPts.df_inCity$INC_CNT[is.na(CrimeObsPts.df_inCity$INC_CNT)] <- 0
CrimeObsPts.df_inCity <- CrimeObsPts.df_inCity[ord,]
# Used for pixel reconstruction (needs grids of rectangular area)
CrimeObsPts.df_full<- as.data.frame(CrimeObsPts.raster,xy=TRUE)
names(CrimeObsPts.df_full) <- c("X_COORD","Y_COORD","INC_CNT")
CrimeObsPts.df_full$INC_CNT[is.na(CrimeObsPts.df_full$INC_CNT)] <- 0
CrimeObsPts.df_full <- CrimeObsPts.df_full[ord.full,]


# Extract City Boundary 
isInCity2_mat <- matrix(isInCity2,nrow=r@ncols,ncol=r@nrows)
Bndy_mat <- DetectEdge(isInCity2_mat)
isOnBndy <- as.vector(Bndy_mat)  

h <- grd.full@grid@cellsize  # kernel bandwidth for KDE
#EnodeSet <- round(2^seq(8.5,11,by=0.5)) # expected number of nodes
EnodeSet <- seq(400,1600,by=200) # expected number of nodes
# EpixSet <- EnodeSet # expected number of pixels
bndyPercent <- 0.1 # the percent of nodes on boundary
#SmParamSet <- 10^seq(-5,-3,by=0.5) # Gibbs prior parameter
SmParamSet <- seq(0.0005,0.0015,by=0.0005) # Gibbs prior parameter
iterMax <- 30 # maximal EM iterations
reltol <- 1e-8 # relevance tolerance for convergence
eps <- 1e-20

NumNode <- rep(0,length(EnodeSet))
# peak SNR and representation error (RMSE%)
Mesh.ML_PSNR <- rep(0,length(EnodeSet))
Mesh.ML_repErr <- rep(0,length(EnodeSet))
Mesh.MAP_PSNR <- matrix(0,nrow=length(EnodeSet),ncol=length(SmParamSet))
Mesh.MAP_repErr <- matrix(0,nrow=length(EnodeSet),ncol=length(SmParamSet))
Pix.ML_PSNR <- rep(0,length(EnodeSet))
Pix.ML_repErr <- rep(0,length(EnodeSet))
Pix.MAP_PSNR <- matrix(0,nrow=length(EnodeSet),ncol=length(SmParamSet))
Pix.MAP_repErr <- matrix(0,nrow=length(EnodeSet),ncol=length(SmParamSet))

## generate density
DenIm <- generateDenIm(CrimeHistPts,r,KDEgrd.full,bw=h)

KDE.Hist_full <- KDEgrd.full[,c("X_COORD","Y_COORD")]
KDE.Hist_full$VALUE <- as.vector(DenIm)
KDE.Hist_inCity <- KDE.Hist_full[isInCity2,]  
KDE.Hist_inCity$DENVAL <- KDE.Hist_inCity$VALUE/sum(KDE.Hist_inCity$VALUE) 

ptm <- proc.time()
for (i in 1:length(EnodeSet)){
  ENode <- EnodeSet[i] 
  BndyNode <- round(bndyPercent*ENode)
  
  ## generate constrained Mesh
  meshList <- generateConstraintMesh2(CrimeHistPts,r,KDEgrd.full,ENode,isOnBndy,BndyNode,bw=h,plot=FALSE)
  MeshNode_raw <- meshList$meshNode
  MeshTri <- meshList$Tri
  
  ## generate pixel-based grid
  PixList <- generatePixGrd(ENode,RegGrd,r)
  PixGrd <- PixList$grd
  PixRaster <- PixList$raster

  ## ML
  # Mesh
  ReconList <- MeshRecon(CrimeObsPts.df_inCity,MeshNode_raw,MeshTri,KDEgrd,iterMax=iterMax,Estimation="ML",reltol=reltol,eps=eps)
  Recon.city <- ReconList$Recon.city
  MeshNode <- ReconList$MeshNode
  Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE) 
  NumNode[i] <- nrow(MeshNode)
  err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL)) 
  Mesh.ML_PSNR[i] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2) 
  Mesh.ML_repErr[i] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0)) 
  
  # Pixel
  Recon.city <- PixRecon(CrimeObsPts.df_full,PixGrd,KDEgrd.full,iterMax=iterMax,Estimation="ML",interpMethod='linear',reltol=reltol,eps=eps)
  Recon.city <- Recon.city[isInCity2,]
  Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE)
  err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL))
  Pix.ML_PSNR[i] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2) 
  Pix.ML_repErr[i] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0)) 

  ## MAP
  for (j in 1:length(SmParamSet)){
    # Mesh
    ReconList <- MeshRecon(CrimeObsPts.df_inCity,MeshNode_raw,MeshTri,KDEgrd,iterMax=iterMax,Estimation="MAP",SmoothParam=SmParamSet[j])
    Recon.city <- ReconList$Recon.city
    MeshNode <- ReconList$MeshNode
    if (any(is.na(Recon.city$VALUE))){
      # MAP-EM alogrithm may yield negative estimates which violates Poisson assumption
      Mesh.MAP_PSNR[i,j] <- NA
      Mesh.MAP_repErr[i,j] <- NA
    }else{
      Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE)
      err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL))
      Mesh.MAP_PSNR[i,j] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2)
      Mesh.MAP_repErr[i,j] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0))
    }

    # Pixel
    Recon.city <- PixRecon(CrimeObsPts.df_full,PixGrd,KDEgrd.full,iterMax=iterMax,Estimation="MAP",SmoothParam=SmParamSet[j],
                           Raster=PixRaster,type="queen",interpMethod='linear',reltol=reltol,eps=eps)
    Recon.city <- Recon.city[isInCity2,]
    if (any(is.na(Recon.city$VALUE))){
      # MAP-EM alogrithm may yield negative estimates which violates Poisson assumption
      Mesh.MAP_PSNR[i,j] <- NA
      Mesh.MAP_repErr[i,j] <- NA
    }else{
      Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE)
      err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL))
      Pix.MAP_PSNR[i,j] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2)
      Pix.MAP_repErr[i,j] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0))
    }

  }
}
proc.time()-ptm

source("PlotRecon_mesh_pix.R")
