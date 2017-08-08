# Evaluate mesh grid and pixel grid reconstruction of crime density map
# Nodal and pixel values are estimated by ML and MAP estimation 

setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/MeshModeling/")

PredTarget <- "ViolentCrime"
source("SetupGrid.R")
source("SetupCrimeData.R")

source("IMfunction.R")
source("PixGenFunction.R")
source("PixReconFunction.R")
source("MeshEvalFunction.R")

startDate.pix <- as.Date("2001-01-01")
endDate.pix <- as.Date("2010-12-31")
CrimeHistPts <- subset(CrimeData,DATEOCC>=startDate.pix & DATEOCC<=endDate.pix,select=c("X_COORD","Y_COORD","INC_CNT"))

startDate.obs <- as.Date("2011-01-01")
endDate.obs <- as.Date("2012-12-31")
CrimeObsPts <- subset(CrimeData,DATEOCC>=startDate.obs & DATEOCC<=endDate.obs,select=c("X_COORD","Y_COORD","INC_CNT"))

CrimeObsPts.raster <- rasterize(CrimeObsPts[,c("X_COORD","Y_COORD")], r, CrimeObsPts$INC_CNT, fun=sum)
CrimeObsPts.df_full<- as.data.frame(CrimeObsPts.raster,xy=TRUE)
names(CrimeObsPts.df_full) <- c("X_COORD","Y_COORD","INC_CNT")
CrimeObsPts.df_full$INC_CNT[is.na(CrimeObsPts.df_full$INC_CNT)] <- 0
CrimeObsPts.df_full <- CrimeObsPts.df_full[ord.full,]

h <- grd.full@grid@cellsize  # kernel bandwidth for KDE
EpixSet <- seq(200,1700,by=300) # expected number of pixels
SmParamSet <- c(0.001,0.005,0.01,0.1,1) # Gibbs prior parameter
iterMax <- 30 # maximal EM iterations
reltol <- 1e-8 # relevance tolerance for convergence
eps <- 1e-20

Pix.ML_PSNR <- rep(0,length(EpixSet))
Pix.ML_repErr <- rep(0,length(EpixSet))
Pix.MAP_PSNR <- matrix(0,nrow=length(EpixSet),ncol=length(SmParamSet))
Pix.MAP_repErr <- matrix(0,nrow=length(EpixSet),ncol=length(SmParamSet))

## generate density
DenIm <- generateDenIm(CrimeHistPts,r,KDEgrd.full,bw=h)

KDE.Hist_full <- KDEgrd.full[,c("X_COORD","Y_COORD")]
KDE.Hist_full$VALUE <- as.vector(DenIm)
KDE.Hist_inCity <- KDE.Hist_full[isInCity2,]  
KDE.Hist_inCity$DENVAL <- KDE.Hist_inCity$VALUE/sum(KDE.Hist_inCity$VALUE) 

ptm <- proc.time()
for (i in 1:length(EpixSet)){
  EPix <- EpixSet[i] 
  
  ## generate pixel-based grid
  PixList <- generatePixGrd(EPix,RegGrd,r)
  PixGrd <- PixList$grd
  PixRaster <- PixList$raster
    
  # ML
  Recon.city <- PixRecon(CrimeObsPts.df_full,PixGrd,KDEgrd.full,iterMax=iterMax,Estimation="ML",interpMethod='linear',reltol=reltol,eps=eps)
  Recon.city <- Recon.city[isInCity2,]
  Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE)

  # Evaluation
  err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL))
  Pix.ML_PSNR[i] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2) # peak SNR
  Pix.ML_repErr[i] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0)) # representation error (percent RMSE)
  
  # MAP
  for (j in 1:length(SmParamSet)){
    Recon.city <- PixRecon(CrimeObsPts.df_full,PixGrd,KDEgrd.full,iterMax=iterMax,Estimation="MAP",SmoothParam=SmParamSet[j],
                          Raster=PixRaster,type="rook",interpMethod='linear',reltol=reltol,eps=eps)
    Recon.city <- Recon.city[isInCity2,]
    Recon.city$DENVAL <- Recon.city$VALUE/sum(Recon.city$VALUE)

    err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,Recon.city$DENVAL))
    Pix.MAP_PSNR[i,j] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2) # peak SNR
    Pix.MAP_repErr[i,j] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0)) # representation error (percent RMSE)
  }
  
}
proc.time()-ptm

source("PlotRecon_pix.R")