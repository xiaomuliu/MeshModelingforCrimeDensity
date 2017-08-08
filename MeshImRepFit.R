# Evaluate Mesh representation of crime density map
# The evaluation is done using MDL and PSNR of least square fit through mesh grids

setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/MeshModeling/")

PredTarget <- "ViolentCrime"
source("SetupGrid.R")
source("SetupCrimeData.R")

## Crime density and its corresponding mesh
source("MeshGenFunction.R")
source("IMfunction.R")
source("MeshFitFunction.R")
source("MeshEvalFunction.R")

startDate.mesh <- as.Date("2001-01-01")
endDate.mesh <- as.Date("2010-12-31")
CrimeHistPts <- subset(CrimeData,DATEOCC>=startDate.mesh & DATEOCC<=endDate.mesh,select=c("X_COORD","Y_COORD","INC_CNT"))

# Constrainted Mesh
# Extract City Boundary 
isInCity2_mat <- matrix(isInCity2,nrow=r@ncols,ncol=r@nrows)
Bndy_mat <- DetectEdge(isInCity2_mat)
isOnBndy <- as.vector(Bndy_mat)  

h <- grd.full@grid@cellsize # kernel bandwidth for KDE
EnodeSet <- round(2^seq(8,11.5,by=0.5)) # expected number of nodes
bndyPercent <- 0.1 # the percent of nodes on boundary
NumNode <- rep(0,length(EnodeSet))
Mesh.MDL <- rep(0,length(EnodeSet))
Mesh.PSNR <- rep(0,length(EnodeSet))

ptm <- proc.time()
for (i in 1:length(EnodeSet)){
  ENode <- EnodeSet[i] 
  BndyNode <- round(bndyPercent*ENode)
  
  ## generate constrained Mesh
  meshList <- generateConstraintMesh2(CrimeHistPts,r,KDEgrd.full,ENode,isOnBndy,BndyNode,bw=h,plot=FALSE)
  MeshNode_raw <- meshList$meshNode
  MeshTri <- meshList$Tri
  DenIm <- meshList$pixelIm
  
  KDE.full <- KDEgrd.full[,c("X_COORD","Y_COORD")]
  KDE.full$VALUE <- as.vector(DenIm)
  KDE.inCity <- KDE.full[isInCity2,]  
  KDE.inCity$DENVAL <- KDE.inCity$VALUE/sum(KDE.inCity$VALUE) 

  # LS fit
  FitList <- MeshLSfit(KDE.inCity,MeshNode_raw,MeshTri,KDEgrd)
  KDE_hat.city <- FitList$KDE_hat.city
  MeshNode <- FitList$MeshNode
  
  NumNode[i] <- nrow(MeshNode)
  # Convert to density values
  # KDE_hat.city$DENVAL <- KDE_hat.city$VALUE/sum(KDE_hat.city$VALUE) 
  
  # Evaluation
  err_sigma <- sqrt(MSE(KDE.inCity$DENVAL,KDE_hat.city$VALUE)) 
  Mesh.MDL[i] <- MDL(nrow(KDE.inCity),NumNode[i],err_sigma)
  Mesh.PSNR[i] <- PSNR(KDE.inCity$DENVAL,err_sigma^2) # peak SNR
}
proc.time()-ptm

source("PlotLSfit.R")
save.image(file="MeshImRepFit.RData")