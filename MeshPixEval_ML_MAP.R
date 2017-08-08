setwd("/Users/xiaomuliu/CrimeProject/SpatioTemporalModeling/MeshModeling/")
source("MeshNodeEval_ML_MAP.R")

# Compare with pixel-based ML reconstruction
# Theoratically, pix-ML estimates are the actual counts of corresponding pixel cells.

# setup X-Y grid ratio
XYratio <- (X_range[2]-X_range[1])/(Y_range[2]-Y_range[1])
NumPix <- NumNode
Pix.ML_PSNR <- rep(0,length(NumPix))
Pix.ML_repErr <- rep(0,length(NumPix))

for (i in 1:length(NumPix)){
  # the ratio of the coarse grid (with number of pixels for evaluation) over the reference (fine) grid
  scaleRatio <- sqrt(NumPix[i]/nrow(RegGrd))
  
  PixRaster <- raster(ncol=round(scaleRatio*grd.full@grid@cells.dim[1]),nrow=round(scaleRatio*grd.full@grid@cells.dim[2]),
                      xmn=grd.full@bbox[1,1],xmx=grd.full@bbox[1,2],ymn=grd.full@bbox[2,1],ymx=grd.full@bbox[2,2],crs=CRS(prj)) 
  CrimeObsPixRaster <- rasterize(CrimeObsPts.df_inCity[,c("X_COORD","Y_COORD")],PixRaster,CrimeObsPts.df_inCity$INC_CNT,fun=sum)
  
  CrimeObsPixProj <- resample(CrimeObsPixRaster,r,method='bilinear')  # disaggregate to resolution of raster r
  # CrimeObsPixProj <- disaggregate(CrimeObsPixRaster, round(1/scaleRatio), method='bilinear')  # disaggregate to resolution of raster r
  
  CrimeObsPix.df_inCity <- as.data.frame(CrimeObsPixProj,xy=TRUE)[isInCity,]
  names(CrimeObsPix.df_inCity) <- c("X_COORD","Y_COORD","VALUE")
  CrimeObsPix.df_inCity$VALUE[is.na(CrimeObsPix.df_inCity$VALUE)] <- 0
  # # The summation of values in disaggregated cells should be equal to the value in the bigger cell
  # # This step is unnecessary as the density is calculated in the end
  # CrimeObsPix.df_inCity$VALUE <- CrimeObsPix.df_inCity$VALUE * scaleRatio
  
  # re-arrange the order according to KDE grid
  CrimeObsPix.df_inCity <- CrimeObsPix.df_inCity[ord,]

  CrimeObsPix.df_inCity$DENVAL <- CrimeObsPix.df_inCity$VALUE/sum(CrimeObsPix.df_inCity$VALUE) 

  # Evaluation
  err_sigma <- sqrt(MSE(KDE.Hist_inCity$DENVAL,CrimeObsPix.df_inCity$DENVAL)) 
  Pix.ML_PSNR[i] <- PSNR(KDE.Hist_inCity$DENVAL,err_sigma^2) # peak SNR
  Pix.ML_repErr[i] <- err_sigma/sqrt(MSE(KDE.Hist_inCity$DENVAL,0)) # representation error (percent RMSE)
}

source("PlotRecon_mesh_pix.R")
