cl <- rainbow(length(SmParamSet)+1)
ptChar <- 13:(13+length(SmParamSet))
yRange <- c(min(min(Mesh.MAP_PSNR,na.rm=T),min(Mesh.ML_PSNR))-0.1*diff(range(Mesh.MAP_PSNR,na.rm=T)),
            max(max(Mesh.MAP_PSNR,na.rm=T),max(Mesh.ML_PSNR))+0.1*diff(range(Mesh.MAP_PSNR,na.rm=T)))
par(pty="s")
plot(NumNode,Mesh.ML_PSNR,type="b",col=cl[1],pch=ptChar[1],lty=1,cex=0.8,xlab="Number of nodes",ylab="PSNR",
     ylim=yRange,axes=FALSE)
for (j in 1:length(SmParamSet)){
  lines(NumNode,Mesh.MAP_PSNR[,j],type="b",col=cl[j+1],pch=ptChar[j+1],lty=1,cex=0.8)
}
axis(side=1,at=NumNode)
axis(2)
box()
grid(nx=NULL, ny=NULL, col="lightgray", lty="dotted",lwd=2)
# par(new=TRUE)
# plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n",xlab="",ylab="")
legend("topright",legend=c("ML",paste0("MAP",as.character(1:length(SmParamSet)))),
       col=cl, pch=ptChar,lty=1,lwd=2,cex=1,inset=c(0,0),xpd=TRUE)


yRange <- 100*c(min(min(Mesh.MAP_repErr,na.rm=T),min(Mesh.ML_repErr))-0.1*diff(range(Mesh.MAP_repErr,na.rm=T)),
                max(max(Mesh.MAP_repErr,na.rm=T),max(Mesh.ML_repErr))+0.1*diff(range(Mesh.MAP_repErr,na.rm=T)))
par(pty="s")
plot(NumNode,100*Mesh.ML_repErr,type="b",col=cl[1],pch=ptChar[1],lty=1,cex=0.8,xlab="Number of nodes",ylab="Representation Error (%)",
     ylim=yRange,axes=FALSE)
for (j in 1:length(SmParamSet)){
  lines(NumNode,100*Mesh.MAP_repErr[,j],type="b",col=cl[j+1],pch=ptChar[j+1],lty=1,cex=0.8)
}
axis(side=1,at=NumNode)
axis(2)
box()
grid(nx=NULL, ny=NULL, col="lightgray", lty="dotted",lwd=2)
# par(new=TRUE)
# plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n",xlab="",ylab="")
legend("bottomright",legend=c("ML",paste0("MAP",as.character(1:length(SmParamSet)))),
       col=cl, pch=ptChar,lty=1,lwd=2,cex=1,inset=c(0,0),xpd=TRUE)


CrimeObsPts.inCity_raster <- rasterize(CrimeObsPts.df_inCity[,c("X_COORD","Y_COORD")],r,CrimeObsPts.df_inCity$INC_CNT,fun=sum)
par(mfrow=c(1,2),oma=c(0,1.5,0,1.5), mar=c(5.1,5,4.1,4))
# plot.new()
plot(CrimeObsPts.inCity_raster,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
plot(city.shp, border="black",add=TRUE)

Recon.inCity_raster <- rasterize(Recon.city[,c("X_COORD","Y_COORD")],r,Recon.city$DENVAL,fun=sum)
plot(Recon.inCity_raster,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
plot(city.shp, border="black",add=TRUE)