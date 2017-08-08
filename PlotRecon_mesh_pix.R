jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# PSNR
cl <- rainbow(length(SmParamSet)+1)
yRange <- c(min(min(Mesh.MAP_PSNR,na.rm=T),min(Pix.MAP_PSNR,na.rm=T),min(Mesh.ML_PSNR),min(Pix.ML_PSNR))
            -0.1*max(diff(range(Mesh.MAP_PSNR,na.rm=T)),diff(range(Pix.MAP_PSNR,na.rm=T))),
            max(max(Mesh.MAP_PSNR,na.rm=T),max(Pix.MAP_PSNR,na.rm=T),max(Mesh.ML_PSNR),max(Pix.ML_PSNR))
            +0.1*max(diff(range(Mesh.MAP_PSNR,na.rm=T)),diff(range(Pix.MAP_PSNR,na.rm=T))))
# Expand right side of clipping rect to make room for the legend
par(pty="m",mar=par()$mar+c(0,0,0,4))

plot(NumNode,Mesh.ML_PSNR,type="b",col=cl[1],pch=16,lty=1,cex=1,xlab="Number of nodes/pixels",ylab="PSNR",
     ylim=yRange,axes=FALSE,cex.lab=1.3)
lines(NumNode,Pix.ML_PSNR,type="b",col=cl[1],pch=17,lty=2,cex=1)
for (j in 1:length(SmParamSet)){
  lines(NumNode,Mesh.MAP_PSNR[,j],type="b",col=cl[j+1],pch=16,lty=1,cex=1)
  lines(NumNode,Pix.MAP_PSNR[,j],type="b",col=cl[j+1],pch=17,lty=2,cex=1)
}

axis(side=1,at=NumNode,lwd.ticks=0.7,cex.axis=0.7)
axis(side=2,lwd.ticks=0.7,cex.axis=0.7)
box()
grid(nx=round(1.1*length(NumNode)), ny=NULL, col="lightgray", lty="dotted",lwd=1)
par(new=TRUE)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n",xlab="",ylab="")
legend("topright",legend=c("Mesh-ML",paste0("Mesh-MAP",as.character(1:length(SmParamSet))),"Pix-ML",paste0("Pix-MAP",as.character(1:length(SmParamSet)))),
       col=rep(cl,2), pch=rep(c(16,17),each=length(cl)),lty=rep(c(1,2),each=length(cl)),lwd=1,cex=0.4,inset=c(-0.2,0),xpd=TRUE)

# ----------RMSE-----------------

# Restore default clipping rect
par(mar=c(5, 4, 4, 2) + 0.1)

yRange <- 100*c(min(min(Mesh.MAP_repErr,na.rm=T),min(Pix.MAP_repErr,na.rm=T),min(Mesh.ML_repErr),min(Pix.ML_repErr))
                -0.1*max(diff(range(Mesh.MAP_repErr,na.rm=T)),diff(range(Pix.MAP_repErr,na.rm=T))),
                max(max(Mesh.MAP_repErr,na.rm=T),max(Pix.MAP_repErr,na.rm=T),max(Mesh.ML_repErr),max(Pix.ML_repErr))
                +0.1*max(diff(range(Mesh.MAP_repErr,na.rm=T)),diff(range(Pix.MAP_repErr,na.rm=T))))
par(pty="m",mar=par()$mar+c(0,0,0,4))
plot(NumNode,100*Mesh.ML_repErr,type="b",col=cl[1],pch=16,lty=1,cex=1,xlab="Number of nodes/pixels",ylab="Representation Error (RMSE%)",
     ylim=yRange,axes=FALSE,cex.lab=1.3)
lines(NumNode,100*Pix.ML_repErr,type="b",col=cl[1],pch=17,lty=2,cex=1)
for (j in 1:length(SmParamSet)){
  lines(NumNode,100*Mesh.MAP_repErr[,j],type="b",col=cl[j+1],pch=16,lty=1,cex=1)
  lines(NumNode,100*Pix.MAP_repErr[,j],type="b",col=cl[j+1],pch=17,lty=2,cex=1)
}

axis(side=1,at=NumNode,lwd.ticks=0.7,cex.axis=0.7)
axis(side=2,lwd.ticks=0.7,cex.axis=0.7)
box()
grid(nx=round(1.1*length(NumNode)), ny=NULL, col="lightgray", lty="dotted",lwd=1)
par(new=TRUE)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n",xlab="",ylab="")
legend("topright",legend=c("Mesh-ML",paste0("Mesh-MAP",as.character(1:length(SmParamSet))),"Pix-ML",paste0("Pix-MAP",as.character(1:length(SmParamSet)))),
       col=rep(cl,2), pch=rep(c(16,17),each=length(cl)),lty=rep(c(1,2),each=length(cl)),lwd=1,cex=0.4,inset=c(-0.2,0),xpd=TRUE)
par(mar=c(5, 4, 4, 2) + 0.1)