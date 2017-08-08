jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

par(pty="s")
plot(NumNode,Mesh.MDL,type="b",col="red",pch=16,lty=1,lwd=1.5,cex=1,xlab="Number of nodes",ylab="MDL",axes=FALSE,cex.lab=1.3)
axis(side=1,at=NumNode)
axis(2)
box()
grid(nx=NULL, ny=NULL, col="lightgray", lty="dotted",lwd=2)

plot(NumNode,Mesh.PSNR,type="b",col="red",pch=16,lty=1,cex=0.8,xlab="Number of nodes",ylab="PSNR (dB)",axes=FALSE)
axis(side=1,at=NumNode)
axis(2)
box()
grid(nx=NULL, ny=NULL, col="lightgray", lty="dotted",lwd=2)

KDE.inCity_raster <- rasterize(KDE.inCity[,c("X_COORD","Y_COORD")],r,KDE.inCity$DENVAL,fun=sum)
par(mfrow=c(1,2),oma=c(0,1.5,0,1.5), mar=c(5.1,5,4.1,4))
# plot.new()
plot(KDE.inCity_raster,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
plot(city.shp, border="black",add=TRUE)

KDE_hat.city_raster <- rasterize(KDE_hat.city[,c("X_COORD","Y_COORD")],r,KDE_hat.city$VALUE,fun=sum)
plot(KDE_hat.city_raster,col=jet.colors(256),xlab="X coordinate",ylab="Y coordinate")
plot(city.shp, border="black",add=TRUE)