MSE <- function(x,x_hat){return(mean((x-x_hat)^2))}

PSNR <- function(x,mse)10*log10(max(x)/mse)

MDL <- function(NumPix,NumNode,err_sigma){return (NumPix/2*log(2*pi*err_sigma^2)+NumPix/2+(NumNode+1)/2*log(NumPix))}
