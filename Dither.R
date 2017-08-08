library(matlab)
errorDiffusion <- function(I, type=1, threshold){
  #  Error diffusion method for image dithering.
  #  I: The input image
  #  threshold: threshold
  #  type: 1 or 2
  #     type = 1: Floyd & Steinberg filter
  #   	type = 2: Stucki filter
  
  if (type == 1) {
    #  Floyd & Steinberg filter
    Filter <- matrix(c(0,3,0,5,7,1),nrow=2,ncol=3)
    Filter <- 1/sum(colSums(Filter)) * Filter
  }
  else if (type == 2){
    # Stucki filter
    Filter <- matrix(c(0,2,1,0,4,2,0,8,4,8,4,2,4,2,1),nrow=3,ncol=5)
    Filter <- 1/sum(colSums(Filter)) * Filter
  }
  else{
    Output <- rep(0,1)
    return(Output)
  }
  
  fy <- nrow(Filter)
  fx <- ncol(Filter)
  #  Output image, will be filled pixel by pixel
  Output <- matrix(0,nrow=nrow(I),ncol=ncol(I))
  
  #  Zero-padding input image, so that we don't have to deal with edges when
  #  applying the filter. The extra pixels are removed in the output image
  I <- padarray(I, c(fy,fx))
  padx <- fx
  pady <- fy
  
  Iy <- nrow(I)
  Ix <- ncol(I)
  
  # For each pixel of the images (disregarding the zero-padded ones)
  for (y in (pady+1):(Iy-pady)){
    for (x in (padx+1):(Ix-padx)){
      # Threshold image and save to output 
      # (taking into account the zero padding introduced)
      oy <- y - pady
      ox <- x - padx
      Output[oy,ox] <- (I[y,x] > threshold)
      
      # Calculate error
      # error = I[y,x] - Output[oy,ox]
      error <- I[y,x] - 2*threshold*Output[oy,ox]
      
      # Get position where the filter should be applied
      rect <- filterPosition(x, y, fx, fy)
      
      # Distribute error according to the filter
      I[rect$ymin:rect$ymax, rect$xmin:rect$xmax] <- I[rect$ymin:rect$ymax, rect$xmin:rect$xmax] + error * Filter
    }
  }
  
  return(Output)
}


filterPosition <- function(x, y, fx, fy){
  #  Rectangle of image where the filter should be applied, based on the
  #  current processed pixel (x,y)
  k <- floor(fx/2)
  
  ymin <- y
  ymax <- y + fy - 1
  xmin <- x - (fx - k - 1)
  xmax <- x + k
  
  return(list(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
}