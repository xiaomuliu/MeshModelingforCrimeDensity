DetectEdge <- function(ImSegment){
  N <- nrow(ImSegment)
  M <- ncol(ImSegment)
  Edge <- matrix(FALSE,nrow=N,ncol=M)
  for (i in 1:N){
    for (j in 1:M){
      if (i==1 | j==1 | i==N | j==M){
        Edge[i,j] <- ImSegment[i,j]==TRUE
      }else{
        Edge[i,j] <- ImSegment[i,j]==TRUE & (ImSegment[i-1,j]==FALSE | ImSegment[i+1,j]==FALSE | 
                                                    ImSegment[i,j-1]==FALSE | ImSegment[i,j+1]==FALSE)
      }                                           
    }
  }
  return(Edge)
}
