h <- grd.full@grid@cellsize  # kernel bandwidth for KDE
#ENode <- 800 # expected number of mesh nodes
Enode <- 1024
bndyPercent <- 0.1 # the percent of nodes on boundary
SmParam <- 10^(-4) # Gibbs prior parameter
iterMax <- 30 # maximal EM iterations
reltol <- 1e-8 # relevance tolerance for convergence
eps <- 1e-20 # small value to replace zeros in Possion rate
meshHistSpan <- 365*2 # time span of observed data for generating mesh
meshObsSpan <- 365*2