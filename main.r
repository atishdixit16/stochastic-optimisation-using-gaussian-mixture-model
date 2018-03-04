source("multiVariateGMM.r" )

rastrigin <- function(x) { #x is a trials-by-dim matrix
        A = 10
        n = ncol(x)
        return ( apply( x, 1, function(i)  A*n + sum(i^2 - A*cos(2*pi*i))  )  )
}


dimN = 4
limMat <- matrix(c(-5,5),dimN,2,byrow=T)
sol <- gmmOpt(func=rastrigin, iter=200, trials=1000, dim = dimN, lim = limMat, display = TRUE, wtUpdateMethod="exp", StdDevWtUpdate='flat') 

print(sol)
