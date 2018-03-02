source("multiVariateGMM.r" )

rastrigin <- function(x) { #x is a trials-by-dim matrix
        A = 10
        n = ncol(x)
        return ( apply( x, 1, function(i)  A*n + sum(i^2 - A*cos(2*pi*i))  )  )
}

limMat <- matrix(c(-5,-5,5,5),2,2)
sol <- gmmOpt(func=rastrigin, iter=100, trials=1000, dim = 2, lim = limMat, display = TRUE, wtUpdateMethod="exp") 

print(sol)
