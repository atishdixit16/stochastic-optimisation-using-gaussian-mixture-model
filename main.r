source("multiVariateGMM.r" )

rastrigin <- function(x) { #x is a trials-by-dim matrix
        A = 10
        n = ncol(x)
        return ( apply( x, 1, function(i)  A*n + sum(i^2 + A*cos(2*pi*i))  )  )
}

sol <- gmmOpt(func=rastrigin, 100, 50, dim = 2, lim = c(-5,5), display = TRUE) 

print(sol)
