source("multiVariateGMM.r" )

rastrigin <- function(x) { #x is a trials-by-dim matrix
        A = 10
        n = ncol(x)
        return ( apply( x, 1, function(i)  A*n + sum(i^2 - A*cos(2*pi*i))  )  )
}

rosenbrock <- function(x) { #x is a trials-by-dim matrix
	n <- ncol(x)
	return ( apply( x , 1 , function(i) sum (100*(i[1:(n-1)]^2 - i[2:n])^2 + (i[1:(n-1)] - 1)^2) ) )
}

dimN = 2
limMat <- matrix(c(-5,5),dimN,2,byrow=T)
sol <- gmmOpt(func=rosenbrock, iter=300, trials=2000, dim = dimN, lim = limMat, display = TRUE, wtUpdateMethod="exp", StdDevWtUpdate='flat') 

print(sol)
