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

linSys2 <- function(x) { #x is a trials-by-2 matrix (answer=4,1)
	stopifnot(ncol(x) == 2 )
	A <- matrix(c(2,1,-1,1),2,2)
	b <- matrix(c(7,5), 2,1)
	b_try <- apply(x,1,function(i) A%*%i )
	ans <- apply(b_try,2, function(i) sqrt(sum((i-b)**2)/2)  )  #root mean square as objective function
	return ( ans )
}

dimN = 2
limMat <- matrix(c(-5,5),dimN,2,byrow=T)
sol <- gmmOpt(func=linSys2, iter=300, trials=2000, dim = dimN, lim = limMat, display = TRUE, wtUpdateMethod="exp", StdDevWtUpdate='flat') 

print(sol)
