gmmOpt <- function(func, iter, trials, dim = 1) {
	lim = c(5,10)
	# samples is a trials-by-dim matrix
	samples <- NULL
	for (i in 1:trials)
		samples <- rbind(samples,  runif(dim,lim[1],lim[2]))
        mus <- NULL
        for (i in 1:iter) {
                funcN <- func(samples)
                maxN <- samples[ which(funcN==max(funcN))[1] ,  ]
                mus <- rbind(mus, maxN)
                samples <- GaussMixture( size=trials, mus, f=func, weigthUpdateMethod="flat")
        }
        return ( samples[ which(funcN==max(funcN))[1] ] )
}

rastrigin <- function(x) { #x is a trials-by-dim matrix
	A = 10
	n = ncol(x)
	return ( apply( x, 1, function(i) A*n + sum(i^2 + A*cos(2*pi*i)) )  )
}

library(MASS)
GaussMixture <- function(size, mus, f, weigthUpdateMethod="flat") {
	fmus <- f(mus)
        fMax <- max(fmus)
        fMin <- min(fmus)
        if (weigthUpdateMethod=="flat")
        	wts <- fmus/sum(fmus)
        if (weigthUpdateMethod=="exp") {
        	wts <- exp(fmus-fMax)
                wts <- wts/sum(wts)
        }
        n = ncol(mus)
        components <- sample(1:nrow(mus),prob=wts,size=size,replace=TRUE)
	samples <- NULL
	for (i in 1:size)
               	samples <- rbind( samples, mvrnorm(n=1,mu=mus[components[i], ],Sigma=diag(n)) )
        samples
}

