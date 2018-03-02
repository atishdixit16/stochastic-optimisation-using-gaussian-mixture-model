gmmOpt <- function(func, iter, trials, dim ,  lim , display, wtUpdateMethod="flat") {
        # samples is a trials-by-dim matrix
        samples <- NULL
        for (i in 1:trials) {
		oneSample <- NULL
		for (j in 1:dim)
			oneSample <- c( oneSample, runif(1,lim[j,1],lim[j,2]) )
                samples <- rbind(samples, oneSample)
	}
        mus <- NULL
        minima = Inf
        epsilon = 1e-5
        for (i in 1:iter) {
		if (display==TRUE) {
			if (dim==2) {
                        	x <- seq(lim[1,1],lim[1,2],0.1)
                        	y <- seq(lim[2,1],lim[2,2],0.1)
	                        lx <- length(x)
	                        ly <- length(y)
        	                mat <- matrix(0,lx, ly)
                	        for (i in 1:lx)
                        	        for (j in 1:ly)
                                	        mat[i,j] = func(t(matrix(c(x[i],y[j]))))
	                        image(x,y,mat)
        	                points(samples[,1],samples[,2],pch='+')
			}
                }

                funcN <- func(samples)
                optN <- samples[ which(funcN==min(funcN))[1] ,  ]
                mus <- rbind(mus, optN)
                samples <- GaussMixture( size=trials, mus, f=func, wtUpdateMethod)
                if ( abs(minima -  min(funcN)) <= epsilon )
                        break
                else
                        minima = min(funcN)
        }
        return ( samples[ which(funcN==min(funcN))[1] ] )
}

		       
library(MASS)
GaussMixture <- function(size, mus, f, weigthUpdateMethod="flat") {
        fmus <- f(mus)
        fMax <- max(fmus)
        fMin <- min(fmus)
        if (weigthUpdateMethod=="flat")
                wts <- sum(fmus)/fmus
        if (weigthUpdateMethod=="exp") 
                wts <- exp(fMin - fmus)
        wts =  wts/sum(wts) 
        n = ncol(mus)
        components <- sample(1:nrow(mus),prob=wts,size=size,replace=TRUE)
        samples <- NULL
        for (i in 1:size)
                samples <- rbind( samples, mvrnorm(n=1,mu=mus[components[i], ],Sigma=diag(n)) )
        samples
}

