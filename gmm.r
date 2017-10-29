sampleN <- function(size, mus=NULL, f, weigthUpdateMethod="exp") {
	if (length(mus) >  0) {
		fmus <- f(mus)
		fMax <- max(fmus)
		fMin <- min(fmus)
		if (weigthUpdateMethod=="flat")
			wts <- fmus/sum(fmus)
		if (weigthUpdateMethod=="exp") {
			wts <- exp(fmus-fMax)
			wts <- wts/sum(wts)
		}
		#if (length(mus) > 1)
		#	sds <- 1 - ((fmus - fMin)/(fMax - fMin)) + 1
		#else
		#	sds <- 2
		sds <- rep(5,length(mus))
		components <- sample(1:length(mus),prob=wts,size=size,replace=TRUE)
		samples <- rnorm(n=size,mean=mus[components],sd=sds[components])
	}
	samples
}

gaussMixFunc <- function(x, wts, mus, sds) {
        l <- length(wts)
        y <- 0
        for (i in 1:l)
                y <- y + wts[i]*dnorm(x, mus[i], sds[i])
        y
}

rastrigin <- function(x) {
        return( (10 + ( x^2 - 10*cos(2*pi*x) ) )*(-1) + 20 )
}

StyblinskiTang <- function(x) {
        return( ( x**4 - 16*x**2 + 5*x )*(-0.5)  )
}

Ackeley <- function(x) {
        return( (-20*exp(-0.2*sqrt(0.5*(x^2))) - exp(0.5*(cos(2*pi*x))) + exp(1) + 20)*(-1) + 20 )
}

Levi <- function(x) {
        return( ( (sin(3*pi*x))**2 + (x-1)^2 + 1 ) *(-1) + 20  )
}

Eggholder <- function(x) {
        return( (-47*sin(sqrt(abs(x/2 + 47))) - x*sin(sqrt(abs(x+47))) )*(-1) )
}

gmmOpt <- function(func, iter, trials, display=FALSE) {
	samples <- runif(trials,5,10)
	mus <- NULL
	for (i in 1:iter) {
		funcN <- func(samples)
		maxN <- samples[ which(funcN==max(funcN))[1] ]
		mus <- c(mus, maxN)
		samples <- sampleN(size=trials, mus, f=func, weigthUpdateMethod="exp")
		if (display==TRUE) {
			x <- seq(-10,10,0.05)
			plot(x, func(x), type='l', ylim=c(0,50), xlim=c(-10,10))
			fmus <- func(mus)
	                fMax <- max(fmus)
        	        fMin <- min(fmus)
                        wts <- exp(fmus-fMax)
                        wts <- wts/sum(wts)
			#if (length(mus) > 1)
                        #	sds <- 1 - ((fmus - fMin)/(fMax - fMin)) + 0.5
                	#else
                        #	sds <- 1.5
			sds <- rep(5,length(mus))
			lines(x, gaussMixFunc(x, wts, mus, sds))
			abline(v=maxN)
			readline()
		}
	}
	return ( samples[ which(funcN==max(funcN))[1] ] )
}

f1 <- StyblinskiTang 
answer <- gmmOpt(func=f1, iter=30, trials=500, display=TRUE)


