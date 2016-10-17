#source("/Users/dgy/Desktop/SDS385/R/Exercise2.r");
#main(X,y,m,beta);



#data preparation
N<- 2000;
P<- 50;
X<- matrix(rnorm(N * P), nrow=N);
beta<- matrix(rnorm(P),ncol=1);
pr<- 1/(1 + exp(-X %*% beta));
m<- matrix(rep(1,N),ncol=1);
y<- matrix(rbinom(N,1,pr),ncol=1);

#data preparation end

#weight is the used in computing the exponentially weighted moving average
main<- function(x,y,m,beta,round = 20000, SGDstep = 0.001, batchsize = 1, weight = 0.01){
	
	SGDResult <- SGD(x,y,m,round,SGDstep, batchsize, weight);
	#jpeg('~/Desktop/betaVSbetahat.jpg')
	#plotBetaVSBetahat(beta,SGDResult$betahat)
	#dev.off()
	
	#jpeg('~/Desktop/nll.jpg')
	plotNLL(round, SGDResult);
    #dev.off()
}

plotBetaVSBetahat <- function(beta,betahat){
	plot(beta, betahat, xlab="estimated beta", ylab="ground truth beta");
}

plotNLL <- function(round, SGDResult){
	
	p.x <- c(1:round);#x coordinates in the plot
	plot(p.x,SGDResult$nll1, col="red", xlab="round", ylab="negative log likelihood", main = "green:running average, red:likelihood based on all data, blue: exponentially weighted running average");

	lines(p.x, SGDResult$nll2, col="green");
	lines(p.x, SGDResult$nll3, col="blue");
}


SGD<-function(x, y, m, round, step, batchsize, weight){
	#initialize b as a vector of 1s
	b = matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );

	result.nll1<-c();
	result.nll2<-c();
	result.nll3<-c();

	for(r in 1:round){
		batch <- sample(1:nrow(x), batchsize);
		batchx <- matrix(x[batch,],nrow=batchsize);
		batchy <- matrix(y[batch,],nrow=batchsize);
		batchm <- matrix(m[batch,],nrow=batchsize);

		g = gradient( batchx, batchy, b, batchm);
		b = b - step*g;
		
		result.nll1<-c( result.nll1, nll(x,y,b,m) );

		currentnll <- nll(batchx, batchy, b, batchm);

		if(r > 1){
			runningAverage <- result.nll2[r-1] * (r-1)/r + currentnll/r;
			result.nll3 <- c(result.nll3, weight * currentnll + (1-weight) * result.nll3[r-1]);
		}
		else{
		    runningAverage <- currentnll;
			result.nll3 <- c(result.nll3, currentnll);
		}
		result.nll2 <- c(result.nll2, runningAverage);

		
	}
	result <- list("nll1" = result.nll1,"nll2" = nrow(x) * result.nll2,"nll3" = nrow(x) * result.nll3, "betahat" = b);
	return(result);
}




#compute negative log likelihood, b is beta, x is a N by P matrix, y and m are N by 1 matrices
nll<-function(x,y,b,m){
	l=0;
	l=l+ sum( log( choose(m,y) ) );
	l=l+ sum( (x %*% b) * y );
	l=l+ sum( m*log( ( 1/( 1+exp(x %*% b) ) ) ) );
	return(-l);
}

#compute gradient
gradient<-function(x,y,b,m){
	g<- t(x) %*% ( m/( 1+exp(-x %*% b) )-y );
	return(g);
}







#Question D
#SGD with Robbins-Monro rule
SGDRMStep<-function(x, y, m, round, C,alpha,tzero, batchsize, weight){
	#initialize b as a vector of 0s
	b = matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );

	result.nll1<-c();
	result.nll2<-c();
	result.nll3<-c();

	for(r in 1:round){
		
		
		batch <- sample(1:nrow(x), batchsize);
		batchx <- matrix(x[batch,],nrow=batchsize);
		batchy <- matrix(y[batch,],nrow=batchsize);
		batchm <- matrix(m[batch,],nrow=batchsize);

		g = gradient( batchx, batchy, b, batchm);
		
		step = (r+tzero)^(-alpha)*C;
		b = b - step*g;
		
		result.nll1<-c( result.nll1, nll(x,y,b,m) );

		currentnll <- nll(batchx, batchy, b, batchm);

		if(r > 1){
			runningAverage <- result.nll2[r-1] * (r-1)/r + currentnll/r;
			result.nll3 <- c(result.nll3, weight * currentnll + (1-weight) * result.nll3[r-1]);
		}
		else{
		    runningAverage <- currentnll;
			result.nll3 <- c(result.nll3, currentnll);
		}
		result.nll2 <- c(result.nll2, runningAverage);

		
	}
	result <- list("nll1" = result.nll1,"nll2" = nrow(x) * result.nll2,"nll3" = nrow(x) * result.nll3, "betahat" = b);
	return(result);
}

#testing different C and alpha, with SGD of fixed step size as reference
main2<- function(x,y,m,beta,round = 20000, tzero=1,C,alpha, batchsize = 1, weight = 0.5){
	
	SGDResult1 <- SGDRMStep(x,y,m,round,C,alpha,tzero, batchsize, weight);
	SGDResult2 <- SGD(x,y,m,round, 0.01, batchsize, weight);
	
	#plotBetaVSBetahat(beta,SGDResult$betahat);
	jpeg('~/Desktop/rmsgd.jpg')
	plot(1:round,SGDResult2$nll1, col="red", xlab="round", ylab="negative log likelihood", main = "green:RMSGD, red:SGD with fixed step size 0.01");
	lines(1:round,SGDResult1$nll1, col="green");
	dev.off()

}
# seems C=1 and alpha = 0.6 and t0 = 1 is a great combination


#Question E
# SGD with Polyak–Ruppert averaging
SGDPR<-function(x, y, m, round, burnIn, C,alpha,tzero, batchsize, weight){
	#initialize b as a vector of 0s
	b = matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );
	result.nll1<-c();
	result.nll2<-c();
	result.nll3<-c();
	
	for(r in 1:burnIn){
		batch <- sample(1:nrow(x), batchsize);
		batchx <- matrix(x[batch,],nrow=batchsize);
		batchy <- matrix(y[batch,],nrow=batchsize);
		batchm <- matrix(m[batch,],nrow=batchsize);

		g = gradient( batchx, batchy, b, batchm);
		
		step = (r+tzero)^(-alpha)*C;
		b = b - step*g;		
		
		result.nll1 <- c(result.nll1, 400);
		result.nll3 <- c(result.nll3, 400);
		result.nll2 <- c(result.nll2, 400);
	}
	
	bAverage = b;

	
	for(r in 1:round){
		
		
		batch <- sample(1:nrow(x), batchsize);
		batchx <- matrix(x[batch,],nrow=batchsize);
		batchy <- matrix(y[batch,],nrow=batchsize);
		batchm <- matrix(m[batch,],nrow=batchsize);

		g = gradient( batchx, batchy, b, batchm);
		
		step = (r+burnIn+tzero)^(-alpha)*C;
		b = b - step*g;
		
		bAverage = b/(r+1) + bAverage*r/(r+1);
		result.nll1<-c( result.nll1, nll(x,y, bAverage,m) );

		currentnll <- nll(batchx, batchy, bAverage, batchm);

		if(r > 1){
			runningAverage <- result.nll2[r+burnIn-1] * (r-1)/r + currentnll/r;
			result.nll3 <- c(result.nll3, weight * currentnll + (1-weight) * result.nll3[r-1]);
		}
		else{
		    runningAverage <- currentnll;
			result.nll3 <- c(result.nll3, currentnll);
		}
		result.nll2 <- c(result.nll2, runningAverage);

		
	}
	result <- list("nll1" = result.nll1,"nll2" = nrow(x) * result.nll2,"nll3" = nrow(x) * result.nll3, "betahat" = bAverage);
	return(result);
}


main3<- function(x,y,m,beta,round = 10000, burnIn=10000, tzero=1,C=1,alpha=0.6, batchsize = 1, weight = 0.5){
	
	SGDResult1 <- SGDPR(x,y,m,round, burnIn, C,alpha,tzero, batchsize, weight);
	SGDResult2 <- SGDRMStep(x,y,m,round + burnIn,C,alpha, tzero, batchsize, weight);
	#plotBetaVSBetahat(beta,SGDResult$betahat);
	
	jpeg('~/Desktop/prsgd.jpg')
	plot(1:(round+burnIn),SGDResult1$nll1, col="red", xlab="round", ylab="negative log likelihood", main = "green:standard SGD(with RM step), red: PR averaging + SGD(with RM step)");
	lines(1:(round+burnIn),SGDResult2$nll1,col="green");
	dev.off()
}
#Polyak–Ruppert averaging makes the likelihood much more stable 
