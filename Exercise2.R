#source("/Users/dgy/Desktop/SDS385/R/Exercise2.r");
#main(X,y,m,beta);



#data preparation
N<- 200;
P<- 5;
X<- matrix(rnorm(N * P), nrow=N);
beta<- matrix(rnorm(P),ncol=1);
pr<- 1/(1 + exp(-X %*% beta));
m<- matrix(rep(1,N),ncol=1);
y<- matrix(rbinom(N,1,pr),ncol=1);

#data preparation end

#weight is the used in computing the exponentially weighted moving average
main<- function(x,y,m,beta,round = 20000, SGDstep = 0.001, batchsize = 1, weight = 0.5){
	
	SGDResult <- SGD(x,y,m,round,SGDstep, batchsize, weight);
	
	plotBetaVSBetahat(beta,SGDResult$betahat);
	
	#plotNLL(round, SGDResult);

}

plotBetaVSBetahat <- function(beta,betahat){
	plot(beta, betahat, xlab="estimated beta", ylab="ground truth beta");
}

plotNLL <- function(round, SGDResult){
	
	p.x <- c(1:round);#x coordinates in the plot
	plot(p.x,SGDResult$nll1, col="red", xlab="round", ylab="negative log likelihood", main = "green:running average likelihood, red:likelihood based on all data, blue: exponentially weighted running average");

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

