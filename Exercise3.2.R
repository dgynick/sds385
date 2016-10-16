#source("/Users/dgy/Desktop/SDS385/R/Exercise3.2.r");

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
main<- function(x, y, m, beta, round = 2000, weight = 0.5){
	
	SGDResult1 <- SGDBacktracking(x, y, m, round, c = 0.8, p = 0.8, initialStep = 1, batchsize = 1, weight);
	SGDResult2 <- QNewton(x, y, m, round, c = 0.8, p = 0.8, initialStep = 1, weight);
	SGDResult3 <- Newton(x, y, m, round, weight);
	
	p.x <- c(1:round);#x coordinates in the plot
	plot(p.x,SGDResult1$nll1, col="red", xlab="round", ylab="negative log likelihood", main = "red: SGD, green: Quasi-Newton, blue: Newton");
	lines(p.x, SGDResult2$nll1, col="green");
	lines(p.x, SGDResult3$nll1, col="blue");
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
	g<- crossprod(x, ( m/( 1+exp(-x %*% b) )-y ));
	return(g);
}

ngradient<-function(x,y,b,m){
	
	temp <- exp(-x %*% b);#temp is a vector
	d <- diag( c( m* temp/(1+temp)^2 ), nrow = nrow(x), ncol = nrow(x) ); 

	hessian <-  t(x) %*% d %*% x;
	return(hessian);
}

updateQngradient <- function(s,y,H){
	rho = 1/crossprod(y,s)[1,1];
	I = diag(nrow(H));

	updatedH = (I - rho * (s %*% t(y))) %*% H %*% (I - rho * s %*% t(y)) + rho * crossprod(t(s));
	return(updatedH);
}

SGDBacktracking<-function(x, y, m, round, c, p, initialStep, batchsize, weight){
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

		g = gradient(batchx, batchy, b, batchm);
		
		step = initialStep;
		while( nll(batchx, batchy, b-step*g, batchm) > nll(batchx, batchy, b, batchm) - c * step * crossprod(g)){
			step = step * p;
		}
		b = b - step * g;
		
		
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

QNewton<-function(x, y, m, round,  c, p, initialStep, weight){
	#initialize b as a vector of 1s
	b = matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );

	result.nll1<-c();
	result.nll2<-c();
	result.nll3<-c();
	
	#initialize the estimate of Hessian as I
	H = diag(ncol(x));
	
	for(r in 1:round){
		
		g = gradient(x, y, b, m);
		direction = H %*% g;
		
		step = initialStep;
		while( nll(x, y, b - step* direction, m) > nll(x, y, b, m) - c * step * crossprod(direction, g)){
			step = step * p;
		}
		
		H = updateQngradient(-step * direction, gradient(x, y, b - step * direction, m) - g, H);
		b = b - step * direction;
		
		
		result.nll1<-c( result.nll1, nll(x,y,b,m) );

		currentnll <- nll(x, y, b, m);

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



Newton<-function(x, y, m, round, weight, step = 1){
	#initialize b as a vector of 1s
	b = matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );

	result.nll1<-c();
	result.nll2<-c();
	result.nll3<-c();
	
	
	for(r in 1:round){

		
		g = gradient(x, y, b, m);
		direction = solve(ngradient(x, y, b, m), g);
		
		b = b - step * direction;
		
		
		result.nll1<-c( result.nll1, nll(x,y,b,m) );

		currentnll <- nll(x, y, b, m);

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