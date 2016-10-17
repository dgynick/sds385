#source("/Users/dgy/Desktop/SDS385/R/Exercise4.1.r");
#main(X,y,m,beta);



#data preparation
N<- 2000;
P<- 50;
sparsity <- 0.8
X<- matrix(rnorm(N * P), nrow=N);
mask <- matrix(rbinom(N*P, 1, 1 - sparsity), nrow=N)
X <- X*mask
beta<- matrix(rnorm(P),ncol=1);
pr<- 1/(1 + exp(-X %*% beta));
m<- matrix(rep(1,N),ncol=1);
y<- matrix(rbinom(N,1,pr),ncol=1);

#data preparation end
main<- function(x,y,m,beta,round = 2000, c = 0.8, p = 0.8, initialStep = 1, batchsize = 1, weight = 0.5){
	
	SGDResult1 <- SGDBacktracking(x,y,m,round,c,p,initialStep, batchsize, weight);
	SGDResult2 <- ImprovedSGDBacktracking(x,y,m,round,c,p,initialStep, 3, 15, weight);
	SGDResult3 <- AdaGrad(x,y,m,round,c,p,initialStep, batchsize, weight);
	
	jpeg("~/Desktop/exercise4.1sparse.jpg")
	p.x <- c(1:round);#x coordinates in the plot
	plot(p.x,SGDResult1$nll1, col="red", xlab="round", ylab="negative log likelihood", main = "red:backtrackingSGD, green: improved SGD, blue: adagrad ");

	lines(p.x, SGDResult2$nll1, col="green");
	lines(p.x, SGDResult3$nll1, col="blue");
	dev.off()
}

gradient<-function(x,y,b,m){
	g<- crossprod(x, ( m/( 1+exp(-x %*% b) )-y ));
	return(g);
}

#compute negative log likelihood, b is beta, x is a N by P matrix, y and m are N by 1 matrices
nll<-function(x,y,b,m){
	l=0;
	l=l+ sum( log( choose(m,y) ) );
	l=l+ sum( (x %*% b) * y );
	l=l+ sum( m*log( ( 1/( 1+exp(x %*% b) ) ) ) );
	return(-l)
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

		g = gradient( batchx, batchy, b, batchm);
		
		step = initialStep;
		while(nll(batchx, batchy, b-step*g, batchm) > nll(batchx, batchy, b, batchm) - c * step * crossprod(g)){
			step = step*p;
		}
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


ImprovedSGDBacktracking<-function(x, y, m, round, c, p, initialStep, subsampleSize, newStepFrequency, weight){
	#initialize b as a vector of 1s
	b = matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );

	result.nll1<-c();
	result.nll2<-c();
	result.nll3<-c();
	
	step = 0;
	for(r in 1:round){
		
		if(r %% newStepFrequency == 1){
			batch <- sample(1:nrow(x), subsampleSize);
			batchx <- matrix(x[batch,], nrow = subsampleSize);
			batchy <- matrix(y[batch,], nrow = subsampleSize);
			batchm <- matrix(m[batch,], nrow = subsampleSize);
			g = gradient(batchx, batchy, b, batchm);
			
			step = initialStep;
			while(nll(batchx, batchy, b - step * g, batchm) > nll(batchx, batchy, b, batchm) - c * step * crossprod(g)){
				step = step*p;
			}
			step = step / subsampleSize;#IS THIS CORRECT???????
		}
		
		
		batch <- sample(1:nrow(x), 1);
		batchx <- matrix(x[batch,],nrow=1);
		batchy <- matrix(y[batch,],nrow=1);
		batchm <- matrix(m[batch,],nrow=1);

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

AdaGrad<-function(x, y, m, round, c, p, initialStep, batchsize, weight){
	#initialize b as a vector of 1s
	b = matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );

	result.nll1<-c();
	result.nll2<-c();
	result.nll3<-c();
	
	
	A = matrix(rep(0.001, ncol(x)), nrow = ncol(x));
	
	for(r in 1:round){
		batch <- sample(1:nrow(x), batchsize);
		batchx <- matrix(x[batch,],nrow=batchsize);
		batchy <- matrix(y[batch,],nrow=batchsize);
		batchm <- matrix(m[batch,],nrow=batchsize);

		g = gradient(batchx, batchy, b, batchm);
		A = A + g ^ 2;
		direction = g / sqrt(A);

		step = initialStep;
		while(nll(batchx, batchy, b - step * direction, batchm) > nll(batchx, batchy, b, batchm) - c * step * crossprod(direction, g)){
			step = step*p;
		}
		b = b - step * direction;
		
		
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

