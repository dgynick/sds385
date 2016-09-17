
library("Matrix");
library("microbenchmark");


#multiply A by a diagnal vector w
multdiag<-function(A,w){
	for(i in 1:length(w)){
		A[,i]=w[i]*A[,i];
	}
	return(A);
}


InversionMethod<- function(A,b){
	temp1<- solve(A);
	return(temp1 %*% b);
}

CholeskyMethod<- function(A,b){

	R= chol(A);
	temp1<- forwardsolve(t(R),b);
	return(backsolve(R,temp1));
}

#N is number of observations, p is the length of each observation
Simulation <- function(N,P){
	X<- matrix(rnorm(N*P),nrow=N);
	y<- rnorm(N);
	w<- rep(1,N);
	temp1<-multdiag(t(X),w);
	A<- temp1 %*% X;
	b<- temp1 %*% y;
	print(microbenchmark(InversionMethod(A,b)));
	print(microbenchmark(CholeskyMethod(A,b)));
}









SparseInversionMethod<- function(A,b){
	temp2<- solve(A);
	return(temp2 %*% b);
}

#This turned out to work not so well even in sparse case
SparseCholeskyMethod2<- function(A,b){

	temp3<- expand(Cholesky(A));
	
	R<- t(temp3$L) %*% (temp3$P);
	temp4<- forwardsolve(t(R),b);
	return(backsolve(R,temp4));
}

SparseCholeskyMethod<- function(A,b){
	R= chol(A);
	temp2<- forwardsolve(t(R),b);
	return(backsolve(R,temp2));
}




#s is the probability of a sparse cell in X
SparseSimulation <- function(N,P,s){
	X<- matrix(rnorm(N*P),nrow=N);
	mask<- matrix(rbinom(N*P,1,s),nrow=N);
	Xs<- Matrix(mask*X,sparse=T);
	X<- Matrix(mask*X);
	y<- Matrix(rnorm(N));
	w<- diag(rep(1,N));

	temp1<- t(Xs)%*%w;
	As<- Matrix(temp1 %*% Xs);
	
	temp2<- t(X)%*%w;
	A<- Matrix(temp2 %*% X);
	
	b<- Matrix(temp2 %*% y);
	
	print(microbenchmark(SparseInversionMethod(A,b)));

	print(microbenchmark(SparseCholeskyMethod(As,b)));

	print(microbenchmark(SparseCholeskyMethod(A,b)));
}
