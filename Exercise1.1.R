#source("/Users/dgy/Desktop/SDS385/R/Exercise1.1.r")
library("Matrix");
library("microbenchmark");
X<- Matrix(1:8,nrow=4,ncol=2)
w<- c(1,1,1,1)
y<- c(1,5,6,9)


InversionMethod<- function(A,b){
	temp1<- solve(A);
	return(temp1 %*% b);
}

QRMethod<- function(A,b){
	temp1<- qr.solve(A,b);
	return(temp1);
}

CholeskyMethod<- function(A,b){

	R= chol(A);
	temp1<- forwardsolve(t(R),b);
	return(backsolve(R,temp1));
}
Simulation <- function(N,P){
	X<- matrix(rnorm(N*P),nrow=N);
	y<- rnorm(N);
	w<- diag(rep(1,N));
	temp1<-crossprod(X,w);
	A<- temp1 %*% X;
	b<- temp1 %*% y;
	print(microbenchmark(InversionMethod(A,b)));
	print(microbenchmark(CholeskyMethod(A,b)));
	print(microbenchmark(QRMethod(A,b)));
}










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
	
	print(microbenchmark(InversionMethod(A,b)));

	print(microbenchmark(SparseCholeskyMethod(As,b)));

	print(microbenchmark(SparseCholeskyMethod(A,b)));
}