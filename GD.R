#source("/Users/dgy/Desktop/SDS385/R/GD.r");

#data preparation
raw<- read.csv("/Users/dgy/Desktop/SDS385/R/wdbc.csv",header=F);
raw_matrix<- as.matrix(raw);
y<- matrix(ifelse(raw_matrix[,2]=="M",1,0),ncol=1);

#use the first ten attributes of the data provided, and scale it
x<- scale(apply(raw_matrix[,3:12],2,as.numeric));

#append dumb 1s to simulate intercept
x<- cbind(x,rep(1,length(y)));

m<- matrix(rep(1,569),nrow=length(y));
#data preparation end

#x,y,m are matrices
main<- function(x,y,m,round=200,GDstep=0.01,Newtonstep=1){

	p.x<-c(1:round);#x coordinates in the plot
	
	GDResult.nll<-GD(x,y,m,round,GDstep);
	NewtonResult.nll<-Newton(x,y,m,round,Newtonstep);
	
	plot(p.x,GDResult.nll,col="red",xlab="round",ylab="negative log likelihood",main="green: newton's method, red: gradient descent");
	lines(p.x,NewtonResult.nll,col="green");

}


GD<-function(x,y,m,round,step){
	#initialize b as a vector of 1s
	b=matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );
	result.nll<-c();
	
	#beta=0.8; for backtracking line search
	for(r in 1:round){
		
		g= gradient( x,y,b,m );
		b= b-step*g;
		result.nll<-c( result.nll,nll(x,y,b,m) );
		
		#print(g);
		#print(nll(x,y,b,m));
		#if(negativeL(x,y,b-g,m)> nl- step/2*sum(g*g)){
		#	step=beta*step;
		#}
	}
	#print(b);
	return(result.nll);
}

Newton<-function(x,y,m,round,step){
	#initialize b as a vector of 1s
	b=matrix( rep( 0,dim(x)[2] ),nrow=dim(x)[2] );
	result.nll<-c();
	
	#beta=0.8; for backtracking line search
	for(r in 1:round){
		
		g= ngradient(x,y,b,m);
		b= b-step*g;
		result.nll<-c(result.nll,nll(x,y,b,m));
		
		#print(g);
		#print(nll(x,y,b,m));
		#if(negativeL(x,y,b-g,m)> nl- step/2*sum(g*g)){
		#	step=beta*step;
		#}
	}
	#print(b);
	return(result.nll);
}

#compute negative log likelihood 
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


#compute newton's gradient, inverse hessian matrix times gradient, for the complete formula, please refer to question c
ngradient<-function(x,y,b,m){
	
	temp<- exp(-x %*% b);#temp is a vector
	d<- diag( c( m* temp/(1+temp)^2 ) ); 
	hessian<- t(x) %*% d %*% x;
	g<- solve(hessian,gradient(x,y,b,m));
	return(g);
	
}


