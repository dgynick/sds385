#source("/Users/dgy/Desktop/SDS385/R/GD.r");
#GD(x,y,m);

#scale
#get around log0 by adding 0.01, or any other way?
#change to matrix speed up the performance

raw<- read.csv("/Users/dgy/Desktop/SDS385/R/wdbc.csv",header=F);
raw_matrix<- as.matrix(raw);

y<- matrix(ifelse(raw_matrix[,2]=="M",1,0),ncol=1);
x<- scale(apply(raw_matrix[,3:12],2,as.numeric));
x<- cbind(x,rep(1,length(y)));


m<- matrix(rep(1,569),nrow=length(y));



#set b=1

GD<- function(x,y,m,round=200,method="g"){
	b=matrix(rep(0,dim(x)[2]),nrow=dim(x)[2]);
	step=0.001;
	beta=0.8;
	
	if(method=="g"){
		for(r in 1:round){
		
			g= gradient(x,y,b,m);
			b= b-step*g;
		
			#print(g);
			print(nll(x,y,b,m));

			#if(negativeL(x,y,b-g,m)> nl- step/2*sum(g*g)){
			#	step=beta*step;
			#}
		}
		print(b);
	}
	
	else if(method=="n"){
		for(r in 1:round){
		
			g= ngradient(x,y,b,m);
			b= b-step*g;
		
			#print(g);
			print(nll(x,y,b,m));

			#if(negativeL(x,y,b-g,m)> nl- step/2*sum(g*g)){
			#	step=beta*step;
			#}
		}
		print(b);
	}
	
}

nll<-function(x,y,b,m){
	l=0;
	l=l+ sum(log(choose(m,y)));
	l=l+ sum((x %*% b)*y);
	l=l+ sum(m*log((1/(1+exp(x%*%b)))));
	return(-l);
}

gradient<-function(x,y,b,m){
	g<- t(x) %*% (m/(1+exp(-x %*% b))-y);
	return(g);
}

ngradient<-function(x,y,b,m){
	temp<- exp(-x%*%b);
	w<- diag(c(m* temp/(1+temp)^2));
	hessian<- t(x)%*%w%*%x;
	g<- solve(hessian) %*% gradient(x,y,b,m);
	return(g);
}


