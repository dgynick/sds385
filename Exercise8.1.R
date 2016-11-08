#source("/Users/dgy/Desktop/SDS385/R/Exercise8.1.r")
require("Matrix")
require(microbenchmark)

makeD2_sparse = function (dim1, dim2)  {

    D1 = bandSparse(dim1 * dim2, m = dim1 * dim2, k = c(0, 1), 
        diagonals = list(rep(-1, dim1 * dim2), rep(1, dim1 * 
            dim2 - 1)))
    D1 = D1[(seq(1, dim1 * dim2)%%dim1) != 0, ]
    D2 = bandSparse(dim1 * dim2 - dim1, m = dim1 * dim2, k = c(0, 
        dim1), diagonals = list(rep(-1, dim1 * dim2), rep(1, 
        dim1 * dim2 - 1)))
    return(rBind(D1, D2))
}

Y = read.csv("/Users/dgy/Desktop/SDS385/fmri/fmri_z.csv", header= FALSE)
Y = as.matrix(Y)
Y = Y[-1,]
Y = data.matrix(apply(Y, 2, as.numeric))

m = nrow(Y)
n = ncol(Y)

D = makeD2_sparse(m, n)
L = crossprod(D)

Y = c(Y)


lambda = 0.1
A = lambda * L + diag(nrow(L))
b = Y


main <- function(){

	print(microbenchmark((result1 = CholeskyMethod(A,b)), times = 1))
	print(microbenchmark((result2 = CG(A,b)), times = 1))
	plot(result1$x, result2$finalx, col = "red", main="Conjugate gradient against cholesky, yellow: after 2 iteration, blue: after 3, red:after 6")
	points(result1$x, result2$finalx,col="red")
	points(result1$x, result2$x[,2],col="yellow")
	points(result1$x, result2$x[,3],col="green")

}

#to be tested
plotHeatMap <- function(){
	#A = lambda * L + diag(nrow(L))
	#b = Y
	
	result2 = CG(A,b)
	image(matrix(result2$finalx, nrow = m))
}




CholeskyMethod<- function(A,b){
	result = list()
	R= chol(A);
	temp1<- forwardsolve(t(R),b);
	result$x = backsolve(R,temp1)
	return(result);
}


test <-function(){
	#for test only
	A = matrix(c(3,4,6,4,6,9,6,9,14), nrow=3)
	b = c(1,2,3)
	
	print(microbenchmark((result1 = CholeskyMethod(A,b)), times = 1))
	print(microbenchmark((result2 = CG(A,b)), times = 1))
	plot(result1$x, result2$finalx, col = "red", main="Conjugate gradient against cholesky, yellow: after 2 iteration, blue: after 3, red:after 6")
	points(result1$x, result2$finalx,col="red")
	points(result1$x, result2$x[,2],col="yellow")
	points(result1$x, result2$x[,3],col="green")
	
}

#conjugate gradient without preconditioning
CG <- function(A,b){
	result = list()
	prevx = rep(0, ncol(A))
	
	#residual Ax-b
	prevr = -b
	
	#the descent direction
	prevp = -prevr
	
	result$x = matrix(prevx, ncol = 1)
	crossprodprevr = crossprod(prevr)
	
	objectives = vector()
	
	for(iteration in 1:5){
		aprevp = A %*% prevp
		
		#step size
		alpha = as.numeric(crossprodprevr/crossprod(prevp, aprevp))

		x = prevx + alpha * prevp
		
		r = prevr + alpha * aprevp
		crossprodr = crossprod(r)
		
		#compute next descent direction
		beta = as.numeric(crossprodr/crossprodprevr)
		p = -r + beta * prevp
        
        prevp = p
        prevr = r
        prevx = x
        crossprodprevr = crossprodr
        
        result$x = cbind(result$x, x)
	    #objectives = c(objectives, as.numeric(crossprod(A %*% x - b)))
	}
	#result$objectives = objectives
	result$finalx = prevx
	return(result)
}
