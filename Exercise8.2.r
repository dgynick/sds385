#source("/Users/dgy/Desktop/SDS385/R/Exercise8.2.r")
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

Y = c(Y)

objective <- function(x, Y, lambda, D){
	return(as.numeric(0.5 * crossprod(x - Y) + lambda * sum(abs(D %*% x))))
}


main <- function(){
	microbenchmark(result1 = ADMM(0.1, 10000), times = 1)
	#plot(1:length(result1$objectives), result1$objectives)
	image(result1$x)
}


ADMM <- function(lambda, max_iteration, absThr = 0.0001, relThr = 0.0001, dynamicRho = TRUE){

	rho = 0.1
	
	x = rep(0, length(Y))
	prevv = rep(0, length(Y))
	prevu = prevv
	
	prevs = rep(0, nrow(D))
	prevt = prevs

	#stopping criterion
	primalThr <- function(x, v, r, s){
		return(((length(Y) + nrow(D)) ^ 0.5) * absThr + relThr * (max(as.numeric(crossprod(c(x, r))), as.numeric(crossprod(c(v, s))) ^ 0.5)))
	}
	dualThr <- function(rho, u, t){
		return(((length(Y) + nrow(D)) ^ 0.5) * absThr + relThr * (as.numeric(crossprod(c(u, t) * rho)) ^ 0.5))
	}
	
	
	temp = crossprod(D)
	diag(temp) = diag(temp) + 1
	
	#prechache factorization
	#R = chol(temp)
	
	result = list("converged" = FALSE)
	objectives = vector()
	
	for(iteration in 1: max_iteration){
		
		#update x
		x = (Y + rho * prevv - rho * prevu)/(rho + 1)

		#update r 
		r = sign(prevs - prevt) * pmax(abs(prevs - prevt) - lambda/rho, 0)
		
		#update v and s

		#the code when using cholesky decomposition
		#temp1 = forwardsolve(t(R), x + prevu + as.vector(crossprod(D, r + prevt)))
		#v = backsolve(R, temp1)
		
		v = as.vector(solve(temp, x + prevu + as.vector(crossprod(D, r + prevt))))
		s = as.vector(D %*% v)
		
		#update u and t
		u = prevu + x - v
		t = prevt + r - s 
		
		objectives = c(objectives, objective(x, Y, lambda, D))
		
		#check convergence
		primalRes = as.numeric(crossprod(c(x - v, r - s)) ^ 0.5)
		dualRes = as.numeric(crossprod(rho * c(prevv - v, prevs - s)) ^ 0.5)
		
		if(primalRes < primalThr(x, v, r ,s) && dualRes < dualThr(rho, u, t)){
			#converged
			result$converged = TRUE
			result$iterations = iteration
			break
		}
		
		
		#rho change, recompute factorization, rescale u
		if(dynamicRho){
			if(primalRes > 5 * dualRes){
				
				#print("increasing rho at iteration")
				#print(iteration)
				#print(" ")
				rho = rho * 2
				u = u/2
				t = t/2

			}
			else if(primalRes * 5 < dualRes){
				#print("decreasing rho at iteration")
				#print(iteration)
				#print(" ")
				rho = rho/2
				u = u * 2
				t = t * 2

			}
		}
		
		prevv = v
		prevs = s
		prevu = u
		prevt = t
		
	}
	result$x = x
	result$objectives = objectives
	return(result)
}