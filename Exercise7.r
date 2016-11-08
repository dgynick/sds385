#source("/Users/dgy/Desktop/SDS385/R/Exercise7.r")
library("Matrix");

X = read.csv("~/Desktop/SDS385/diabetes/diabetesX.csv", header= FALSE)
X = as.matrix(X)
X = X[-1,]
X = data.matrix(apply(X, 2, as.numeric))


Y = as.matrix(read.csv("~/Desktop/SDS385/diabetes/diabetesY.csv", header= FALSE))
Y = data.matrix(apply(Y, 2, as.numeric))


#objective function without intercept
objective <- function(X, beta, Y, lambda){
	return(crossprod(Y - X %*% beta)/2/nrow(X) + lambda * sum(abs(beta)))
}

ADMM <- function(X,Y,lambda, max_iteration, absThr = 0.0001, relThr = 0.0001){
	X = scale(X)
	Y = scale(Y)
	
	n = nrow(X)
	p = ncol(X)
	
	XtX = crossprod(X)/n
	XtY = crossprod(X,Y)/n
	

	rho = 0.1
	
	
	prevu = rep(0, p)
	prevz = rep(0, p)
	beta = rep(0, p)

	#stopping criterion
	primalThr <- function(beta, z){
		(p^0.5) * absThr + relThr * (max(crossprod(beta), crossprod(z)) ^ 0.5)
	}
	dualThr <- function(rho,u){
		(p^0.5) * absThr + relThr * (crossprod(u * rho) ^ 0.5)
	}
	
	#prechache factorization
	R = chol(XtX + rho * diag(p))
	
	result = list("converged" = FALSE)
	objectives = vector()
	
	for(iteration in 1: max_iteration){
		#solve beta
		b = XtY + rho * prevz - rho * prevu
		temp1 = forwardsolve(t(R), b)
		beta = backsolve(R, temp1)
		objectives = c(objectives, objective(X, beta, Y, lambda))
		
		#solve z
		z = sign(beta + prevu) * pmax(abs(beta + prevu) - lambda/rho, 0)
		
		#update dual var u
		u = prevu + (beta - z)

		
		#check convergence
		primalRes = crossprod(beta - z) ^ 0.5
		dualRes = crossprod(rho * (prevz - z)) ^ 0.5
		
		if(primalRes < primalThr(beta, z) && dualRes < dualThr(rho,u)){
			#converged
			result$converged = TRUE
			result$iterations = iteration
			break
		}
		
		
		#rho change, recompute factorization, rescale u
		if(TRUE){
			if(primalRes > 10 * dualRes){
				
				#print("increasing rho at iteration")
				#print(iteration)
				#print(" ")
				rho = rho * 2
				u = u/2
				R = chol(XtX + rho*diag(p))
			}
			else if(primalRes * 10 < dualRes){
				#print("decreasing rho at iteration")
				#print(iteration)
				#print(" ")
				rho = rho/2
				u = u * 2
				R = chol(XtX + rho*diag(p))
			}
		}
		
		prevu = u
		prevz = z

	}

	result$beta = beta
	result$z = z
	result$objectives = objectives
	result$dualRes = dualRes
	result$primalRes = primalRes
	return(result)
}

compareConvergence<- function(lambda){
	source("/Users/dgy/Desktop/SDS385/R/Exercise6.r")
	result = ADMM(X,Y,lambda,30000)
	result1 = proximalG(X,Y,0.001, lambda, result$iterations)
	result2 = fastProximalG(X,Y,0.001, lambda, result$iterations)
	plot(1:result$iterations, result$objectives, col="blue", main = "green: proximal gradient red: accelerated proximal gradient blue: ADMM")
	lines(1:result$iterations, result1$objectives[1:result$iterations], col= "green")
	lines(1:result$iterations, result2$objectives[1:result$iterations], col= "red")
}

compareBeta<-  function(lambda){
	result = ADMM(X,Y,lambda,100000)
	library(glmnet)
	X = scale(X)
	Y = scale(Y)
	fit = glmnet(X,Y, lambda = c(lambda), intercept = FALSE, alpha = 1)
	plot(fit$beta, result$beta)
}
