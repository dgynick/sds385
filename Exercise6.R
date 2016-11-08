
#source("/Users/dgy/Desktop/SDS385/R/Exercise6.r")
#proximalG(X,Y,0.001,0.01,30000)
#proximalGwithIntercept(X,Y,0.05,0.01,30000)
#fastProximalG(X,Y,0.001,0.01,30000)

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

proximalG <- function(X,Y,step,lambda, max_iteration){
	X = scale(X)
	Y = scale(Y)
	
	n = nrow(X)
	p = ncol(X)
	
	
	beta = rep(0, p)
	
	XtX = crossprod(X)
	XtY = crossprod(X,Y)
	
	
	##Proximal gradient
	objectives = vector()
	for(i in 1:max_iteration){
		
		#compute u
		u = as.vector(beta - step * (XtX %*% beta - XtY)/n)

 		#update beta without intercept
 		beta = sign(u) * pmax(abs(u) - lambda * step, 0)	

 		objectives = c(objectives, objective(X,beta,Y,lambda))
	}
	#track convergence
	#plot(1:max_iteration, objectives)
	
	#print(beta)
	
	#compare the predicted Y of the first ten observations
	#print(X[1:10,] %*% beta)
	library(glmnet)
	fit = glmnet(X,Y, lambda = c(lambda), intercept = FALSE, alpha = 1)
	#print(predict(fit, newx = X[1:10,] , s = lambda, type= "response"))
	
	
	#compare beta with ground truth
	ground_truth_beta = fit$beta
	#plot(beta, ground_truth_beta)
	
	Result = list("objectives"= objectives, "beta" = beta)
	return(Result)
}

proximalGwithIntercept <- function(X,Y,step,lambda, max_iteration){
	
	n = nrow(X)
	#include intercept
	X = cbind(X, rep(1, n))

	p = ncol(X)
	
	
	beta = rep(0, p)
	
	XtX = crossprod(X)
	XtY = crossprod(X,Y)


	objectives = vector()
	for(i in 1:max_iteration){

		#compute u
		u = as.vector(beta - step * (XtX %*% beta - XtY)/n)
		
 		#update beta with intercept
		beta[-p] = (sign(u) * pmax(abs(u) - lambda * step, 0))[-p]
		beta[p] = u[p]
		objectives = c(objectives, objective(X,beta,Y,lambda))
	}
	#plot(1:max_iteration, objectives)
	
	print(X[1:10,] %*% beta)

	library(glmnet)
	fit = glmnet(X[,-p], Y, lambda = c(lambda), alpha = 1)
	print(predict(fit, newx = X[1:10,-p] , s = lambda, type= "response"))
	
	
	ground_truth_beta = fit$beta
	plot(beta[1:(p-1)], ground_truth_beta)
}

fastProximalG <- function(X,Y,step,lambda, max_iteration){
	X = scale(X)
	Y = scale(Y)
	
	n = nrow(X)
	p = ncol(X)
	
	
	beta = rep(0, p)
	
	XtX = crossprod(X)
	XtY = crossprod(X,Y)
	
	
	Z = beta
	objectives = vector()
	Sprev = 0
	for(i in 1:max_iteration){
		
		#compute u
		u = as.vector(Z - step * (XtX %*% Z - XtY)/n)
		betaNext = sign(u) * pmax(abs(u) - lambda * step, 0)		
		objectives = c(objectives, objective(X, betaNext,Y,lambda))
		
	    Snext = ((1 + (1 + 4 * Sprev^2)) ^ 0.5) / 2
		Z = betaNext + (Sprev-1)/Snext * (betaNext - beta)
		
 		Sprev = Snext
 		beta = betaNext
	}
	#track convergence
	#plot(1:max_iteration, objectives)
	
	
	#print(X[1:10,] %*% beta)

	library(glmnet)
	fit = glmnet(X,Y, lambda = c(lambda/2/n), intercept = FALSE)
	#print(predict(fit, newx = X[1:10,] , s = lambda, type= "response"))
	
	#compare beta with ground truth
	ground_truth_beta = fit$beta
	#plot(beta, ground_truth_beta)
	
	Result = list("objectives"= objectives)
	return(Result)
}

compareAccelerate <- function(iteration){
	result1 = proximalG(X,Y,0.001,0.01, iteration)
	result2 = fastProximalG(X,Y,0.001,0.01, iteration)
	plot(1: iteration, result1$objectives, col="green", main = "green: proximal gradient red: accelerated proximal gradient")
	lines(1: iteration, result2$objectives, col= "red")
}