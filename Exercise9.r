#source("/Users/dgy/Desktop/SDS385/R/Exercise9.r")

main <- function(n,p){
	X = rnorm(n*p)
	X = matrix(X, nrow = n)
	result = singleFactorPMDlasso(X, n^0.5/2, p^0.5/2)
	difference = X - result$Xhat
	print(X)
	print(difference)
}


softThresh <- function(a, lambda, lower, upper, length, criteria = 0.0001){
	relativeCriteria = criteria * length
	
	delta = (lower + upper)/2
	result = sign(a) * pmax(abs(a) - delta, 0)

	result = result/(sum(result ^ 2) ^ 0.5)
	norm1 = sum(abs(result))
	
	if(norm1 > lambda + relativeCriteria){
		return(softThresh(a, lambda, delta, upper, length))
	}
	else if(norm1 < lambda - relativeCriteria){
		return(softThresh(a, lambda, lower, delta, length))
	}
	return(result)
}

#c1 c2 should be between 0 and 1
singleFactorPMDlasso <- function(X, c1, c2, max_iteration = 10000, convergenceCriteria = 0.0001){
	n = nrow(X)
	p = ncol(X)
	
	lambda1 = max(1, c1 * n)
	lambda2 = max(1, c2 * p)
	
	stop1 = convergenceCriteria * n
	stop2 = convergenceCriteria * p
	
	result = list()
	result$converged = FALSE
	prevv = as.vector((eigen(crossprod(X), symmetric = TRUE))$vectors[,1])
	prevu = rep(0, n)
	
	for(iteration in 1: max_iteration){
		#update u
		xv = as.vector(X %*% prevv)
		u = xv/(sum(xv ^ 2))^0.5
		if(sum(abs(u)) > lambda1){
			u = softThresh(xv, lambda1, 0, max(abs(xv)), n)
		}
		
		#update v
		xtu = as.vector(crossprod(X, u))
		v = xtu/(sum(xtu ^2))^0.5
		if(sum(abs(v)) > lambda2){
			v = softThresh(xtu, lambda2, 0, max(abs(xtu)), p)
		}
		
		if(sum((u - prevu)^2) < stop1^2 && sum((v - prevv)^2) < stop2^2){
			result$converged = TRUE
			result$iterations = iteration
			prevv = v
			prevu = u
			break
		}
		prevv = v
		prevu = u
	}
	d = as.numeric(crossprod(prevu, X) %*% v)
	
	result$u = prevu
	result$v = prevv
	result$d = d
	result$Xhat = crossprod(t(u),v)*d
	result$Fdistance = sum((result$Xhat - X) ^ 2)
	return(result)
}

MultiFactorPMDlasso <- function(X, k, c1, c2, max_iteration = 10000, convergenceCriteria = 0.0001){
	if(k == 1){
		return(singleFactorPMDlasso(X, c1, c2))
	}
	result = list()
	d = vector()
	Xleft = X
	for(i in 1:k){
		temp = singleFactorPMDlasso(Xleft, c1, c2)
		d = c(d, temp$d)
		Xleft = Xleft - temp$Xhat
		result$u = cbind(result$u, temp$u)
		result$v = cbind(result$v, temp$v)
	}
	result$d = diag(d)
	result$Xhat  = result$u %*% result$d %*% t(result$v)
	result$Fdistance = sum((result$Xhat - X) ^ 2)
	return(result)
}

TFIDF <- function(X){
	docfreq = colSums(ifelse(X>0, 1, 0))
	docfreq = docfreq/nrow(X)
	processrow <- function(row){
		return((log(pmax(1,row)) + ifelse(row > 0, 1, 0))/docfreq)
	}
	return(t(apply(X, 1, processrow)))
}


socialMarketing <- function(k){
	
	X = read.csv("/Users/dgy/Desktop/SDS385/social_marketing.csv")
	X = X[,-1]
	#scale X: 1.scale each row to sum =1 2. scale each row with max = 1  3. scale each row to mean = 0 and same standard deviation 4. tfidf 
	X = TFIDF(X)
	
	result = MultiFactorPMDlasso(X, k, 0.7, 0.7)
	
	return(result)	
}

#result = socialMarketing(k)
#u = result$u[, ?]
#X[u>0.01,]
# 1st u: people with spam
# 4th u: people with cooking, personal fittness, outdoors, health nutrition
# 5th u: people with school, family, food, sports_fandom, religion, parenting

