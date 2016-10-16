X = read.csv("~/Desktop/SDS385/diabetes/diabetesX.csv", header= FALSE)
X = as.matrix(X)
X = X[-1,]
X = data.matrix(apply(X, 2, as.numeric))

Y = as.matrix(read.csv("~/Desktop/SDS385/diabetes/diabetesY.csv", header= FALSE))
Y = data.matrix(apply(Y, 2, as.numeric))

#Q(A)
library(glmnet)
fit = glmnet(X,Y)
plot(fit)

df = fit$df
dev = fit$dev.ratio
lamda = fit$lambda

SSE = var(Y)*(nrow(Y)-1)/nrow(Y)
insampleErr = SSE * (1-dev)
print(insampleErr)


#Q(B)
crossvalidation <- function(lamda, fold){
	order = sample.int(nrow(X))
	Moose = 0
	
	testsize = floor(nrow(X)/fold)
	
	for(i in 1:fold){
		start = (i-1)* testsize + 1
		end = i*testsize
		if(i==fold){
			end = nrow(X)
		}
		
		trainX = X[-order[start:end],]
		trainY = Y[-order[start:end],]
		
		testX = X[order[start:end],]
		testY = Y[order[start:end],]
		
		fit = glmnet(trainX, trainY, lambda = c(lamda))
		Yhat = predict(fit,newx=testX,s=c(lamda))[,1]
		
		Moose = Moose + sum((testY-Yhat)^2)/(end-start+1)/fold
	}
	return(Moose)
}

range = -4:4
myCVMoose =  vector()
for(loglamda in range){
	myCVMoose = c(myCVMoose, crossvalidation(exp(loglamda),5))
}

#sanity check
cvfit = cv.glmnet(X, Y)

plot(cvfit, main="red represents the moose estimate by cv.glmnet, green represents my estimate of moose, blue represents the in sample error")
lines(range, myCVMoose, col="green")
lines(log(lamda), insampleErr, col = "blue")

#Q(C)
Yhat = predict(fit, newx=X,s=c(0))
sigmasquare = sum((Y-Yhat)^2)/(nrow(X)-ncol(X))

penalty = 2*df/nrow(X) *sigmasquare
Cp = insampleErr + penalty

plot(cvfit, main="brown represents Cp, red represents the moose estimate by cv.glmnet, green represents my estimate of moose, blue represents the in sample error")
lines(range, myCVMoose, col="green")
lines(log(lamda), insampleErr, col = "blue")
lines(log(lamda), Cp, col = "brown")
