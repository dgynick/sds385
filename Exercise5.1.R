


Toy <- function(samplesize, sparsity, lamdalist){
	Theta = rep(3, samplesize)
	mask = rbinom(samplesize, 1, 1 - sparsity)
	
	mean = Theta*mask
	data = rnorm(samplesize, mean = mean)
	for(lamda in lamdalist){
		print(lamda)
		estimate = function(x){
			if(x > lamda){
				return (x-lamda)
			}
			if(x < -lamda){
				return (x+lamda)
			}
			return(0)
		}
		estimation = unlist(lapply(data,estimate))
		#plot(mean, estimation)
	
		MSE = sum((estimation - mean)^2)/samplesize 
		print(MSE)
	}
}

Toy(10000,0.1,c(0,0.05,0.1,0.15,0.2)) #optimal around 0.05 and 0.1

Toy(10000,0.3,c(0,0.1,0.2,0.3,0.4,0.5)) #optimal around 0.3

Toy(10000,0.5,c(0,0.2,0.4,0.6,0.8,1)) #optimal around 0.4

Toy(10000,0.7,c(0,0.3,0.6,0.9,1.2,1.5)) #optimal around 0.6

