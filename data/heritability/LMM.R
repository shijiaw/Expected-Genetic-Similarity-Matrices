set.seed(100)
library(MASS)

Y_0 <- c(175, 160, 166, 131, 172.5, 150, 170, 180)

Y <- (Y_0 - mean(Y_0))/sd(Y_0)

X <- matrix(rep(0, 8), nc = 1)
site <- 1
sigmab2 <- 1

nrow4K <- 8
K_read <- read.csv("K.csv", header = FALSE)
K <- matrix(1, nr = nrow4K, nc = nrow4K)
for(i in 1:nrow4K){
  for(j in 1:nrow4K){
    K[i, j] <- as.numeric(K_read[i,j])
    if(i == j){
      K[i, j] <- K[i, j] + 10^(-7)
    }
  }
}

#eigen <- eigen(K)
#lambda <- eigen$values 
#vec <- eigen$vectors

#K_inv <- eigen$vectors%*%diag(1/lambda)%*%t(eigen$vectors)
#K_inv%*%K

## simulate random effects ###
b <- rep(0, 8)

b <- mvrnorm(n = 1, mu = rep(0, 8), Sigma = 0.25*K)

Gm <- (X[,site]-mean(X[,site]))
#1/sum(Gm*Gm)*sum(Gm*Y)

#1/sum(Gm*Gm)*sum(Gm*Y[,1])

Niter <- 5000
betastore <- rep(NA, Niter)
bstore <- rnorm(length(b),0,0.4)
sigmae2store <- rep(NA, Niter)
sigmab2store <- rep(NA, Niter)
betastore[1] <- 0.3
sigmae2store[1] <- 0.5
sigmab2store[1] <- 1
mub <- matrix(NA, nr = 8, nc = 1)
#bKb <- rep(NA, nrow(Y)/nrow4K)
#yGb <- rep(NA, nrow(Y)/nrow4K)
bKb <- 1
yGb <- 1
aostar <- 2+length(Y)/2
gostar <- 2+length(Y)/2
##Conduct MCMC sampling ###
for(iter in 2:Niter){
  #sigmabeta2 <- 1/(sum(Gm*Gm)/sigmae2store[iter-1]+1)
  #mubeta <- sigmabeta2/sigmae2store[iter-1]*sum(Gm*(Y-b))
  #betastore[iter] <- rnorm(1, mean = mubeta, sd = sqrt(sigmabeta2))
  betastore[iter] <- 0
    
  Sigmabm <- solve(diag(1/sigmae2store[iter-1], nrow4K) + solve(K)/sigmab2store[iter-1])
  for(index in 1:1){
    mub[,index] <- Sigmabm%*%((Y-Gm*betastore[iter])/sigmae2store[iter-1])
    bstore <- mvrnorm(1, mu = mub[,index], Sigma = Sigmabm)
    bKb[index] <- bstore%*%solve(K)%*%bstore/2
    yGbtemp <- Y - bstore-Gm*betastore[iter]
    yGb[index] <- sum(yGbtemp*yGbtemp)/2
  }
  
  bostar <- 1+sum(bKb)
  
  sigmab2store[iter] <- 1/rgamma(1, shape = aostar, scale = 1/bostar)
  
  hostar <- 1+sum(yGb)
  
  sigmae2store[iter] <- 1/rgamma(1, shape = gostar, scale = 1/hostar)
  
}

gname = c("traceplots.eps",sep="")  
postscript(gname,width=6,height=4.8,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
plot(sigmae2store, type = 'l')
plot(sigmab2store, type = 'l')
dev.off()

h2 <- mean(sigmab2store[-(1:2000)])/(mean(sigmab2store[-(1:2000)])+mean(sigmae2store[-(1:2000)]))



