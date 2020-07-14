set.seed(100)
library(MASS)
library(scales)

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

Niter <- 10000
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

library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)

data_sigmag2 <- data.frame(sigmag2 = sigmab2store, Iteration = as.vector(1:Niter) )
data_sigmae2 <- data.frame(sigmae2 = sigmae2store, Iteration = as.vector(1:Niter) )
h1 <- mean(sigmab2store[-(1:5000)])
h2 <- mean(sigmae2store[-(1:5000)])

g1 <- ggplot(data_sigmag2, aes(x=Iteration, y=sigmag2)) + geom_line() + theme_minimal() + geom_hline(yintercept=h1, linetype="dotted", size=1, color = "red")+ xlab(expression(sigma[g]^2)) + ylab("") + scale_color_manual(values = c('#595959', 'blue'))
ggsave(file = 'Figure6topa.pdf', plot = ggplot2::last_plot(), height = 3, width = 3, units = "in", dpi = 900)


g2 <- ggplot(data_sigmae2, aes(x=Iteration, y=sigmae2)) + geom_line() + theme_minimal() + geom_hline(yintercept=h2, linetype="dotted", size=1, color = "red")+ xlab(expression(sigma[e]^2)) + ylab("") + scale_color_manual(values = c('#595959', 'blue'))
ggsave(file = 'Figure6topb.pdf', plot = ggplot2::last_plot(), height = 3, width = 3, units = "in", dpi = 900)


gname = c("MCMC.eps",sep="")  
postscript(gname,width=8,height=3.2,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(g1,
          g2,
          ncol = 2, nrow = 1, common.legend = TRUE)
dev.off()

df_p <- data.frame(
  p11 = sigmab2store[-(1:5000)], p21 = sigmae2store[-(1:5000)], p31 = sigmab2store[-(1:5000)]/(sigmab2store[-(1:5000)]+sigmae2store[-(1:5000)])
)


mean_p11 <- mean(sigmab2store[-(1:5000)])
mean_p21 <- mean(sigmae2store[-(1:5000)])
mean_p31 <- mean(sigmab2store[-(1:5000)]/(sigmab2store[-(1:5000)]+sigmae2store[-(1:5000)]))
p11_25 <- quantile(sigmab2store[-(1:5000)], 0.025)
p11_975 <- quantile(sigmab2store[-(1:5000)], 0.975)
p21_25 <- quantile(sigmae2store[-(1:5000)], 0.025)
p21_975 <- quantile(sigmae2store[-(1:5000)], 0.975)
p31_25 <- quantile(sigmab2store[-(1:5000)]/(sigmab2store[-(1:5000)]+sigmae2store[-(1:5000)]), 0.025)
p31_975 <- quantile(sigmab2store[-(1:5000)]/(sigmab2store[-(1:5000)]+sigmae2store[-(1:5000)]), 0.975)



hist_p11 <- ggplot(data=df_p, aes(p11)) + 
  geom_histogram(aes(y=..count../sum(..count..)), position="identity", color="black", fill="black", alpha=0.4, binwidth = 0.01)+ xlab(expression(sigma[g]^2)) + 
  theme_minimal()+ ylab("Density") + geom_vline(xintercept=mean_p11, color = "red", size=1)+ xlim(0, 2) + 
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02), limits = c(0, 0.025)) +
  geom_vline(xintercept=p11_25, linetype="dotted", color = "blue", size=1)+ geom_vline(xintercept=p11_975, linetype="dotted", color = "blue", size=1)

ggsave(file = 'Figure6a.pdf', plot = ggplot2::last_plot(), height = 2, width = 2, units = "in", dpi = 900)

hist_p21 <- ggplot(data=df_p, aes(p21)) + 
  geom_histogram(aes(y=..count../sum(..count..)), position="identity", color="black", fill="black", alpha=0.4, binwidth = 0.01)+ xlab(expression(sigma[e]^2)) + 
  theme_minimal()+ ylab("") + geom_vline(xintercept=mean_p21, color = "red", size=1)+xlim(0, 2) +   
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02) , limits = c(0, 0.025)) +
  geom_vline(xintercept=p21_25, linetype="dotted", color = "blue", size=1)+ geom_vline(xintercept=p21_975, linetype="dotted", color = "blue", size=1)

ggsave(file = 'Figure6b.pdf', plot = ggplot2::last_plot(), height = 2, width = 2, units = "in", dpi = 900)

hist_p31 <- ggplot(data=df_p, aes(p31)) + 
  geom_histogram(aes(y=..count../sum(..count..)), position="identity", color="black", fill="white", alpha=0.4, binwidth = 0.01)+ xlab(expression(h^2)) + 
  theme_minimal()+ ylab("") + geom_vline(xintercept=mean_p31, color = "red", size=1)+  
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02), limits = c(0,0.025)) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), limits = c(0,1)) +
  geom_vline(xintercept=p31_25, linetype="dotted", color = "blue", size=1)+ geom_vline(xintercept=p31_975, linetype="dotted", color = "blue", size=1)

ggsave(file = 'Figure6c.pdf', plot = ggplot2::last_plot(), height = 2, width = 2, units = "in", dpi = 900)

gname = c("sigma_hist.eps",sep="")  
postscript(gname,width=10,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(hist_p11,
          hist_p21,
          hist_p31,
          ncol = 3, nrow = 1, common.legend = TRUE)

dev.off()


h2 <- mean(sigmab2store[-(1:5000)])/(mean(sigmab2store[-(1:5000)])+mean(sigmae2store[-(1:5000)]))



