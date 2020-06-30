library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)
setwd("/Users/oudomame/Desktop/sfu_vault/msdir")

to.upper<-function(X) X[upper.tri(X,diag=TRUE)]

K_seq <- to.upper(read.csv("output_seq/scenario4/K_1.csv", sep = " ", header = FALSE))
K_tree <- to.upper(read.csv("output_tree/scenario4/K_1.csv", sep = ",", header = FALSE))


for(i in 2:100){
  name_i_seq <- paste("output_seq/scenario4/K_",i,".csv", sep = "")
  name_i_tree <- paste("output_tree/scenario4/K_",i,".csv", sep = "")
  Ki_seq <- to.upper(read.csv(name_i_seq, sep = " ", header = FALSE))
  Ki_tree <- to.upper(read.csv(name_i_tree, sep = ",", header = FALSE))
  K_seq <- c(K_seq, Ki_seq)
  K_tree <- c(K_tree, Ki_tree)
}



g1 <- ggplot(data.frame(sequence = unlist(K_seq), tree = unlist(K_tree)), aes(x=sequence, y=tree)) + geom_point(color = 'blue', size = 0.5) + theme_bw()+ xlim(-1, 4.5)+ ylim(-1, 4.5)+ geom_abline(slope = 1, intercept= 0,color="red", linetype = 2)+ ylab("") + xlab("N = 20, M = 1000, theta = 0.3")
gname = c("sim/scenario4/s1000m20mu03.eps",sep="")  
postscript(gname,width=3,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
g1
dev.off()



