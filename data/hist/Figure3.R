library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)

to.upper<-function(X) X[upper.tri(X,diag=TRUE)]

K_seq <- to.upper(read.csv("G/K_1.csv", sep = " ", header = FALSE))
K_tree <- to.upper(read.csv("T/K_1.csv", sep = ",", header = FALSE))
K_S <- to.upper(read.csv("spatialK/K_1.csv", sep = ",", header = FALSE))
K_MDS <- to.upper(read.csv("MDS/K_1.csv", sep = ",", header = TRUE,row.names = 1))

for(i in 2:100){
  name_i_seq <- paste("G/K_",i,".csv", sep = "")
  name_i_tree <- paste("T/K_",i,".csv", sep = "")
  name_i_S <- paste("spatialK/K_",i,".csv", sep = "")
  name_i_MDS <- paste("MDS/K_",i,".csv", sep = "")
  Ki_seq <- to.upper(read.csv(name_i_seq, sep = " ", header = FALSE))
  Ki_tree <- to.upper(read.csv(name_i_tree, sep = ",", header = FALSE))
  Ki_S <- to.upper(read.csv(name_i_S, sep = ",", header = FALSE))
  Ki_MDS <- to.upper(read.csv(name_i_MDS, sep = ",", header = TRUE,row.names = 1))
  K_seq <- c(K_seq, Ki_seq)
  K_tree <- c(K_tree, Ki_tree)
  K_S <- c(K_S, Ki_S)
  K_MDS <- c(K_MDS, Ki_MDS)
}

dif <- c(unlist(K_seq-K_tree), unlist(K_seq-K_S), unlist(K_seq-K_MDS))
Method <- c(rep('Expected Genetic Similarity Matrix', length(unlist(K_seq-K_tree))), rep('Kernel Method', length(unlist(K_seq-K_tree))), rep('Multidimensional Scaling', length(unlist(K_seq-K_tree))))
df = data.frame(dif = dif, Method = Method)
g1 <- ggplot(df, aes(factor(Method), dif))
#g1 <- ggplot(df, aes(x=dif, fill = Difference, color = Difference)) +
#  geom_histogram(aes(y = (..count..)/sum(..count..)), position="identity", binwidth=0.01)+ xlab("") + ylab("proportion")+ theme_bw()+ xlim(-1.5, 0.5) +
#  theme(legend.position="top")


#g1 <- ggplot(data.frame(sequence = unlist(K_seq), tree = unlist(K_tree)), aes(x=sequence, y=tree)) + geom_point(color = 'blue', size = 0.5) + theme_bw()+ xlim(-1, 4.5)+ ylim(-1, 4.5)+ geom_abline(slope = 1, intercept= 0,color="red", linetype = 2)+ ylab(expression(K[ij]^T))+ xlab(expression(K[ij]^G))
gname = c("violin.eps",sep="")  
postscript(gname,width=6,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
g1 + geom_violin(aes(fill = Method))+ ylim(-1.5,0.5)+ theme_minimal()+xlab("")+ylab("pairwise difference (Method - Ground_truth)")+ 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
       strip.background = element_rect(fill="lightblue"),
        legend.position = "top")+coord_flip()
dev.off()

cor(unlist(K_seq), unlist(K_S), method = c("pearson", "kendall", "spearman"))
cor.test(unlist(K_seq), unlist(K_S), method=c("pearson", "kendall", "spearman"))

cor(unlist(K_seq), unlist(K_MDS), method = c("pearson", "kendall", "spearman"))
cor.test(unlist(K_seq), unlist(K_MDS), method=c("pearson", "kendall", "spearman"))

