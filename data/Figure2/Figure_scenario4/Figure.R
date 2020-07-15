library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)
#setwd("/Users/oudomame/Desktop/sfu_vault/msdir")

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



g1 <- ggplot(data.frame(sequence = unlist(K_seq), tree = unlist(K_tree)), aes(x=sequence, y=tree)) + geom_point(color = 'blue', size = 0.5) + theme_minimal()+ xlim(-1, 4.5)+ ylim(-1, 4.5)+ geom_abline(slope = 1, intercept= 0,color="red", linetype = 2)+ ylab("Expected Genetic Similarity Matrix")+ xlab("Ground_truth")+ ggtitle("N = 20, M = 1000, theta = 0.3")
gname = c("s1000m20mu03.eps",sep="")  
ggsave(file = 'Figure2d.png', plot = ggplot2::last_plot(), height = 3, width = 3, units = "in", dpi = 900)




