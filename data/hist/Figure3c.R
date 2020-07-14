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
Difference <- c(rep('G-T', length(unlist(K_seq-K_tree))), rep('G-S', length(unlist(K_seq-K_tree))), rep('G-MDS', length(unlist(K_seq-K_tree))))

df = data.frame(dif = dif[Difference == 'G-MDS'], Difference = Difference[Difference == 'G-MDS'])

g1 <- ggplot(df, aes(x=dif, fill = Difference, color = Difference)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), position="identity", binwidth=0.01, colour="#5e9bfc", fill="#5e9bfc")+ xlab("Pairwise Difference") + ylab("Proportion")+ theme_minimal()+ xlim(-1.5, 0.5) + ylim(0,0.25) +
  theme(legend.position="none") + ggtitle('Multidimensional Scaling -vs- Ground Truth')

ggsave(file = 'Figure3c.pdf', plot = ggplot2::last_plot(), height = 2,
  width = 7, units = "in", dpi = 900)

