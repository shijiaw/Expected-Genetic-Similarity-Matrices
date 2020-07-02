library(xtable)
library(Rmisc) 
library(stringr)
library(reshape2)
#library(hrbrthemes)
library(ggplot2)
setwd("/Users/oudomame/Desktop/sfu_vault/msdir/sim/sim3")

#to.upper<-function(X) X[upper.tri(X,diag=TRUE)]

K_tree <- read.csv("K.csv", sep = ",", header = FALSE)
var1 <- rep(c('H. habilis', 'H. rudolfensis', 'Georgian H. erectus', 	
                      'African H. erectus s.s.', 'Asian H. erectus s.s.', 'H. antecessor', 'H. neanderthalensis', 	
                      'H. sapiens'), 8)
var2 <- rep(c('H. habilis', 'H. rudolfensis', 'Georgian H. erectus', 	
                      'African H. erectus s.s.', 'Asian H. erectus s.s.', 'H. antecessor', 'H. neanderthalensis', 	
                      'H. sapiens'), each = 8)
xtable(K_tree,digits = 3)


data <- data.frame(var1 = var1, var2 = var2, Similarity = melt(K_tree)$value)
g1 <- ggplot(data = data, aes(x=var1, y=var2, fill=Similarity)) + 
   geom_tile()+ scale_fill_gradient(low="grey", high="black")+ xlab("")+ ylab("")+
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.line.x = element_line(size = 3, colour = "white"), panel.background = element_rect(fill = "white"))+
  labs(title = "Genetic similarity for ancient hominins")


gname = c("gsm4real.eps",sep="")  
postscript(gname,width=6,height=4.8,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
g1
dev.off()
