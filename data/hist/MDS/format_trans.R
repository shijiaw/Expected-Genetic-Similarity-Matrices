library(lessR)

for(index in 1:100){
  inputname <- paste('MDS_K', index,'.tsv', sep = '')
  mycor <- read.table(file = inputname, sep = '\t', header = TRUE, row.names = 1)
  newMat <- matrix(0, nr = nrow(mycor), nc = ncol(mycor))
  for(i in 1:nrow(mycor)){
    for(j in 1:ncol(mycor)){
      newMat[i,j] <- as.numeric(mycor[i,j])
    }
  }
  
  colnames(newMat) <- colnames(mycor)
  rownames(newMat) <- colnames(newMat)
  
  #names <- rep(NA, 20)
  #for(i in 1:20){
  #  names[i] <- paste("leaf_", i, sep = '')
  #}
  cor <- corReorder(R = newMat, var=c(leaf_1, leaf_2, leaf_3, leaf_4, leaf_5, leaf_6, leaf_7, leaf_8, leaf_9, leaf_10, leaf_11, leaf_12, leaf_13, leaf_14, leaf_15, leaf_16, leaf_17, leaf_18, leaf_19, leaf_20))
  output_name <- paste("K_", index, ".csv", sep = '')
  write.csv(cor,output_name, row.names = TRUE)
}



