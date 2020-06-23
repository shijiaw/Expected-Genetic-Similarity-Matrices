# create a .sh file
write("sudo python scripts/phylogeny_distance.py --lmm tree1.newick > MDS_K1.tsv", file="k3.sh")
for(i in 2:100){
  command <- paste("sudo python scripts/phylogeny_distance.py --lmm tree",i,".newick > MDS_K",i,".tsv", sep = "")
  write(command, file="k3.sh", append=T)
}
#write("NewFile2", file="k3.sh", append=T)
