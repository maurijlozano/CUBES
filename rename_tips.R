#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R
library(ape)
tree <- read.tree(file = "outtree",)
names <- read.table("rename.txt", sep=" ")
tree$tip.label<-names[[2]][match(tree$tip.label, names[[1]])]
write.tree(tree,file = "tRNA_NJ.tree")
