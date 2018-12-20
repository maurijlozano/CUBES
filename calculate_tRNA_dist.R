#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("A tRNA content matrix must be supplied", call.=FALSE)
} 

sptRNAs <- read.table(args[1], header=TRUE, sep="\t",dec=",")
len <- length (sptRNAs[,1])
ncod <- length(sptRNAs[1,3:length(sptRNAs[1,])])

dist = matrix(nrow=len,ncol=len)

for ( i in 1:len){
    dist[i,i] = 0
}

#calculate distance func.
calculate_tRNA_dist <- function(i,j){
    sim_AB  = 0

    A <- sptRNAs[i, 3:length(sptRNAs[1,])]
    B <- sptRNAs[j, 3:length(sptRNAs[1,])]

    for (cod in 1:ncod){
        if (A[,cod]==0 && B[,cod]==0) {sim_AB = sim_AB +1}
        else {
            sim_AB = sim_AB + ((2*min(A[,cod],B[,cod]))/(A[,cod]+B[,cod]))
        }
    }
    dist_AB <- 1 - (sim_AB/ncod)
    return(dist_AB)
}


for ( i in 1:(len-1)){
    for ( j in (i+1):len){
        dist[i,j] = calculate_tRNA_dist(i,j)
        dist[j,i] =dist[i,j]
    }
}

library('stringr')
rows <- sptRNAs[,2]
rows <- str_pad(rows, 10, side = "right", pad = " ")

write.table(len, file = "len", sep = " ", dec = ".", row.names = F, col.names = F)
write.table(rows, file = "row.mat", sep = " ", dec = ".", row.names = F, col.names = F)
write.table(dist, file = "dist.mat", sep = " ", dec = ".", row.names = F, col.names = F)
write.table(sptRNAs[,c(2,1)], file = "rename.txt", sep = " ", dec = ".", row.names = F, col.names = F)
