#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R
library("tidyr")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

setwd(args[1])
name <- args[2]

#load data files

hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")

data = read.table("fmdata.txt", header=F, sep=" ") 
data <- cbind(data, rep(1, length(data[,1])), rep(1,length(data[,1])))
data <- data[,-2]
colnames(data) <- c(as.character(unlist(hdata[,1:60])), "AUG", "UGG")

corder2 <- c("GCU","GCC","GCA","GCG","CGU","CGC","CGA","CGG","AGA","AGG","AAU","AAC","GAU","GAC","UGU","UGC","CAA","CAG","GAA","GAG","GGU","GGC","GGA","GGG","CAU","CAC","AUU","AUC","AUA","UUA","UUG","CUU","CUC","CUA","CUG","AAA","AAG","AUG","UUU","UUC","CCU","CCC","CCA","CCG","UCU","UCC","UCA","UCG","AGU","AGC","ACU","ACC","ACA","ACG","UGG","UAU","UAC","GUU","GUC","GUA","GUG")
t<-t(data)
t <- t[match(corder2,row.names(t)),]
fdata <- t(t)
fdata <- cbind(fdata[,1:37],fdata[,39:54],fdata[,56:61])
fdata <- data.frame(fdata)
rownames(fdata) <- data[,1]


fdata1 <- fdata[ grep("C.{1,2}$",rownames(fdata)), ]
fdata1 <- data.frame(fdata1)
phe <- fdata[ grep("^PHE",rownames(fdata)), ]

for (i in 1:ncol(fdata1)){fdata1[,i]<-as.numeric(as.character(fdata1[,i]))}
for (i in 1:ncol(phe)){phe[,i]<-as.numeric(as.character(phe[,i]))}

delta <- fdata1[length(fdata1$GUG),]-fdata1[1,]
deltaphe <- phe[1,]-fdata1[length(fdata1$GUG),] 

sopt_DCBS <- read.table("s_opts_DCBS_GNM2.txt", header=TRUE, sep="\t", dec=".")
sopt_DCBS <- sopt_DCBS[order(sopt_DCBS[,11], decreasing=T),]
Sopt <- as.numeric(sopt_DCBS[1,1:10])

hdata3 <- read.table('trna.txt', header=T, row.names = NULL, sep=" ")
chead <- read.table("../c_head.txt", header=FALSE, sep=" ", stringsAsFactors = FALSE)
chead <- rbind(chead,c("W","UGG"),c("M","AUG"),c("TER","UAA"),c("TER2","UAG"),c("TER3","UGA"))
hdata3 <- hdata3[match(chead$V2,hdata3$codon),]

trnas <- cbind(chead,trna=hdata3[,2])
corder = read.table('../codOrder.txt', header=F, sep=" ")
trnas <- trnas[match(corder$V2,trnas$V2),]
trnas$trna <- as.numeric(as.character(trnas$trna))

rownames(trnas)=NULL
rownames(delta)="DeltaC"

aminoacids <- c("C","D","E","F","H","K","N","Q","Y","I","A","G","P","T","V","L","R","S")
aminoacids3letter <- c("Cys","Asp","Glu","Phe","His","Lys","Asn","Gln","Tyr","Ile","Ala","Gly","Pro","Trp","Val","Leu","Arg","Ser")
aancod <- c(2,2,2,2,2,2,2,2,2,3,4,4,4,4,4,6,6,6)
aa.table <- as.data.frame(cbind(aminoacids,aminoacids3letter,aancod))
aa.table$aancod <- as.numeric(as.character(aa.table$aancod))

aatable <- gather(delta, Codon, Freq)
aatable2 <- gather(deltaphe, Codon, Freq)
aatable <- cbind(aatable, Dphe=aatable2[,2])

custom.sort3 <- function(x){
	data.frame(x)
	thdbase <- substring(x, 3)
	scndbase <- substring(x, 2,2)
	fstbase <- substring(x, 1,1)
	x <- cbind(x,fstbase,scndbase,thdbase)
	ord <- c("U","C","A","G")
	factor(x[,2],levels=ord)
	factor(x[,3],levels=ord)
	x[,4] <- factor(x[,4],levels=ord)
	x[order(x[,2],x[,4]),]
}

mydata <- NULL
for (i in aminoacids){
	codons <- chead[which(chead[,1] == i),2]
	codons <- custom.sort3(codons)[,1]
	for (j in codons){
		mycodons <- aatable[which(aatable$Codon == codons[which(codons == j)]),]
		mycodons <- cbind(mycodons, AA = i, AA3letter=aa.table[which(aa.table$aminoacids == i),2])
		mydata <- rbind(mydata,mycodons)	
	}
}
rownames(mydata)=NULL



get.ws3 <- function (tRNA, s = NULL, sking) 
{
    if (is.null(s)) 
        s <- c(0, 0, 0, 0, 0.41, 0.28, 0.9999, 0.68, 0.89, 0.5)
    p = 1 - s
    W = NULL
	box4 <- c(5,17,21,29,33,37,49,53,61)
	box2 <- c(1,9,13,25,41,45,57)
	codonWorder <- c(1,2,3,4,9,10,11,12,13,14,15,16,25,26,27,28,41,42,43,44,45,46,47,48,57,58,59,60,5,6,7,8,17,18,19,20,21,22,23,24,29,30,31,32,33,34,35,36,37,38,39,40,49,50,51,52,53,54,55,56,61,62,63,64)
   
	for (i in box2) W = c(W, p[1] * tRNA[i] + p[5] * tRNA[i + 1], p[2] * tRNA[i + 1] + p[6] * tRNA[i], p[3] * tRNA[i + 2] + p[7] * tRNA[i], p[4] * tRNA[i + 3] + p[8] * tRNA[i + 2])
	for (i in box4) W = c(W, p[1] * tRNA[i] + p[5] * tRNA[i + 1] + p[10] * tRNA[i + 2], p[2] * tRNA[i + 1] + p[6] * tRNA[i], p[3] * tRNA[i + 2] + p[7] * tRNA[i], p[4] * tRNA[i + 3] + p[8] * tRNA[i + 2])
	
	W <- W[match(1:64,codonWorder)]
	W[36] = p[4] * tRNA[36]

    if (sking == 1) 
        W[35] = p[9]
    W = W[-c(11, 12, 15, 36)]
    return(W)
}

if (trnas$trna[35] == 0){sk <- 1 } else {sk <- 0}
ws3 <- get.ws3(tRNA=trnas$trna,s=Sopt, sking=sk)
trnas <- cbind(trnas[-c(11,12,15,36),], ws3)

trnas.norm <- NULL
for (i in aminoacids){
	wsnorm <- trnas[which(trnas[,1] == i),]
	wsmax <- max(wsnorm[,4])
	wsnorm[,4] <- wsnorm[,4]/wsmax
	trnas.norm <- rbind(trnas.norm,wsnorm)	
}

trnas.norm <- as.data.frame(trnas.norm)


heatmap.table <- cbind(trnas.norm[match(mydata[,1],trnas.norm[,2]),],mydata[,c(2,3,5)])
rownames(heatmap.table)=paste(heatmap.table[,7],"|",heatmap.table[,2], sep="")
colnames(heatmap.table)[1]<-"AA"
colnames(heatmap.table)[2]<-"Codon"
colnames(heatmap.table)[3]<-"tRNAs"
colnames(heatmap.table)[4]<-"w"
colnames(heatmap.table)[5]<-"\u0394Cn-C1"
colnames(heatmap.table)[6]<-"\u0394PHE-Cn"
heatmap.tablebk <- heatmap.table

library('pheatmap')

heatmap.table <- heatmap.table[,4:(ncol(heatmap.table)-1)]
heatmap.table[,1] <- heatmap.table[,1]+1


paletteLength <- 50
myColor <- c(colorRampPalette(c("firebrick3", "#D13838"))(2),             colorRampPalette(c("#D54A4A", "#DD6E6E"))(2), colorRampPalette(c("#E18080", "#F2C8C8"))(8), colorRampPalette(c("white"))(1), colorRampPalette(c("#D0D0E7", "#7373B9"))(8), colorRampPalette(c("#5C5CAE", "#2E2E97"))(2), colorRampPalette(c("#17178B", "#000080"))(2), colorRampPalette(c("white", "navy"))(25))

myBreaks <- c(seq(-1, -0.7, length.out=2), seq(-0.65, -0.4, length.out=2), seq(-0.35, -0.05, length.out=8), seq(-0.025, 0.025, length.out=2), seq(0.05, 0.35, length.out=8), seq(0.4, 0.65, length.out=2), seq(0.7, 1, length.out=2), seq(0.9999,2, length.out=ceiling(paletteLength/2)))



svg(filename="delta-heatmap_U.svg", width=1, height=10, pointsize=10)
p<-pheatmap(heatmap.table, cellwidth = 5, cellheight = 10, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = F,color=myColor, breaks=myBreaks,fontsize = 6, gaps_row =c(2,4,6,8,10,12,14,16,18,21,25,29,33,37,41,47,53))
print(p)
dev.off()
