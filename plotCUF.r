#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

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

library("tidyr")

fdata1 <- fdata[ grep("C.{1,2}$|^PHE",rownames(fdata)), ]
fdata1 <- data.frame(fdata1)
Dist <- read.table("dist.txt", header=TRUE, sep="\t", dec=",", stringsAsFactors = FALSE)
Dist[,2] <- apply(as.matrix(Dist[,2]),2,function(x) x-x[1])
phedist <- length(Dist[,1])+1
Clast <- as.numeric(as.character(Dist[phedist-1,2]))
Dist <- as.data.frame(rbind(Dist, c("PHE",Dist[phedist-1,2]+Dist[phedist-1,2]/10)))

fdata1 <- data.frame(cbind(SET=rownames(fdata1),X=Dist[,2],fdata1))
aatable <- gather(fdata1, Codon, Freq, -SET, -X)

aminoacids <- c("C","D","E","F","H","K","N","Q","Y","I","A","G","P","T","V","L","R","S")

chead <- read.table("../c_head.txt", header=FALSE, sep=" ", stringsAsFactors = FALSE)
chead <- rbind(chead,c("W","UGG"),c("M","AUG"),c("TER","UAA"),c("TER2","UAG"),c("TER3","UGA"))

mydata <- NULL
for (i in aminoacids){
	codons <- chead[which(chead[,1] == i),2]
	for (j in codons){
		mycodons <- aatable[which(aatable$Codon == codons[which(codons == j)]),]
		mycodons <- cbind(mycodons, AA=rep(i,length(mycodons[,1])))
		mydata <- rbind(mydata,mycodons)	
	}
}

#library(directlabels)
library(ggpubr)
library('ggplot2')
library('ggrepel')
library('ggthemes')
#library(extrafont)


custom.sort <- function(x){
	thdbase <- substring(x[,3], 3)
	fstbase <- substring(x[,3], 1,1)
	x <- cbind(x,fstbase,thdbase)
	ord <- c("U","C","A","G")
	factor(x[,6],levels=ord)
	x[,7] <- factor(x[,7],levels=ord)
	x[order(x[,6],x[,7]),]
}


sopt_DCBS <- read.table("s_opts_DCBS_GNM.txt", header=TRUE, sep="\t", dec=".")
sopt_DCBS <- sopt_DCBS[order(sopt_DCBS[,10], decreasing=T),]
Sopt <- as.numeric(sopt_DCBS[1,1:9])

hdata3 <- read.table('trna.txt', header=T, row.names = NULL, sep=" ")
hdata3 <- hdata3[match(chead$V2,hdata3$codon),]
trnas <- cbind(chead,trna=hdata3[,2])
corder = read.table('../codOrder.txt', header=F, sep=" ")
trnas <- trnas[match(corder$V2,trnas$V2),]
trnas$trna <- as.numeric(as.character(trnas$trna))

rownames(trnas)=NULL

get.ws3 <- function (tRNA, s = NULL, sking) 
{
    if (is.null(s)) 
        s <- c(0, 0, 0, 0, 0.41, 0.28, 0.9999, 0.68, 0.89)
    p = 1 - s
    W = NULL
    for (i in seq(1, 61, by = 4)) W = c(W, p[1] * tRNA[i] + p[5] * 
        tRNA[i + 1], p[2] * tRNA[i + 1] + p[6] * tRNA[i], p[3] * 
        tRNA[i + 2] + p[7] * tRNA[i], p[4] * tRNA[i + 3] + p[8] * 
        tRNA[i + 2])
    W[36] = p[4] * tRNA[36]
    if (sking == 1) 
        W[35] = p[9]
    W = W[-c(11, 12, 15, 36)]
    return(W)
}

custom.sort3 <- function(x){
	thdbase <- substring(x[,2], 3)
	fstbase <- substring(x[,2], 1,1)
	x <- cbind(x,fstbase,thdbase)
	ord <- c("U","C","A","G")
	factor(x[,5],levels=ord)
	x[,6] <- factor(x[,6],levels=ord)
	x[order(x[,5],x[,6]),]
}

if (trnas$trna[35] == 0){sk <- 1 } else {sk <- 0}
ws3 <- get.ws3(tRNA=trnas$trna,s=Sopt, sking=sk)
trnas <- cbind(trnas[-c(11,12,15,36),], ws3)


for (i in aminoacids) {
	myplot <- mydata[which(mydata$AA==i),]
	mytrnas <- trnas[which(trnas$V1==i),]
	myplot <- custom.sort(myplot)
	mytrnas <- custom.sort3(mytrnas)

	mytrnas$trna <- as.numeric(mytrnas$trna)
	mytrnas$ws3 <- as.numeric(mytrnas$ws3)
	if (i == "C" | i == "D" | i == "F" | i == "H" | i == "N" | i == "Y") {mycolors <- c("#0f4c8aff","#8a0000ff") } else if (i == "E" | i == "K" | i == "Q") {mycolors <- c("#808000ff","#cf2e7dff") } else if (i == "I") {mycolors <- c("#0f4c8aff","#8a0000ff","#008000ff") } else if (i == "A" | i == "G" | i == "P" | i == "T" | i == "V") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff") } else if (i == "L") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff","#808000ff", "#cf2e7dff") }  else if (i == "R") {mycolors <- c("#808000ff","#cf2e7dff","#61b8ffff", "#cc0000ff","#008000ff", "#800080ff") } else if (i == "S") {mycolors <- c("#0f4c8aff","#8a0000ff","#61b8ffff", "#cc0000ff","#808000ff", "#800080ff") }

	#mycolors <- c('dodgerblue4','red4','dodgerblue3','red3',  'steelblue1','tomato1')
	myplot[,3] <- factor(myplot[,3],unique(myplot[,3]) , ordered = TRUE)
	mytrnas$V2 <- factor(mytrnas$V2,unique(mytrnas$V2) , ordered = TRUE)
	vlineint <-	(Clast + Clast/30)
	vv <- NULL
	tt <- NULL
	for (j in myplot$Codon) {
		vv <- c(vv, mytrnas[which(mytrnas$V2==j),4])
	}
	for (k in myplot$Codon) {
		tt <- c(tt, mytrnas[which(mytrnas$V2==k),3])
	}
		
#	myplot$Codon <- paste(myplot$Codon,"[",tt,"]","(",round(vv,2),")", sep="")	
	myplot$Codon <- paste(myplot$Codon,"(",round(vv,2),")", sep="")	
	
	myplot$Freq <- as.numeric(myplot$Freq)
	myplot$X <- as.numeric(as.character(myplot$X))
	#order factors
	myplot$Codon = factor(myplot$Codon)
	myplot$Codon = factor(myplot$Codon,levels(myplot$Codon)[unique(myplot$Codon)])

	p1 <- ggplot(myplot, aes(x=X,y=Freq, group=Codon, color=Codon)) +
		#geom_point(aes(color=Codon),show.legend=F)  +
		geom_line(aes(color=Codon), size=1)+
		#geom_text_repel(aes(label=Codon, fontface=2), size = 2, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
		xlab("") +
		ylab("f") + 
		scale_color_manual(values=mycolors) +
		geom_vline(xintercept = vlineint , linetype="dashed", color = "black", size=0.5) +
		#geom_dl(aes(label = Codon), method = list(dl.trans(x = x - .4, y = y - .5), "last.points")) +
		#geom_dl(aes(label = Codon), method = list(dl.trans(x = x - .4, y = y - .5), "first.points")) +
		theme_calc()+
		#scale_color_manual(values=c("dodgerblue3"))+
		guides(colour = guide_legend(override.aes = list(size=3))) +
		theme(axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.position = c(0.6, 0.5),
		#legend.position="bottom",
		legend.title=element_blank(),
		legend.text = element_text(colour="black", size=18, face="bold", family="Arial Narrow"),
	#	legend.box = "horizontal",
		legend.background = element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank())
    nam <- paste("p", i, sep = "")
    assign(nam, p1)
}

for (i in aminoacids) {
	if (i == "C" | i == "D" | i == "F" | i == "H" | i == "N" | i == "Y") {mycolors <- c("#0f4c8aff","#8a0000ff") } else if (i == "E" | i == "K" | i == "Q") {mycolors <- c("#808000ff","#cf2e7dff") } else if (i == "I") {mycolors <- c("#0f4c8aff","#8a0000ff","#008000ff") } else if (i == "A" | i == "G" | i == "P" | i == "T" | i == "V") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff") } else if (i == "L") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff","#808000ff", "#cf2e7dff") }  else if (i == "R") {mycolors <- c("#808000ff","#cf2e7dff","#61b8ffff", "#cc0000ff","#008000ff", "#800080ff") } else if (i == "S") {mycolors <- c("#0f4c8aff","#8a0000ff","#61b8ffff", "#cc0000ff","#808000ff", "#800080ff") }

	#mycolors <- c('dodgerblue4','red4','dodgerblue3','red3',  'steelblue1','tomato1')
	mytrnas <- trnas[which(trnas$V1==i),]
	mytrnas <- custom.sort3(mytrnas)
	#order factors
	mytrnas$V2 = factor(mytrnas$V2)
	mytrnas$V2 = factor(mytrnas$V2,levels(mytrnas$V2)[unique(mytrnas$V2)])


	p2 <- ggplot(mytrnas, aes(x=V2,y=trna, fill=V2)) +
    geom_bar(stat="identity",show.legend=F)+
	geom_text(aes(label=trna,fontface=2), vjust="inward", size=6,show.legend=F) +
#	geom_label(aes(label=trna), color = "white", vjust="inward", size=3,show.legend=F) +
	scale_fill_manual(values=mycolors) +
	theme_calc()+
    #scale_color_manual(values=c("dodgerblue3"))+
	theme(axis.text.x = element_blank(),
	axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
	legend.position = c(0.5, 0.5),
    #legend.position="bottom",
	legend.title=element_blank(),
	legend.background =element_blank(),
	legend.box = "horizontal",
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank())
	nam <- paste("pt", i, sep = "")
    assign(nam, p2)
}

a <- ggarrange(pC,pD,pE,pF,pH,pK,pN,pQ,pY,pI,pA,pG,pP,pT,pV,pL,pR,pS, ncol = 18, align = "hv",labels = aminoacids, hjust = -2.5, vjust = 2, font.label = list(size = 20, color = "black", face = "bold"))
b <- ggarrange(ptC,ptD,ptE,ptF,ptH,ptK,ptN,ptQ,ptY,ptI,ptA,ptG,ptP,ptT,ptV,ptL,ptR,ptS, ncol = 18, align = "hv",hjust = -2.5, vjust = 2)

name2 <- paste(name,"AAfreq_trnas_wi.svg",sep="_")
svg(filename=name2, width=24, height=3.5, pointsize=16)
ggarrange(a,b, nrow = 2, align = "hv",hjust = -2.5, vjust = 2, heights = c(2.5,0.8))
dev.off()

