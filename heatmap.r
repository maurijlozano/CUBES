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
heatmapdata <- (as.numeric(as.character(t(fdata1[length(fdata1[,1]),]))) - as.numeric(as.character(t(fdata1[1,]))))
heatmapdata <- data.frame(heatmapdata)


hdata.tb <- rbind(c(heatmapdata[1:4,1],0,0),c(heatmapdata[5:10,1]),c(heatmapdata[11:12,1],0,0,0,0),c(heatmapdata[13:14,1],0,0,0,0),c(heatmapdata[15:16,1],0,0,0,0),c(heatmapdata[17:18,1],0,0,0,0),c(heatmapdata[19:20,1],0,0,0,0),c(heatmapdata[21:24,1],0,0),c(heatmapdata[25:26,1],0,0,0,0),c(heatmapdata[27:29,1],0,0,0),c(heatmapdata[30:35,1]),c(heatmapdata[36:37,1],0,0,0,0),c(heatmapdata[38:39,1],0,0,0,0),c(heatmapdata[40:43,1],0,0),c(heatmapdata[44:49,1]),c(heatmapdata[50:53,1],0,0),c(heatmapdata[54:55,1],0,0,0,0),c(heatmapdata[56:59,1],0,0))


#esto es rscu phe -> voy por acá!!!!
fdata2 <- data.frame(t(fdata[ grep("^PHE$",rownames(fdata)), ]))
fdata2$PHE <- as.numeric(as.character(fdata2$PHE))
hdata.tb2 <- rbind(c(fdata2[1:4,1],0,0),c(fdata2[5:10,1]),c(fdata2[11:12,1],0,0,0,0),c(fdata2[13:14,1],0,0,0,0),c(fdata2[15:16,1],0,0,0,0),c(fdata2[17:18,1],0,0,0,0),c(fdata2[19:20,1],0,0,0,0),c(fdata2[21:24,1],0,0),c(fdata2[25:26,1],0,0,0,0),c(fdata2[27:29,1],0,0,0),c(fdata2[30:35,1]),c(fdata2[36:37,1],0,0,0,0),c(fdata2[38:39,1],0,0,0,0),c(fdata2[40:43,1],0,0),c(fdata2[44:49,1]),c(fdata2[50:53,1],0,0),c(fdata2[54:55,1],0,0,0,0),c(fdata2[56:59,1],0,0))

#ntRNA
hdata3 <- read.table('trna.txt', header=T, row.names = NULL, sep=" ")
hdata3 <- hdata3[match(corder2,hdata3$codon),]

hdata.tb3 <- rbind(c(hdata3[1:4,2],0,0),c(hdata3[5:10,2]),c(hdata3[11:12,2],0,0,0,0),c(hdata3[13:14,2],0,0,0,0),c(hdata3[15:16,2],0,0,0,0),c(hdata3[17:18,2],0,0,0,0),c(hdata3[19:20,2],0,0,0,0),c(hdata3[21:24,2],0,0),c(hdata3[25:26,2],0,0,0,0),c(hdata3[27:29,2],0,0,0),c(hdata3[30:35,2]),c(hdata3[36:37,2],0,0,0,0),c(hdata3[39:40,2],0,0,0,0),c(hdata3[41:44,2],0,0),c(hdata3[45:50,2]),c(hdata3[51:54,2],0,0),c(hdata3[56:57,2],0,0,0,0),c(hdata3[58:61,2],0,0))

heatmap_data <- NULL
for ( i in 1:6){
	heatmap_data <- cbind(heatmap_data, DeltaC=hdata.tb[,i], PHE=hdata.tb2[,i], tRNA=hdata.tb3[,i])

}

aas <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Phe","Pro","Ser","Thr","Tyr","Val")
rownames(heatmap_data) <- aas
maxDeltaC <- max(c(heatmap_data[,1],heatmap_data[,4],heatmap_data[,7],heatmap_data[,10],heatmap_data[,13],heatmap_data[,16]))


hdata.tb4 <- rbind(c(as.character(hdata3[1:4,1]),"",""),c(as.character(hdata3[5:10,1])),c(as.character(hdata3[11:12,1]),"","","",""),c(as.character(hdata3[13:14,1]),"","","",""),c(as.character(hdata3[15:16,1]),"","","",""),c(as.character(hdata3[17:18,1]),"","","",""),c(as.character(hdata3[19:20,1]),"","","",""),c(as.character(hdata3[21:24,1]),"",""),c(as.character(hdata3[25:26,1]),"","","",""),c(as.character(hdata3[27:29,1]),"","",""),c(as.character(hdata3[30:35,1])),c(as.character(hdata3[36:37,1]),"","","",""),c(as.character(hdata3[39:40,1]),"","","",""),c(as.character(hdata3[41:44,1]),"",""),c(as.character(hdata3[45:50,1])),c(as.character(hdata3[51:54,1]),"",""),c(as.character(hdata3[56:57,1]),"","","",""),c(as.character(hdata3[58:61,1]),"",""))

nullvec <- rep("",18)

heatmap_data2 <- NULL
for ( i in 1:6){
	heatmap_data2 <- cbind(heatmap_data2, Codon=hdata.tb4[,i], nullvec, tRNA=hdata.tb3[,i])

}



heatmap_data[,c(2,5,8,11,14,17)] <- heatmap_data[,c(2,5,8,11,14,17)]+1
heatmap_data[,c(3,6,9,12,15,18)] <- heatmap_data[,c(3,6,9,12,15,18)]+2


#heatmaps
library('pheatmap')

paletteLength <- 50
myColor <- c(colorRampPalette(c("firebrick3", "white"))(ceiling(paletteLength/4)), colorRampPalette(c("white", "navy"))((floor(paletteLength/4))), colorRampPalette(c("white", "navy"))(ceiling(paletteLength/4)),colorRampPalette(c("white", "navy"))(floor(paletteLength/4)))
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(heatmap_data), -0.00000001, length.out=ceiling(paletteLength/4)), 0, seq(0.00000001, maxDeltaC, length.out=floor(paletteLength/4)), seq(1,1.9999, length.out=ceiling(paletteLength/4)),seq(2,max(heatmap_data), length.out=floor(paletteLength/4)))

name1 <- paste(name,".pdf",sep="")

pheatmap(heatmap_data,filename = name1, cellwidth = 25, cellheight = 15, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = FALSE, color=myColor, breaks=myBreaks, display_numbers = matrix(ifelse(heatmap_data2 > 0.99, heatmap_data2, ""), nrow(heatmap_data)), number_color="black", gaps_col=c(3,6,9,12,15))

#matrix(ifelse(heatmap_data > 0.99, heatmap_data, ""), nrow(heatmap_data))

#pheatmap(hdata.tb,filename = "deltaCheatmap.pdf", cellwidth = 25, cellheight = 15, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = FALSE, color=myColor, breaks=myBreaks)
#pheatmap(hdata.tb,filename = "pheheatmap.pdf", cellwidth = 25, cellheight = 15, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = FALSE, color=colorRampPalette(c("white", "navy"))(50))
#pheatmap(hdata.tb3,filename = "trnaheatmap.pdf", cellwidth = 25, cellheight = 15, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = FALSE, color=colorRampPalette(c("white", "navy"))(50))


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
chead <- read.table("../c_head.txt", header=FALSE, sep=" ")

aminoacids <- c("C","D","E","F","H","K","N","Q","Y","I","A","G","P","T","V","L","R","S")

mydata <- NULL
for (i in aminoacids){
	codons <- chead[which(chead[,1] == i),2]
	for (j in codons){
		mycodons <- aatable[which(aatable$Codon == codons[which(codons == j)]),]
		mycodons <- cbind(mycodons, AA=rep(i,length(mycodons[,1])))
		mydata <- rbind(mydata,mycodons)	
	}
}

hdata3 <- hdata3[match(chead$V2,hdata3$codon),]
trnas <- cbind(chead,trna=hdata3[,2])

custom.sort <- function(x){
	thdbase <- substring(x[,3], 3)
	fstbase <- substring(x[,3], 1,1)
	x <- cbind(x,fstbase,thdbase)
	ord <- c("U","C","A","G")
	factor(x[,6],levels=ord)
	x[,7] <- factor(x[,7],levels=ord)
	x[order(x[,6],x[,7]),]
}


#library(directlabels)
library(ggpubr)
library('ggplot2')
library('ggrepel')
library('ggthemes')

for (i in aminoacids) {
	myplot <- mydata[which(mydata$AA==i),]
	myplot <- custom.sort(myplot)
	myplot$Freq <- as.numeric(myplot$Freq)
	myplot$X <- as.numeric(as.character(myplot$X))
	if (i == "C" | i == "D" | i == "F" | i == "H" | i == "N" | i == "Y") {mycolors <- c("#0f4c8aff","#8a0000ff") } else if (i == "E" | i == "K" | i == "Q") {mycolors <- c("#808000ff","#cf2e7dff") } else if (i == "I") {mycolors <- c("#0f4c8aff","#8a0000ff","#008000ff") } else if (i == "A" | i == "G" | i == "P" | i == "T" | i == "V") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff") } else if (i == "L") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff","#808000ff", "#cf2e7dff") }  else if (i == "R") {mycolors <- c("#808000ff","#cf2e7dff","#61b8ffff", "#cc0000ff","#008000ff", "#800080ff") } else if (i == "S") {mycolors <- c("#0f4c8aff","#8a0000ff","#61b8ffff", "#cc0000ff","#808000ff", "#800080ff") }
	myplot[,3] <- factor(myplot[,3],unique(myplot[,3]) , ordered = TRUE)
	vlineint <-	(Clast + Clast/30)
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
	legend.position = c(0.5, 0.5),
    #legend.position="bottom",
	legend.title=element_blank(),
	legend.text = element_text(colour="black", size=16, face="bold"),
	legend.background =element_rect(fill=alpha("white",0.75), size=0.5, linetype="solid"),
	legend.box = "horizontal",
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank())
    nam <- paste("p", i, sep = "")
    assign(nam, p1)
}

custom.sort3 <- function(x){
	thdbase <- substring(x[,2], 3)
	fstbase <- substring(x[,2], 1,1)
	x <- cbind(x,fstbase,thdbase)
	ord <- c("U","C","A","G")
	factor(x[,4],levels=ord)
	x[,5] <- factor(x[,5],levels=ord)
	x[order(x[,4],x[,5]),]
}




for (i in aminoacids) {
	if (i == "C" | i == "D" | i == "F" | i == "H" | i == "N" | i == "Y") {mycolors <- c("#0f4c8aff","#8a0000ff") } else if (i == "E" | i == "K" | i == "Q") {mycolors <- c("#808000ff","#cf2e7dff") } else if (i == "I") {mycolors <- c("#0f4c8aff","#8a0000ff","#008000ff") } else if (i == "A" | i == "G" | i == "P" | i == "T" | i == "V") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff") } else if (i == "L") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff","#808000ff", "#cf2e7dff") }  else if (i == "R") {mycolors <- c("#808000ff","#cf2e7dff","#61b8ffff", "#cc0000ff","#008000ff", "#800080ff") } else if (i == "S") {mycolors <- c("#0f4c8aff","#8a0000ff","#61b8ffff", "#cc0000ff","#808000ff", "#800080ff") }
	mytrnas <- trnas[which(trnas$V1==i),]
	mytrnas <- custom.sort3(mytrnas)
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


#axis.text.x = element_text(face="bold", size=12, angle=90)
#matrix(ifelse(mytrna > 0.99, mytrnas, "")
a <- ggarrange(pC,pD,pE,pF,pH,pK,pN,pQ,pY,pI,pA,pG,pP,pT,pV,pL,pR,pS, ncol = 18, align = "hv",labels = aminoacids, hjust = -2.5, vjust = 2, font.label = list(size = 20, color = "black", face = "bold"))
b <- ggarrange(ptC,ptD,ptE,ptF,ptH,ptK,ptN,ptQ,ptY,ptI,ptA,ptG,ptP,ptT,ptV,ptL,ptR,ptS, ncol = 18, align = "hv",hjust = -2.5, vjust = 2)

name2 <- paste(name,"AAfreq_trnas.svg",sep="_")
svg(filename=name2, width=24, height=3.5, pointsize=16)
ggarrange(a,b, nrow = 2, align = "hv",hjust = -2.5, vjust = 2, heights = c(2.5,0.8))
dev.off()


if (length(grep("C.{1,2}[_ ]woPHE$",rownames(fdata))) != 0) {

	fdatasp <- fdata[ grep("C.{1,2}[_ ]woPHE$",rownames(fdata)), ]
	fdatasp <- data.frame(fdatasp)
	fdatasp <- data.frame(cbind(SET=rownames(fdatasp),X=1:length(fdatasp[,1]),fdatasp))
	aatablesp <- gather(fdatasp, Codon, Freq, -SET, -X)


	mydatasp <- NULL
	for (i in aminoacids){
		codons <- chead[which(chead[,1] == i),2]
		for (j in codons){
			mycodonssp <- aatablesp[which(aatablesp$Codon == codons[which(codons == j)]),]
			mycodonssp <- cbind(mycodonssp, AA=rep(i,length(mycodonssp[,1])))
			mydatasp <- rbind(mydatasp,mycodonssp)	
		}
	}


	for (i in aminoacids) {
		myplot <- mydatasp[which(mydatasp$AA==i),]
		myplot <- custom.sort(myplot)
		myplot$Freq <- as.numeric(myplot$Freq)
		myplot$X <- as.numeric(as.character(myplot$X))
		if (i == "C" | i == "D" | i == "F" | i == "H" | i == "N" | i == "Y") {mycolors <- c("#0f4c8aff","#8a0000ff") } else if (i == "E" | i == "K" | i == "Q") {mycolors <- c("#808000ff","#cf2e7dff") } else if (i == "I") {mycolors <- c("#0f4c8aff","#8a0000ff","#008000ff") } else if (i == "A" | i == "G" | i == "P" | i == "T" | i == "V") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff") } else if (i == "L") {mycolors <- c("#61b8ffff","#cc0000ff","#008000ff", "#800080ff","#808000ff", "#cf2e7dff") }  else if (i == "R") {mycolors <- c("#808000ff","#cf2e7dff","#61b8ffff", "#cc0000ff","#008000ff", "#800080ff") } else if (i == "S") {mycolors <- c("#0f4c8aff","#8a0000ff","#61b8ffff", "#cc0000ff","#808000ff", "#800080ff") }
		myplot[,3] <- factor(myplot[,3],unique(myplot[,3]) , ordered = TRUE)
		vlineint <-	(Clast + Clast/30)

		p1 <- ggplot(myplot, aes(x=X,y=Freq, group=Codon, color=Codon)) +
		#geom_point(aes(color=Codon),show.legend=F)  +
		geom_line(aes(color=Codon), size=1)+

		#geom_text_repel(aes(label=Codon, fontface=2), size = 2, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
		xlab("") +
		ylab("f") + 
		#geom_dl(aes(label = Codon), method = list(dl.trans(x = x - .4, y = y - .5), "last.points")) +
		#geom_dl(aes(label = Codon), method = list(dl.trans(x = x - .4, y = y - .5), "first.points")) +
		scale_color_manual(values=mycolors) +
		theme_calc()+
		#scale_color_manual(values=c("dodgerblue3"))+
		guides(colour = guide_legend(override.aes = list(size=3))) +
		theme(axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.position = c(0.5, 0.5),
		#legend.position="bottom",
		legend.title=element_blank(),
		legend.text = element_text(colour="black", size=16, face="bold"),
		legend.background =element_rect(fill=alpha("white",0.75), size=0.5, linetype="solid"),
		legend.box = "horizontal",
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank())
		nam <- paste("psp", i, sep = "")
		assign(nam, p1)
	}

	woPHE <- ggarrange(pspC,pspD,pspE,pspF,pspH,pspK,pspN,pspQ,pspY,pspI,pspA,pspG,pspP,pspT,pspV,pspL,pspR,pspS, ncol = 18, align = "hv",labels = aminoacids, hjust = -2.5, vjust = 2)

	name3 <- paste(name,"AAfreq_trnas_woPHE.svg",sep="_")
	p <- ggarrange(a,woPHE,b, nrow = 3, align = "hv",hjust = -2.5, vjust = 2, heights = c(1.5,1.5,0.6))
	svg(filename=name3, width=30, height=5, pointsize=16)
	print(p)
	dev.off()

}



