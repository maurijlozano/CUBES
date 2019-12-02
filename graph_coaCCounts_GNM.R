#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

#If coreGenes are not correctly show is because the locus tags doesn't coincide between PHE and the other cores.

setwd(args[1])
hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")
hdata <- cbind(hdata, "AUG","UGG") 
data = read.table("GNM.counts", header=F, sep=" ") 
colnames(data) <- c(as.character(unlist(hdata[,1:62])))
rownames(data)<-data[,1]
#remove met and trp
data <- data[,c(-1,-61,-62)]
coreGenes <- read.table("coreGenes.txt", header=T, sep="\t") 
n <- dim(coreGenes)[2]
colnames(coreGenes)[grep("[Ss][iI][nN][gG][lL]",colnames(coreGenes))] <- "Singletons"
colfunc <- colorRampPalette(c("dodgerblue3", "firebrick1"))
colores<-c(colfunc(n-2),"red","blue")

#libraries
library('ggplot2')
library('ggrepel')
library('ggthemes')
library('gtools')
#CA
library("FactoMineR")
library("factoextra")

res.ca <- CA(data, graph = FALSE)

#variance percent
varpc <- get_eigenvalue(res.ca)[1:2,2]
xlabtext <- paste("C1 (",round(varpc[1],2)," %)", sep="")
ylabtext <- paste("C2 (",round(varpc[2],2)," %)", sep="")

#get genes and codons...
genes <- data.frame(get_ca_row(res.ca)$coord[,1:2])
genes <- cbind(genes, SET=rep("GNM",dim(genes)[1]))
codons <- data.frame(get_ca_col(res.ca)$coord[,1:2])

genes1 <- genes

#fviz_ca_biplot(res.ca, repel = TRUE, label=F)

for(i in 1:dim(coreGenes)[2]){
	genesi <- na.omit(genes[match(coreGenes[,i],rownames(data)),1:2])
	genesi <- cbind(genesi, SET=rep(colnames(coreGenes)[i],dim(genesi)[1]))
	genes <- rbind(genes, genesi)
}



genes2 <- NULL
for(i in 1:dim(coreGenes)[2]){
	genesi <- na.omit(genes[match(coreGenes[,i],rownames(data)),1:2])
	genesi <- cbind(genesi, SET=rep(colnames(coreGenes)[i],dim(genesi)[1]))
	genes2 <- rbind(genes2, genesi)
}

#order.factor.levels
genes2$SET <- as.character(genes2$SET)
genes2$SET <- factor(genes2$SET, levels = unique(mixedsort(genes2$SET)))

# centroid.. also with -> aggregate(genes2[-3], by=list(genes2$SET),mean)
centroids <- aggregate(cbind(genes2$Dim.1,genes2$Dim.2)~SET,genes2,mean)
colnames(centroids) <- c("SET","Dim.1","Dim.2")

modal.labels <- rbind(centroids[1,],centroids[round((dim(centroids)[1]-2)/2),],centroids[(dim(centroids)[1]-2):dim(centroids)[1],])

##lims
#correctionX <- (max(centroids$Dim.1) - min(centroids$Dim.1))/10
#correctionY <- (max(centroids$Dim.2) - min(centroids$Dim.2))/10
#minx <- min(centroids$Dim.1) - correctionX
#maxx <- max(centroids$Dim.1) + correctionX
#miny <- min(centroids$Dim.2) - correctionY
#maxy <- max(centroids$Dim.2) + correctionY


#plots
svg(filename="CA_genesCCounts_GNM.svg", width=5, height=4, pointsize=10)

	ggplot(genes, aes(y = Dim.2, x = Dim.1)) +
	geom_point(data=genes1, color="grey", alpha = 0.8, size=2) +
#	geom_point(data=genes2, aes(color=SET), size=2, show.legend=F, alpha=0.3) +
#	stat_ellipse(aes(color=SET), show.legend=F) +
	geom_point(data=centroids, aes(fill=SET), size=4, shape=21, color="black", show.legend=F) +
	geom_text_repel(data=modal.labels, aes(label=SET, fontface=2), size = 4, arrow = NULL, point.padding=0.2, min.segment.length=1, segment.alpha=0.5) +
#	geom_label_repel(data=modal.labels, aes(label=SET, fontface=2), size = 2.5, seed=1234, fill="white", alpha=0.6, arrow = NULL) +
#	geom_label_repel(data=modal.labels, aes(label=SET, fontface=2), size = 2.5, seed=1234, fill=NA, arrow = NULL) +
	scale_fill_manual(values=colores) +
	scale_color_manual(values=colores) +
#	xlim(minx,maxx) + ylim(miny,maxy) +
	xlab(xlabtext) + ylab(ylabtext) +
#	labs(fill = "Core-SET") +
	theme_calc() +
	theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=16, face="bold"),
	  axis.title.y = element_text(size=16, face="bold"),
	  legend.title = element_text(size=13, face="bold"),
	  legend.text = element_text(size=12)
	) 

dev.off()



#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
wophe <- as.numeric(args[2])

if (wophe == 0) {
	genes <- genes1
	genes2 <- NULL
	corewophe <- coreGenes[,which(colnames(coreGenes) != "PHE")]
	corephe <- coreGenes[,which(colnames(coreGenes) == "PHE")]
	corephe <- corephe[grep("^$",corephe,invert=T)]
	for(i in 1:dim(corewophe)[2]){
		#setdiff(a$x, b$y)  elements in a not in b
		genesi <- setdiff(corewophe[,i],corephe)
		genesi <- genesi[grep("^$",genesi,invert=T)]
		genesi <- na.omit(genes[match(genesi,rownames(data)),1:2])
		genesi <- cbind(genesi, SET=rep(colnames(corewophe)[i],dim(genesi)[1]))
		genes <- rbind(genes, genesi)
		genes2 <- rbind(genes2, genesi)
	}

	#order.factor.levels
	genes2$SET <- as.character(genes2$SET)
	genes2$SET <- factor(genes2$SET, levels = unique(mixedsort(genes2$SET)))

	# centroid.. also with -> aggregate(genes2[-3], by=list(genes2$SET),mean)
	centroidsPHE <- centroids[grep("PHE",centroids$SET),]
	centroids <- aggregate(cbind(genes2$Dim.1,genes2$Dim.2)~SET,genes2,mean)
	colnames(centroids) <- c("SET","Dim.1","Dim.2")
	centroids <- rbind(centroids,centroidsPHE)
	colnames(centroids) <- c("SET","Dim.1","Dim.2")
	colores<-c(colfunc(n-2),"blue","red")

	modal.labels <- rbind(centroids[1,],centroids[round((dim(centroids)[1]-2)/2),],centroids[(dim(centroids)[1]-2):dim(centroids)[1],])

	#plots
		p <- ggplot(genes, aes(y = Dim.2, x = Dim.1)) +
		geom_point(data=genes1, color="grey", alpha = 0.8, size=2) +
	#	geom_point(data=genes2, aes(color=SET), size=2, show.legend=F, alpha=0.3) +
	#	stat_ellipse(aes(color=SET), show.legend=F) +
		geom_point(data=centroids, aes(fill=SET), size=4, shape=21, color="black", show.legend=F) +
		geom_text_repel(data=modal.labels, aes(label=SET, fontface=2), size = 4, arrow = NULL, point.padding=0.2, min.segment.length=1, segment.alpha=0.5) +
	#	geom_label_repel(data=modal.labels, aes(label=SET, fontface=2), size = 2.5, seed=1234, fill="white", alpha=0.6, arrow = NULL) +
	#	geom_label_repel(data=modal.labels, aes(label=SET, fontface=2), size = 2.5, seed=1234, fill=NA, arrow = NULL) +
		scale_fill_manual(values=colores) +
		scale_color_manual(values=colores) +
	#	xlim(minx,maxx) + ylim(miny,maxy) +
		xlab(xlabtext) + ylab(ylabtext) +
	#	labs(fill = "Core-SET") +
		theme_calc() +
		theme(axis.text.x = element_text(size=14),
		  axis.text.y = element_text(size=14),
		  axis.title.x = element_text(size=16, face="bold"),
		  axis.title.y = element_text(size=16, face="bold"),
		  legend.title = element_text(size=13, face="bold"),
		  legend.text = element_text(size=12)
		) 

	svg(filename="CA_genesCCounts_GNM_woPHE.svg", width=5, height=4, pointsize=10)
		print(p)
	dev.off()

}



#***************************************************************************************************#***************************************************************************************************
#***************************************************************************************************
#codons plot
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

fdata1 <- fdata[ grep("C.{1,2}$|^PHE|^[Ss][Ii][Nn][Gg][Ll][Ee]",rownames(fdata)), ]
fdata1 <- data.frame(fdata1)
for (i in 1:ncol(fdata1)){fdata1[,i]<-as.numeric(as.character(fdata1[,i]))}

#it is phe-c1 but I didn't change the variable name..
cn_c1 <- fdata1[grep("[Pp][Hh][Ee]",rownames(fdata1)),] - fdata1[1,]

aminoacids <- c("C","D","E","F","H","K","N","Q","Y","I","A","G","P","T","V","L","R","S")
chead <- read.table("../c_head.txt", header=FALSE, sep=" ", stringsAsFactors = FALSE)
colnames(chead) <- c("AA","Codon")

maxChangeCodon <- NULL
for (i in aminoacids){
	maxC <- 0
	for (j in chead[grep(i,chead[,1]),]$Codon){
		maxCN <- cn_c1[,grep(j,colnames(cn_c1))]
		if (maxCN > maxC){
			maxC <- maxCN
			maxCcod <- j
		}
	}
	maxChangeCodon <- c(maxChangeCodon, maxCcod)
}

codonsM <- codons[grep(paste(maxChangeCodon, collapse='|'), rownames(codons)),] 
dataU <- subset(codonsM, grepl("UCU|CUU|CCU|CGU|ACU|GUU|GCU|GGU",rownames(codonsM)))
dataC <- subset(codonsM, grepl("UUC|UAC|UGC|CAC|AUC|AAC|GAC",rownames(codonsM)))


svg(filename="CA_codonsCCounts_GNM.svg", width=5, height=4, pointsize=10)

ggplot (codons, aes( y = Dim.2, x = Dim.1)) +
#draw outer most circle - C and U bias
geom_point(color = "grey" , size=3) +
geom_point(data=codonsM, fill = "tan3", size=4, colour="black", shape=21, stroke = .5) +
geom_point(data=dataU, fill = "#61b8ffff", size=4, colour="black", shape=21, stroke = .5) +
geom_point(data=dataC, fill = "#8a0000ff", size=4, colour="black", shape=21, stroke = .5) +
xlab(xlabtext) + ylab(ylabtext) + labs(fill = "Legend") +
geom_text_repel(data=codonsM, aes(label=rownames(codonsM), fontface=2), size = 4, arrow = NULL) +
labs(fill = "W (Inner circle)") +
labs(color = "Trna copies (Outer circle)") +
theme_calc() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=16, face="bold"),
	  axis.title.y = element_text(size=16, face="bold")
) 

dev.off()


