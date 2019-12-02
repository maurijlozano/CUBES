#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 


setwd(args[1])
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

#plots
library('ggplot2')
library('ggrepel')
library('ggthemes')


data = read.table("genes_GNM.coa", header=TRUE, sep=" ") #check field separator
cores <- data[ grep("C.{1,2}$",data$class), ]

data = read.table(args[2], header=TRUE, row.names=1, sep=" ")[,1:2] #check field separator
if (cores[1,3] > cores[length(cores$class),3]){data$Axis1 <- -data$Axis1}
if (cores[1,4] > cores[length(cores$class),4]){data$Axis2 <- -data$Axis2}

codonsM <- data[grep(paste(maxChangeCodon, collapse='|'), rownames(data)),] 
dataU <- subset(codonsM, grepl("UCU|CUU|CCU|CGU|ACU|GUU|GCU|GGU",rownames(codonsM)))
dataC <- subset(codonsM, grepl("UUC|UAC|UGC|CAC|AUC|AAC|GAC",rownames(codonsM)))

varpc = read.table("coa_var_axis1_2_GNM.txt", header=TRUE, sep=" ") #check field separator
varpc <- varpc*100
xlabtext <- paste("C1 (",varpc[1,1]," %)", sep="")
ylabtext <- paste("C2 (",varpc[1,2]," %)", sep="")

svg(filename="CA_codons_GNM_rv.svg", width=5, height=4, pointsize=10)

ggplot (data, aes( y = Axis2, x = Axis1)) +
#draw outer most circle - C and U bias
geom_point(color = "grey" , size=3) +
geom_point(data=codonsM, fill = "tan3", size=3, colour="black", shape=21, stroke = .5) +
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

