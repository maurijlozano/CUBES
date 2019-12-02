#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## program...
library('ggplot2')
library('ggrepel')
library('ggthemes')
library('gtools')

setwd(args[1])
data = read.table(args[2], header=TRUE, sep=" ") #check field separator


cores <- data[ grep("C.{1,2}$",data$class), ]
if (cores[1,3] > cores[length(cores$class),3]){data$Axis1 <- -data$Axis1}
if (cores[1,4] > cores[length(cores$class),4]){data$Axis2 <- -data$Axis2}

genes <- data[ which(data$class=='genes'), ]
modal <- data[ grep("C.{1,2}$|^Single|^PHE",data$class), ]
modal$class <- as.character(modal$class)
modal$class <- factor(modal$class, levels = mixedsort(modal$class))
single_phe <- data[ grep("^Single|^PHE",data$class), ]
single_phe$class <- as.character(single_phe$class)
single_phe$class[2]<-"Singletons"

modal.labels <- rbind(modal[1,],modal[round((dim(modal)[1]-2)/2),],modal[(dim(modal)[1]-2),],single_phe)


varpc = read.table("coa_var_axis1_2_GNM.txt", header=TRUE, sep=" ") #check field separator
varpc <- varpc*100
xlabtext <- paste("C1 (",varpc[1,1]," %)", sep="")
ylabtext <- paste("C2 (",varpc[1,2]," %)", sep="")

ncores <- dim(modal)[1] - 2
colfunc <- colorRampPalette(c("dodgerblue3", "firebrick1"))
colores<-c(colfunc(ncores),"red","blue")


svg(filename="CA_genes_GNM_rv.svg", width=5, height=4, pointsize=10)

ggplot (data, aes( y = Axis2, x = Axis1)) +
geom_point(data = genes, color='gray', size = 2, alpha = 0.8) +
geom_point(data = modal, aes(fill = class), size=4, colour="black", shape=21, stroke = 0.5, show.legend=F) +
xlab(xlabtext) + ylab(ylabtext) + #labs(fill = "Core-genome Set")  +
geom_text_repel(data = modal.labels, aes(label=class, fontface=2), size = 4, arrow = NULL,point.padding=0.2, min.segment.length=1, segment.alpha=0.5) +
scale_fill_manual(values=colores) +
theme_calc() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=16, face="bold"),
	  axis.title.y = element_text(size=16, face="bold"),
	  legend.title = element_text(size=13, face="bold"),
	  legend.text = element_text(size=12)
) 

dev.off()


