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

setwd(args[1])
data = read.table(args[2], header=TRUE, sep=" ") #check field separator

cores <- data[ grep("C.{1,2}$",data$class), ]
if (cores[1,3] > cores[length(cores$class),3]){data$Axis1 <- -data$Axis1}
if (cores[1,4] > cores[length(cores$class),4]){data$Axis2 <- -data$Axis2}
genes <- data[ which(data$class=='genes'), ]
modal <- data[ grep("C.{1,2}$|^PHE|^Single",data$class), ]
woPHE <- data[ grep("woPHE",data$class), ]

varpc = read.table("coa_var_axis1_2_GNM.txt", header=TRUE, sep=" ") #check field separator
varpc <- varpc*100
xlabtext <- paste("C1 (",varpc[1,1]," %)", sep="")
ylabtext <- paste("C2 (",varpc[1,2]," %)", sep="")


png(filename = "CA_genes_GNM.png", width = 1000, height = 800, units = "px", pointsize = 9,
    bg = "white", res = 230, type = "cairo-png")

ggplot (data, aes( y = Axis2, x = Axis1)) +
geom_point(data = genes, color='gray', size = 1, alpha = I(0.3)) +
geom_point(data = modal, fill = "dodgerblue3", size=3, colour="black", shape=21, stroke = 0.5) +
xlab(xlabtext) + ylab(ylabtext) + labs(fill = "Legend")  +
geom_text_repel(data = modal, aes(label=class, fontface=2), size = 3, point.padding = 1 , arrow = arrow(length = unit(0.01, 'npc'))) +
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 

dev.off()

