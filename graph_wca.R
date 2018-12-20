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
library('ade4')

setwd(args[1])
hdata <- read.table("../c_head.txt", header=FALSE, sep=" ")
data = read.table(args[2], header=F, sep=" ") #check field separator
data = data[,1:61]
colnames(data) <- c("ID", "SET", as.character(unlist(hdata[,2])))
row.names(data) <- data[,1]

wca <- dudi.coa(t(data[,3:61]), scan = FALSE, nf = 2)
facaa <- as.factor(hdata[,1])
scua <- wca(wca, facaa, scan = FALSE, nf = 2)

#scua$co -> genes in columns, transposed matrix


#*****************************************************************
#****************      coa of counts    **************************
#*****************************************************************
#
#coagenes <- cbind(wca$co, data[,2]) 
#names(coagenes) <- c("C1","C2","Class") 
#
#genes <- coagenes[ which(coagenes$Class=='genes'), ]
#modal <- coagenes[ grep("C.{1,2}$|^PHE",row.names(coagenes)), ]
#woPHE <- coagenes[ grep("woPHE",row.names(coagenes)), ]
#
#svg(filename="CA_genes.svg", width=5, height=4, pointsize=10)
#
#ggplot (coagenes, aes( y = C2, x = C1)) +
#geom_point(data = genes, color='gray', size = 1, alpha = I(0.3)) +
#geom_point(data = modal, fill = "dodgerblue3", size=2, colour="black", shape=21, stroke = 0.5) +
#ggtitle("Correspondance Analysis RSCU") + xlab("C1") + ylab("C2") + labs(fill = "Legend")  +
#geom_text_repel(data = modal, aes(label=row.names(modal), fontface=2), size = 2, point.padding = 1 , arrow = arrow(length = unit(0.01, 'npc'))) +
#theme_calc()
#
#dev.off()

wcagenes <- cbind(scua$co, data[,2])
names(wcagenes) <- c("C1","C2","Class")

genes <- wcagenes[ which(wcagenes$Class=='genes'), ]
modals <- wcagenes[ which(wcagenes$Class=='modal'), ]
modal <- modals[ grep("C.{1,2}$|^PHE",row.names(modals)), ]
woPHE <- modals[ grep("woPHE",row.names(modals)), ]

svg(filename="WCA_genes.svg", width=5, height=4, pointsize=10)

ggplot (wcagenes, aes( y = C2, x = C1)) +
geom_point(data = genes, color='gray', size = 1, alpha = I(0.3)) +
geom_point(data = modal, fill = "dodgerblue3", size=3, colour="black", shape=21, stroke = 0.5) +
ggtitle("Correspondance Analysis RSCU") + xlab("C1") + ylab("C2") + labs(fill = "Legend")  +
geom_text_repel(data = modal, aes(label=row.names(modal), fontface=2), size = 3, point.padding = 1 , arrow = arrow(length = unit(0.01, 'npc'))) +
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 

dev.off()

svg(filename="WCA_genes2.svg", width=5, height=4, pointsize=10)

ggplot (wcagenes, aes( y = C2, x = C1)) +
geom_point(data = modal, fill = "dodgerblue3", size=3, colour="black", shape=21, stroke = 0.5) +
ggtitle("Correspondance Analysis RSCU") + xlab("C1") + ylab("C2") + labs(fill = "Legend")  +
geom_text_repel(data = modal, aes(label=row.names(modal), fontface=2), size = 3, point.padding = 1 , arrow = arrow(length = unit(0.01, 'npc'))) +
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 

dev.off()


#*************************************************************
#****************        woPHE       **************************
#*************************************************************

if (length(woPHE$Class) != 0) {

    svg(filename="WCA_genes_woPHE.svg", width=5, height=4, pointsize=10)

    p <- ggplot (wcagenes, aes( y = C2, x = C1)) +
    geom_point(data = genes, color='gray', size = 1, alpha = I(0.3)) +
    geom_point(data = modal, fill = "dodgerblue3", size=3, colour="black", shape=21, stroke = 0.5) +
    geom_point(data = woPHE, fill = "firebrick3", size=3, colour="black", shape=21, stroke = 0.5) +
    ggtitle("Correspondance Analysis RSCU") + xlab("C1") + ylab("C2") + labs(fill = "Legend")  +
    geom_text_repel(data = rbind(modal,woPHE), aes(label=row.names(rbind(modal,woPHE)), fontface=2), size = 3, point.padding = 1 , arrow = arrow(length = unit(0.01, 'npc'))) +
    theme_calc() +
	theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 
    print(p)
    dev.off()

    svg(filename="WCA_genes2_woPHE.svg", width=5, height=4, pointsize=10)

    p <-ggplot (wcagenes, aes( y = C2, x = C1)) +
    geom_point(data = modal, fill = "dodgerblue3", size=3, colour="black", shape=21, stroke = 0.5) +
    geom_point(data = woPHE, fill = "firebrick3", size=3, colour="black", shape=21, stroke = 0.5) +
    ggtitle("Correspondance Analysis RSCU") + xlab("C1") + ylab("C2") + labs(fill = "Legend")  +
    geom_text_repel(data = rbind(modal,woPHE), aes(label=row.names(rbind(modal,woPHE)), fontface=2), size = 3, point.padding = 1 , arrow = arrow(length = unit(0.01, 'npc'))) +
    theme_calc() +
	theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
	) 
    print(p)
    dev.off()
}



