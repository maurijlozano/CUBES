#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

#Calcular el indice en R


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
## program...
library('ggplot2')
library('ggrepel')
library('ggthemes')
library('data.table')

hdata <- read.table("freq_head.txt", header=FALSE, sep=" ")

setwd(args[1])

#data: FUC vs Cores
#test if file exist
destfile="fmdata.txt"    
if (!file.exists(destfile)) {
stop("Frequence data file does not exist", call.=FALSE)
}

destfile="rgf.txt"    
if (!file.exists(destfile)) {
    stop("RGF data file does not exist", call.=FALSE)
}

destfile="dist.txt"    
if (!file.exists(destfile)) {
    stop("Distance data file does not exist", call.=FALSE)
}


fdata <- read.table("fmdata.txt", header=FALSE, sep=" ")
rgf <- read.table("rgf.txt", header=TRUE, sep="\t", dec=",")
Dist <- read.table("dist.txt", header=TRUE, sep="\t", dec=",")
Dist[,2] <- apply(as.matrix(Dist[,2]),2,function(x) x-x[1])

fdata <- cbind(fdata[,1], fdata[,3:61])
colnames(fdata) <- as.character(unlist(hdata[1,]))


#producto de matrices para calcular el indice Suma.prod... 
fdata = cbind(fdata, as.matrix(fdata[,2:60])%*%as.matrix(rgf$RGF))
colnames(fdata)[61] <- "Index"

fdata.C <- cbind(fdata[fdata$SET %like% "C[0-9]*$|^PHE", ], Class=c(rep("C", length(fdata[fdata$SET %like% "C[0-9]*$", 1])),rep("PHE", length(fdata[fdata$SET %like% "^PHE", 1]))))
fdata.woPHE <- cbind(fdata[fdata$SET %like% "C[0-9]*_woPHE$|^PHE", ], Class=c(rep("C_woPHE", length(fdata[fdata$SET %like% "C[0-9]*_woPHE$", 1])),rep("PHE", length(fdata[fdata$SET %like% "^PHE", 1]))))
C = cbind(fdata.C[fdata.C$SET %like% "C",], Dist = Dist[,2])
PHE = fdata[fdata$SET %like% "^PHE",]


svg(filename="IvsD_modal.svg", width=8, height=6, pointsize=10)


ggplot(C, aes(x=Dist, y=Index, group=Class, color=Class)) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
geom_hline(yintercept = PHE$Index[1], linetype="dashed", color = "red", , size=1) +
annotate("text", min(Dist$Dist), PHE$Index[1], vjust = 1.5, label = "PHE") +
geom_point(show.legend = FALSE) +
geom_text_repel(aes(label=SET, fontface=2), size = 4, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
ylim(min(fdata$Index)-0.25, max(fdata$Index)+0.25) +
ggtitle("Adaptation Index vs Distance") + 
xlab("Evolutionary distance") +
ylab("I") + 
theme_calc() +
theme(axis.text.x = element_text(size=13),
	  axis.text.y = element_text(size=13),
	  axis.title.x = element_text(size=15, face="bold"),
	  axis.title.y = element_text(size=15, face="bold")
) 

dev.off()

if (length(fdata.woPHE$SET) > 1){
    C.woPHE = cbind(fdata.woPHE[fdata.woPHE$SET %like% "C",], Dist = Dist[,2])
    CwoPHE = rbind(C,C.woPHE)

    svg(filename="IvsD_modal_woPHE.svg", width=8, height=6, pointsize=10)

    p <- ggplot(CwoPHE, aes(x=Dist, y=Index, group=Class, color=Class)) +
    geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
    geom_hline(yintercept = PHE$Index[1], linetype="dashed", color = "red", , size=1) +
    annotate("text", min(Dist$Dist), PHE$Index[1], vjust = 1.5, label = "PHE") +
    geom_point(show.legend = FALSE) +
    geom_text_repel(aes(label=SET, fontface=2), size = 3, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
    scale_color_manual(values=c('black','firebrick'))+
    ylim(min(CwoPHE$Index,PHE$Index)-0.25, max(CwoPHE$Index,PHE$Index)+0.25) +
    ggtitle("Adaptation Index vs Distance") + 
    xlab("Evolutionary distance") +
    ylab("I") + 
    theme_calc() +
	theme(axis.text.x = element_text(size=13),
	  axis.text.y = element_text(size=13),
	  axis.title.x = element_text(size=15, face="bold"),
	  axis.title.y = element_text(size=15, face="bold")
	) 
    print(p)
    dev.off()
}


#Writes table with all the data

write.table(fdata, file = "results_modal.txt", sep = "\t", dec = ".", row.names = TRUE)

#write.table(data, file = "data.txt", append = FALSE, quote = TRUE, sep = "\t",
#                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#                 col.names = TRUE, qmethod = c("escape", "double"),
#                 fileEncoding = "")
