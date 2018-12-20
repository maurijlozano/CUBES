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
destfile="fdata.txt"    
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



fdata <- read.table("fdata.txt", header=FALSE, sep=" ")
rgf <- read.table("rgf.txt", header=TRUE, sep="\t", dec=",")
Dist <- read.table("dist.txt", header=TRUE, sep="\t", dec=",")
Dist <- Dist[order(Dist$SET),]
Dist[,2] <- apply(as.matrix(Dist[,2]),2,function(x) x-x[1])

colnames(fdata) <- as.character(unlist(hdata[1,]))
fdata[,1] <- as.character(fdata[,1])
fdata[fdata[,1] %like% "single",1] <- rep("Singletons",length(fdata[fdata[,1] %like% "single",1]))


#Matrix product to calculate the adaptation Index... 
fdata = cbind(fdata, as.matrix(fdata[,2:60])%*%as.matrix(rgf$RGF))
colnames(fdata)[61] <- "Index"


results <- aggregate(fdata[,2:61], by=list(fdata$SET), function(x) c(mean = mean(x), sd = sd(x), SE = sd(x)/sqrt(length(x))))
colnames(results)[1] <- "SET"

C.res <- cbind(results[results$SET %like% "C[0-9]*$|^PHE", grep("SET|Index",names(results))], Class=c(rep("C", length(results[results$SET %like% "C[0-9]*$", 1])),rep("PHE", length(results[results$SET %like% "^PHE", 1]))))

woPHE.res <- cbind(results[results$SET %like% "C[0-9]*.woPHE$|^PHE", grep("SET|Index",names(results))], Class=c(rep("C_woPHE", length(results[results$SET %like% "C[0-9]*.woPHE$", 1])),rep("PHE", length(results[results$SET %like% "^PHE", 1]))))

C = cbind(C.res[C.res$SET %like% "C",], Dist = Dist[,2])
PHE = results[results$SET %like% "^PHE",]


svg(filename="IvsD1.svg", width=8, height=6, pointsize=10)


ggplot(C, aes(x=Dist, y=Index[,1], group=Class, color=Class)) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
geom_hline(yintercept = PHE$Index[1], linetype="dashed", color = "red", , size=1) +
annotate("text", min(Dist$Dist), PHE$Index[1], vjust = 1.5, label = "PHE") +
geom_point(show.legend = FALSE) +
geom_text_repel(aes(label=SET, fontface=2), size = 4, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
ylim(min(results$Index[,1])-0.25, max(results$Index[,1])+0.25) +
ggtitle("Adaptation Index vs Distance") + 
xlab("Evolutionary distance") +
ylab("I") + 
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=13, face="bold"),
	  axis.title.y = element_text(size=13, face="bold")
) 

dev.off()


svg(filename="IvsD2_SD.svg", width=8, height=6, pointsize=10)

ggplot(C, aes(x=Dist, y=Index[,1], group=Class, color=Class)) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
geom_point(show.legend = FALSE) +
geom_errorbar(aes(ymin=C$Index[,1]-C$Index[,2], ymax=C$Index[,1]+C$Index[,2]), width=(max(Dist$Dist)-min(Dist$Dist))/100) +
geom_hline(yintercept = PHE$Index[1], linetype="dashed", color = "red", , size=1) +
annotate("text", min(Dist$Dist), PHE$Index[1], vjust = 1.5, label = "PHE") +
geom_text_repel(aes(label=SET, fontface=2), size = 4, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
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


svg(filename="IvsD3_SE.svg", width=8, height=6, pointsize=10)

ggplot(C, aes(x=Dist, y=Index[,1], group=Class, color=Class)) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
geom_point(show.legend = FALSE) +
geom_errorbar(aes(ymin=C$Index[,1]-C$Index[,3], ymax=C$Index[,1]+C$Index[,3]), width=(max(Dist$Dist)-min(Dist$Dist))/100) +
geom_hline(yintercept = PHE$Index[1], linetype="dashed", color = "red", , size=1) +
annotate("text", min(Dist$Dist), PHE$Index[1], vjust = 1.5, label = "PHE") +
geom_text_repel(aes(label=SET, fontface=2), size = 4, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
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

if (length(woPHE.res$SET) > 1) {
    woPHE = cbind(woPHE.res[woPHE.res$SET %like% "C",], Dist = Dist[,2])
    CwoPHE = rbind(C,woPHE)

    svg(filename="IvsD_woPHE.svg", width=8, height=6, pointsize=10)

    p <- ggplot(CwoPHE, aes(x=Dist, y=Index[,1], group=Class, color=Class)) +
    geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
    geom_hline(yintercept = PHE$Index[1], linetype="dashed", color = "red", , size=1) +
    annotate("text", min(Dist$Dist), PHE$Index[1], vjust = 1.5, label = "PHE") +
    geom_point(show.legend = FALSE) +
    geom_text_repel(aes(label=SET, fontface=2), size = 3, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
    scale_color_manual(values=c('black','firebrick3'))+
    ylim(min(CwoPHE$Index[,1],PHE$Index[,1])-0.25, max(CwoPHE$Index[,1],PHE$Index[,1])+0.25) +
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

write.table(results, file = "results.txt", sep = "\t", dec = ".", row.names = TRUE)


#write.table(data, file = "data.txt", append = FALSE, quote = TRUE, sep = "\t",
#                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#                 col.names = TRUE, qmethod = c("escape", "double"),
#                 fileEncoding = "")
