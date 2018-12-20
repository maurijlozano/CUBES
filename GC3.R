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
destfile="fcdata.txt"    
if (!file.exists(destfile)) {
stop("Frequence data file does not exist", call.=FALSE)
}



destfile="dist.txt"    
if (!file.exists(destfile)) {
    stop("Distance data file does not exist", call.=FALSE)
}



Dist <- read.table("dist.txt", header=TRUE, sep="\t", dec=",")
Dist <- Dist[order(Dist$SET),]
Dist[,2] <- apply(as.matrix(Dist[,2]),2,function(x) x-x[1])


#*******************************************
#Calculate %GC on the third condon base: NN@

fcdata <- read.table("fcdata.txt", header=FALSE, sep=" ")
colnames(fcdata) <- as.character(unlist(hdata[1,]))
fcdata[,1] <- as.character(fcdata[,1])
fcdata[fcdata[,1] %like% "[Ss]ingle",1] <- rep("Singletons",length(fcdata[fcdata[,1] %like% "single",1]))


GC3 <- fcdata[grep("..G|..C",names(fcdata))]
fcdata.sum <- cbind.data.frame(SET=fcdata[,"SET"], total = apply(as.matrix(fcdata[,-1]),1, sum))
GC3.sum <- cbind.data.frame(SET=fcdata[,"SET"], GC3 = apply(as.matrix(GC3),1,sum))

GC3 <- cbind.data.frame(SET=fcdata[,"SET"], GC3 = ((GC3.sum[,-1])/(fcdata.sum[,-1]))*100 )

GC3.res <- aggregate(GC3[,2], by=list(GC3[,1]), function(x) c(mean = mean(x), sd = sd(x), SE = sd(x)/sqrt(length(x))))
colnames(GC3.res)[1] <- "SET"


GC3.c <- cbind.data.frame(GC3.res[GC3.res$SET %like% "C[0-9]*$|^PHE",], Class=c(rep("C", length(GC3.res[GC3.res$SET %like% "C[0-9]*$", 1])),rep("PHE", length(GC3.res[GC3.res$SET %like% "^PHE", 1]))))
GC3.woPHE <- cbind.data.frame(GC3.res[GC3.res$SET %like% "C[0-9]*_woPHE$|^PHE",], Class=c(rep("C_woPHE", length(GC3.res[GC3.res$SET %like% "C[0-9]*_woPHE$", 1])),rep("PHE", length(GC3.res[GC3.res$SET %like% "^PHE", 1]))))
GC3.SINGLE = GC3.res[GC3.res$SET %like% "^Singletons$",]

GC3.C = cbind(GC3.c[GC3.c$SET %like% "C[0-9]*$",], Dist = Dist[,2])
GC3.PHE = GC3.c[GC3.c$SET %like% "^PHE$",]


svg(filename="GC3vsD.svg", width=8, height=6, pointsize=10)

ggplot(GC3.C, aes(x=Dist, y=GC3.C$x[,1], group=Class, color=Class)) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
geom_hline(yintercept = GC3.PHE$x[,1], linetype="dashed", color = "red", size=1) +
annotate("text", min(Dist$Dist), GC3.PHE$x[,1], vjust = 1.5, label = "PHE", size=6) +
geom_point(show.legend = FALSE, size=3) +
geom_text_repel(aes(label=SET, fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
ylim(min(GC3.c$x[,1])-1, max(GC3.c$x[,1])+1) +
ggtitle("% GC3 vs Distance") + 
xlab("Evolutionary distance") +
ylab("% GC3") + 
theme_calc() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()

svg(filename="GC3vsD_SD.svg", width=8, height=6, pointsize=10)

ggplot(GC3.C, aes(x=Dist, y=GC3.C$x[,1], group=Class, color=Class)) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
geom_hline(yintercept = GC3.PHE$x[,1], linetype="dashed", color = "red", , size=1) +
annotate("text", min(Dist$Dist), GC3.PHE$x[,1], vjust = 1.5, label = "PHE", size=6) +
geom_point(show.legend = FALSE, size=3) +
geom_errorbar(aes(ymin=GC3.C$x[,1]-GC3.C$x[,2], ymax=GC3.C$x[,1]+GC3.C$x[,2]), width=(max(Dist$Dist)-min(Dist$Dist))/100, show.legend = FALSE) +
geom_text_repel(aes(label=SET, fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
ggtitle("% GC3 vs Distance") + 
xlab("Evolutionary distance") +
ylab("% GC3") + 
theme_calc() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()

svg(filename="GC3vsD_SE.svg", width=8, height=6, pointsize=10)

ggplot(GC3.C, aes(x=Dist, y=GC3.C$x[,1], group=Class, color=Class)) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
geom_hline(yintercept = GC3.PHE$x[,1], linetype="dashed", color = "red", , size=1) +
annotate("text", min(Dist$Dist), GC3.PHE$x[,1], vjust = 1.5, label = "PHE", size=6) +
geom_point(show.legend = FALSE, size=3) +
geom_errorbar(aes(ymin=GC3.C$x[,1]-GC3.C$x[,3], ymax=GC3.C$x[,1]+GC3.C$x[,3]), width=(max(Dist$Dist)-min(Dist$Dist))/100, show.legend = FALSE) +
geom_text_repel(aes(label=SET, fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
ggtitle("% GC3 vs Distance") + 
xlab("Evolutionary distance") +
ylab("% GC3") + 
theme_calc() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()


if (length(GC3.woPHE$SET) > 1) {

	GC3.woPHE = cbind.data.frame(GC3.woPHE[GC3.woPHE$SET %like% "C[0-9]*",], Dist = Dist[,2])
	GC3CwoPHE <- rbind.data.frame(GC3.C,GC3.woPHE)
	GC3CwoPHE$Dist <- as.numeric(as.character(GC3CwoPHE$Dist))
	GC3CwoPHE$x[,1] <- as.numeric(as.character(GC3CwoPHE$x[,1]))
	class.select <- GC3.woPHE[1,3]
    svg(filename="GC3vsD_woPHE.svg", width=8, height=6, pointsize=10)

    p <- ggplot(GC3CwoPHE, aes(x=Dist, y=GC3CwoPHE$x[,1], group=Class, color=Class)) +
    geom_line(aes(linetype=Class, color=Class, alpha=Class), show.legend = FALSE, size=1) +
    geom_hline(yintercept = GC3.PHE$x[,1], linetype="dashed", color = "red", , size=1) +
    annotate("text", min(Dist$Dist), GC3.PHE$x[,1], vjust = 1.5, label = "PHE", size=6) +
    geom_point(show.legend = FALSE, size=3, aes(alpha=Class)) +
    geom_text_repel(aes(label=ifelse(Class == as.character(class.select), NA, as.character(GC3CwoPHE$SET)), fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
    scale_color_manual(values=c('black','firebrick'))+
	scale_alpha_manual(values=c(1,0.5))+
    ylim(min(GC3CwoPHE$x[,1],GC3.PHE$x[,1])-0.25, max(GC3CwoPHE$x[,1],GC3.PHE$x[,1])+0.25) +
    ggtitle("% GC3 vs Distance") + 
    xlab("Evolutionary distance") +
    ylab("% GC3") + 
    theme_calc() +
	theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
	) 
    print(p)
    dev.off()
}




if (length(GC3.SINGLE$SET) == 1) {

    svg(filename="GC3vsD_ws.svg", width=8, height=6, pointsize=10)
	p <- ggplot(GC3.C, aes(x=Dist, y=GC3.C$x[,1], group=Class, color=Class)) +
	geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
	geom_hline(yintercept = GC3.PHE$x[,1], linetype="dashed", color = "red", , size=1) +
	geom_hline(yintercept = GC3.SINGLE$x[,1], linetype="dashed", color = "blue", , size=1) +
	annotate("text", min(Dist$Dist), GC3.SINGLE$x[,1], vjust = 1.5, hjust = .2, label = "Singletons", size=6) +
	annotate("text", min(Dist$Dist), GC3.PHE$x[,1], vjust = 1.5, hjust = .2, label = "PHE", size=6) +
	geom_point(show.legend = FALSE, size=3) +
	geom_text_repel(aes(label=SET, fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
	scale_color_manual(values=c('black'))+
	ylim(min(GC3.c$x[,1],GC3.SINGLE$x[,1])-1, max(GC3.c$x[,1],GC3.SINGLE$x[,1])+1) +
	ggtitle("% GC3 vs Distance") + 
	xlab("Evolutionary distance") +
	ylab("% GC3") + 
	theme_calc() +
	theme(axis.text.x = element_text(size=14),
		  axis.text.y = element_text(size=14),
		  axis.title.x = element_text(size=18, face="bold"),
		  axis.title.y = element_text(size=18, face="bold")
		) 
    print(p)
	dev.off()
}


#Writes table with all the data
write.table(GC3.res, file = "resultsGC3.txt", sep = "\t", dec = ".", row.names = TRUE)

#write.table(data, file = "data.txt", append = FALSE, quote = TRUE, sep = "\t",
#                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#                 col.names = TRUE, qmethod = c("escape", "double"),
#                 fileEncoding = "")
