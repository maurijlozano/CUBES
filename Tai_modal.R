#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

#Calculate tAI vs dist on R


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

setwd(args[1])

#data: FUC vs Cores
#test if file exist
destfile="modal.counts"    
if (!file.exists(destfile)) {
stop("Counts data file does not exist", call.=FALSE)
}

destfile="s_opts.txt"    
if (!file.exists(destfile)) {
    stop("S_NcAj data file does not exist", call.=FALSE)
}

destfile="s_opts_DCBS.txt"    
if (!file.exists(destfile)) {
    stop("S_BCBS data file does not exist", call.=FALSE)
}


destfile="dist.txt"    
if (!file.exists(destfile)) {
    stop("Distance data file does not exist", call.=FALSE)
}

#load data
hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")
sopt <- read.table("s_opts.txt", header=TRUE, sep="\t", dec=".")
sopt_DCBS <- read.table("s_opts_DCBS.txt", header=TRUE, sep="\t", dec=".")
sopt_DCBS <- sopt_DCBS[order(sopt_DCBS[,10], decreasing=T),]
Dist <- read.table("dist.txt", header=TRUE, sep="\t", dec=",")
Dist[,2] <- apply(as.matrix(Dist[,2]),2,function(x) x-x[1])
corder = read.table('../codOrder.txt', header=F, sep=" ")

#data for tai calcualtion - modals


f.counts <- read.table("modal.counts", header=FALSE, row.names = 1, sep=" ")
colnames(f.counts) <- c(as.character(unlist(hdata[,2:60])), "AUG", "UGG")
corder2 <- corder[-c(11,12,15),]
t<-t(f.counts)
t <- t[match(corder2$V2,row.names(t)),]
fdata <- t(t)


#trnas
corder = read.table('../codOrder.txt', header=F, sep=" ")
trna = read.table('trna.txt', header=T, row.names = NULL, sep=" ")
trna <- trna[match(corder$V2,trna$codon),]

corder2 <- corder[-c(11,12,15),]
t<-t(fdata)
t <- t[match(corder2$V2,row.names(t)),]
f.o <- t(t)


library('tAI')
Sopt <- as.numeric(sopt[1,1:9])
Sopt_DCBS <- as.numeric(sopt_DCBS[1,1:9])
if (trna$trna[35] == 0){sk <- 1 } else {sk <- 0}
ws <- get.ws(tRNA=trna$trna,s=Sopt, sking=sk)
tai <- get.tai(f.o[,-33], ws)

ws_DCBS <- get.ws(tRNA=trna$trna,s=Sopt_DCBS, sking=sk)
tai_DCBS <- get.tai(f.o[,-33], ws_DCBS)

ws_data <- cbind(round(ws,2),round(ws_DCBS,2),trna[-c(11,12,15,33),])
#fdata <- cbind(tAI_NAj=tai, tAI_DCBS=tai_DCBS)
#plot(fdata)

fdata <- as.data.frame(cbind(SET=as.character(rownames(f.counts)),tai_DCBS))



#plots

fdata.C <- cbind(fdata[rownames(fdata) %like% "C[0-9]*$|^PHE", ], Class=c(rep("C", length(fdata[rownames(fdata) %like% "C[0-9]*$", 1])),rep("PHE", length(fdata[rownames(fdata) %like% "^PHE", 1]))))
fdata.woPHE <- cbind(fdata[rownames(fdata) %like% "C[0-9]*_woPHE$|^PHE", ], Class=c(rep("C_woPHE", length(fdata[rownames(fdata) %like% "C[0-9]*_woPHE$", 1])),rep("PHE", length(fdata[rownames(fdata) %like% "^PHE", 1]))))
fdata.C <- as.data.frame(fdata.C)
C = cbind(fdata.C[rownames(fdata.C) %like% "C",], Dist = Dist[,2])
PHE = fdata[rownames(fdata) %like% "^PHE",]
PHE$tai_DCBS <- as.numeric(as.character(PHE$tai_DCBS))
fdata$tai_DCBS <- as.numeric(as.character(fdata$tai_DCBS))
C$tai_DCBS <- as.numeric(as.character(C$tai_DCBS))

yminimo <- min(C$tai_DCBS,PHE$tai_DCBS[1] ) - (.1 * (max(C$tai_DCBS,PHE$tai_DCBS[1]) -min(C$tai_DCBS,PHE$tai_DCBS[1])))
ymaximo <- max(C$tai_DCBS,PHE$tai_DCBS[1]) + (.1 * (max(C$tai_DCBS,PHE$tai_DCBS[1]) -min(C$tai_DCBS,PHE$tai_DCBS[1])))


svg(filename="tAi_modal.svg", width=8, height=6, pointsize=10)


ggplot(C, aes(x=Dist, y=tai_DCBS, group=Class, color=Class)) +
geom_point(show.legend = FALSE, size=3) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
annotate("text", min(Dist$Dist), PHE$tai_DCBS[1], vjust = 1.5, label = "PHE", size=6) +
geom_hline(yintercept = PHE$tai_DCBS[1]  , linetype="dashed", color = "red", size=1) +
geom_text_repel(aes(label=SET, fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
ylim(yminimo,ymaximo) +
xlab("Evolutionary distance") +
ylab("m-tAI") + 
theme_calc() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()

if (length(fdata.woPHE$SET) > 1){
    C.woPHE = cbind(fdata.woPHE[fdata.woPHE$SET %like% "C",], Dist = Dist[,2])
    CwoPHE = rbind(C,C.woPHE)
    CwoPHE$tai_DCBS <- as.numeric(as.character(CwoPHE$tai_DCBS))

	yminimo <- min(CwoPHE$tai_DCBS,PHE$tai_DCBS) - (.1 * (max(CwoPHE$tai_DCBS,PHE$tai_DCBS) - min(CwoPHE$tai_DCBS,PHE$tai_DCBS)))
	ymaximo <- max(CwoPHE$tai_DCBS,PHE$tai_DCBS) + (0.1 * (max(CwoPHE$tai_DCBS,PHE$tai_DCBS)-min(CwoPHE$tai_DCBS,PHE$tai_DCBS)))
	class.select <- C.woPHE[1,3]

    svg(filename="tAi_modal_woPHE.svg", width=8, height=6, pointsize=10)

    p <- ggplot(CwoPHE, aes(x=Dist, y=tai_DCBS, group=Class, color=Class)) +
    geom_line(aes(linetype=Class, color=Class, alpha=Class), show.legend = FALSE, size=1) +
    geom_hline(yintercept = PHE$tai_DCBS[1], linetype="dashed", color = "red", , size=1) +
    annotate("text", min(Dist$Dist), PHE$tai_DCBS[1], vjust = 1.5, label = "PHE", size=6) +
    geom_point(aes(alpha=Class), show.legend = FALSE, size=3) +
    geom_text_repel(aes(label=ifelse(Class == as.character(class.select), NA, as.character(CwoPHE$SET)), fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
    scale_color_manual(values=c('black','firebrick'))+
    scale_alpha_manual(values=c(1,0.5))+
    ylim(yminimo, ymaximo) +
    xlab("Evolutionary distance") +
    ylab("m-tAI") + 
    theme_calc()+
	theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
) 
    print(p)

    dev.off()

}


#write.table(fdata, file = "results_modal.txt", sep = "\t", dec = ".", row.names = TRUE)

#write.table(data, file = "data.txt", append = FALSE, quote = TRUE, sep = "\t",
#                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#                 col.names = TRUE, qmethod = c("escape", "double"),
#                 fileEncoding = "")
