#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

#Calculate tAI vs dist on R using complete genome and adds U:U interaction


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
destfile="modal_ws.counts"    
if (!file.exists(destfile)) {
stop("Counts data file does not exist", call.=FALSE)
}

destfile="s_opts_DCBS_GNM2.txt"    
if (!file.exists(destfile)) {
    stop("S_BCBS data file does not exist", call.=FALSE)
}


destfile="dist.txt"    
if (!file.exists(destfile)) {
    stop("Distance data file does not exist", call.=FALSE)
}

#load data
hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")
sopt_DCBS <- read.table("s_opts_DCBS_GNM2.txt", header=TRUE, sep="\t", dec=".")
sopt_DCBS <- sopt_DCBS[order(sopt_DCBS[,11], decreasing=T),]
Dist <- read.table("dist.txt", header=TRUE, sep="\t", dec=",")
Dist[,2] <- apply(as.matrix(Dist[,2]),2,function(x) x-x[1])
corder = read.table('../codOrder.txt', header=F, sep=" ")

#data for tai calcualtion - modals


f.counts <- read.table("modal_ws.counts", header=FALSE, row.names = 1, sep=" ")
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

Sopt_DCBS <- as.numeric(sopt_DCBS[1,1:10])


get.ws3 <- function (tRNA, s = NULL, sking) 
{
    if (is.null(s)) 
        s <- c(0, 0, 0, 0, 0.41, 0.28, 0.9999, 0.68, 0.89, 0.5)
    p = 1 - s
    W = NULL
	box4 <- c(5,17,21,29,33,37,49,53,61)
	box2 <- c(1,9,13,25,41,45,57)
	codonWorder <- c(1,2,3,4,9,10,11,12,13,14,15,16,25,26,27,28,41,42,43,44,45,46,47,48,57,58,59,60,5,6,7,8,17,18,19,20,21,22,23,24,29,30,31,32,33,34,35,36,37,38,39,40,49,50,51,52,53,54,55,56,61,62,63,64)
   
	for (i in box2) W = c(W, p[1] * tRNA[i] + p[5] * tRNA[i + 1], p[2] * tRNA[i + 1] + p[6] * tRNA[i], p[3] * tRNA[i + 2] + p[7] * tRNA[i], p[4] * tRNA[i + 3] + p[8] * tRNA[i + 2])
	for (i in box4) W = c(W, p[1] * tRNA[i] + p[5] * tRNA[i + 1] + p[10] * tRNA[i + 2], p[2] * tRNA[i + 1] + p[6] * tRNA[i], p[3] * tRNA[i + 2] + p[7] * tRNA[i], p[4] * tRNA[i + 3] + p[8] * tRNA[i + 2])
	
	W <- W[match(1:64,codonWorder)]
	W[36] = p[4] * tRNA[36]

    if (sking == 1) 
        W[35] = p[9]
    W = W[-c(11, 12, 15, 36)] 
	w = W/max(W)
    if (sum(w == 0) > 0) {
        ws <- w[w != 0]
        gm <- exp(sum(log(ws))/length(ws))
        w[w == 0] = gm
    }

    return(w)
}

if (trna$trna[35] == 0){sk <- 1 } else {sk <- 0}

ws_DCBS <- get.ws3(tRNA=trna$trna,s=Sopt_DCBS, sking=sk)
tai_DCBS <- get.tai(f.o[,-33], ws_DCBS)


fdata <- as.data.frame(cbind(SET=as.character(rownames(f.counts)),tai_DCBS))
rownames(fdata)[rownames(fdata) %like% "^single"] = "Single"


#plots

fdata.C <- cbind(fdata[rownames(fdata) %like% "C[0-9]*$|^PHE|^Single", ], Class=c(rep("C", length(fdata[rownames(fdata) %like% "C[0-9]*$", 1])),"PHE","Single"))
fdata.woPHE <- cbind(fdata[rownames(fdata) %like% "C[0-9]*_*woPHE$|^PHE", ], Class=c(rep("C_woPHE", length(fdata[rownames(fdata) %like% "C[0-9]*_woPHE$", 1])),rep("PHE", length(fdata[rownames(fdata) %like% "^PHE", 1]))))
fdata.C <- as.data.frame(fdata.C)
C = cbind(fdata.C[rownames(fdata.C) %like% "C",], Dist = Dist[,2])
PHE = fdata[rownames(fdata) %like% "^PHE",]
Single = fdata[rownames(fdata) %like% "^Single",]
PHE$tai_DCBS <- as.numeric(as.character(PHE$tai_DCBS))
Single$tai_DCBS <- as.numeric(as.character(Single$tai_DCBS))
fdata$tai_DCBS <- as.numeric(as.character(fdata$tai_DCBS))
C$tai_DCBS <- as.numeric(as.character(C$tai_DCBS))


yminimo <- min(fdata$tai_DCBS,PHE$tai_DCBS[1] ) - (.1 * (max(fdata$tai_DCBS,PHE$tai_DCBS[1]) -min(fdata$tai_DCBS,PHE$tai_DCBS[1])))
ymaximo <- max(fdata$tai_DCBS,PHE$tai_DCBS[1]) + (.1 * (max(fdata$tai_DCBS,PHE$tai_DCBS[1]) -min(fdata$tai_DCBS,PHE$tai_DCBS[1])))

svg(filename="tAi_modal_GNM2.svg", width=8, height=6, pointsize=10)


ggplot(C, aes(x=Dist, y=tai_DCBS, group=Class, color=Class)) +
geom_point(show.legend = FALSE, size=3) +
geom_line(aes(linetype=Class, color=Class), show.legend = FALSE, size=1) +
annotate("text", min(Dist$Dist), PHE$tai_DCBS[1], vjust = 1.5, label = "PHE", size=6) +
geom_hline(yintercept = PHE$tai_DCBS[1]  , linetype="dashed", color = "red", size=1) +
annotate("text", min(Dist$Dist), Single$tai_DCBS[1], vjust = 1.5, hjust = -0.1, label = "Singletons", size=6) +
geom_hline(yintercept = Single$tai_DCBS[1]  , linetype="dashed", color = "Blue", size=1) +
geom_text_repel(aes(label=SET, fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
scale_color_manual(values=c('black'))+
ylim(yminimo, ymaximo) +
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

	class.select <- C.woPHE[1,3]

    svg(filename="tAi_modal_woPHE_GNM2.svg", width=8, height=6, pointsize=10)

    p <- ggplot(CwoPHE, aes(x=Dist, y=tai_DCBS, group=Class, color=Class)) +
    geom_line(aes(linetype=Class, color=Class, alpha=Class), show.legend = FALSE, size=1) +
    geom_hline(yintercept = PHE$tai_DCBS[1], linetype="dashed", color = "red", , size=1) +
    annotate("text", min(Dist$Dist), PHE$tai_DCBS[1], vjust = 1.5, label = "PHE", size=6) +
    annotate("text", min(Dist$Dist), Single$tai_DCBS[1], vjust = 1.5, hjust = -0.1, label = "Singletons", size=6) +
    geom_hline(yintercept = Single$tai_DCBS[1]  , linetype="dashed", color = "Blue", size=1) +
    geom_point(aes(alpha=Class), show.legend = FALSE, size=3) +
    geom_text_repel(aes(label=ifelse(Class == as.character(class.select), NA, as.character(CwoPHE$SET)), fontface=2), size = 6, point.padding = 1 , show.legend = FALSE, arrow = arrow(length = unit(0.01, 'npc'))) +
    scale_color_manual(values=c('black','firebrick'))+
    scale_alpha_manual(values=c(1,0.5))+
    ylim(yminimo, ymaximo) +
    xlab("Evolutionary distance") +
    ylab("m-tAI") + 
    theme_calc() +
	theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
	) 
    print(p)

    dev.off()

}

